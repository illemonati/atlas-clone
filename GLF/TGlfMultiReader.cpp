/*
 * TGlfMultireader.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#include "TGlfMultiReader.h"
#include <algorithm>
#include <cstdint>
#include <exception>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include "TBgzWriter.h"
#include "coretools/Files/TGzWriter.h"
#include "coretools/Files/TWriter.h"
#include "coretools/Main/TError.h"
#include "genometools/GenotypeTypes.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/strongTypes.h"
#include "coretools/Types/weakTypes.h"

namespace GLF {
using coretools::instances::logfile;
using coretools::instances::parameters;
using genometools::HighPrecisionPhredIntProbability;

using coretools::str::toString;

namespace impl {

TGlfReader * nextChr(const std::vector<TGlfReader *>& ps, bool onlyData) {
	if (onlyData) {
		return *std::min_element(ps.begin(), ps.end(), [](auto p1, auto p2) {
			if (p1->refId() < p2->refId()) return true;
			if (p1->refId() == p2->refId())
				return p1->position() < p2->position();
			else
				return false;
		});
	} else {
		return *std::min_element(ps.begin(), ps.end(), [](auto p1, auto p2) { return p1->refId() < p2->refId(); });
	}
}

void _checkChromosomeInfo(const TGlfChromosome & _curChr, const std::vector<TGlfReader *>& ps) {
	// check that all files share the same chromosome info
	for (const auto p: ps) {
		TGlfChromosome *chr = p->pointerToChr(_curChr.refId());
		if (chr) {
			if (chr->name() != _curChr.name())
				UERROR("Chrosomome names differ between files '", ps[0]->name(), "' and '",
					   p->name(), "': '", _curChr.name(), "' != '", chr->name(), + "'!");
			if (chr->length() != _curChr.length())
				UERROR("Chrosomome lengths differ between files '", ps[0]->name(), "' and '",
					   p->name(), "': '", _curChr.length(), "' != '", chr->length(), "'!");
		} else {
			logfile().list(p->name(), " does not have contig ", _curChr.name(), ", considering all data empty.");
		}
	}
};

} // namespace impl

TMultiGLFDataOneAllelicCombination fill(coretools::TConstView<GLF::TMultiGLFDataSample> samples,
		  genometools::AllelicCombination alleleicCombination) {
	TMultiGLFDataOneAllelicCombination storage;
	storage.reserve(samples.size());
	for (const auto &s : samples) {
		if (!s.hasData()) {
			continue;
		}
		if (s.isHaploid())
			storage.emplace_back(s[first(alleleicCombination)],
								 s[second(alleleicCombination)]);
		else
			storage.emplace_back(s[homoFirst(alleleicCombination)],
								 s[het(alleleicCombination)],
								 s[homoSecond(alleleicCombination)]);
	}
	return storage;
};

uint32_t totalDepth(const TMultiGLFData &samples) noexcept {
	return std::accumulate(samples.begin(), samples.end(), uint32_t{}, [](auto tot, auto s) {return tot + s.depth();});
};

//----------------------------------------------------
// TGlfMultiReader
//----------------------------------------------------
TGlfMultiReader::TGlfMultiReader() {
	_minDepth = parameters().getParameterWithDefault<size_t>("minDepth", 0);
	if (_minDepth > 0) logfile().list("Will only keep sites with depth >= " + toString(_minDepth) + ".");

	_windowSize = parameters().getParameterWithDefault<size_t>("window", 100000);
	if (_windowSize == 0) UERROR("Window size must be at least 1!");
};

void TGlfMultiReader::_openGLFs() {
	_numGLFs     = _GLFNames.size();
	_GLFIsActive.resize(_numGLFs, false);

	// open files
	_GLFs.clear();
	_GLFs.reserve(_numGLFs);
	_readersOpened = true;
	logfile().startIndent("Opening " + toString(_numGLFs) + " GLF files:");
	for (const auto &name : _GLFNames) {
		logfile().listFlush("Opening GLF '" + name + "' ...");
		_GLFs.emplace_back(name);
		logfile().done();
	}
	logfile().endIndent();
	_setAllInactive();
};

void TGlfMultiReader::openGLFs(const std::vector<std::string> &FileNames) {
	_GLFNames = FileNames;
	_openGLFs();
};

void TGlfMultiReader::openGLFs() {
	using namespace coretools::str;
	const auto parameter = parameters().getParameter<std::string>("glf");
	// assume that GLF file names are given in a file if string does not contain ".gz"
	if (!stringContains(parameter, ".gz")) {
		logfile().list("Reading glf input names from file '" + parameter + "'");
		std::ifstream in(parameter.c_str());
		if (!in) UERROR("Failed to open file '", parameter, "'!");

		// read file
		while (in.good() && !in.eof()) {
			std::string line;
			std::getline(in, line);
			std::vector<std::string> vec;
			fillContainerFromStringWhiteSpace(line, vec, true);
			// skip empty lines
			if (vec.size() > 0) { _GLFNames.push_back(vec[0]); }
		}
		in.close();
	} else {
		parameters().fillParameterIntoContainer("glf", _GLFNames, ',');
	}
	_openGLFs();
};

void TGlfMultiReader::addReference(const std::string &FastaFile) {
	fastaReader.open(FastaFile);
};

//-------------------------------------
// set active / inactive
//-------------------------------------
int TGlfMultiReader::_getGLFIndexFromName(const std::string &name) const {
	const auto res = std::find(_GLFNames.cbegin(), _GLFNames.cend(), name);
	if (res == _GLFNames.cend()) UERROR("GLF with name '", name, "' not in TGlfMultiReader!");
	return std::distance(_GLFNames.begin(), res);
};

void TGlfMultiReader::_setActive(size_t index) {
	if (index >= _numGLFs) UERROR("Index out of range in TGlfMultiReader::setActive(const int index)!");
	if (!_GLFIsActive[index]) {
		_GLFIsActive[index] = true;
		_activeGLFs.push_back(&_GLFs[index]);
	}
};

void TGlfMultiReader::_setAllInactive() {
	_GLFIsActive.assign(_numGLFs, false);
	_activeGLFs.clear();
};

void TGlfMultiReader::_prepareParsing() {
	for (TGlfReader *it : _activeGLFs) {
		it->rewind();
		it->readNext();
	}
	// where to start?
	_jumpToNextPosition();
};

bool TGlfMultiReader::_jumpToNextPosition() {
	auto min      = impl::nextChr(_activeGLFs, _minSamplesWithData);
	_curChr       = *min->curChr();
	_curRefId     = min->refId();
	_windowStart  = _minSamplesWithData > 0 ? min->position() : 0; // 0 if not jump to position
	impl::_checkChromosomeInfo(_curChr, _activeGLFs);

	return !min->eof();
};

void TGlfMultiReader::setActive(int index) {
	_setAllInactive();
	_setActive(index);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string &name) {
	setActive(_getGLFIndexFromName(name));
};

void TGlfMultiReader::setActive(int index1, int index2) {
	_setAllInactive();
	_setActive(index1);
	_setActive(index2);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string &name1, const std::string &name2) {
	setActive(_getGLFIndexFromName(name1), _getGLFIndexFromName(name2));
};

void TGlfMultiReader::setActive(const std::vector<int> &indexes) {
	_setAllInactive();
	for (const auto i : indexes) _setActive(i);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::vector<std::string> &names) {
	_setAllInactive();
	for (const auto &name : names) {
		_setActive(_getGLFIndexFromName(name));
	}
	_prepareParsing();
};

void TGlfMultiReader::setAllActive() {
	_setAllInactive();
	for (size_t i = 0; i < _numGLFs; ++i) _setActive(i);
	_prepareParsing();
};

//-------------------------------------
// Looping over active files
//-------------------------------------
bool TGlfMultiReader::_moveToNextChromosome() {
	// increment chromosome ref_char id
	++_curRefId;

	// advance all active files behind
	bool allFilesReachedEnd = true;
	for (TGlfReader *it : _activeGLFs) {
		while (!it->eof() && it->refId() < _curRefId) it->jumpToNextChr();
		if (!it->eof()) allFilesReachedEnd = false;
	}

	// check if we reached end of all files
	if (allFilesReachedEnd) return false;

	// get name and length from first active file not at end
	_jumpToNextPosition();
	return true;
};

std::vector<size_t> TGlfMultiReader::readWindow() {
	_windowStart += _dataWindow.size();
	if (_windowStart >= _curChr.length()) {
		if(!_moveToNextChromosome()) return {};
	}
	const size_t N    = std::min(_windowSize, _curChr.length() - _windowStart);
	_dataWindow.assign(N, {});
	_numActive.assign(N, 0);

	bool allEOF=true;
	for (auto reader : _activeGLFs) {
		// find first data in window
		while (!reader->eof() && (reader->refId() == _curRefId) && (reader->position() < _windowStart)) {
			reader->readNext();
		}
		if (!reader->eof()) allEOF = false;

		// fill everything as noData
		for (size_t iW = 0; iW < N; ++iW) {
			if (!reader->eof() && (reader->refId() == _curRefId) && reader->position() - _windowStart == iW) {
				_dataWindow[iW].emplace_back(reader->genotypeLikelihoodsGLF(), reader->depth());
				++_numActive[iW];
				reader->readNext();
			} else {
				_dataWindow[iW].emplace_back(reader->pointerToChr(_curRefId)->isHaploid());
			}
		}
	}
	if (allEOF) return {};
	std::vector<size_t> ids;
	ids.reserve(_dataWindow.size());
	for (size_t i = 0; i < _dataWindow.size(); ++i) {
		if (_numActive[i] >= _minSamplesWithData) {
			ids.push_back(i);
		}
	}
	if (ids.empty()) return readWindow();
	return ids;
}

std::vector<std::string> TGlfMultiReader::namesOfActiveFiles() const {
	std::vector<std::string> vec;

	// sample names are file names without glf ending
	for (TGlfReader *it : _activeGLFs) { vec.emplace_back(it->name()); }
	return vec;
};

std::vector<std::string> TGlfMultiReader::sampleNamesOfActiveFiles() const {
	std::vector<std::string> vec;

	// sample names are file names without glf ending
	for (TGlfReader *it : _activeGLFs) { vec.emplace_back(coretools::str::readBeforeLast(it->name(), ".glf")); }
	return vec;
};

}; // end namespace GLF
