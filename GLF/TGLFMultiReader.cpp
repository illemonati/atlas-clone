/*
 * TGlfMultireader.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#include "TGLFMultiReader.h"
#include "TAlleles.h"

#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringProperties.h"

namespace GLF {
using coretools::instances::logfile;
using coretools::instances::parameters;

using coretools::str::toString;

namespace impl {

TGLFReader *nextChr(const std::vector<TGLFReader *> &ps, bool onlyData) {
	if (onlyData) {
		return *std::min_element(ps.begin(), ps.end(), [](auto p1, auto p2) {
			if (p1->curChromosome().refID() < p2->curChromosome().refID()) return true;
			if (p1->curChromosome().refID() == p2->curChromosome().refID())
				return p1->position() < p2->position();
			else
				return false;
		});
	} else {
		return *std::min_element(ps.begin(), ps.end(), [](auto p1, auto p2) {
			return p1->curChromosome().refID() < p2->curChromosome().refID();
		});
	}
}

} // namespace impl

//--------------------------------------------
// TGlfVector
//--------------------------------------------

void TGLFVector::_openFiles(const std::vector<std::string>& FileNames) {
	// open GLF files
	const auto _numGLFs = FileNames.size();
	_GLFs.clear();
	_GLFs.reserve(_numGLFs);

	logfile().startIndent("Opening " + toString(_numGLFs) + " GLF files:");

	for (size_t i = 0; i < FileNames.size(); ++i) {
		if (_sampleNamesProvided) {
			logfile().listFlush("Opening GLF file '" + FileNames[i] + "' of sample '", _sampleNames[i], "' ...");
		} else {
			logfile().listFlush("Opening GLF file '" + FileNames[i] + "' ...");
		}
		_GLFs.emplace_back(FileNames[i]);
		logfile().done();
	}
	logfile().endIndent();

	// check that they contain the same chromosomes
	if (_GLFs.size() > 1) {
		for (auto g = _GLFs.begin() + 1; g != _GLFs.end(); ++g) {
			if (!_GLFs[0].index().hasSameChromosomes(g->index())) {
				UERROR("GLF files '", _GLFs[0].name(), "' and '", g->name(), "' do not contain the same chromosomes!");
			}
		}
	}
}

TGLFVector::TGLFVector() {
	using namespace coretools::str;
	_GLFs.clear();
	_sampleNamesProvided = false;

	logfile().startIndent("Opening GLF files (parameter 'glf'):");

	std::vector<std::string> GLFFileNames;
	// read file if provided, else parse command line
	const auto FileNames = parameters().get<std::string>("glf");
	if (stringContains(FileNames, ",") || FileNames.substr(FileNames.size() - 6, 6) == "glf.gz") {
		logfile().list("Parsing GLF file names from comma separated list.");
		parameters().fill("glf", GLFFileNames);
	} else {
		logfile().list("Reading GLF info from file '", FileNames, "'.");

		// read file
		coretools::TInputFile in(FileNames, coretools::FileType::NoHeader);
		if (in.numCols() > 1) {
			for (; !in.empty(); in.popFront()) {
				GLFFileNames.emplace_back(in.get(0));
				_sampleNames.emplace_back(in.get(1));
			}
			_sampleNamesProvided = true;
		} else {
			for (; !in.empty(); in.popFront()) { GLFFileNames.emplace_back(in.get(0)); }
		}

		if (_sampleNamesProvided) {
			logfile().conclude("Read ", GLFFileNames.size(), " file names and corresponding sample names.");
		} else {
			logfile().conclude("Read ", GLFFileNames.size(), " file names.");
		}
	}

	// read sample names unless already provided
	if (!_sampleNamesProvided) {
		// read sample names if provided, else use GLF file names
		if (parameters().exists("sampleNames")) {
			const auto sampleNames = parameters().get<std::string>("sampleNames");

			if (!stringContains(sampleNames, ",") && std::filesystem::exists(sampleNames)) {
				logfile().list("Reading sample neames from file '", sampleNames, "'. (parameter 'sampleNames')");
				coretools::TInputFile in(sampleNames, coretools::FileType::NoHeader);
				for (; in.empty(); in.popFront()) { _sampleNames.emplace_back(in.get(0)); }
			} else {
				parameters().fill("sampleNames", _sampleNames);
			}

			if (_sampleNames.size() != GLFFileNames.size()) {
				UERROR("Number of provided sample names does not match number of GLF files!");
			}
			_sampleNamesProvided = true;
		} else {
			logfile().list(
			    "Will deduce sample names from GLF file names. (use 'sampleNames' to provide alternative names)");
			_sampleNames.reserve(GLFFileNames.size());
			for (auto &f : GLFFileNames) { _sampleNames.emplace_back(coretools::str::readBeforeLast(f, ".glf")); }
		}

		// if there are duplicates, add suffix
		bool foundDuplicates = false;
		for (size_t i = 0; i < _sampleNames.size(); ++i) {
			for (size_t j = i + 1; j < _sampleNames.size(); ++j) {
				int counter = 1;
				if (_sampleNames[i] == _sampleNames[j]) {
					_sampleNames[j] += "." + coretools::str::toString(counter++);
					if (!foundDuplicates) {
						logfile().startIndent("Duplicate samples will be rename as follows:");
						foundDuplicates = true;
					}
					logfile().list(_sampleNames[i], " -> ", _sampleNames[j]);
				}
			}
		}
		if (foundDuplicates) { logfile().endIndent(); }
	}

	// open GLF files
	_openFiles(GLFFileNames);
	logfile().endIndent();
}

//----------------------------------------------------
// TGlfMultiReader
//----------------------------------------------------
TGLFMultiReader::TGLFMultiReader() {
	_minDepth = parameters().get<size_t>("minDepth", 0);
	if (_minDepth > 0) logfile().list("Will only keep sites with depth >= " + toString(_minDepth) + ".");

	_windowSize = parameters().get<size_t>("window", 100000);
	if (_windowSize == 0) UERROR("Window size must be at least 1!");

	_GLFIsActive.resize(_GLFs.size(), false);
	setAllActive();
}

void TGLFMultiReader::addReference(const std::string &FastaFile) { _fastaReader.open(FastaFile); }

//-------------------------------------
// set active / inactive
//-------------------------------------

void TGLFMultiReader::_setActive(size_t index) {
	if (index >= _GLFs.size()) UERROR("Index out of range in TGlfMultiReader::setActive(const int index)!");
	if (!_GLFIsActive[index]) {
		_GLFIsActive[index] = true;
		_activeGLFs.push_back(&_GLFs[index]);
	}
}

void TGLFMultiReader::_prepareParsing() {
	for (TGLFReader *it : _activeGLFs) {
		it->rewind();
	}
	// where to start?
	_jumpToNextPosition();
}

void TGLFMultiReader::_jumpToNextPosition() {
	auto min = impl::nextChr(_activeGLFs, _minSamplesWithData);

	if (min->empty()) { return; }

	_curRefId = min->curChromosome().refID();
	if (_minSamplesWithData > 0) {
		_curWindow.move(min->position(), _windowSize);
	} else {
		_curWindow.move(curChr().from(), _windowSize);
	}
	_dataWindow.clear();
}

void TGLFMultiReader::setAllActive() {
	for (size_t i = 0; i < _GLFs.size(); ++i) _setActive(i);
	_prepareParsing();
}

//-------------------------------------
// Looping over active files
//-------------------------------------
bool TGLFMultiReader::_moveToNextChromosome() {
	// increment chromosome ref_char id
	++_curRefId;
	if (_curRefId > _GLFs[0].lastRefID()) return false;

	// advance all active files behind
	bool allFilesReachedEnd = true;
	for (TGLFReader *it : _activeGLFs) {
		it->jumpToChr(_curRefId);
		if (!it->empty()) allFilesReachedEnd = false;
	}

	// check if we reached end of all files
	if (allFilesReachedEnd) return false;

	// get name and length from first active file not at end
	_jumpToNextPosition();
	return true;
}

std::vector<size_t> TGLFMultiReader::readWindow() {
	if (!_dataWindow.empty()) _curWindow.move(_curWindow.to(), _windowSize);

	while (_curWindow.from() >= curChr().to()) {
		if (!_moveToNextChromosome()) return {};
	}

	if (curChr().to() < _curWindow.to()) { _curWindow.move(_curWindow.from(), curChr().to()); }
	const size_t N = _curWindow.size();
	if (N == 0) readWindow();

	_dataWindow.assign(N, {});
	static std::vector<size_t> numActive;
	numActive.assign(N, 0);

	bool allEOF = true;

	for (auto reader : _activeGLFs) {
		// find first data in window
		reader->jumpToPositionOrBeyond(_curWindow.from());
		if (!reader->empty()) allEOF = false;

		// fill everything as noData
		for (size_t iW = 0; iW < N; ++iW) {
			if (!reader->empty() && reader->position() == _curWindow.from() + iW) {
				_dataWindow[iW].emplace_back(reader->front(), reader->depth());
				++numActive[iW];
				reader->popFront();
			} else {
				_dataWindow[iW].emplace_back(reader->curChromosome().isHaploid());
			}
		}
	}
	if (allEOF) return {};

	std::vector<size_t> ids;
	ids.reserve(_dataWindow.size());

	for (size_t i = 0; i < _dataWindow.size(); ++i) {
		if (numActive[i] >= _minSamplesWithData) { ids.push_back(i); }
	}

	if (ids.empty()) { return readWindow(); }
	return ids;
}

std::vector<size_t> TGLFMultiReader::readWindow(const SiteSubset::TAlleles &Alleles) {
	if (!Alleles) return readWindow();

	// Find chromosome in allele-list
	while (Alleles.empty(_curWindow.refID())) {
		if (!_moveToNextChromosome()) return {};
	}

	if (_dataWindow.empty()) { // start of chromosome
		_curWindow.move(Alleles.begin(_curWindow.refID())->position, _windowSize);
	} else {
		_curWindow.move(_curWindow.to(), _windowSize);
	}

	if (_curWindow.from() >= curChr().to()) {
		if (!_moveToNextChromosome()) return {};
	}

	const auto begin = Alleles.begin(_curWindow);
	if (begin == Alleles.end()) return readWindow(Alleles);

	auto latest = begin;
	while ((latest + 1) != Alleles.end() && _curWindow.within((latest + 1)->position)) ++latest;

	if (curChr().to() < _curWindow.to()) { _curWindow.move(_curWindow.from(), curChr().to()); }
	const size_t N = latest->position - _curWindow.from() + 1;
	if (N == 0) readWindow(Alleles);

	_dataWindow.assign(N, {});
	static std::vector<size_t> numActive;
	numActive.assign(N, 0);

	bool allEOF = true;

	for (auto reader : _activeGLFs) {
		// find first data in window
		auto it = begin;
		reader->jumpToPositionOrBeyond(it->position);
		if (!reader->empty()) allEOF = false;

		// fill everything as noData
		for (size_t iW = 0; iW < N; ++iW) {
			if (!reader->empty() && reader->position() == _curWindow.from() + iW) {
				_dataWindow[iW].emplace_back(reader->front(), reader->depth());
				++numActive[iW];
				reader->popFront();
			} else {
				_dataWindow[iW].emplace_back(reader->curChromosome().isHaploid());
			}
		}
	}

	if (allEOF) return {};
	std::vector<size_t> ids;
	ids.reserve(N);

	for (auto it = begin; it != latest + 1; ++it) {
		const auto iW = it->position - _curWindow.from();
		if (numActive[iW] >= _minSamplesWithData) { ids.push_back(iW); }
	}

	if (ids.empty()) { return readWindow(Alleles); }
	return ids;
}

std::vector<std::string> TGLFMultiReader::namesOfActiveFiles() const {
	std::vector<std::string> vec;

	// sample names are file names without glf ending
	for (TGLFReader *it : _activeGLFs) { vec.emplace_back(it->name()); }
	return vec;
}

std::vector<std::string> TGLFMultiReader::sampleNamesOfActiveFiles() const {
	std::vector<std::string> vec;

	for (size_t i = 0; i < _GLFs.size(); ++i) {
		if (_GLFIsActive[i]) { vec.emplace_back(_GLFs.sampleName(i)); }
	}
	return vec;
}

}; // end namespace GLF
