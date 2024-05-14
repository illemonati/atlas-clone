/*
 * TGlfMultireader.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#include "TGlfMultiReader.h"

#include "coretools/Main/TParameters.h"
#include "coretools/Strings/fillContainer.h"
#include "coretools/Strings/stringProperties.h"
#include "coretools/Files/TInputFile.h"

namespace GLF {
using coretools::instances::logfile;
using coretools::instances::parameters;

using coretools::str::toString;

namespace impl {

TGlfReader * nextChr(const std::vector<TGlfReader *>& ps, bool onlyData) {
	if (onlyData) {
		return *std::min_element(ps.begin(), ps.end(), [](auto p1, auto p2) {
			if (p1->curChromosome().refID() < p2->curChromosome().refID()) return true;
			if (p1->curChromosome().refID() == p2->curChromosome().refID())
				return p1->position() < p2->position();
			else
				return false;
		});
	} else {
		return *std::min_element(ps.begin(), ps.end(), [](auto p1, auto p2) { return p1->curChromosome().refID() < p2->curChromosome().refID(); });
	}
}

} // namespace impl

//--------------------------------------------
// TGlfVector
//--------------------------------------------
void TGlfVector::open(std::string_view FileNames) {
	using namespace coretools::str;
	_GLFs.clear();
	_GLFNames.clear();

	// assume that GLF file names are given in a file if string does not contain ".gz"
	if (!stringContains(FileNames, ".gz")) {
		logfile().list("Reading glf input names from file '", FileNames, "'");	
		
		// read file
		coretools::TInputFile in(FileNames, coretools::FileType::NoHeader);
		for(; in.empty(); in.popFront()){			
			_GLFNames.emplace_back(in.get(0));
		}		
	} else {
		parameters().fill("glf", _GLFNames);
	}
	
	// open GLF files
	size_t _numGLFs     = _GLFNames.size();	
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

	//check that they contain the same chromosomes
	if(_GLFs.size() > 1){
		for(auto g = _GLFs.begin() + 1; g != _GLFs.end(); ++g){
			if(!_GLFs[0].index().hasSameChromosomes(g->index())){
				UERROR("GLF files '", _GLFs[0].name(), "' and '", g->name(), "' do not contain the same chromosomes!");
			}
		}
	}
}

void TGlfVector::openFromParameters(){
	const auto parameter = parameters().get<std::string>("glf");
	open(parameter);
}

int TGlfVector::_getGLFIndexFromName(const std::string &name) const {
	const auto res = std::find(_GLFNames.cbegin(), _GLFNames.cend(), name);
	if (res == _GLFNames.cend()) UERROR("GLF with name '", name, "' not in TGlfMultiReader!");
	return std::distance(_GLFNames.begin(), res);
};


//----------------------------------------------------
// TGlfMultiReader
//----------------------------------------------------
TGlfMultiReader::TGlfMultiReader() {
	_minDepth = parameters().get<size_t>("minDepth", 0);
	if (_minDepth > 0) logfile().list("Will only keep sites with depth >= " + toString(_minDepth) + ".");

	_windowSize = parameters().get<size_t>("window", 100000);
	if (_windowSize == 0) UERROR("Window size must be at least 1!");
};

void TGlfMultiReader::openGLFs() {		
	_GLFs.openFromParameters();
	
	_GLFIsActive.resize(_GLFs.size(), false);
	_setAllInactive();
};

void TGlfMultiReader::addReference(const std::string &FastaFile) {
	fastaReader.open(FastaFile);
};

const genometools::TChromosomes& TGlfMultiReader::chromosomes(){
	return _GLFs[0].chromosomes();
}

//-------------------------------------
// set active / inactive
//-------------------------------------

void TGlfMultiReader::_setActive(size_t index) {
	if (index >= _GLFs.size()) UERROR("Index out of range in TGlfMultiReader::setActive(const int index)!");
	if (!_GLFIsActive[index]) {
		_GLFIsActive[index] = true;
		_activeGLFs.push_back(&_GLFs[index]);
	}
};

void TGlfMultiReader::_setAllInactive() {	
	_GLFIsActive.assign(_GLFs.size(), false);
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
	_curChr = &min->curChromosome();	
	_curRefId = _curChr->refID();	
	_windowStart  = _minSamplesWithData > 0 ? min->position() : 0; // 0 if not jump to position

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
	for (size_t i = 0; i < _GLFs.size(); ++i) _setActive(i);
	_prepareParsing();
};

//-------------------------------------
// Looping over active files
//-------------------------------------
bool TGlfMultiReader::_moveToNextChromosome() {
	// increment chromosome ref_char id
	++_curRefId;
	if(_curRefId > _GLFs[0].chromosomes().size()) return false; 	

	// advance all active files behind
	bool allFilesReachedEnd = true;
	for (TGlfReader *it : _activeGLFs) {
		it->jumpToChr(_curRefId);		
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
	if (_windowStart >= _curChr->length()) {
		if(!_moveToNextChromosome()) return {};
	}
	const size_t N    = std::min(_windowSize, _curChr->length() - _windowStart);
	_dataWindow.assign(N, {});
	_numActive.assign(N, 0);

	bool allEOF=true;
	for (auto reader : _activeGLFs) {
		// find first data in window
		while (!reader->eof() && (reader->curChromosome().refID() == _curRefId) && (reader->position() < _windowStart)) {
			reader->readNext();
		}
		if (!reader->eof()) allEOF = false;

		// fill everything as noData
		for (size_t iW = 0; iW < N; ++iW) {
			if (!reader->eof() && (reader->curChromosome().refID() == _curRefId) && reader->position() - _windowStart == iW) {
				_dataWindow[iW].emplace_back(reader->genotypeLikelihoodsGLF(), reader->depth());
				++_numActive[iW];
				reader->readNext();
			} else {
				_dataWindow[iW].emplace_back(reader->curChromosome().isHaploid());
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
