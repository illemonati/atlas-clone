/*
 * TGlfMultireader.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#include "TGlfMultiReader.h"
#include "TAlleles.h"

#include "coretools/Main/TParameters.h"
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

void TGlfVector::_openFiles(){
	// open GLF files
	size_t _numGLFs     = _GLFFileNames.size();	
	_GLFs.clear();
	_GLFs.reserve(_numGLFs);

	logfile().startIndent("Opening " + toString(_numGLFs) + " GLF files:");
	for(size_t i=0; i < _GLFFileNames.size(); ++i){
		if(_sampleNamesProvided){
			logfile().listFlush("Opening GLF file '" + _GLFFileNames[i] + "' of sample '", _sampleNames[i], "' ...");
		} else {
			logfile().listFlush("Opening GLF file '" + _GLFFileNames[i] + "' ...");
		}
		_GLFs.emplace_back(_GLFFileNames[i]);
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

void TGlfVector::openFromParameters() {
	using namespace coretools::str;	
	_GLFs.clear();
	_GLFFileNames.clear();
	_sampleNamesProvided = false;

	logfile().startIndent("Opening GLF files (parameter 'glf'):");

	// read file if provided, else parse command line
	const auto FileNames = parameters().get<std::string>("glf");
	if(stringContains(FileNames, ",") || FileNames.substr(FileNames.size()-6, 6) == "glf.gz"){
		logfile().list("Parsing GLF file names from comma separated list.");	
		parameters().fill("glf", _GLFFileNames);		
	} else {	
		logfile().list("Reading GLF info from file '", FileNames, "'.");	
		
		// read file
		coretools::TInputFile in(FileNames, coretools::FileType::NoHeader);
		if(in.numCols() > 1){
			for(; !in.empty(); in.popFront()){			
				_GLFFileNames.emplace_back(in.get(0));
				_sampleNames.emplace_back(in.get(1));
			}		
			_sampleNamesProvided = true;
		} else {
			for(; !in.empty(); in.popFront()){			
				_GLFFileNames.emplace_back(in.get(0));
			}		
		}

		if(_sampleNamesProvided){
			logfile().conclude("Read ", _GLFFileNames.size(), " file names and corresponding sample names.");
		} else {
			logfile().conclude("Read ", _GLFFileNames.size(), " file names.");
		}
	}
	
	// read sample names unless already provided
	if(!_sampleNamesProvided){
		// read sample names if provided, else use GLF file names
		if (parameters().exists("sampleNames")) {								
			const auto sampleNames = parameters().get<std::string>("sampleNames");

			if (!stringContains(sampleNames, ",") && std::filesystem::exists(sampleNames)) {		
				logfile().list("Reading sample neames from file '", sampleNames, "'. (parameter 'sampleNames')");
				coretools::TInputFile in(sampleNames, coretools::FileType::NoHeader);
				for(; in.empty(); in.popFront()){								
					_sampleNames.emplace_back(in.get(0));
				}
			} else {
				parameters().fill("sampleNames", _sampleNames);
			}			

			if (_sampleNames.size() != _GLFFileNames.size()) {
				UERROR("Number of provided sample names does not match number of GLF files!");
			}			
			_sampleNamesProvided = true;
		} else {
			logfile().list("Will deduce sample names from GLF file names. (use 'sampleNames' to provide alternative names)");			
			_sampleNames.reserve(_GLFFileNames.size());
			for (auto& f : _GLFFileNames) {
				_sampleNames.emplace_back(coretools::str::readBeforeLast(f, ".glf"));
			}
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
	_openFiles();	
	logfile().endIndent();
}

size_t TGlfVector::index(const std::string &name) const {
	const auto res = std::find(_GLFFileNames.cbegin(), _GLFFileNames.cend(), name);
	if (res == _GLFFileNames.cend()) UERROR("GLF with name '", name, "' not in TGlfMultiReader!");
	return std::distance(_GLFFileNames.begin(), res);
}

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

	if(min->eof()){ return false;  }

	_curChr = &min->curChromosome();	
	_curRefId = _curChr->refID();	
	if(_minSamplesWithData > 0){
		_curWindow.move(min->position(), _windowSize);
	} else {
		_curWindow.move(_curChr->from(), _windowSize);
	}
	_dataWindow.clear();

	return true;
};

void TGlfMultiReader::setActive(int index) {
	_setAllInactive();
	_setActive(index);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string &name) {
	setActive(_GLFs.index(name));
};

void TGlfMultiReader::setActive(int index1, int index2) {
	_setAllInactive();
	_setActive(index1);
	_setActive(index2);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string &name1, const std::string &name2) {
	setActive(_GLFs.index(name1), _GLFs.index(name2));
};

void TGlfMultiReader::setActive(const std::vector<int> &indexes) {
	_setAllInactive();
	for (const auto i : indexes) _setActive(i);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::vector<std::string> &names) {
	_setAllInactive();
	for (const auto &name : names) {
		_setActive(_GLFs.index(name));
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
	if(_curRefId > _GLFs[0].lastRefID()) return false; 	

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

std::vector<size_t> TGlfMultiReader::readWindow(const SiteSubset::TAlleles& Alleles) {
	if (!_dataWindow.empty()) _curWindow.move(_curWindow.to(), _windowSize);

	if (_curWindow.from() >= _curChr->to()) {
		if(!_moveToNextChromosome()) return {};
	}

	if (Alleles) {
		auto begin = Alleles.begin(_curWindow);
		if (begin == Alleles.end()) return readWindow(Alleles);

		_curWindow.move(begin->position, _curWindow.size());
	} else if (_curChr->to() < _curWindow.to()){
		_curWindow.move(_curWindow.from(), _curChr->to());
	}
	const size_t N = _curWindow.size();

	_dataWindow.assign(N, {});
	_numActive.assign(N, 0);

	bool allEOF=true;

	for (auto reader : _activeGLFs) {
		// find first data in window
		reader->jumpToAtOrFirstAfterPosition(_curWindow.from());
		if (!reader->eof()) allEOF = false;

		// fill everything as noData
		for (size_t iW = 0; iW < N; ++iW) {
			if (!reader->eof() && reader->position() == _curWindow.from() + iW) {
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

	if (ids.empty()){		
		return readWindow(Alleles);
	} 
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
	
	for(size_t i = 0; i < _GLFs.size(); ++i){	
		if(_GLFIsActive[i]){
		 	vec.emplace_back(_GLFs.sampleName(i));
		}
	}
	return vec;
};


}; // end namespace GLF
