/*
 * TSiteSubset.cpp
 *
 *  Created on: Nov 25, 2019
 *      Author: vivian
 */

#include "TSiteSubset.h"

namespace GenotypeLikelihoods{

//-------------------------------------------------
// TSiteSubsetSite
//-------------------------------------------------
TSiteSubsetSite::TSiteSubsetSite(const uint32_t Position, const Base Ref, const Base Alt){
	position = Position;
	ref = Ref;
	alt = Alt;
};

void TSiteSubsetSite::print() const{
	std::cout << " " << position << "(" << ref << "," << alt << ")";
};


//-------------------------------------------------
// TSiteSubsetWindow
//-------------------------------------------------
TSiteSubsetWindow::TSiteSubsetWindow(uint32_t Start, uint32_t End){
	_start = Start;
	_end = End;
};
TSiteSubsetWindow::~TSiteSubsetWindow(){};

void TSiteSubsetWindow::addPosition(const uint32_t pos, const Base ref, const Base alt){
	_positions.emplace(pos, ref, alt);
};

void TSiteSubsetWindow::print(){
	std::cout << "[" << _start << ", " << _end << "]:";
	for(auto& p : _positions){
		p.print();
	}
	std::cout << std::endl;
};

size_t TSiteSubsetWindow::size(){
	return _positions.size();
};

//-------------------------------------------------
// TSiteSubsetChr
//-------------------------------------------------
TSiteSubsetChr::TSiteSubsetChr(uint32_t & WindowSize){
	windowSize = WindowSize;
};

std::map<uint32_t, TSiteSubsetWindow>::iterator TSiteSubsetChr::_findWindow(const unsigned int & pos){
	int w = (double) pos / (double) windowSize;
	return windows.find(w);
}

std::map<uint32_t, TSiteSubsetWindow>::iterator TSiteSubsetChr::_findOrCreateWindow(const unsigned int & pos){
	auto it = _findWindow(pos);
	if(it == windows.end()){
		//insert window
		int w = (double) pos / (double) windowSize;
		windows.emplace(std::piecewise_construct, std::forward_as_tuple(w), std::forward_as_tuple(w*windowSize, (w+1)*windowSize - 1));
		it = _findWindow(pos);
	}
	return it;
};

void TSiteSubsetChr::addPosition(const uint32_t Position, const Base Ref, const Base Alt){
	auto it = _findOrCreateWindow(Position);
	it->second.addPosition(Position, Ref, Alt);
};

void TSiteSubsetChr::print(const std::string chrName){
	std::cout << "Chromosome '" << chrName << "':" << std::endl;
	for(auto& w : windows){
		w.second.print();
	}
};

bool TSiteSubsetChr::hasPositionsInWindow(const unsigned int & windowStart){
	return _findWindow(windowStart)!= windows.end();
};

std::vector<TSiteSubsetSite>& TSiteSubsetChr::getPositionInWindow(const unsigned int & windowStart){
	auto w = _findWindow(windowStart);
	if(w == windows.end()){
		throw "TSiteSubset Error: window '" + toString(windowStart) + "' does not exist!";
	}
	return w->second._positions;
};

size_t TSiteSubsetChr::size(){
	size_t s = 0;
	for(auto& w : windows){
		s += w.second.size();
	}
	return s;
};

//-------------------------------------------------
// TSiteSubset
//-------------------------------------------------
std::map<std::string, TSiteSubsetChr>::iterator TSiteSubset::_findOrCreateChromosome(const std::string chrName){
	auto c = _chromosomes.find(chrName);
	if(c == _chromosomes.end()){
		_chromosomes.emplace(std::piecewise_construct, std::forward_as_tuple(chrName), std::forward_as_tuple(chrName, _windowSize));
		c = _chromosomes.find(chrName);
	}
	return c;
};

void TSiteSubset::_checkAlleles(const std::string & chr, const uint32_t & pos, const Base & ref, const Base & alt, const std::string & refAllele, const std::string & altAllele){
	if(ref == N){
			throw "Unknown reference allele '" + refAllele + "' on chr " + chr + " at " + toString(pos+1) + "!";
		}
		if(alt == N){
			throw "Unknown alternative allele '" + altAllele + "' on chr " + chr + " at " + toString(pos+1) + "!";
		}

		if(ref == alt && !_storesInvariantSites){
			throw "Site on chr " + chr + " at " + toString(pos+1) + " is invariant!";
		}
		if(ref != alt && _storesInvariantSites){
			throw "Site on chr " + chr + " at " + toString(pos+1) + " is polymorphic!";
		}
};

void TSiteSubset::_readFile(TLog* logfile){
	logfile->listFlushTime("Reading sites to be used from '" + _filename + "' ...");
	//open file
	std::ifstream sitesFile(_filename.c_str());
	if(!sitesFile) throw "Failed to open sites file '" + _filename + "'!";

	//tmp variables
	TGenotypeMap genoMap;
	long lineNum = 0;
	std::vector<std::string> vec;

	//read file
	while(sitesFile.good() && !sitesFile.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(sitesFile, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 4)
				throw "Wrong number of columns in sites file '" + _filename + "' on line " + toString(lineNum) + "! Expected 4 but read " + toString(vec.size()) + ".";

			//get chromosome
			auto c = _findOrCreateChromosome(vec[0]);

			//extract positions
			uint32_t pos = stringToInt(vec[1]) - 1; //make 0-based
			Base ref = genoMap.getBase(vec[2][0]);
			Base alt = genoMap.getBase(vec[3][0]);

			//check alleles
			_checkAlleles(c->first, pos, ref, alt, vec[2], vec[3]);

			//add position
			c->second.addPosition(pos, ref, alt);
		}
	}

	//close file
	sitesFile.close();

	//report
	logfile->doneTime();
	logfile->conclude("Parsed " + toString(size()) + " sites on " + toString(_chromosomes.size()) + " chromosomes.");
};

void TSiteSubset::_readFile(BAM::TFastaBuffer & reference, BAM::TChromosomes & Chromosomes, TLog* logfile){ //version that checks witth fasta reference
	logfile->listFlushTime("Reading sites to be used from '" + _filename + "' ...");
	//open file
	std::ifstream sitesFile(_filename.c_str());
	if(!sitesFile) throw "Failed to open sites file '" + _filename + "'!";

	//tmp variables
	TGenotypeMap genoMap;
	long lineNum = 0;
	std::vector<std::string> vec;

	//report conflicts with fasta file
	bool conflictsFound = false;
	TOutputFile conflicts;

	//read file
	while(sitesFile.good() && !sitesFile.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(sitesFile, vec);
		//skip empty lines
		if(vec.size() > 0){
			throw "Wrong number of columns in sites file '" + _filename + "' on line " + toString(lineNum) + "! Expected 4 but read " + toString(vec.size()) + ".";

			//get chromosome
			auto c = _findOrCreateChromosome(vec[0]);

			//extract positions
			uint32_t pos = stringToInt(vec[1]) - 1; //make 0-based
			Base ref = genoMap.getBase(vec[2][0]);
			Base alt = genoMap.getBase(vec[3][0]);

			//check alleles
			_checkAlleles(c->first, pos, ref, alt, vec[2], vec[3]);

			//add position
			c->second.addPosition(pos, ref, alt);

			//check with reference
			Base trueRef = genoMap.getBase(reference.refAt(Chromosomes.getIndexFromName(c->first), pos));
			if(trueRef != ref && trueRef != alt){
				//conflict with fasta
				if(!conflictsFound){
					conflicts.open(_filename + ".conflicts.gz");
					conflicts.writeHeader({"chr", "position", "reference", "allele1", "allele2"});
				}
				conflicts << c->first << pos+1 << genoMap.getBaseAsChar(trueRef) << genoMap.getBaseAsChar(ref) << genoMap.getBaseAsChar(alt) << std::endl;
			}

			//add position
			c->second.addPosition(vec, c->first, _storesInvariantSites);
		}
	}
	//report
	logfile->doneTime();
	logfile->conclude("Parsed " + toString(size()) + " sites on " + toString(_chromosomes.size()) + " chromosomes.");

	//write conflicts, if any
	if(conflictsFound){
		logfile->conclude("Reference conflicted with provided alleles at " + toString(conflicts.lineNumber() - 1) + " positions!");
		logfile->conclude("These positions were written to " + conflicts.name());
	}

	//close file
	sitesFile.close();
};

TSiteSubset::TSiteSubset(std::string Filename, uint32_t & WindowSize, TLog* logfile, bool InvariantSites){
	_filename = Filename;
	_windowSize = WindowSize;
	_storesInvariantSites = InvariantSites;
	_readFile(logfile);
	curChrIt = _chromosomes.end();
};

TSiteSubset::TSiteSubset(std::string Filename, uint32_t WindowSize, TLog* logfile, bool InvariantSites, BAM::TFastaBuffer & Reference, BAM::TChromosomes & Chromosomes){
	_filename = Filename;
	_windowSize = WindowSize;
	_storesInvariantSites = InvariantSites;
	_readFile(Reference, Chromosomes, logfile);
	curChrIt = _chromosomes.end();
};

void TSiteSubset::setChr(const std::string chr){
	curChrIt = _chromosomes.find(chr);
};

void TSiteSubset::print(){
	std::cout << "Sites File '" << _filename << "':" << std::endl;
	for(auto& c : _chromosomes){
		c.second.print(c.first);
	}
};

bool TSiteSubset::hasPositionsInWindow(const unsigned int & windowStart){
	if(curChrIt == _chromosomes.end()){
		return false;
	}

	return curChrIt->second.hasPositionsInWindow(windowStart);
}

std::vector<TSiteSubsetSite>& TSiteSubset::getPositionInWindow(const unsigned int & windowStart){
	if(!hasPositionsInWindow(windowStart)){
		return empty;
	}

	return curChrIt->second.getPositionInWindow(windowStart);
};

size_t TSiteSubset::size(){
	size_t size = 0;
	for(auto& c : _chromosomes){
		size += c.second.size();
	}
	return size;
};

}; //end namespace
