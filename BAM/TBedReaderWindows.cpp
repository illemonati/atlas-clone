/*
 * TBedReader.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: linkv
 */

#include <TBedReaderWindows.h>
#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <utility>
#include "genometools/GenomePositions/TChromosomes.h"
#include "coretools/Main/TLog.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM{

using coretools::str::toString;
using coretools::str::fromString;
using coretools::instances::logfile;

//-----------------------
// TBedReaderChromosome
//-----------------------

auto TBedReaderChromosome::findWindow(size_t pos) const {
	int w = (double) pos / (double) windowSize;
	return windows.find(w);
}

auto TBedReaderChromosome::findWindow(size_t pos) {
	int w = (double) pos / (double) windowSize;
	return windows.find(w);
}

auto TBedReaderChromosome::findOrCreateWindow(size_t pos){
	auto windowIt = findWindow(pos);
	if(windowIt == windows.end()){
		//insert window
		int w = (double) pos / (double) windowSize;
		windows.emplace(std::pair<int, TBedReaderWindow>(w, {w*windowSize, (w+1)*windowSize - 1}));
		return findWindow(pos);
	}
	return windowIt;
}

void TBedReaderChromosome::addPosition(std::vector<std::string> & tmp, size_t & numPositionsAdded, size_t siteLimit){
	uint64_t start = fromString<uint64_t>(tmp[1]);
	uint64_t end = fromString<uint64_t>(tmp[2]);

	//add to counter
	numPositionsAdded += end - start;
	if (siteLimit > 0)
		if(numPositionsAdded > siteLimit)
			end -= numPositionsAdded - siteLimit;

	//identify window
	auto windowIt = findOrCreateWindow(start);

	//add position to that window
	//Note BED is already 0 indexed
	for(size_t i=start; i<end; ++i){
		if(i >= windowIt->second.end) findOrCreateWindow(i);
		windowIt->second.addPosition(i);
	}
};

bool TBedReaderChromosome::hasPositionsInWindow(size_t windowStart) const{
	auto windowIt = findWindow(windowStart);
	if(windowIt == windows.end()) return false;
	return true;
};

std::vector<size_t>& TBedReaderChromosome::getPositionInWindow(size_t windowStart){
	auto windowIt = findWindow(windowStart);
	if(windowIt == windows.end()) UERROR("TBedReader Error: window '", windowStart, "' does not exist!");
	return windowIt->second.positions;
};

size_t TBedReaderChromosome::size() const{
	size_t s = 0;
	for (const auto& w: windows) {
		s += w.second.size();
	}
	return s;
};

//-----------------------
// TBedReader
//-----------------------
void TBedReaderWindows::readFile(const genometools::TChromosomes & chromosomeList, size_t siteLimit, bool adaptRegions){
	//open file
	std::istream* myStream = NULL;
	if(filename.find(".gz")) myStream = new gz::igzstream(filename.c_str());
	else myStream = new std::ifstream(filename.c_str());
	if(!*myStream) UERROR("Failed to open BED file '", filename, "'!");

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	curChr = "";

	//read file
	while((*myStream).good() && !(*myStream).eof() && (numPositionsAdded < siteLimit || siteLimit == 0)){
		++lineNum;
		std::string line;
		std::getline(*myStream, line);

		coretools::str::fillContainerFromStringWhiteSpace(line, vec, true);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 3) UERROR("Less than three columns in bed file '", filename, "' on line ", lineNum, + "!");
			if(fromString<int>(vec[1]) < 0 || fromString<int>(vec[2]) < 0) UERROR("Negative value in file '", filename, "' on line ", lineNum, "!");
			if(fromString<int>(vec[2]) <= fromString<int>(vec[1])) UERROR("Error: End position <= start position ('", filename, "', line ", lineNum, ")");
			//get chromosome
			if(!chromosomeList.exists(vec[0])) logfile().warning("Chromosome '" + vec[0] + "' from BED file is not present in the BAM header!");
			if(vec[0] != curChr){
				chrIt = chromosomes.find(vec[0]);
				if(chrIt == chromosomes.end()){
					chromosomes.emplace(std::pair<std::string, TBedReaderChromosome>(vec[0], {vec[0], windowSize}));
					chrIt = chromosomes.find(vec[0]);
				}
				curChr = vec[0];
			}
			if(adaptRegions){
				if(fromString<uint32_t>(vec[1]) < chromosomeList.getChromosome(vec[0]).end().position() && fromString<uint32_t>(vec[2]) > chromosomeList.getChromosome(vec[0]).start().position()){
					if(fromString<uint32_t>(vec[2]) > chromosomeList.getChromosome(vec[0]).end().position()) 
						vec[2] = toString(chromosomeList.getChromosome(vec[0]).end().position());
					if(fromString<uint32_t>(vec[1]) < chromosomeList.getChromosome(vec[0]).start().position())
						vec[1] = toString(chromosomeList.getChromosome(vec[0]).start().position());
					//add positions
					chrIt->second.addPosition(vec, numPositionsAdded, siteLimit);
				}
			} else {
				if(fromString<uint32_t>(vec[1]) > chromosomeList.getChromosome(vec[0]).end().position() || fromString<uint32_t>(vec[1]) < chromosomeList.getChromosome(vec[0]).start().position()) 
					UERROR("Start position for chromosome ", vec[0], " in file '", filename, "' is outside of this chromosome.");
				if(fromString<uint32_t>(vec[2]) > chromosomeList.getChromosome(vec[0]).end().position() || fromString<uint32_t>(vec[2]) < chromosomeList.getChromosome(vec[0]).start().position()) 
					UERROR("End position for chromosome ", vec[0], " in file '", filename, "' is outside of this chromosome.");
				//add positions
				chrIt->second.addPosition(vec, numPositionsAdded, siteLimit);
			}

		}
	}

	//close file
	delete myStream;
};

void TBedReaderWindows::setChr(std::string_view chr){
	curChr = chr;
};

bool TBedReaderWindows::hasPositionsInWindow(size_t windowStart) const{
	auto it = chromosomes.find(curChr);
	if(it == chromosomes.end()) return false;
	else return it->second.hasPositionsInWindow(windowStart);
}

std::vector<size_t>& TBedReaderWindows::getPositionInWindow(uint32_t & windowStart){
	//find chromosome
	chrIt = chromosomes.find(curChr);
	if(chrIt == chromosomes.end()) UERROR("TBedReader Error: chromosome '", curChr, "' does not exist!");
	return chrIt->second.getPositionInWindow(windowStart);
};

size_t TBedReaderWindows::size() const{
	size_t s=0;
	for (const auto& chr : chromosomes) {
		s += chr.second.size();
	}
	return s;
};

size_t TBedReaderWindows::getNumChromosomes() const{
	return chromosomes.size();
};

bool TBedReaderWindows::containsChromosome(std::string chrName) const{
	return chromosomes.count(chrName);
}

TBedReaderChromosome* TBedReaderWindows::findChromosome(std::string chrName) {
	return &(chromosomes.find(chrName)->second);
}

void TBedReaderWindows::listInitializedChromosomes(std::vector<std::string> &initializedChromosomes) {
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
		if (chrIt->second.size() > 0){
			initializedChromosomes.push_back(chrIt->second.name);
		}
	}
}

}; //end namespace


