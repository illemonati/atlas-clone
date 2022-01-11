/*
 * TBedReader.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: linkv
 */

#include <TBedReaderWindows.h>

namespace BAM{

using coretools::str::toString;
using coretools::str::convertString;

//-----------------------
// TBedReaderWindow
//-----------------------

TBedReaderWindow::TBedReaderWindow(uint32_t Start, uint32_t End){
	hasData = false;
	start = Start;
	end = End;
};
TBedReaderWindow::~TBedReaderWindow(){};
void TBedReaderWindow::addPosition(uint32_t & pos){
	positions.push_back(pos);
};

void TBedReaderWindow::print(){
	std::cout << "[" << start+1 << ", " << end+1 << "]:";
	for(std::vector<uint32_t>::iterator it=positions.begin(); it!=positions.end(); ++it) std::cout << " " << *it + 1;
	std::cout << std::endl;
};

uint32_t TBedReaderWindow::size(){
	return positions.size();
}

//-----------------------
// TBedReaderChromosome
//-----------------------

TBedReaderChromosome::TBedReaderChromosome(std::string & Name, uint32_t & WindowSize){
	name = Name;
	windowSize = WindowSize;
};

TBedReaderChromosome::~TBedReaderChromosome(){
	//delete all windows
	for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt){
		delete windowIt->second;
	}
	windows.clear();
};

void TBedReaderChromosome::findWindow(uint32_t pos){
	int w = (double) pos / (double) windowSize;
	windowIt = windows.find(w);
}

void TBedReaderChromosome::findOrCreateWindow(uint32_t pos){
	findWindow(pos);
	if(windowIt == windows.end()){
		//insert window
		int w = (double) pos / (double) windowSize;
		windows.insert(std::pair<int, TBedReaderWindow*>(w, new TBedReaderWindow(w*windowSize, (w+1)*windowSize - 1)));
		findWindow(pos);
	}
}

void TBedReaderChromosome::addPosition(std::vector<std::string> & tmp, uint32_t & numPositionsAdded){
	uint64_t start = convertString<uint64_t>(tmp[1]);
	uint64_t end = convertString<uint64_t>(tmp[2]);

	//add to counter
	numPositionsAdded += end - start;

	//identify window
	findOrCreateWindow(start);

	//add position to that window
	//Note BED is already 0 indexed
	for(uint32_t i=start; i<end; ++i){
		if(i >= windowIt->second->end) findOrCreateWindow(i);
		windowIt->second->addPosition(i);
	}
};

void TBedReaderChromosome::print(){
	for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt) windowIt->second->print();
};

bool TBedReaderChromosome::hasPositionsInWindow(uint32_t windowStart){
	findWindow(windowStart);
	if(windowIt == windows.end()) return false;
	return true;
};

std::vector<uint32_t>& TBedReaderChromosome::getPositionInWindow(const uint32_t windowStart){
	findWindow(windowStart);
	if(windowIt == windows.end()) throw "TBedReader Error: window '" + toString(windowStart) + "' does not exist!";
	return windowIt->second->positions;
};

uint32_t TBedReaderChromosome::size(){
	uint32_t s = 0;
	for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt)
		s += windowIt->second->size();
	return s;
};

//-----------------------
// TBedReader
//-----------------------
void TBedReaderWindows::readFile(const TChromosomes & chromosomeList, uint32_t siteLimit, TLog* logfile){
	//open file
	std::istream* myStream = NULL;
	if(filename.find(".gz")) myStream = new gz::igzstream(filename.c_str());
	else myStream = new std::ifstream(filename.c_str());
	if(!*myStream) throw "Failed to open BED file '" + filename + "'!";

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
			if(vec.size() < 3) throw "Less than three columns in bed file '" + filename + "' on line " + toString(lineNum) + "!";
			//get chromosome
			if(!chromosomeList.exists(vec[0])) logfile->warning("Chromosome '" + vec[0] + "' from BED file is not present in the BAM header!");
			if(vec[0] != curChr){
				chrIt = chromosomes.find(vec[0]);
				if(chrIt == chromosomes.end()){
					chromosomes.insert(std::pair<std::string, TBedReaderChromosome*>(vec[0], new TBedReaderChromosome(vec[0], windowSize)));
					chrIt = chromosomes.find(vec[0]);
				}
				curChr = vec[0];
			}

			//add positions
			chrIt->second->addPosition(vec, numPositionsAdded);
		}
	}

	//close file
	delete myStream;
};

TBedReaderWindows::TBedReaderWindows(std::string Filename, uint32_t WindowSize, const TChromosomes & chromosomeList, uint32_t siteLimit, TLog* logfile){
	filename = Filename;
	windowSize = WindowSize;
	numPositionsAdded = 0;
	curChr = "";
	readFile(chromosomeList, siteLimit, logfile);
};

TBedReaderWindows::~TBedReaderWindows(){
	//delete all chromosomes
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
		delete chrIt->second;
	}
	chromosomes.clear();
};

void TBedReaderWindows::setChr(const std::string & chr){
	curChr = chr;
};

void TBedReaderWindows::print(){
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt) chrIt->second->print();
};

bool TBedReaderWindows::hasPositionsInWindow(uint32_t windowStart){
	chrIt = chromosomes.find(curChr);
	if(chrIt == chromosomes.end()) return false;
	else return chrIt->second->hasPositionsInWindow(windowStart);
}

std::vector<uint32_t>& TBedReaderWindows::getPositionInWindow(uint32_t & windowStart){
	//find chromosome
	chrIt = chromosomes.find(curChr);
	if(chrIt == chromosomes.end()) throw "TBedReader Error: chromosome '" + curChr + "' does not exist!";
	return chrIt->second->getPositionInWindow(windowStart);
};

uint32_t TBedReaderWindows::size(){
	uint32_t s=0;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt)
		s += chrIt->second->size();
	return s;
};

uint32_t TBedReaderWindows::getNumChromosomes(){
	return chromosomes.size();
};

}; //end namesapce


