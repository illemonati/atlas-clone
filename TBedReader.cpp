/*
 * TBedReader.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: linkv
 */

#include "TBedReader.h"

//-----------------------
// TBedReaderWindow
//-----------------------

TBedReaderWindow::TBedReaderWindow(unsigned int Start, unsigned int End){
	hasData = false;
	start = Start;
	end = End;
};
TBedReaderWindow::~TBedReaderWindow(){};
void TBedReaderWindow::addPosition(unsigned int & pos){
	positions.push_back(pos);
};

void TBedReaderWindow::print(){
	std::cout << "[" << start+1 << ", " << end+1 << "]:";
	for(std::vector<unsigned int>::iterator it=positions.begin(); it!=positions.end(); ++it) std::cout << " " << *it + 1;
	std::cout << std::endl;
};

unsigned int TBedReaderWindow::size(){
	return positions.size();
}

//-----------------------
// TBedReaderChromosome
//-----------------------

TBedReaderChromosome::TBedReaderChromosome(std::string & Name, unsigned int & WindowSize){
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

void TBedReaderChromosome::findWindow(const unsigned int & pos){
	int w = (double) pos / (double) windowSize;
	windowIt = windows.find(w);
}

void TBedReaderChromosome::findOrCreateWindow(const unsigned int & pos){
	findWindow(pos);
	if(windowIt == windows.end()){
		//insert window
		int w = (double) pos / (double) windowSize;
		windows.insert(std::pair<int, TBedReaderWindow*>(w, new TBedReaderWindow(w*windowSize, (w+1)*windowSize - 1)));
		findWindow(pos);
	}
}

void TBedReaderChromosome::addPosition(std::vector<std::string> & tmp, unsigned int & numPositionsAdded){
	long start = stringToLong(tmp[1]);
	long end = stringToLong(tmp[2]);

	//add to counter
	numPositionsAdded += end - start;

	//identify window
	findOrCreateWindow(start);

	//add position to that window
	//Note BED is already 0 indexed
	for(unsigned int i=start; i<end; ++i){
		if(i >= windowIt->second->end) findOrCreateWindow(i);
		windowIt->second->addPosition(i);
	}
};

void TBedReaderChromosome::print(){
	for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt) windowIt->second->print();
};

bool TBedReaderChromosome::hasPositionsInWindow(const unsigned int & windowStart){
	findWindow(windowStart);
	if(windowIt == windows.end()) return false;
	return true;
};

std::vector<unsigned int>& TBedReaderChromosome::getPositionInWindow(const unsigned int windowStart){
	findWindow(windowStart);
	if(windowIt == windows.end()) throw "TBedReader Error: window '" + toString(windowStart) + "' does not exist!";
	return windowIt->second->positions;
};

unsigned int TBedReaderChromosome::size(){
	unsigned int s = 0;
	for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt)
		s += windowIt->second->size();
	return s;
};

//-----------------------
// TBedReader
//-----------------------

void TBedReader::readFile(BamTools::SamSequenceDictionary & Sequences, unsigned int siteLimit, TLog* logfile){
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
	while((*myStream).good() && !(*myStream).eof() && (numPositionsAdded < siteLimit || siteLimit == -1)){
		++lineNum;
		std::string line;
		std::getline(*myStream, line);

		fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 3) throw "Less than three columns in bed file '" + filename + "' on line " + toString(lineNum) + "!";
			//get chromosome
			if(!Sequences.Contains(vec[0])) logfile->warning("Chromosome '" + vec[0] + "' from BED file is not present in the BAM header!");
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



TBedReader::TBedReader(std::string Filename, const unsigned int & WindowSize, BamTools::SamSequenceDictionary & Sequences, unsigned int siteLimit, TLog* logfile){
	filename = Filename;
	windowSize = WindowSize;
	numPositionsAdded = 0;
	curChr = "";
	readFile(Sequences, siteLimit, logfile);
};

TBedReader::~TBedReader(){
	//delete all chromosomes
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
		delete chrIt->second;
	}
	chromosomes.clear();
};

void TBedReader::setChr(const std::string & chr){
	curChr = chr;
};

void TBedReader::print(){
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt) chrIt->second->print();
};

bool TBedReader::hasPositionsInWindow(const unsigned int & windowStart){
	chrIt = chromosomes.find(curChr);
	if(chrIt == chromosomes.end()) return false;
	else return chrIt->second->hasPositionsInWindow(windowStart);
}

std::vector<unsigned int>& TBedReader::getPositionInWindow(unsigned int & windowStart){
	//find chromosome
	chrIt = chromosomes.find(curChr);
	if(chrIt == chromosomes.end()) throw "TBedReader Error: chromosome '" + curChr + "' does not exist!";
	return chrIt->second->getPositionInWindow(windowStart);
};

unsigned int TBedReader::size(){
	unsigned int s=0;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt)
		s += chrIt->second->size();
	return s;
};

unsigned int TBedReader::getNumChromosomes(){
	return chromosomes.size();
};




