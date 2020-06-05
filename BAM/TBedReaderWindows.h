/*
 * TBedReader.h
 *
 *  Created on: Oct 6, 2015
 *      Author: wegmannd
 */

#ifndef TBEDREADER_H_
#define TBEDREADER_H_

#include <fstream>
#include <vector>
#include <map>
#include "stringFunctions.h"
#include "TLog.h"
#include "gzstream.h"
#include "TChromosomes.h"

//read sorted bed files window by window
//store all data in chr / window combinations using vectors
//Store all positions 0-based, as in TWindow

namespace BAM{

class TBedReaderWindow{
public:
	bool hasData;
	uint32_t start, end;
	std::vector<uint32_t> positions;

	TBedReaderWindow(uint32_t Start, uint32_t End);
	~TBedReaderWindow();
	void addPosition(uint32_t & pos);
	void print();
	uint32_t size();
};

class TBedReaderChromosome{
public:
	std::string name;
	std::map<uint32_t, TBedReaderWindow*> windows;
	std::map<uint32_t, TBedReaderWindow*>::iterator windowIt;
	int windowSize;

	TBedReaderChromosome(std::string & Name, uint32_t & WindowSize);
	~TBedReaderChromosome();
	void findWindow(const uint32_t & pos);
	void findOrCreateWindow(const uint32_t & pos);
	void addPosition(std::vector<std::string> & tmp, uint32_t & numPositionsAdded);
	void print();
	bool hasPositionsInWindow(const uint32_t & windowStart);
	std::vector<uint32_t>& getPositionInWindow(uint32_t windowStart);
	uint32_t size();
};

class TBedReaderWindows{
private:
	std::map<std::string, TBedReaderChromosome*> chromosomes;
	std::map<std::string, TBedReaderChromosome*>::iterator chrIt;
	uint32_t windowSize;
	std::string curChr;
	uint32_t numPositionsAdded;

	void readFile(const TChromosomes & chromosomeList, uint32_t siteLimit, TLog* logfile);

public:
	std::string filename;
	TBedReaderWindows(std::string Filename, const uint32_t & WindowSize, const TChromosomes & chromosomes, uint32_t siteLimit, TLog* logfile);
	~TBedReaderWindows();
	void setChr(const std::string & chr);
	void print();
	bool hasPositionsInWindow(const uint32_t & windowStart);
	std::vector<uint32_t>& getPositionInWindow(uint32_t & windowStart);
	uint32_t size();
	uint32_t getNumChromosomes();

};

}; //end namesapce

#endif /* TBEDREADER_H_ */
