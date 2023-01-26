/*
 * TBedReader.h
 *
 *  Created on: Oct 6, 2015
 *      Author: wegmannd
 */

#ifndef TBEDREADER_H_
#define TBEDREADER_H_

#include <stdint.h>
#include <map>
#include <string>
#include <vector>
#include "coretools/Math/counters.h"
namespace genometools { class TChromosomes; }

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
	coretools::TCountDistribution<> distPerSites;

	TBedReaderChromosome(std::string & Name, uint32_t & WindowSize);
	~TBedReaderChromosome();
	void findWindow(uint32_t pos);
	void findOrCreateWindow(uint32_t pos);
	void addPosition(std::vector<std::string> & tmp, uint32_t & numPositionsAdded, uint32_t siteLimit);
	void print();
	bool hasPositionsInWindow(uint32_t windowStart);
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

	void readFile(const genometools::TChromosomes & chromosomeList, uint32_t siteLimit, bool adaptRegions);

public:
	std::string filename;
	TBedReaderWindows(std::string Filename, uint32_t WindowSize, const genometools::TChromosomes & chromosomes, uint32_t siteLimit, bool adaptRegions);
	TBedReaderWindows();
	~TBedReaderWindows();
	void setChr(const std::string & chr);
	void print();
	bool hasPositionsInWindow(uint32_t windowStart);
	std::vector<uint32_t>& getPositionInWindow(uint32_t & windowStart);
	uint32_t size();
	uint32_t getNumChromosomes();
	bool containsChromosome(std::string chrName) const;
	TBedReaderChromosome* findChromosome(std::string chrName) const;
	void listInitializedChromosomes(std::vector<std::string> &initializedChromosomes);
};

}; //end namesapce

#endif /* TBEDREADER_H_ */
