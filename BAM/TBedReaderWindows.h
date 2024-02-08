/*
 * TBedReader.h
 *
 *  Created on: Oct 6, 2015
 *      Author: wegmannd
 */

#ifndef TBEDREADER_H_
#define TBEDREADER_H_

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "coretools/Counters/TCountDistribution.h"

namespace genometools { class TChromosomes; }

//read sorted bed files window by window
//store all data in chr / window combinations using vectors
//Store all positions 0-based, as in TWindow

namespace BAM{

class TBedReaderWindow{
public:
	bool hasData = false;
	size_t start, end;
	std::vector<size_t> positions;

	TBedReaderWindow(size_t Start, size_t End) : start(Start), end(End) {};
	void addPosition(size_t pos) {positions.push_back(pos);};
	uint32_t size() const noexcept {return positions.size();};
};

class TBedReaderChromosome{
public:
	std::string name;
	std::map<size_t, TBedReaderWindow> windows;
	size_t windowSize;
	coretools::TCountDistribution<> distPerSites;

	TBedReaderChromosome(std::string_view Name, size_t WindowSize): name(Name), windowSize(WindowSize) {};
	auto findWindow(size_t pos) const;
	auto findWindow(size_t pos);
	auto findOrCreateWindow(size_t pos);
	void addPosition(std::vector<std::string> & tmp, size_t & numPositionsAdded, size_t siteLimit);
	bool hasPositionsInWindow(size_t windowStart) const;
	std::vector<size_t>& getPositionInWindow(size_t windowStart);
	size_t size() const;
};

class TBedReaderWindows{
private:
	std::string filename;
	std::map<std::string, TBedReaderChromosome> chromosomes;
	std::map<std::string, TBedReaderChromosome>::iterator chrIt;
	size_t windowSize;
	std::string curChr = "";
	size_t numPositionsAdded = 0;

	void readFile(const genometools::TChromosomes & chromosomeList, size_t siteLimit, bool adaptRegions);

public:
	TBedReaderWindows(std::string_view Filename, size_t WindowSize, const genometools::TChromosomes & chromosomes, size_t siteLimit, bool adaptRegions) : filename(Filename), windowSize(WindowSize) {readFile(chromosomes, siteLimit, adaptRegions);};
	void setChr(std::string_view chr);
	bool hasPositionsInWindow(size_t windowStart) const;
	std::vector<size_t>& getPositionInWindow(uint32_t & windowStart);
	size_t size() const;
	size_t getNumChromosomes() const;
	bool containsChromosome(const std::string& chrName) const;
	TBedReaderChromosome* findChromosome(const std::string& chrName);
	void listInitializedChromosomes(std::vector<std::string> &initializedChromosomes);
};

}; //end namesapce

#endif /* TBEDREADER_H_ */
