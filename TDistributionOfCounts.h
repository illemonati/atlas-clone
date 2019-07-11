/*
 * TDepthCounts.h
 *
 *  Created on: Mar 18, 2019
 *      Author: phaentu
 */

#ifndef TDISTRIBUTIONOFCOUNTS_H_
#define TDISTRIBUTIONOFCOUNTS_H_

#include "TFile.h"
#include <string>

class TDistributionOfCounts{
private:
	int maxValue;
	int numBins;
	int lastBin;
	std::string whatIsCounted;

	unsigned int* counts;

public:
    TDistributionOfCounts(int MaxValue, std::string nameOfCounts);
    ~TDistributionOfCounts();

    void add(int depth);
    void add(int depth, int number);
    unsigned int totNumSites();
    void writeCounts(std::string filename);
    void writeNormalizedCounts(std::string filename);
    void writeCumulativeCounts(std::string filename);
    void writeNormalizedCumulativeCounts(std::string filename);
    void writeQuantiles(std::string filename);
};


#endif /* TDISTRIBUTIONOFCOUNTS_H_ */
