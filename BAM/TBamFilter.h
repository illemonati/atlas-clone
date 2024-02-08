/*
 * TBamFilter.h
 *
 *  Created on: May 23, 2020
 *      Author: phaentu
 */

#ifndef BAM_TBAMFILTER_H_
#define BAM_TBAMFILTER_H_

#include <string>
#include "coretools/Counters/TCountDistributionVector.h"

namespace BAM{

//-----------------------------------------------------
//TBamFileFilter
//-----------------------------------------------------

class TBamFilter{
	bool _enabled = false;
	coretools::TCountDistributionVector<> _counter;
	std::string _reason; //used for reporting

public:
	void enable(std::string_view Reason, size_t numRG, size_t numChrom) {
		_enabled = true;
		_reason  = Reason;
		_counter.resize(numRG);
		_counter.resizeDistributions(numChrom);
	}
	void clone(const TBamFilter& Filter, size_t numRG, size_t numChrom) {
		_enabled = Filter._enabled;
		_reason  = Filter._reason;

		if (_enabled) {
			_counter.resize(numRG);
			_counter.resizeDistributions(numChrom);
		} else {
			_counter.resize(0);
		}
	}

	void disable() noexcept { _enabled = false; }
	operator bool() const noexcept {return _enabled;}
	void filterOut(std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID, coretools::TOutputFile& Log);
	void summary(size_t total, size_t readGroup) const;
	void fillHeader(std::vector<std::string> &header) const; 
	void printCounts(coretools::TOutputFile &out, size_t rg_ID) const;
	void printCombinedCounts(coretools::TOutputFile &out) const;
	size_t getCountsPerChromosome(size_t ref_ID) const;
	size_t getCountsAtReadGroupAndChromosome(size_t rg_ID, size_t ref_ID) const;
	size_t getCombinedCounts() const;
	size_t getCounts(size_t rg_ID) const;
};

}; //end namespace

#endif /* BAM_TBAMFILTER_H_ */
