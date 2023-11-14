/*
 * TBamFilter.h
 *
 *  Created on: May 23, 2020
 *      Author: phaentu
 */

#ifndef BAM_TBAMFILTER_H_
#define BAM_TBAMFILTER_H_

#include <string>

#include "coretools/Math/TNumericRange.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Math/counters.h"

#include "TSequencedBase.h"


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
	void disable() noexcept { _enabled = false; }
	bool filters() const{ return _enabled; };
	void filterOut(std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID, coretools::TOutputFile* const log);
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
