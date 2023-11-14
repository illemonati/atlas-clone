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
	const coretools::TCountDistributionVector<>& numFiltered() const { return _counter; }
protected:
	bool _enabled = false;
	coretools::TCountDistributionVector<> _counter;
	std::string _reason; //used for reporting
	coretools::TOutputFile* _log = nullptr;

public:
	TBamFilter() = default;
	virtual ~TBamFilter() = default;
	void enable(std::string_view Reason, size_t numRG, size_t numChrom) {
		_enabled = true;
		_reason  = Reason;
		resizeCounter(numRG, numChrom);
	}
	void disable() noexcept {_enabled = false;}
	void resizeCounter(size_t numRG, size_t numChrom);
	bool filters() const{ return _enabled; };
	void setReason(std::string_view reason);
	void setLog(coretools::TOutputFile& Log);
	void filterOut(std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID);
	void summary(size_t total, size_t readGroup) const;
	const std::string& getReason() const { return _reason; }
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
