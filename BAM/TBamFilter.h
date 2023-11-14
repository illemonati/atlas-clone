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
private:
	const coretools::TCountDistributionVector<>& numFiltered() const { return _counter; }
protected:
	bool _keep = true;
	coretools::TCountDistributionVector<> _counter;
	std::string _reason; //used for reporting
	coretools::TOutputFile* _log = nullptr;

public:
	TBamFilter() = default;
	virtual ~TBamFilter() = default;
	void keep();
	void resizeCounter(size_t numRG, size_t numChrom);
	bool filters() const{ return !_keep; };
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

//Filters out if pass(false), keeps if pass(true)
class TBamFilterBool:public TBamFilter{
public:
	TBamFilterBool(){};
	void filter(std::string_view Reason, size_t numRG, size_t numChrom);
	bool pass(bool state, std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID);
};

class TBamFilterRange:public TBamFilter {
private:
	coretools::TNumericRange<size_t> _range;
public:
	TBamFilterRange() = default;
	TBamFilterRange(size_t max): _range(0, true, max, true) {}

	void filter(const coretools::TNumericRange<size_t> & Range, std::string_view Reason, size_t numRG, size_t numChrom){
		_keep = false;
		_range = Range;
		_reason = Reason;
		resizeCounter(numRG, numChrom);
	}

	const coretools::TNumericRange<size_t> range() const {
		return _range;
	}

	std::string rangeString() const{
		return _range.rangeString();
	}

	bool pass(size_t value, std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID) {
		if(!_keep && !_range.within(value)){
			filterOut(alignmentName, isSecondMate, readGroup, chromosomeID);
			return false;
		}
		return true;
	}
};

}; //end namespace

#endif /* BAM_TBAMFILTER_H_ */
