/*
 * TBamFilter.h
 *
 *  Created on: May 23, 2020
 *      Author: phaentu
 */

#ifndef BAM_TBAMFILTER_H_
#define BAM_TBAMFILTER_H_

#include <array>
#include <memory>
#include <set>
#include <stdint.h>
#include <string>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "coretools/Files/TFile.h"
#include "coretools/Math/TNumericRange.h"
#include "TSequencedBase.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Math/counters.h"


namespace BAM{

//-----------------------------------------------------
//TBamFileFilter
//-----------------------------------------------------
class TBamFileFilter{
protected:
	bool _keep = true;
	coretools::TCountDistributionVector<> _counter;
	std::string _reason; //used for reporting
	coretools::TOutputFile* _log = nullptr;

public:
	TBamFileFilter() = default;
	virtual ~TBamFileFilter() = default;
	void keep();
	void resizeCounter(size_t numRG, size_t numChrom);
	bool filters() const{ return !_keep; };
	void setReason(std::string_view reason);
	void setLog(coretools::TOutputFile& Log);
	void filterOut(std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID);
	void summary(size_t total, size_t readGroup) const;
	const coretools::TCountDistributionVector<>& numFiltered() const { return _counter; }
	const std::string& getReason() const { return _reason; }
	void fillHeader(std::vector<std::string> &header) const; 
	void printCounts(coretools::TOutputFile &out, size_t rg_ID) const;
	void printCountsPerChromosome(coretools::TOutputFile &out, size_t ref_ID) const;
	void printCombinedCounts(coretools::TOutputFile &out) const;
	size_t getCounts(size_t rg_ID) const;
	size_t getCountsPerChromosome(size_t ref_ID) const;
	size_t getCountsAtReadGroupAndChromosome(size_t rg_ID, size_t ref_ID) const;
	size_t getCombinedCounts() const;
};

//Filters out if pass(false), keeps if pass(true)
class TBamFileFilterBool:public TBamFileFilter{
public:
	TBamFileFilterBool(){};
	void filter(std::string_view Reason, size_t numRG, size_t numChrom);
	bool pass(bool state, std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID);
};

template <typename T = size_t>
class TBamFileFilterRange:public TBamFileFilter {
private:
	coretools::TNumericRange<T> _range;
public:
	TBamFileFilterRange() = default;
	TBamFileFilterRange(T max): _range(0, true, max, true) {}

	void filter(const coretools::TNumericRange<T> & Range, std::string_view Reason, size_t numRG, size_t numChrom){
		_keep = false;
		_range = Range;
		_reason = Reason;
		resizeCounter(numRG, numChrom);
	};

	const coretools::TNumericRange<T> range() const {
		return _range;
	};

	std::string rangeString() const{
		return _range.rangeString();
	};

	bool pass(const T & value, std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID) {
		if(!_keep && !_range.within(value)){
			filterOut(alignmentName, isSecondMate, readGroup, chromosomeID);
			return false;
		}
		return true;
	};
};

//-------------------------------------
// TBaseFilter
//-------------------------------------
class TBaseFilter{
protected:
	bool _filter;

public:
	explicit constexpr TBaseFilter() : _filter(false) {};
	virtual ~TBaseFilter() = default;

	constexpr operator bool() const{
		return _filter;
	};

	virtual bool pass(const TSequencedBase & base) const = 0;
};

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter : public TBaseFilter{
private:
	coretools::TNumericRange<genometools::PhredIntProbability> _range;

public:
	TQualityFilter();

	bool pass(const TSequencedBase & base) const override{
		return _range.within(base.originalQuality_phredInt);
	};
};

//-------------------------------------
// TContextFilter
//-------------------------------------
class TContextFilter : public TBaseFilter{
private:
	//std::array<bool, static_cast<uint8_t>(genometools::BaseContext::max) + 1> _keptContexts;
	coretools::TStrongArray<bool, genometools::BaseContext, coretools::index(genometools::BaseContext::max) + 1> _keptContexts;

public:
	explicit TContextFilter();
	bool pass(const TSequencedBase & base) const override;
};


//-----------------------------------------------------
//TAlignmentBlacklist
//-----------------------------------------------------
class TAlignmentList{
private:
	std::set<std::string> _list;

public:
	TAlignmentList() = default;
	TAlignmentList(std::string_view filename);
	void addFromFile(std::string_view filename);
	void add(std::string_view alignment);
	void remove(const std::string& alignment);
	void clear();
	bool isInBlacklist(const std::string& alignment) const;
};

}; //end namespace

#endif /* BAM_TBAMFILTER_H_ */
