
#ifndef BAM_TBAMFILEFILTERS_H_
#define BAM_TBAMFILEFILTERS_H_

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/Math/counters.h"
#include "coretools/enum.h"
#include <cstddef>
#include <cstdint>
#include <string_view>
#include "TBamFilter.h"
namespace BAM{

enum class FilterType : size_t {
	min,
	ReadLength = min,
	MappedLength,
	MappingQuality,
	FragmentLength,
	maxRange,
	Duplicate = maxRange,
	SoftClippedRation,
	ImproperPairs,
	Unmapped,
	FailedQC,
	Secondary,
	Supplementary,
	LongerThanFragment,
	ReadGroup,
	FwdStrand,
	RevStrand,
	FirstMate,
	SecondMate,
	Blacklist,
	External,
	max
};

class TBamFilters {
	coretools::TStrongArray<coretools::TNumericRange<size_t>, FilterType, coretools::index(FilterType::maxRange)> _ranges;
	coretools::TStrongArray<TBamFilter, FilterType> _counters;
	coretools::TOutputFile* _log = nullptr;
	bool _enabled = false;

public:

	bool enabled() const noexcept {return _enabled;}

	void filterOut(FilterType Filter, std::string_view alignmentName, bool isSecondMate, size_t readGroup,
			  int64_t chromosomeID) {
		auto &filter = _counters[Filter];
		if (filter.filters()) {
			filter.filterOut(alignmentName, isSecondMate, readGroup, chromosomeID, _log);
		}
	}

	bool pass(FilterType Filter, bool Pass, std::string_view alignmentName, bool isSecondMate, size_t readGroup,
			  int64_t chromosomeID) {
		auto &filter = _counters[Filter];
		if (filter.filters() && !Pass) {
			filter.filterOut(alignmentName, isSecondMate, readGroup, chromosomeID, _log);
			return false;
		}
		return true;
	}
	bool pass(FilterType Filter, size_t Value, std::string_view alignmentName, bool isSecondMate, size_t readGroup,
	          int64_t chromosomeID) {
		assert(Filter < FilterType::maxRange);
		auto &filter = _counters[Filter];
		if (filter.filters() && !_ranges[Filter].within(Value)) {
			filter.filterOut(alignmentName, isSecondMate, readGroup, chromosomeID, _log);
			return false;
		}
		return true;
	}

	void enable(FilterType Filter, std::string_view Reason, size_t numRG, size_t numChrom) {
		_counters[Filter].enable(Reason, numRG, numChrom);
		_enabled = true;
		
	}
	void enable(FilterType Filter, const coretools::TNumericRange<size_t> & Range, std::string_view Reason, size_t numRG, size_t numChrom) {
		if (Filter >= FilterType::maxRange) DEVERROR("Cannot Do Rangefilter on Type ", coretools::index(Filter), "!");
		_ranges[Filter] = Range;
		enable(Filter, Reason, numRG, numChrom);
	}

	void disable(FilterType Filter) {
		_counters[Filter].disable();
	}

	const TBamFilter& operator[](FilterType t) const noexcept {return _counters[t];}
	TBamFilter& operator[](FilterType t) noexcept {return _counters[t];}

	const coretools::TNumericRange<size_t>& range(FilterType Filter) const noexcept {
		assert(Filter < FilterType::maxRange);
		return _ranges[Filter];
	}

	void setLog(coretools::TOutputFile& Log) {
		_log = &Log;
	}

	void fillHeader(std::vector<std::string> &Header) const {
		for (auto& f: _counters) f.fillHeader(Header);
	}

	void printCombinedCounts(coretools::TOutputFile &Out) const {
		for (auto& f: _counters) f.printCombinedCounts(Out);
	}

	void printCounts(coretools::TOutputFile &Out, size_t rg_ID) const {
		for (auto& f: _counters) f.printCounts(Out, rg_ID);
	}

	void summary(size_t Total, size_t ReadGroup) const {
		for (auto& f: _counters) f.summary(Total, ReadGroup);
	}
};
}

#endif
