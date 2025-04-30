#ifndef BAM_TBAMFILEFILTERS_H_
#define BAM_TBAMFILEFILTERS_H_

#include <string_view>

#include "TAlignmentList.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TNumericRange.h"

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
	// Use clone-function
	coretools::TStrongArray<TBamFilter, FilterType> _filters;

	// Do not clone!
	coretools::TOutputFile _log;
	size_t _numRG         = 0;
	size_t _numChrom      = 0;

	// Clone
	coretools::TStrongArray<coretools::TNumericRange<size_t>, FilterType, coretools::index(FilterType::maxRange)> _ranges;
	TAlignmentList _blacklist;
	double _softClipRatio = 1.;
	bool _enabled         = false;

public:
	TBamFilters(bool Enable = false);

	void clone(const TBamFilters& Filters);
	bool enabled() const noexcept {return _enabled;}

	double softClipRation() const noexcept {return _softClipRatio;}
	const TAlignmentList& blacklist() const noexcept {return _blacklist;}

	void resize(size_t numRG, size_t numChrom, std::string_view Filename);

	void filterOut(FilterType Filter, std::string_view alignmentName, bool isSecondMate, size_t readGroup,
			  int64_t chromosomeID) {
		assert(_numRG > 0);
		auto &filter = _filters[Filter];
		if (filter) {
			filter.filterOut(alignmentName, isSecondMate, readGroup, chromosomeID, _log);
		}
	}

	bool pass(FilterType Filter, bool Pass, std::string_view alignmentName, bool isSecondMate, size_t readGroup,
			  size_t chromosomeID) {
		assert(_numRG > 0);
		auto &filter = _filters[Filter];
		if (filter && !Pass) {
			filter.filterOut(alignmentName, isSecondMate, readGroup, chromosomeID, _log);
			return false;
		}
		return true;
	}
	bool pass(FilterType Filter, size_t Value, std::string_view alignmentName, bool isSecondMate, size_t readGroup,
	          size_t chromosomeID) {
		assert(_numRG > 0);
		assert(Filter < FilterType::maxRange);
		auto &filter = _filters[Filter];
		if (filter && !_ranges[Filter].within(Value)) {
			filter.filterOut(alignmentName, isSecondMate, readGroup, chromosomeID, _log);
			return false;
		}
		return true;
	}

	void enable(FilterType Filter, std::string_view Reason) {
		_filters[Filter].enable(Reason, _numRG, _numChrom);
		_enabled = true;
		
	}
	void enable(FilterType Filter, const coretools::TNumericRange<size_t> & Range, std::string_view Reason) {
		if (Filter >= FilterType::maxRange) DEVERROR("Cannot Do Rangefilter on Type ", coretools::index(Filter), "!");
		_ranges[Filter] = Range;
		enable(Filter, Reason);
	}

	void disable(FilterType Filter) {
		_filters[Filter].disable();
	}

	const TBamFilter& operator[](FilterType t) const noexcept {return _filters[t];}
	TBamFilter& operator[](FilterType t) noexcept {return _filters[t];}

	const coretools::TNumericRange<size_t>& range(FilterType Filter) const noexcept {
		assert(Filter < FilterType::maxRange);
		return _ranges[Filter];
	}

	void fillHeader(std::vector<std::string> &Header) const {
		for (auto& f: _filters) if (f) f.fillHeader(Header);
	}

	void printCombinedCounts(coretools::TOutputFile &Out) const {
		for (auto& f: _filters) if (f) f.printCombinedCounts(Out);
	}

	void printCounts(coretools::TOutputFile &Out, size_t rg_ID) const {
		for (auto& f: _filters) if (f) f.printCounts(Out, rg_ID);
	}

	bool summary(size_t Total, size_t ReadGroup) const {
		bool summaryWritten = false;
		for (auto& f: _filters){
			if (f){
				summaryWritten = summaryWritten || f.summary(Total, ReadGroup);
			} 
		} 
		return summaryWritten; //return true if at least one summary was written
	}
};
} // namespace BAM

#endif
