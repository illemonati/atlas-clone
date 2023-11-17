
#ifndef TBASEFILTER_H
#define TBASEFILTER_H

#include "TSequencedBase.h"
#include "coretools/Math/TNumericRange.h"
namespace GenomeTasks {

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter {
private:
	coretools::TNumericRange<genometools::PhredIntProbability> _range;

public:
	TQualityFilter();

	bool pass(const BAM::TSequencedBase & base) const {
		return _range.within(base.originalQuality_phredInt);
	};
};

//-------------------------------------
// TContextFilter
//-------------------------------------
class TContextFilter {
private:
	//std::array<bool, static_cast<uint8_t>(genometools::BaseContext::max) + 1> _keptContexts;
	bool _filter;
	coretools::TStrongArray<bool, genometools::BaseContext, coretools::index(genometools::BaseContext::max) + 1> _keptContexts;

public:
	explicit TContextFilter();
	bool pass(const BAM::TSequencedBase & base) const noexcept {
		return _keptContexts[base.context()];
	}
	constexpr operator bool() const noexcept {
		return _filter;
	}
};

}

#endif
