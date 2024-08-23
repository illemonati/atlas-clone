#ifndef TBASEFILTER_H
#define TBASEFILTER_H

#include "coretools/Math/TNumericRange.h"

#include "TSequencedBase.h"

namespace GenomeTasks {

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter {
private:
	coretools::TNumericRange<coretools::PhredInt> _range;

public:
	TQualityFilter();

	bool pass(const BAM::TSequencedBase & base) const {
		return _range.within(base.originalQuality);
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
