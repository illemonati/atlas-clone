#ifndef TBASEFILTER_H
#define TBASEFILTER_H

#include "coretools/Math/TNumericRange.h"

#include "TSequencedData.h"

namespace BAM {

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter {
private:
	coretools::TNumericRange<coretools::PhredInt> _range;

public:
	TQualityFilter();

	bool pass(const BAM::TSequencedData & data) const {
		return _range.within(data.originalQuality);
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
	bool pass(const BAM::TSequencedData & data) const noexcept {
		return _keptContexts[data.context()];
	}
	constexpr operator bool() const noexcept {
		return _filter;
	}
};

}

#endif
