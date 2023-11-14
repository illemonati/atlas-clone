
#ifndef TBASEFILTER_H
#define TBASEFILTER_H

#include "TSequencedBase.h"
#include "coretools/Math/TNumericRange.h"
namespace GenomeTasks {

//-------------------------------------
// TBaseFilter
//-------------------------------------
class TBaseFilter{
protected:
	bool _filter;

public:
	explicit constexpr TBaseFilter() : _filter(false) {}
	virtual ~TBaseFilter() = default;

	constexpr operator bool() const{
		return _filter;
	}

	virtual bool pass(const BAM::TSequencedBase & base) const = 0;
};

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter : public TBaseFilter{
private:
	coretools::TNumericRange<genometools::PhredIntProbability> _range;

public:
	TQualityFilter();

	bool pass(const BAM::TSequencedBase & base) const override{
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
	bool pass(const BAM::TSequencedBase & base) const override;
};

}

#endif
