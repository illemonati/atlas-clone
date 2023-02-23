#include "TPMDType.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods {

using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using genometools::Base;
using namespace coretools::str;

namespace impl {
TBaseLikelihoods baseLikelihoods(coretools::Probability pCT, coretools::Probability pGA,
								 const TBaseLikelihoods &baseLikelihoodsNoPMD) {
	TBaseLikelihoods baseLikelihoods(baseLikelihoodsNoPMD);
	baseLikelihoods[Base::C] =
		(1.0 - pCT) * baseLikelihoodsNoPMD[Base::C] + pCT * baseLikelihoodsNoPMD[Base::T];
	baseLikelihoods[Base::G] =
		(1.0 - pGA) * baseLikelihoodsNoPMD[Base::G] + pGA * baseLikelihoodsNoPMD[Base::A];
	return baseLikelihoods;
}

TBaseProbabilities massFunction(Base b, coretools::Probability pCT, coretools::Probability pGA,
								const TBaseLikelihoods &baseLikelihoodsNoPMD) {
	switch (b) {
	case Base::A: return TBaseProbabilities::normalize({1., 0., 0., 0.});
	case Base::C: {
		return TBaseProbabilities::normalize(
			{0., (1. - pCT) * baseLikelihoodsNoPMD[Base::C], 0., pCT * baseLikelihoodsNoPMD[Base::T]});
	}
	case Base::G: {
		return TBaseProbabilities::normalize(
			{pGA * baseLikelihoodsNoPMD[Base::A], 0., (1. - pGA) * baseLikelihoodsNoPMD[Base::G], 0.});
	}
	default: return TBaseProbabilities::normalize({0., 0., 0., 1.}); // case Base::T
	}
}
Base simulate(Base base, coretools::Probability pCT, coretools::Probability pGA) {
	if (base == Base::C) {
		if (randomGenerator().getRand() < pCT) return Base::T;
	}
	else if (base == Base::G) {
		if (randomGenerator().getRand() < pGA) return Base::A;
	}
	// else
	return base;
}
}

TPMDTypeDoubleStrand::TPMDTypeDoubleStrand(const std::vector<std::string> &Details) {
	// expect three elements: type, pmdCT, pmdGA
	constexpr size_t nDetails = 3;
	if (Details.size() != nDetails) {
		UERROR("Cannot initialize PMD type ", name, ": expect ", nDetails, " entries but found ", Details.size(), "!",
			   "\nProvided string: '", concatenateString(Details, ':'), "'.", "\nExpect string of the form '", name,
			   "':functionCT:functionGA'.");
	}
	_pmdCT.reset(makeFunction(Details[1]));
	_pmdGA.reset(makeFunction(Details[2]));
}

void TPMDTypeDoubleStrand::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters){
	_pmdCT->parseEstimationParameters(EstimationParameters);
	_pmdGA->parseEstimationParameters(EstimationParameters);
}

void TPMDTypeDoubleStrand::estimate(const PMDTable_RG &PMDTable,
				    const TPMDEstimationParameters &EstimationParameters) {
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: C->T pattern is the same for forward and reverse reads from their respective 5-prime ends.
	TPMDTable from5(PMDTable[ReadEnd::forward5]);
	from5.add(PMDTable[ReadEnd::reverse5]);
	logfile().startIndent("Learning C-T pattern:");
	_pmdCT->learn(from5, Base::C, Base::T, EstimationParameters);
	logfile().endIndent();

	// Assumption: G->A pattern is the same for forward and reverse reads from their respective 3-prime ends.
	TPMDTable from3(PMDTable[ReadEnd::forward3]);
	from3.add(PMDTable[ReadEnd::reverse3]);
	logfile().startIndent("Learning G-A pattern:");
	_pmdGA->learn(from3, Base::G, Base::A, EstimationParameters);
	logfile().endIndent();
}


TBaseLikelihoods TPMDTypeDoubleStrand::baseLikelihoods(const BAM::TSequencedBase &data,
					       const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::baseLikelihoods(_probCT(data), _probGA(data), baseLikelihoodsNoPMD);
}

TBaseProbabilities TPMDTypeDoubleStrand::massFunction(Base b, const BAM::TSequencedBase &data,
														 const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::massFunction(b, _probCT(data), _probGA(data), baseLikelihoodsNoPMD);
}

void TPMDTypeDoubleStrand::simulate(BAM::TSequencedBase &data) const {
	data.base = impl::simulate(data.base, _probCT(data), _probGA(data));
}
//------------------------------------------------------
// TPMDSingleStrand
//------------------------------------------------------

TPMDTypeSingleStrand::TPMDTypeSingleStrand(const std::vector<std::string> &Details) {
	// expect 2 elements: type, pmdCT
	constexpr size_t nDetails = 3;
	if (Details.size() != nDetails) {
		UERROR("Cannot initialize PMD type ", name, ": expect ", nDetails, " entries but found ", Details.size(), "!",
			   "\nProvided string: '", concatenateString(Details, ':'), "'.", "\nExpect string of the form '", name,
			   "':functionCT:functionGA'.");
	}
	_pmdCT3.reset(makeFunction(Details[1]));
	_pmdCT5.reset(makeFunction(Details[2]));
}

void TPMDTypeSingleStrand::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) {
	_pmdCT3->parseEstimationParameters(EstimationParameters);
	_pmdCT5->parseEstimationParameters(EstimationParameters);
}

void TPMDTypeSingleStrand::estimate(const PMDTable_RG &PMDTable,
				    const TPMDEstimationParameters &EstimationParameters) {
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: 5-prime C->T pattern is the same for forward and reverse reads from their respective
	// 5-prime ends.
	TPMDTable from5(PMDTable[ReadEnd::forward5]);
	from5.add(PMDTable[ReadEnd::reverse5]);

	logfile().startIndent("Learning 5' C-T pattern:");
	_pmdCT5->learn(from5, Base::C, Base::T, EstimationParameters);
	logfile().endIndent();

	// Assumption: 3-prime C->T pattern is the same for forward and reverse reads from their
	// respective 3-prime ends.
	TPMDTable from3(PMDTable[ReadEnd::forward3]);
	from3.add(PMDTable[ReadEnd::reverse3]);

	logfile().startIndent("Learning 3' C-T pattern:");
	_pmdCT3->learn(from3, Base::C, Base::T, EstimationParameters);
	logfile().endIndent();
}

TBaseLikelihoods TPMDTypeSingleStrand::baseLikelihoods(const BAM::TSequencedBase &data,
					       const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::baseLikelihoods(_probCT(data), _probGA(data), baseLikelihoodsNoPMD);
}

TBaseProbabilities TPMDTypeSingleStrand::massFunction(Base b, const BAM::TSequencedBase &data,
														 const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::massFunction(b, _probCT(data), _probGA(data), baseLikelihoodsNoPMD);
}

void TPMDTypeSingleStrand::simulate(BAM::TSequencedBase &data) const {
	data.base = impl::simulate(data.base, _probCT(data), _probGA(data));
}

TPMDType *makeType(std::string_view pmdString) {
	// split by ':'
	std::vector<std::string> details;
	fillContainerFromString(pmdString, details, ":");

	// switch type
	if (details[0] == TPMDTypeNone::name) return new TPMDTypeNone;
	if (details[0] == TPMDTypeSingleStrand::name) return new TPMDTypeSingleStrand(details);
	if (details[0] == TPMDTypeDoubleStrand::name) return new TPMDTypeDoubleStrand(details);

	UERROR("Cannot initialize PMD: unknown PMD type '", details[0], "'!\nUse ", TPMDTypeNone::name, " or ",
		   TPMDTypeSingleStrand::name, " or ", TPMDTypeDoubleStrand::name, ".");
}
}
