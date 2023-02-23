#include "TPMDType.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"

namespace GenotypeLikelihoods {

using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using genometools::Base;

using namespace coretools::str;
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
	// Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	// no PMD for A and C
	TBaseLikelihoods baseLikelihoods(baseLikelihoodsNoPMD);

	// get relevant PMD probabilities
	const auto pmdProb_CT = _probCT(data);
	const auto pmdProb_GA = _probGA(data);

	// add PMD
	baseLikelihoods[Base::C] =
		(1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[Base::C] + pmdProb_CT * baseLikelihoodsNoPMD[Base::T];
	baseLikelihoods[Base::G] =
		(1.0 - pmdProb_GA) * baseLikelihoodsNoPMD[Base::G] + pmdProb_GA * baseLikelihoodsNoPMD[Base::A];
	return baseLikelihoods;
}

TBaseProbabilities TPMDTypeDoubleStrand::massFunction(Base b, const BAM::TSequencedBase &data,
														 const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	switch (b) {
	case Base::A: return TBaseProbabilities::normalize({1., 0., 0., 0.});
	case Base::C: {
		const auto pmdProb_CT = _probCT(data);
		return TBaseProbabilities::normalize(
			{0., (1. - pmdProb_CT) * baseLikelihoodsNoPMD[Base::C], 0., pmdProb_CT * baseLikelihoodsNoPMD[Base::T]});
	}
	case Base::G: {
		const auto pmdProb_GA = _probGA(data);
		return TBaseProbabilities::normalize(
			{pmdProb_GA * baseLikelihoodsNoPMD[Base::A], 0., (1. - pmdProb_GA) * baseLikelihoodsNoPMD[Base::G], 0.});
	}
	default: return TBaseProbabilities::normalize({0., 0., 0., 1.}); // case Base::T
	}
}

void TPMDTypeDoubleStrand::simulate(BAM::TSequencedBase &data) const {
	if (data.base == Base::C) {
		if (randomGenerator().getRand() < _probCT(data)) data.base = Base::T;
	} else if (data.base == Base::G) {
		if (randomGenerator().getRand() < _probGA(data)) data.base = Base::A;
	}
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
	// Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	// no PMD for A, C and G
	TBaseLikelihoods baseLikelihoods(baseLikelihoodsNoPMD);

	const double pmdProb_CT  = _probCT(data);
	baseLikelihoods[Base::C] =
		(1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[Base::C] + pmdProb_CT * baseLikelihoodsNoPMD[Base::T];

	return baseLikelihoods;
}

TBaseProbabilities TPMDTypeSingleStrand::massFunction(Base b, const BAM::TSequencedBase &data,
														 const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	switch (b) {
	case Base::A: return TBaseProbabilities::normalize({1., 0., 0., 0.});
	case Base::C: {
		const double pmdProb_CT = _probCT(data);
		return TBaseProbabilities::normalize(
			{0., (1. - pmdProb_CT) * baseLikelihoodsNoPMD[Base::C], 0., pmdProb_CT * baseLikelihoodsNoPMD[Base::T]});
	}
	case Base::G: return TBaseProbabilities::normalize({0., 0., 1., 0.});
	default: return TBaseProbabilities::normalize({0., 0., 0., 1.}); // Base::T
	}
}

void TPMDTypeSingleStrand::simulate(BAM::TSequencedBase &data) const {
	if (data.base == Base::C) {
		if (randomGenerator().getRand() < _probCT(data)) data.base = Base::T;
	}
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
