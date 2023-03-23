#include "TPMDType.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods {

using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using genometools::Base;
using genometools::Genotype;
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

TBaseProbabilities massFunction(Genotype g, coretools::Probability pCT, coretools::Probability pGA,
								const TBaseLikelihoods &baseLikelihoodsNoPMD) {
	using namespace genometools;
	std::array<Base, 2> bases{first(g), second(g)};
	TBaseData mf{0};
	for (const auto a : bases) {
		switch (a) {
		case Base::A: mf[Base::A] += baseLikelihoodsNoPMD[Base::A]; break;
		case Base::C: {
			mf[Base::C] += (1. - pCT) * baseLikelihoodsNoPMD[Base::C];
			mf[Base::T] += pCT * baseLikelihoodsNoPMD[Base::T];
		} break;
		case Base::G: {
			mf[Base::A] += pGA * baseLikelihoodsNoPMD[Base::A];
			mf[Base::G] += (1. - pGA) * baseLikelihoodsNoPMD[Base::G];
		} break;
		default: mf[Base::T] += baseLikelihoodsNoPMD[Base::T];
		}
	}
	return TBaseProbabilities::normalize(mf);
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

TBaseProbabilities TPMDTypeNone::massFunction(genometools::Genotype g, const TBaseLikelihoods &baseLikelihoodsNoPMD) {
		using namespace genometools;
		const Base a = first(g);
		const Base b = second(g);
		TBaseLikelihoods mf{0};
		mf[a] = baseLikelihoodsNoPMD[a];
		mf[b] = baseLikelihoodsNoPMD[b];
		return TBaseProbabilities::normalize(mf);
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

void TPMDTypeDoubleStrand::estimate(const PMDTable_RG &PMDTable) {
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: C->T pattern is the same for forward and reverse reads from their respective 5-prime ends.
	const auto N = PMDTable.front().size();
	std::vector<double> C_T;
	C_T.reserve(N);
	std::vector<double> T_C;
	C_T.reserve(N);
	std::vector<double> G_A;
	G_A.reserve(N);
	std::vector<double> A_G;
	A_G.reserve(N);
	
	for (size_t i = 0; i < N; ++i) {
		// CT
		C_T.push_back(PMDTable[ReadEnd::forward5][Base::C][Base::T][i]);
		C_T.back() += PMDTable[ReadEnd::reverse5][Base::C][Base::T][i];
		if (C_T.back() < 100) {
			if (i < 10) {
				UERROR("Not sufficient C_T data to estimate PMD model at position ", i, ": ", C_T.back(),
					   ", the first 10 positions must have > 100 data points!\nConsider pooling read groups (parameter "
					   "poolReadGroups).");
			}
			C_T.pop_back();
			break;
		}
		C_T.back() /= (PMDTable[ReadEnd::forward5].sums(Base::C)[i] + PMDTable[ReadEnd::reverse5].sums(Base::C)[i]);

		T_C.push_back(PMDTable[ReadEnd::forward5][Base::T][Base::C][i]);
		T_C.back() += PMDTable[ReadEnd::reverse5][Base::T][Base::C][i];
		T_C.back() /= (PMDTable[ReadEnd::forward5].sums(Base::T)[i]
					   + PMDTable[ReadEnd::reverse5].sums(Base::T)[i]);

		// GA 
		G_A.push_back(PMDTable[ReadEnd::forward3][Base::G][Base::A][i]);
		G_A.back() += PMDTable[ReadEnd::reverse3][Base::G][Base::A][i];
		if (G_A.back() < 100) {
			if (i < 10) {
				UERROR("Not sufficient G_A data to estimate PMD model at position ", i, ": ", G_A.back(),
					   ", the first 10 positions must have > 100 data points!\nConsider pooling read groups (parameter "
					   "poolReadGroups).");
			}
			G_A.pop_back();
			break;
		}
		G_A.back() /= (PMDTable[ReadEnd::forward3].sums(Base::G)[i] + PMDTable[ReadEnd::reverse3].sums(Base::G)[i]);

		A_G.push_back(PMDTable[ReadEnd::forward3][Base::A][Base::G][i]);
		A_G.back() += PMDTable[ReadEnd::reverse3][Base::A][Base::G][i];
		A_G.back() /= (PMDTable[ReadEnd::forward3].sums(Base::A)[i]
					   + PMDTable[ReadEnd::reverse3].sums(Base::A)[i]);
	}
	logfile().startIndent("Learning C-T pattern:");
	_pmdCT->learn(C_T, T_C);
	logfile().endIndent();
	logfile().startIndent("Learning G-A pattern:");
	_pmdGA->learn(G_A, A_G);
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

TBaseProbabilities TPMDTypeDoubleStrand::massFunction(genometools::Genotype g, const BAM::TSequencedBase &data,
														 const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::massFunction(g, _probCT(data), _probGA(data), baseLikelihoodsNoPMD);
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
	_pmdCT5.reset(makeFunction(Details[1]));
	_pmdCT3.reset(makeFunction(Details[2]));
}

void TPMDTypeSingleStrand::estimate(const PMDTable_RG &PMDTable) {
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: 5-prime C->T pattern is the same for forward and reverse reads from their respective
	// 5-prime ends.
	/*
	TPMDTable from5(PMDTable[ReadEnd::forward5]);
	from5.add(PMDTable[ReadEnd::reverse5]);

	logfile().startIndent("Learning 5' C-T pattern:");
	_pmdCT5->learn(from5, Base::C, Base::T);
	logfile().endIndent();

	// Assumption: 3-prime C->T pattern is the same for forward and reverse reads from their
	// respective 3-prime ends.
	TPMDTable from3(PMDTable[ReadEnd::forward3]);
	from3.add(PMDTable[ReadEnd::reverse3]);

	logfile().startIndent("Learning 3' C-T pattern:");
	_pmdCT3->learn(from3, Base::C, Base::T);
	logfile().endIndent();
	*/


}

TBaseLikelihoods TPMDTypeSingleStrand::baseLikelihoods(const BAM::TSequencedBase &data,
					       const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::baseLikelihoods(_probCT(data), _probGA(data), baseLikelihoodsNoPMD);
}

TBaseProbabilities TPMDTypeSingleStrand::massFunction(Base b, const BAM::TSequencedBase &data,
														 const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::massFunction(b, _probCT(data), _probGA(data), baseLikelihoodsNoPMD);
}

TBaseProbabilities TPMDTypeSingleStrand::massFunction(genometools::Genotype g, const BAM::TSequencedBase &data,
														 const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return impl::massFunction(g, _probCT(data), _probGA(data), baseLikelihoodsNoPMD);
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
