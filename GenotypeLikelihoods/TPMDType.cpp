#include "TPMDType.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include <utility>

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

template<size_t End>
constexpr ReadEnd makeForward() {
	if constexpr (End == 3) return ReadEnd::forward3;
	if constexpr (End == 5) return ReadEnd::forward5;
	else static_assert(End == 3);
}

template<size_t End>
constexpr ReadEnd makeReverse() {
	if constexpr (End == 3) return ReadEnd::reverse3;
	if constexpr (End == 5) return ReadEnd::reverse5;
	else static_assert(End == 3);
}

template<size_t End, Base From, Base To>
auto makeFromTo(const PMDTable_RG &PMDTable) {
	// Assumption: From->To pattern is the same for forward and reverse reads from their respective Ends
	const auto N = PMDTable.front().size();
	std::vector<double> from_to;
	from_to.reserve(N);
	std::vector<double> to_from;
	to_from.reserve(N);

	constexpr auto forward = makeForward<End>();
	constexpr auto reverse = makeReverse<End>();

	for (size_t i = 0; i < N; ++i) {
		// CT
		from_to.push_back(PMDTable[forward][i][From][To]);
		from_to.back() += PMDTable[reverse][i][From][To];
		if (from_to.back() < 100) {
			if (i < 10) {
				UERROR("Not sufficient ", From, "-", To, " data to estimate PMD model at position ", i, ": ", from_to.back(),
					   ", the first 10 positions must have > 100 data points!\nConsider pooling read groups (parameter "
					   "poolReadGroups).");
			}
			from_to.pop_back();
			break;
		}
		double s_from = 0;
		double s_to   = 0;
		for (auto b = Base::min; b < Base::max; ++b) {
			s_from += PMDTable[forward][i][From][b];
			s_from += PMDTable[reverse][i][From][b];
			s_to += PMDTable[forward][i][To][b];
			s_to += PMDTable[reverse][i][To][b];
		}
		from_to.back() /= s_from;

		to_from.push_back(PMDTable[forward][i][To][From]);
		to_from.back() += PMDTable[reverse][i][To][From];
		to_from.back() /= s_to;
	}
	return std::make_pair(from_to, to_from);
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
	logfile().startIndent("Learning C-T pattern:");
	const auto [C_T, T_C] = impl::makeFromTo<5, Base::C, Base::T>(PMDTable);
	_pmdCT->learn(C_T, T_C);
	logfile().endIndent();

	logfile().startIndent("Learning G-A pattern:");
	const auto [G_A, A_G] = impl::makeFromTo<3, Base::G, Base::A>(PMDTable);
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
	logfile().startIndent("Learning 5' C-T pattern:");
	const auto [C_T5, T_C5] = impl::makeFromTo<5, Base::C, Base::T>(PMDTable);
	_pmdCT5->learn(C_T5, T_C5);
	logfile().endIndent();

	logfile().startIndent("Learning 3' C-T pattern:");
	const auto [C_T3, T_C3] = impl::makeFromTo<3, Base::C, Base::T>(PMDTable);
	_pmdCT3->learn(C_T3, T_C3);
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
