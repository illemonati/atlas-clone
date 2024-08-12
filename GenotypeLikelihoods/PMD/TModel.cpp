#include "TModel.h"
#include "TSequencedBase.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods::PMD {

using genometools::Base;
using genometools::Genotype;
using BAM::TSequencedBase;
using coretools::P;
using coretools::Probability;

TBaseProbabilities TNoPMD::P_bbar(Base b, const TSequencedBase &,
								  const TBaseLikelihoods &) const noexcept {
	static constexpr TBaseMassFunctions P_bbars{
	    {TBaseProbabilities::normalize({1., 0., 0., 0.}), TBaseProbabilities::normalize({0., 1., 0., 0.}),
	     TBaseProbabilities::normalize({0., 0., 1., 0.}), TBaseProbabilities::normalize({0., 0., 0., 1.})}};
	return P_bbars[b];
}
TBaseProbabilities TNoPMD::P_bbar(Genotype g, const TSequencedBase &,
								  const TBaseLikelihoods &P_dij_bbar) const noexcept {
	const Base a = first(g);
	const Base b = second(g);
	TBaseLikelihoods Psum{P(0.)};
	Psum[a] = P_dij_bbar[a];
	Psum[b] = P_dij_bbar[b];
	return TBaseProbabilities::normalize(Psum);
}

TBaseBaseProbabilities TNoPMD::P_b_bbar(genometools::Genotype g, const BAM::TSequencedBase &,
								const TBaseLikelihoods &P_dij_bbar) const noexcept {
	TBaseBaseProbabilities bbProbs{};
	double sum = 0.;
	for (const auto b: {genometools::first(g), genometools::second(g)}) {
		const auto bbar  = b; // only this is possible
		bbProbs[b][bbar] = P_dij_bbar[bbar];
		sum             += bbProbs[b][bbar];
	}
	for (auto& bProbs: bbProbs) for (auto& prob: bProbs) prob.scale(sum);
	return bbProbs;
}

TBaseProbabilities TWithPMD::P_bbar(Base b, const TSequencedBase &data,
									const TBaseLikelihoods &P_dij_bbar) const noexcept {
	// Assuming genotype g = bb
	// This is also, due to normalization, P_b_bbar with genotype g = bb
	switch (b) {
	case Base::A: return TBaseProbabilities::normalize({1., 0., 0., 0.});
	case Base::C: {
		const auto pCT = _psi.prob<Type::CT>(data);
		return TBaseProbabilities::normalize({0., (1. - pCT) * P_dij_bbar[Base::C], 0., pCT * P_dij_bbar[Base::T]});
	}
	case Base::G: {
		const auto pGA = _psi.prob<Type::GA>(data);
		return TBaseProbabilities::normalize({pGA * P_dij_bbar[Base::A], 0., (1. - pGA) * P_dij_bbar[Base::G], 0.});
	}
	default: return TBaseProbabilities::normalize({0., 0., 0., 1.}); // case Base::T
	}
}

namespace impl {
double bbProbs(Base b, TBaseBaseProbabilities &probs, Probability pCT, Probability pGA,
			   const TBaseLikelihoods &P_dij_bbar) {
	double sum = 0.;
	switch (b) {
	case Base::A: {
		sum += probs[Base::A][Base::A] = P_dij_bbar[Base::A];
		break;
	}
	case Base::C: {
		sum += probs[Base::C][Base::C] = pCT.complement() * P_dij_bbar[Base::C];
		sum += probs[Base::C][Base::T] = pCT * P_dij_bbar[Base::T];
	} break;
	case Base::G: {
		sum += probs[Base::G][Base::G] = pGA.complement() * P_dij_bbar[Base::G];
		sum += probs[Base::G][Base::A] = pGA * P_dij_bbar[Base::A];
	} break;
	default: {
		sum += probs[Base::T][Base::T] = P_dij_bbar[Base::T];
		break;
	}
	}
	return sum;
}
} // namespace impl

TBaseProbabilities TWithPMD::P_bbar(Genotype g, const TSequencedBase &data,
									const TBaseLikelihoods &P_dij_bbar) const noexcept {
	const auto pCT = _psi.prob<Type::CT>(data);
	const auto pGA = _psi.prob<Type::GA>(data);

	TBaseData Psum{0};
	for (const auto b : {first(g), second(g)}) {
		switch (b) {
		case Base::A: Psum[Base::A] += P_dij_bbar[Base::A]; break;
		case Base::C: {
			Psum[Base::C] += (1. - pCT) * P_dij_bbar[Base::C];
			Psum[Base::T] += pCT * P_dij_bbar[Base::T];
		} break;
		case Base::G: {
			Psum[Base::A] += pGA * P_dij_bbar[Base::A];
			Psum[Base::G] += (1. - pGA) * P_dij_bbar[Base::G];
		} break;
		default: Psum[Base::T] += P_dij_bbar[Base::T];
		}
	}
	return TBaseProbabilities::normalize(Psum);
}

TBaseBaseProbabilities TWithPMD::P_b_bbar(Genotype g, const BAM::TSequencedBase &data,
								const TBaseLikelihoods &P_dij_bbar) const noexcept {
	const auto pCT = _psi.prob<Type::CT>(data);
	const auto pGA = _psi.prob<Type::GA>(data);

	TBaseBaseProbabilities bbProbs{};
	const double sum = genometools::isHomozygous(g)
						   ? impl::bbProbs(genometools::first(g), bbProbs, pCT, pGA, P_dij_bbar)
						   : impl::bbProbs(genometools::first(g), bbProbs, pCT, pGA, P_dij_bbar) +
								 impl::bbProbs(genometools::second(g), bbProbs, pCT, pGA, P_dij_bbar);
	for (auto& bProbs: bbProbs) for (auto& prob: bProbs) prob.scale(sum);
	return bbProbs;
}

TBaseLikelihoods TWithPMD::P_dij(const TSequencedBase &data, const TBaseLikelihoods &P_dij_bbar) const noexcept {
	const auto pCT = _psi.prob<Type::CT>(data);
	const auto pGA = _psi.prob<Type::GA>(data);

	TBaseLikelihoods Ps(P_dij_bbar);
	Ps[Base::C] = P((1.0 - pCT) * P_dij_bbar[Base::C] + pCT * P_dij_bbar[Base::T]);
	Ps[Base::G] = P((1.0 - pGA) * P_dij_bbar[Base::G] + pGA * P_dij_bbar[Base::A]);
	return Ps;
}

void TWithPMD::simulate(BAM::TAlignment &aln) const {
	using coretools::instances::randomGenerator;
	using genometools::Base;
	for (auto &d : aln) {
			auto &base = d.base;

			if (base == Base::C) {
				const auto pCT = _psi.prob<Type::CT>(d);
				if (randomGenerator().getRand() < pCT) base = Base::T;
			} else if (base == Base::G) {
				const auto pGA = _psi.prob<Type::GA>(d);
				if (randomGenerator().getRand() < pGA) base = Base::A;
			}
	}
}

} // namespace GenotypeLikelihoods::PMD
