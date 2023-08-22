#include "TModel.h"
#include "TSequencedBase.h"
#include "coretools/devtools.h"

namespace GenotypeLikelihoods::PMD {
TBaseProbabilities TNoPMD::P_bbar(genometools::Base b, const BAM::TSequencedBase &,
								  const TBaseLikelihoods &) const noexcept {
	static constexpr TBaseMassFunctions P_bbars{
	    {TBaseProbabilities::normalize({1., 0., 0., 0.}), TBaseProbabilities::normalize({0., 1., 0., 0.}),
	     TBaseProbabilities::normalize({0., 0., 1., 0.}), TBaseProbabilities::normalize({0., 0., 0., 1.})}};
	return P_bbars[b];
}
TBaseProbabilities TNoPMD::P_bbar(genometools::Genotype g, const BAM::TSequencedBase &,
								  const TBaseLikelihoods &P_dij_bbar) const noexcept {
	using namespace genometools;
	const Base a = first(g);
	const Base b = second(g);
	TBaseLikelihoods mf{0};
	mf[a] = P_dij_bbar[a];
	mf[b] = P_dij_bbar[b];
	return TBaseProbabilities::normalize(mf);
}

TBaseProbabilities TWithPMD::P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept {
	OUT(b, data.base, P_dij_bbar);
	return TBaseProbabilities::normalize({1., 0., 0., 0.});
}
TBaseProbabilities TWithPMD::P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept {
	OUT(g, data.base, P_dij_bbar);
	return TBaseProbabilities::normalize({1., 0., 0., 0.});
}

TBaseLikelihoods TWithPMD::P_dij(const BAM::TSequencedBase &data, const TBaseLikelihoods &P_dij_bbar) const noexcept {
	using genometools::Base;
	const auto pCT = _psi.prob<Type::CT>(data);
	const auto pGA = _psi.prob<Type::GA>(data);
	TBaseLikelihoods ret(P_dij_bbar);

	ret[Base::C] = (1.0 - pCT) * P_dij_bbar[Base::C] + pCT * P_dij_bbar[Base::T];
	ret[Base::G] = (1.0 - pGA) * P_dij_bbar[Base::G] + pGA * P_dij_bbar[Base::A];
	return ret;
}

} // namespace GenotypeLikelihoods::PMD
