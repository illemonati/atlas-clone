#include "TModel.h"

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
}
