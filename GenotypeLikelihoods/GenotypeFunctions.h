#ifndef GENOTYPEFUNCTIONS_H_
#define GENOTYPEFUNCTIONS_H_

#include <numeric>

#include "genometools/Genotypes/Containers.h"

namespace GenotypeLikelihoods{

inline genometools::TGenotypeProbabilities posterior(const genometools::TGenotypeLikelihoods &likelihoods, const genometools::TGenotypeProbabilities &prior) {
	return genometools::TGenotypeProbabilities::normalize(likelihoods, prior, std::multiplies<>());
}

//template<typename Index, size_t N=Index::max, typename Likelihood=Probability>
//TStrongArraySum1<Probability, Index, N> posterior(TStrongArray<Likelihood, Index, N> Likelihoods, TStrongArraySum1<Probability, Index, N> prior)

inline genometools::TBaseLikelihoods fromError(genometools::Base trueBase, coretools::Probability error) {
	using coretools::P;
	genometools::TBaseLikelihoods ret;
	ret.fill(P(error/3.));
	ret[trueBase] = error.complement();
	return ret;
}


constexpr coretools::Probability homozygous(const genometools::TGenotypeProbabilities &ps) {
	using GT = genometools::Genotype;
	using coretools::P;
	return P(ps[GT::AA] + ps[GT::CC] + ps[GT::GG] + ps[GT::TT]);
}
constexpr coretools::Probability heterozygous(const genometools::TGenotypeProbabilities &ps) {
	return homozygous(ps).complement();
}

template<typename Container>
auto frequencies(const Container &vs) {
	const auto tot = std::accumulate(vs.begin(), vs.end(), typename Container::value_type{});
	coretools::TStrongArray<coretools::Probability, typename Container::index_type, Container::capacity> ret;
	std::transform(vs.cbegin(), vs.cend(), ret.begin(), [tot](auto v){return v / tot;});
	return ret;
}
}

#endif
