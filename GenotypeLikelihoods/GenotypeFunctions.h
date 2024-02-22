#ifndef GENOTYPEFUNCTIONS_H_
#define GENOTYPEFUNCTIONS_H_

#include "GenotypeData.h"
#include <cstddef>

namespace GenotypeLikelihoods{

template<template<typename...> typename Container, typename... Args>
TGenotypeLikelihoods getGLH(const Container<TBaseLikelihoods, Args...> &bases, const size_t size) {
	using genometools::Base;
	using GT = genometools::Genotype;
	using genometools::genotype;
	using coretools::P;
	static_assert(std::is_enum_v<Base>);
	// allows for vector to be longer than what is to be used
	// do in log if depth is high
	if (bases.size() > 50) {
		TGenotypeData tmp{0.};
		for (size_t i = 0; i < size; ++i) {
			for (auto b1 = Base::min; b1 < Base::max; ++b1) {
				tmp[genotype(b1, b1)] += log(bases[i][b1]);
				for (auto b2 = coretools::next(b1); b2 < Base::max; ++b2) {
					tmp[genotype(b1, b2)] += log(0.5 * ((double)bases[i][b1] + (double)bases[i][b2]));
				}
			}
		}

		// standardize and de-log
		const auto max = *std::max_element(tmp.begin(), tmp.end());
		TGenotypeLikelihoods ret;
		for (auto i = GT::min; i < GT::max; ++i) ret[i] = P(exp(tmp[i] - max));
		return ret;
	} else { // on natural scale
		TGenotypeLikelihoods ret{P(1.)};
		for (size_t i = 0; i < size; ++i) {
			for (auto b1 = Base::min; b1 < Base::max; ++b1) {
				ret[genotype(b1, b1)] *= bases[i][b1];
				for (auto b2 = coretools::next(b1); b2 < Base::max; ++b2) {
					ret[genotype(b1, b2)] *= coretools::average(bases[i][b1], bases[i][b2]);
				}
			}
		}
		return ret;
	}
}

template<template<typename...> typename Container, typename... Args>
TGenotypeLikelihoods getGLH(const Container<TBaseLikelihoods, Args...> &bases) {
	return getGLH(bases, bases.size());
}

inline TGenotypeProbabilities posterior(const TGenotypeLikelihoods &likelihoods, const TGenotypeProbabilities &prior) {
	return TGenotypeProbabilities::normalize(likelihoods, prior, std::multiplies<>());
}

//template<typename Index, size_t N=Index::max, typename Likelihood=Probability>
//TStrongArraySum1<Probability, Index, N> posterior(TStrongArray<Likelihood, Index, N> Likelihoods, TStrongArraySum1<Probability, Index, N> prior)

inline TBaseLikelihoods fromError(genometools::Base trueBase, coretools::Probability error) {
	using coretools::P;
	TBaseLikelihoods ret;
	ret.fill(P(error/3.));
	ret[trueBase] = error.complement();
	return ret;
}


constexpr coretools::Probability homozygous(const TGenotypeProbabilities &ps) {
	using GT = genometools::Genotype;
	using coretools::P;
	return P(ps[GT::AA] + ps[GT::CC] + ps[GT::GG] + ps[GT::TT]);
}
constexpr coretools::Probability heterozygous(const TGenotypeProbabilities &ps) {
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
