/*
 * TGenotypePrior.h
 *
 *  Created on: Nov 20, 2018
 *      Author: phaentu
 */

#ifndef TGENOTYPEDATA_H_
#define TGENOTYPEDATA_H_

#include <algorithm>
#include <functional>
#include <numeric>
#include <stddef.h>
#include <utility>
#include <cassert>

#include "GenotypeTypes.h"
#include "probability.h"
#include "TStrongArray.h"
#include "TMassFunction.h"
#include "algorithms.h"

namespace GenotypeLikelihoods{

using TBaseProbabilities     = coretools::TStrongMassFunction<coretools::Probability, genometools::Base, 4>;
using TBaseLikelihoods       = coretools::TStrongArray<coretools::Probability, genometools::Base, 4>;
using TBaseData              = coretools::TStrongArray<double, genometools::Base, 4>;
using TBaseCounts            = coretools::TStrongArray<uint32_t, genometools::Base, 5>;
using TGenotypeLikelihoods   = coretools::TStrongArray<coretools::Probability, genometools::Genotype, 10>;
using TGenotypeProbabilities = coretools::TStrongMassFunction<coretools::Probability, genometools::Genotype, 10>;
using TGenotypeData          = coretools::TStrongArray<double, genometools::Genotype, 10>;

template<template<typename...> typename Container, typename... Args>
TGenotypeLikelihoods getGLH(const Container<TBaseLikelihoods, Args...> &bases, const size_t size) {
	using genometools::Base;
	using GT = genometools::Genotype;
	using genometools::genotype;
	// allows for vector to be longer than what is to be used
	// do in log if depth is high
	if (bases.size() > 50) {
		TGenotypeData tmp{0.};
		for (size_t i = 0; i < size; ++i) {
			for (auto b1 = Base::min; b1 < Base::max; ++b1) {
				tmp[genotype(b1, b1)] += log(bases[i][b1]);
				for (auto b2 = genometools::next(b1); b2 < Base::max; ++b2) {
					tmp[genotype(b1, b2)] += log(0.5 * ((double)bases[i][b1] + (double)bases[i][b2]));
				}
			}
		}

		// standardize and de-log
		const auto max = *std::max_element(tmp.begin(), tmp.end());
		TGenotypeLikelihoods ret;
		for (auto i = GT::min; i < GT::max; ++i) ret[i] = exp(tmp[i] - max);
		return ret;
	} else { // on natural scale
		TGenotypeLikelihoods ret{1.};
		for (size_t i = 0; i < size; ++i) {
			for (auto b1 = Base::min; b1 < Base::max; ++b1) {
				ret[genotype(b1, b1)] *= bases[i][b1];
				for (auto b2 = genometools::next(b1); b2 < Base::max; ++b2) {
					ret[genotype(b1, b2)] *= 0.5 * (bases[i][b1] + bases[i][b2]);
				}
			}
		}
		return ret;
	}
}

template<template<typename...> typename Container, typename... Args>
TGenotypeLikelihoods fillGLH(const Container<TBaseLikelihoods, Args...> &bases) {
	return getGLH(bases, bases.size());
}

inline TGenotypeProbabilities posterior(const TGenotypeLikelihoods &likelihoods, const TGenotypeProbabilities &prior) {
	return TGenotypeProbabilities(likelihoods, prior, std::multiplies<>());
}

//template<typename Index, size_t N=Index::max, typename Likelihood=Probability>
//TStrongArraySum1<Probability, Index, N> posterior(TStrongArray<Likelihood, Index, N> Likelihoods, TStrongArraySum1<Probability, Index, N> prior)

inline TBaseLikelihoods fromError(genometools::Base trueBase, coretools::Probability error) {
	TBaseLikelihoods ret;
	ret.fill(error/3.);
	ret[trueBase] = error.complement();
	return ret;
}


constexpr coretools::Probability homozygous(const TGenotypeProbabilities &ps) {
	using GT = genometools::Genotype;
	return ps[GT::AA] + ps[GT::CC] + ps[GT::GG] + ps[GT::TT];
}
constexpr coretools::Probability heterozygous(const TGenotypeProbabilities &ps) {
	return homozygous(ps).complement();
}

template<typename Container>
constexpr auto frequencies(const Container &vs) {
	const auto tot = std::accumulate(vs.begin(), vs.end(), typename Container::value_type{});
	coretools::TStrongArray<coretools::Probability, typename Container::index_type, Container::capacity> ret;
	std::transform(vs.cbegin(), vs.cend(), ret.begin(), [tot](auto v){return v / tot;});
	return ret;
}

} // namespace GenotypeLikelihoods

#endif /* TGENOTYPEDATA_H_ */
