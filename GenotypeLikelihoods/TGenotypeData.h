/*
 * TGenotypePrior.h
 *
 *  Created on: Nov 20, 2018
 *      Author: phaentu
 */

#ifndef TGENOTYPEDATA_H_
#define TGENOTYPEDATA_H_

#include <numeric>
#include <stddef.h>

#include "GenotypeTypes.h"
#include "probability.h"
#include "TStrongArray.h"

namespace GenotypeLikelihoods{

using Likelihood = coretools::Probability;

using TBaseProbabilities     = coretools::TStrongArray<coretools::Probability, genometools::Base, 4>;
using TBaseLikelihoods       = coretools::TStrongArray<Likelihood, genometools::Base, 4>;
using TBaseData              = coretools::TStrongArray<double, genometools::Base, 4>;
using TBaseCounts            = coretools::TStrongArray<uint32_t, genometools::Base, 5>;
using TGenotypeLikelihoods   = coretools::TStrongArray<Likelihood, genometools::Genotype, 10>;
using TGenotypeProbabilities = coretools::TStrongArray<coretools::Probability, genometools::Genotype, 10>;
using TGenotypeData          = coretools::TStrongArray<double, genometools::Genotype, 10>;

template<template<typename... Args> typename Container, typename Type, typename... Args>
void reset(Container<Type, Args...> & c) {
	c.fill(Type{});
}

inline TGenotypeLikelihoods fillGLH(const std::vector<TBaseLikelihoods> &bases, const size_t size) {
	using genometools::Base;
	using GT = genometools::Genotype;
	using genometools::index;
	// allows for vector to be longer than what is to be used
	// do in log if depth is high
	if (bases.size() > 50) {
		// initialize tmp to zero
		TGenotypeData tmp{};
		// add to log genotype data
		for (size_t i = 0; i < size; ++i) {
			for (auto b1 = Base::min; b1 < Base::max; ++b1) {
				for (auto b2 = b1; b2 < Base::max; ++b2) {
					const auto gt = genometools::genotype(b1, b2);
					tmp[gt]      += log(0.5 * (bases[i][b1] + bases[i][b2]));
				}
			}
		}

		// standardize and de-log
		double max           = *std::max_element(tmp.begin(), tmp.end());
		TGenotypeLikelihoods ret;
		for (auto i = GT::min; i < GT::max; ++i) ret[i] = exp(tmp[i] - max);
		return ret;
	} else { // on natural scale
		// initialize tmp to 1.0
		TGenotypeLikelihoods ret{};
		for (size_t i = 0; i < size; ++i) {
			for (auto b1 = Base::min; b1 < Base::max; ++b1) {
				for (auto b2 = b1; b2 < Base::max; ++b2) {
					const auto gt = genometools::genotype(b1, b2);
					ret[gt]      *= 0.5 * (bases[i][b1] + bases[i][b2]);
				}
			}
		}
		return ret;
	}
}

template<typename Container>
void normalize(Container& c) {
	using T = typename Container::value_type;
	const auto tot = std::accumulate(c.begin(), c.end(), T{});
	for (auto & v: c) v /= tot;
};

inline TGenotypeProbabilities posterior(const TGenotypeLikelihoods &likelihoods, const TGenotypeProbabilities &prior) {
	using GT = genometools::Genotype;
	TGenotypeProbabilities ret;
	for (auto gt = GT::min; gt < GT::max; ++gt) ret[gt] = likelihoods[gt] * prior[gt];
	const auto tot = std::accumulate(ret.begin(), ret.end(), coretools::Probability{});
	for (auto & v: ret) v /= tot;
	return ret;
}

inline TBaseLikelihoods fromError(genometools::Base trueBase, coretools::Probability error) {
	TBaseLikelihoods ret;
	ret.fill(error/3.);
	ret[trueBase] = error.complement();
	return ret;
}

} // namespace GenotypeLikelihoods

#endif /* TGENOTYPEDATA_H_ */
