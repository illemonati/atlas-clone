/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"
#include "coretools/Strings/toString.h"
#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "coretools/algorithms.h"
#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods {
using genometools::Base;
using genometools::Genotype;

TGenotypeLikelihoods THaploidDistribution::getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const {
	return TGenotypeLikelihoods({baseLikelihoods[Base::A], 0., 0., 0.,
			                     baseLikelihoods[Base::C], 0., 0.,
								 baseLikelihoods[Base::G], 0.,
								 baseLikelihoods[Base::T]});
}
coretools::Probability THaploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// first == second if Haploid
	return baseLikelihoods[genometools::first(genotype)];
}

double THaploidDistribution::normalize_add(TGenotypeLikelihoods &likelihoods) {
	double sum = 0.;
	// only four
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g = genometools::genotype(b, b);
		likelihoods[g] *= _pi[b];
		sum += likelihoods[g];
	}
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g = genometools::genotype(b, b);
		_piSum[b] += likelihoods[g].scale(sum);
	}
	return sum;
}

void THaploidDistribution::estimate() {
	_pi = TBaseProbabilities::normalize(_piSum);
	_piSum.fill(0.);
}

std::string THaploidDistribution::definition() const noexcept {
	using coretools::str::toString;
	return std::string("AA: ").append(toString(_pi[Base::A]))
		.append(", CC: ").append(toString(_pi[Base::C]))
		.append(", GG: ").append(toString(_pi[Base::G]))
		.append(", TT: ").append(toString(_pi[Base::T]));
}

TGenotypeLikelihoods TDiploidDistribution::getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const {
	return TGenotypeLikelihoods({baseLikelihoods[Base::A], 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::C]),
								 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::G]),
								 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::T]), baseLikelihoods[Base::C],
								 0.5 * (baseLikelihoods[Base::C] + baseLikelihoods[Base::G]),
								 0.5 * (baseLikelihoods[Base::C] + baseLikelihoods[Base::T]), baseLikelihoods[Base::G],
								 0.5 * (baseLikelihoods[Base::G] + baseLikelihoods[Base::T]),
								 baseLikelihoods[Base::T]});
}

coretools::Probability TDiploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// if first == second, then 0.5*first + 0.5*first = first
	return 0.5 * (baseLikelihoods[genometools::first(genotype)] + baseLikelihoods[genometools::second(genotype)]);
}

double TDiploidDistribution::normalize_add(TGenotypeLikelihoods &likelihoods) {
	double sum = 0;
	// all 10
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		likelihoods[g] *= _pi[g];
		sum += likelihoods[g];
	}
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		_piSum[g] += likelihoods[g].scale(sum);
	}
	return sum;
}

void TDiploidDistribution::estimate() {
	// Only update if values make sense, else set to initial value
	const auto sumHomo = _piSum[Genotype::AA] + _piSum[Genotype::CC] + _piSum[Genotype::GG] + _piSum[Genotype::TT];
	const auto sumHet  = _piSum[Genotype::AC] + _piSum[Genotype::AG] + _piSum[Genotype::AT] + _piSum[Genotype::CG] +
						 _piSum[Genotype::CT] + _piSum[Genotype::GT];

	if (sumHomo > sumHet) {
		// estimate looks good enough
		_pi = TGenotypeProbabilities::normalize(_piSum);
	} else {
		// the estimate is definitely wrong. Let's use the default
		_pi = _pi_init;
	}

	_piSum.fill(0.);
}

std::string TDiploidDistribution::definition() const noexcept {
	using coretools::str::toString;
	std::string ret;
	for (auto g = Genotype::min; g < Genotype::max; ++g) {
		ret.append(toString(g)).append(": ").append(toString(_pi[g])).append(", ");
	}
	ret.resize(ret.size() - 2);
	return ret;
}

}; // namespace GenotypeLikelihoods
