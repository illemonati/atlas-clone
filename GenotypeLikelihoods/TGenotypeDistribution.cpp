/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"
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

double THaploidDistribution::weightedSum(const TGenotypeLikelihoods &likelihoods) const {
	return likelihoods[Genotype::AA] * _weights[Genotype::AA] + likelihoods[Genotype::CC] * _weights[Genotype::CC] +
		   likelihoods[Genotype::GG] * _weights[Genotype::GG] + likelihoods[Genotype::TT] * _weights[Genotype::TT];
}

double THaploidDistribution::normalize(TGenotypeLikelihoods &likelihoods) const {
	double sum = 0.;
	// only four
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g = genometools::genotype(b, b);
		likelihoods[g] *= _weights[g];
		sum += likelihoods[g];
	}
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g = genometools::genotype(b, b);
		likelihoods[g].scale(sum);
	}
	return sum;
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
	// if first == second, then 0.5*first + 0.5*first = first -> no if needed
	return 0.5 * (baseLikelihoods[genometools::first(genotype)] + baseLikelihoods[genometools::second(genotype)]);
}

double TDiploidDistribution::weightedSum(const TGenotypeLikelihoods &likelihoods) const {
	return coretools::weightedSum(likelihoods, _weights);
}

double TDiploidDistribution::normalize(TGenotypeLikelihoods &likelihoods) const {
	double sum = 0;
	// all 10
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		likelihoods[g] *= _weights[g];
		sum += likelihoods[g];
	}
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		likelihoods[g].scale(sum);
	}
	return sum;
}

}; // namespace GenotypeLikelihoods
