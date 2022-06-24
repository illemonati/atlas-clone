/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"
#include "GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "algorithms.h"
#include "probability.h"

namespace GenotypeLikelihoods {
using coretools::Probability;
using genometools::Base;
using genometools::Genotype;

TGenotypeLikelihoods THaploidDistribution::getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const {
	return TGenotypeLikelihoods({baseLikelihoods[Base::A], 0., 0., 0., baseLikelihoods[Base::C], 0., 0.,
								 baseLikelihoods[Base::G], 0., baseLikelihoods[Base::G]});
}
coretools::Probability THaploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// first == second if Haploid
	return baseLikelihoods[genometools::first(genotype)];
}

coretools::Probability THaploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods) const {
	return baseLikelihoods[Base::A] * _weights[Genotype::AA] + baseLikelihoods[Base::C] * _weights[Genotype::CC] +
		   baseLikelihoods[Base::G] * _weights[Genotype::GG] + baseLikelihoods[Base::T] * _weights[Genotype::TT];
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

coretools::Probability TDiploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods) const {
	const auto Ls = getGenotypeLikelihoods(baseLikelihoods);
	return coretools::weightedSum(Ls, _weights);
	
}

}; // namespace GenotypeLikelihoods
