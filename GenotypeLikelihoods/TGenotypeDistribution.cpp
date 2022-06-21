/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"
#include "GenotypeTypes.h"
#include "probability.h"

namespace GenotypeLikelihoods{

//-------------------------------------------
// TGenotypeDistribution
// Base class.
//-------------------------------------------
TBaseProbabilities getBaseFrequences(const genometools::Genotype genotype){
	using genometools::Base;
	using genometools::Genotype;
	switch(genotype) {
	case Genotype::AA: return TBaseProbabilities({1., 0., 0., 0.});
	case Genotype::AC: return TBaseProbabilities({.5, .5, 0., 0.});
	case Genotype::AG: return TBaseProbabilities({.5, 0., .5, 0.});
	case Genotype::AT: return TBaseProbabilities({.5, 0., 0., .5});

	case Genotype::CC: return TBaseProbabilities({0., 1., 0., 0.});
	case Genotype::CG: return TBaseProbabilities({0., .5, .5, 0.});
	case Genotype::CT: return TBaseProbabilities({0., .5, 0., .5});

	case Genotype::GG: return TBaseProbabilities({0., 0., 1., 0.});
	case Genotype::GT: return TBaseProbabilities({0., 0., .5, .5});

	case Genotype::TT: return TBaseProbabilities({0., 0., 0., 1.});
	default: return TBaseProbabilities({0.25, 0.25, 0.25, 0.25}); // NN
	}
}

//-------------------------------------------
// TGenotypeDistribution_haploid
// Only haploid genotypes (A,C,G,T) with different frequencies
//-------------------------------------------
TGenotypeDistribution_haploid::TGenotypeDistribution_haploid(){
	_genotypeFrequencies = TGenotypeProbabilities({0.25, 0., 0., 0., 0.25, 0., 0., 0.25, 0., 0.25});
};


}; //end namespace

