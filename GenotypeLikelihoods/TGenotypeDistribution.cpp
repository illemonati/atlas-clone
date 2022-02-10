/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"
#include "GenotypeTypes.h"

namespace GenotypeLikelihoods{

//-------------------------------------------
// TGenotypeDistribution
// Base class.
//-------------------------------------------
void TGenotypeDistribution::fillBaseFrequences(TBaseProbabilities & baseFreq, const genometools::Genotype genotype){
	using genometools::Base;
	using genometools::Genotype;
	if(genotype == Genotype::AA){
		baseFreq[Base::A] = 1.0;
		baseFreq[Base::C] = 0.0;
		baseFreq[Base::G] = 0.0;
		baseFreq[Base::T] = 0.0;
	} else if(genotype == Genotype::AC){
		baseFreq[Base::A] = 0.5;
		baseFreq[Base::C] = 0.5;
		baseFreq[Base::G] = 0.0;
		baseFreq[Base::T] = 0.0;
	} else if(genotype == Genotype::AG){
		baseFreq[Base::A] = 0.5;
		baseFreq[Base::C] = 0.0;
		baseFreq[Base::G] = 0.5;
		baseFreq[Base::T] = 0.0;
	} else if(genotype == Genotype::AT){
		baseFreq[Base::A] = 0.5;
		baseFreq[Base::C] = 0.0;
		baseFreq[Base::G] = 0.0;
		baseFreq[Base::T] = 0.5;
	} else if(genotype == Genotype::CC){
		baseFreq[Base::A] = 0.0;
		baseFreq[Base::C] = 1.0;
		baseFreq[Base::G] = 0.0;
		baseFreq[Base::T] = 0.0;
	} else if(genotype == Genotype::CG){
		baseFreq[Base::A] = 0.0;
		baseFreq[Base::C] = 0.5;
		baseFreq[Base::G] = 0.5;
		baseFreq[Base::T] = 0.0;
	} else if(genotype == Genotype::CT){
		baseFreq[Base::A] = 0.0;
		baseFreq[Base::C] = 0.5;
		baseFreq[Base::G] = 0.0;
		baseFreq[Base::T] = 0.5;
	} else if(genotype == Genotype::GG){
		baseFreq[Base::A] = 0.0;
		baseFreq[Base::C] = 0.0;
		baseFreq[Base::G] = 1.0;
		baseFreq[Base::T] = 0.0;
	} else if(genotype == Genotype::GT){
		baseFreq[Base::A] = 0.0;
		baseFreq[Base::C] = 0.0;
		baseFreq[Base::G] = 0.5;
		baseFreq[Base::T] = 0.5;
	} else if(genotype == Genotype::TT){
		baseFreq[Base::A] = 0.0;
		baseFreq[Base::C] = 0.0;
		baseFreq[Base::G] = 0.0;
		baseFreq[Base::T] = 1.0;
	}
};

//-------------------------------------------
// TGenotypeDistribution_haploid
// Only haploid genotypes (A,C,G,T) with different frequencies
//-------------------------------------------
TGenotypeDistribution_haploid::TGenotypeDistribution_haploid(){
	reset();
};

void TGenotypeDistribution_haploid::reset(){
	using genometools::Base;
	using genometools::Genotype;
	_baseFrequencies [Base::A] = 0.25;
	_baseFrequencies [Base::C] = 0.25;
	_baseFrequencies [Base::G] = 0.25;
	_baseFrequencies [Base::T] = 0.25;

	_genotypeFrequencies[Genotype::AA] = 0.25;
	_genotypeFrequencies[Genotype::AC] = 0.0;
	_genotypeFrequencies[Genotype::AG] = 0.0;
	_genotypeFrequencies[Genotype::AT] = 0.0;
	_genotypeFrequencies[Genotype::CC] = 0.25;
	_genotypeFrequencies[Genotype::CG] = 0.0;
	_genotypeFrequencies[Genotype::CT] = 0.0;
	_genotypeFrequencies[Genotype::GG] = 0.25;
	_genotypeFrequencies[Genotype::GT] = 0.0;
	_genotypeFrequencies[Genotype::TT] = 0.25;
};


}; //end namespace

