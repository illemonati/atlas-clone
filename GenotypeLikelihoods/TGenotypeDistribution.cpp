/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"

namespace GenotypeLikelihoods{

//-------------------------------------------
// TGenotypeDistribution
// Base class.
//-------------------------------------------
void TGenotypeDistribution::fillBaseFrequences(TBaseProbabilities & baseFreq, const genometools::Genotype genotype){
	if(genotype == genometools::AA){
		baseFreq[genometools::A] = 1.0;
		baseFreq[genometools::C] = 0.0;
		baseFreq[genometools::G] = 0.0;
		baseFreq[genometools::T] = 0.0;
	} else if(genotype == genometools::AC){
		baseFreq[genometools::A] = 0.5;
		baseFreq[genometools::C] = 0.5;
		baseFreq[genometools::G] = 0.0;
		baseFreq[genometools::T] = 0.0;
	} else if(genotype == genometools::AG){
		baseFreq[genometools::A] = 0.5;
		baseFreq[genometools::C] = 0.0;
		baseFreq[genometools::G] = 0.5;
		baseFreq[genometools::T] = 0.0;
	} else if(genotype == genometools::AT){
		baseFreq[genometools::A] = 0.5;
		baseFreq[genometools::C] = 0.0;
		baseFreq[genometools::G] = 0.0;
		baseFreq[genometools::T] = 0.5;
	} else if(genotype == genometools::CC){
		baseFreq[genometools::A] = 0.0;
		baseFreq[genometools::C] = 1.0;
		baseFreq[genometools::G] = 0.0;
		baseFreq[genometools::T] = 0.0;
	} else if(genotype == genometools::CG){
		baseFreq[genometools::A] = 0.0;
		baseFreq[genometools::C] = 0.5;
		baseFreq[genometools::G] = 0.5;
		baseFreq[genometools::T] = 0.0;
	} else if(genotype == genometools::CT){
		baseFreq[genometools::A] = 0.0;
		baseFreq[genometools::C] = 0.5;
		baseFreq[genometools::G] = 0.0;
		baseFreq[genometools::T] = 0.5;
	} else if(genotype == genometools::GG){
		baseFreq[genometools::A] = 0.0;
		baseFreq[genometools::C] = 0.0;
		baseFreq[genometools::G] = 1.0;
		baseFreq[genometools::T] = 0.0;
	} else if(genotype == genometools::GT){
		baseFreq[genometools::A] = 0.0;
		baseFreq[genometools::C] = 0.0;
		baseFreq[genometools::G] = 0.5;
		baseFreq[genometools::T] = 0.5;
	} else if(genotype == genometools::TT){
		baseFreq[genometools::A] = 0.0;
		baseFreq[genometools::C] = 0.0;
		baseFreq[genometools::G] = 0.0;
		baseFreq[genometools::T] = 1.0;
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
	_baseFrequencies [genometools::A] = 0.25;
	_baseFrequencies [genometools::C] = 0.25;
	_baseFrequencies [genometools::G] = 0.25;
	_baseFrequencies [genometools::T] = 0.25;

	_genotypeFrequencies[genometools::AA] = 0.25;
	_genotypeFrequencies[genometools::AC] = 0.0;
	_genotypeFrequencies[genometools::AG] = 0.0;
	_genotypeFrequencies[genometools::AT] = 0.0;
	_genotypeFrequencies[genometools::CC] = 0.25;
	_genotypeFrequencies[genometools::CG] = 0.0;
	_genotypeFrequencies[genometools::CT] = 0.0;
	_genotypeFrequencies[genometools::GG] = 0.25;
	_genotypeFrequencies[genometools::GT] = 0.0;
	_genotypeFrequencies[genometools::TT] = 0.25;
};


}; //end namespace

