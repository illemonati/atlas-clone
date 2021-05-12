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
void TGenotypeDistribution::fillBaseFrequences(TBaseData & baseFreq, const BAM::Genotype genotype){
	if(genotype == BAM::AA){
		baseFreq[BAM::A] = 1.0;
		baseFreq[BAM::C] = 0.0;
		baseFreq[BAM::G] = 0.0;
		baseFreq[BAM::T] = 0.0;
	} else if(genotype == BAM::AC){
		baseFreq[BAM::A] = 0.5;
		baseFreq[BAM::C] = 0.5;
		baseFreq[BAM::G] = 0.0;
		baseFreq[BAM::T] = 0.0;
	} else if(genotype == BAM::AG){
		baseFreq[BAM::A] = 0.5;
		baseFreq[BAM::C] = 0.0;
		baseFreq[BAM::G] = 0.5;
		baseFreq[BAM::T] = 0.0;
	} else if(genotype == BAM::AT){
		baseFreq[BAM::A] = 0.5;
		baseFreq[BAM::C] = 0.0;
		baseFreq[BAM::G] = 0.0;
		baseFreq[BAM::T] = 0.5;
	} else if(genotype == BAM::CC){
		baseFreq[BAM::A] = 0.0;
		baseFreq[BAM::C] = 1.0;
		baseFreq[BAM::G] = 0.0;
		baseFreq[BAM::T] = 0.0;
	} else if(genotype == BAM::CG){
		baseFreq[BAM::A] = 0.0;
		baseFreq[BAM::C] = 0.5;
		baseFreq[BAM::G] = 0.5;
		baseFreq[BAM::T] = 0.0;
	} else if(genotype == BAM::CT){
		baseFreq[BAM::A] = 0.0;
		baseFreq[BAM::C] = 0.5;
		baseFreq[BAM::G] = 0.0;
		baseFreq[BAM::T] = 0.5;
	} else if(genotype == BAM::GG){
		baseFreq[BAM::A] = 0.0;
		baseFreq[BAM::C] = 0.0;
		baseFreq[BAM::G] = 1.0;
		baseFreq[BAM::T] = 0.0;
	} else if(genotype == BAM::GT){
		baseFreq[BAM::A] = 0.0;
		baseFreq[BAM::C] = 0.0;
		baseFreq[BAM::G] = 0.5;
		baseFreq[BAM::T] = 0.5;
	} else if(genotype == BAM::TT){
		baseFreq[BAM::A] = 0.0;
		baseFreq[BAM::C] = 0.0;
		baseFreq[BAM::G] = 0.0;
		baseFreq[BAM::T] = 1.0;
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
	_baseFrequencies [BAM::A] = 0.25;
	_baseFrequencies [BAM::C] = 0.25;
	_baseFrequencies [BAM::G] = 0.25;
	_baseFrequencies [BAM::T] = 0.25;

	_genotypeFrequencies[BAM::AA] = 0.25;
	_genotypeFrequencies[BAM::AC] = 0.0;
	_genotypeFrequencies[BAM::AG] = 0.0;
	_genotypeFrequencies[BAM::AT] = 0.0;
	_genotypeFrequencies[BAM::CC] = 0.25;
	_genotypeFrequencies[BAM::CG] = 0.0;
	_genotypeFrequencies[BAM::CT] = 0.0;
	_genotypeFrequencies[BAM::GG] = 0.25;
	_genotypeFrequencies[BAM::GT] = 0.0;
	_genotypeFrequencies[BAM::TT] = 0.25;
};


}; //end namespace

