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
void TGenotypeDistribution::fillBaseFrequences(TBaseData & baseFreq, const Genotype genotype){
	if(genotype == AA){
		baseFreq[A] = 1.0;
		baseFreq[C] = 0.0;
		baseFreq[G] = 0.0;
		baseFreq[T] = 0.0;
	} else if(genotype == AC){
		baseFreq[A] = 0.5;
		baseFreq[C] = 0.5;
		baseFreq[G] = 0.0;
		baseFreq[T] = 0.0;
	} else if(genotype == AG){
		baseFreq[A] = 0.5;
		baseFreq[C] = 0.0;
		baseFreq[G] = 0.5;
		baseFreq[T] = 0.0;
	} else if(genotype == AT){
		baseFreq[A] = 0.5;
		baseFreq[C] = 0.0;
		baseFreq[G] = 0.0;
		baseFreq[T] = 0.5;
	} else if(genotype == CC){
		baseFreq[A] = 0.0;
		baseFreq[C] = 1.0;
		baseFreq[G] = 0.0;
		baseFreq[T] = 0.0;
	} else if(genotype == CG){
		baseFreq[A] = 0.0;
		baseFreq[C] = 0.5;
		baseFreq[G] = 0.5;
		baseFreq[T] = 0.0;
	} else if(genotype == CT){
		baseFreq[A] = 0.0;
		baseFreq[C] = 0.5;
		baseFreq[G] = 0.0;
		baseFreq[T] = 0.5;
	} else if(genotype == GG){
		baseFreq[A] = 0.0;
		baseFreq[C] = 0.0;
		baseFreq[G] = 1.0;
		baseFreq[T] = 0.0;
	} else if(genotype == GT){
		baseFreq[A] = 0.0;
		baseFreq[C] = 0.0;
		baseFreq[G] = 0.5;
		baseFreq[T] = 0.5;
	} else if(genotype == TT){
		baseFreq[A] = 0.0;
		baseFreq[C] = 0.0;
		baseFreq[G] = 0.0;
		baseFreq[T] = 1.0;
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
	_baseFrequencies [A] = 0.25;
	_baseFrequencies [C] = 0.25;
	_baseFrequencies [G] = 0.25;
	_baseFrequencies [T] = 0.25;

	_genotypeFrequencies[AA] = 0.25;
	_genotypeFrequencies[AC] = 0.0;
	_genotypeFrequencies[AG] = 0.0;
	_genotypeFrequencies[AT] = 0.0;
	_genotypeFrequencies[CC] = 0.25;
	_genotypeFrequencies[CG] = 0.0;
	_genotypeFrequencies[CT] = 0.0;
	_genotypeFrequencies[GG] = 0.25;
	_genotypeFrequencies[GT] = 0.0;
	_genotypeFrequencies[TT] = 0.25;
};


}; //end namespace

