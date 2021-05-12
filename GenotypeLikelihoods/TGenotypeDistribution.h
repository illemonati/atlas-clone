/*
 * TGenotypeModel.h
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_

#include "TGenotypeData.h"

namespace GenotypeLikelihoods{

//-------------------------------------------
// TGenotypeDistribution
// Base class.
//-------------------------------------------

class TGenotypeDistribution{
protected:

	TBaseData _baseFrequencies; //reflects expected base frequencies under the model
	TGenotypeProbabilities _genotypeFrequencies;

public:
	TGenotypeDistribution(){};
	virtual ~TGenotypeDistribution(){};

	virtual void reset(){};

	const TBaseData& baseFrequencies(){ return _baseFrequencies; };
	const TGenotypeProbabilities& genotypeFrequencies(){ return _genotypeFrequencies; };

	void fillBaseFrequences(TBaseData & baseFreq, const BAM::Genotype genotype);
};

//-------------------------------------------
// TGenotypeDistribution_haploid
// Only haploid genotypes (A,C,G,T) with different frequencies
//-------------------------------------------

class TGenotypeDistribution_haploid:public TGenotypeDistribution{
private:

public:
	TGenotypeDistribution_haploid();

	void reset() override;
};

}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
