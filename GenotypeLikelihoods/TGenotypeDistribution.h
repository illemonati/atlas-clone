/*
 * TGenotypeModel.h
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_

#include "auxiliaryTools.h"

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

	const TBaseData& baseFrequencies(){ return _baseFrequencies; };
	const TGenotypeProbabilities& genotypeFrequencies(){ return _genotypeFrequencies; };


};

//-------------------------------------------
// TGenotypeDistribution_haploid
// Only haploid genotypes (A,C,G,T) with different frequencies
//-------------------------------------------

class TGenotypeDistribution_haploid:public TGenotypeDistribution{
private:

public:
	TGenotypeDistribution_haploid();




};

}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
