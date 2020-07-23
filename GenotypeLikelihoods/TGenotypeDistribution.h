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

	TBaseData baseFrequencies;

public:
	TGenotypeDistribution(){};
	virtual ~TGenotypeDistribution(){};

	virtual const TBaseData& baseLikelihoods(){
		throw std::runtime_error("void TGenotypeDistribution::const TBaseData& baseLikelihoods(): not implemented for base class!");
	};

};

//-------------------------------------------
// TGenotypeDistribution_haploid
// Only haploid genotypes (A,C,G,T) with different frequencies
//-------------------------------------------

class TGenotypeDistribution_haploid:public TGenotypeDistribution{
private:

public:
	TGenotypeDistribution_haploid();

	const TBaseData& baseLikelihoods(){ return baseFrequencies; };


};

}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
