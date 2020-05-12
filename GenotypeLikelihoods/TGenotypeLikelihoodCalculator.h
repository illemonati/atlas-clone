/*
 * TgenotypeLikelihoods.h
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_
#define GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_

#include "TParameters.h"
#include "TGenotypeDistribution.h"
#include "TPostMortemDamage.h"
#include "TSequencingErrorModel.h"

namespace GenotypeLikelihoods{


class TGenotypeLikelihoodCalculator{
private:
	TLog* logfile;
	TReadGroups* readGroups;
	TReadGroupMap* readGroupMap; //TODO: find way to only initialize in sequencing error models

	TGenotypeDistribution genotypeDistribution;
	TPostMortemDamage pmd;
	TSequencingErrorModels* errorModels; //TODO: find a way not to use a pointer

	//temporary storage
	TGenotypeLikelihoods genotypeLikelihoods;
	std::vector<TBaseLikelihoods> baseLikelihoods;
	TBaseLikelihoods baseLikelihoodsNoPMD;

public:

	TGenotypeLikelihoodCalculator(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile);

	double getErrorRate(const TBaseData & base);
	void calculateGenotypeLikelihoods(const std::vector<TBase*> bases, TGenotypeLikelihoods & genotypeLikelihoods);
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
