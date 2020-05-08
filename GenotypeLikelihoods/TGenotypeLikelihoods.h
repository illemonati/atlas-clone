/*
 * TgenotypeLikelihoods.h
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODS_H_
#define GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODS_H_

#include "TPostMortemDamage.h"
#include "TParameters.h"
#include "TSequencingErrorModel.h"

namespace GenotypeLikelihoods{

class TGenotypeLikelihoods{
private:
	TLog* logfile;
	TReadGroups* readGroups;
	TReadGroupMap* readGroupMap; //TODO: find way to only initialize in sequencing error models

	TPostMortemDamage pmd;
	TSequencingErrorModels* errorModels; //TODO: find a way not to use a pointer

public:

	TGenotypeLikelihoods(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile);

	double getErrorRate(const TBaseData & base);
	void fillGenotypeLikelihoods(const std::vector<TBase> bases, double* genotypeLikelihoods);


};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODS_H_ */
