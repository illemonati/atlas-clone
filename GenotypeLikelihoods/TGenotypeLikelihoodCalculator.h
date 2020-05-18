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
#include "TSequencingErrorModels.h"

namespace GenotypeLikelihoods{


class TGenotypeLikelihoodCalculator{
private:
	bool initialized;
	TLog* logfile;
	TReadGroups* readGroups;
	TReadGroupMap* readGroupMap; //TODO: find way to only initialize in sequencing error models

	TGenotypeDistribution genotypeDistribution;
	TPostMortemDamage pmd;
	TSequencingErrorModels sequencingErrorModels; //TODO: find a way not to use a pointer

	//temporary storage
	TGenotypeLikelihoods genotypeLikelihoods;
	std::vector<TBaseLikelihoods> baseLikelihoods;
	TBaseLikelihoods baseLikelihoodsNoPMD;

public:
	TGenotypeLikelihoodCalculator();
	TGenotypeLikelihoodCalculator(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile);
	~TGenotypeLikelihoodCalculator();

	void init(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile);
	TSequencingErrorModels& getSequencingErrorModels(){ return sequencingErrorModels; };

	double getErrorRate(const TBaseData & base);
	double getErrorWithPMD(const TBaseData & base);
	uint8_t getPhredInt(const TBaseData & base);
	uint8_t getPhredIntWithPMD(const TBaseData & base);
	void recalibrate(TBaseData & base);
	void recalibrateWithPMD(TBaseData & base);
	void recalibrate(TBase* bases, const uint16_t  length);
	void recalibrateWithPMD(TBase* bases, const uint16_t  length);

	void calculateGenotypeLikelihoods(const std::vector<TBase*> bases, TGenotypeLikelihoods & genotypeLikelihoods);
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
