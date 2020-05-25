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
#include "TAlignment.h"

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
	std::vector<TBaseData> baseLikelihoods;
	TBaseData baseLikelihoodsNoPMD;
	TBaseData tmpBaseData;

public:
	TGenotypeLikelihoodCalculator();
	TGenotypeLikelihoodCalculator(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile);
	~TGenotypeLikelihoodCalculator();

	void init(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile);
	TSequencingErrorModels& getSequencingErrorModels(){ return sequencingErrorModels; };

	double getErrorRate(const TBase & base);
	double getErrorWithPMD(const TBase & base);
	uint8_t getPhredInt(const TBase & base);
	uint8_t getPhredIntWithPMD(const TBase & base);
	void recalibrate(TBase & base);
	void recalibrateWithPMD(TBase & base);
	void recalibrate(std::vector<TBase> & bases);
	void recalibrateWithPMD(std::vector<TBase> & bases);

	void calculateGenotypeLikelihoods(const std::vector<TBase*> bases, TGenotypeLikelihoods & genotypeLikelihoods);
	double calculatePMDS(const TBase & base, const Base ref, const double pi);
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
