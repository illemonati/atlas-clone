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
	bool _initialized;
	TLog* _logfile;
	BAM::TReadGroups* _readGroups;
	BAM::TReadGroupMap* _readGroupMap; //TODO: find way to only initialize in sequencing error models

	TGenotypeDistribution _genotypeDistribution;
	TPostMortemDamage _pmd;
	TSequencingErrorModels _sequencingErrorModels; //TODO: find a way not to use a pointer

	//temporary storage
	TGenotypeLikelihoods _genotypeLikelihoods;
	std::vector<TBaseData> _baseLikelihoods;
	TBaseData _baseLikelihoodsNoPMD;
	TBaseData _tmpBaseData;

public:
	TGenotypeLikelihoodCalculator();
	TGenotypeLikelihoodCalculator(TParameters & params, BAM::TReadGroups* ReadGroups, TLog* Logfile);
	~TGenotypeLikelihoodCalculator();

	void init(TParameters & params, BAM::TReadGroups* ReadGroups, TLog* Logfile);
	const TSequencingErrorModels& getSequencingErrorModels() const{ return _sequencingErrorModels; };

	bool hasPMD() const;
	bool recalibrationChangesQualities() const;

	double getErrorRate(const TBase & base) const;
	double getErrorWithPMD(const TBase & base) const;
	uint8_t getPhredInt(const TBase & base) const;
	uint8_t getPhredIntWithPMD(const TBase & base) const;
	void recalibrate(TBase & base) const;
	void recalibrateWithPMD(TBase & base) const;
	void recalibrate(std::vector<TBase> & bases) const;
	void recalibrateWithPMD(std::vector<TBase> & bases) const;

	void calculateGenotypeLikelihoods(const std::vector<TBase*> bases, TGenotypeLikelihoods & genotypeLikelihoods) const;
	double calculatePMDS(const TBase & base, const Base ref, const double pi) const;
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
