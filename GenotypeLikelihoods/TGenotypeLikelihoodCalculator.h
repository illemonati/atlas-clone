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
	bool _initialized;
	TLog* _logfile;
	BAM::TReadGroups* _readGroups;
	BAM::TReadGroupMap* _readGroupMap; //TODO: find way to only initialize in sequencing error models

	TGenotypeDistribution _genotypeDistribution;
	TPostMortemDamage _pmd;
	TSequencingErrorModels _sequencingErrorModels; //TODO: find a way not to use a pointer

	//temporary storage: all are mutable
	mutable TGenotypeLikelihoods _genotypeLikelihoods;
	mutable std::vector<TBaseData> _baseLikelihoods;
	mutable TBaseData _baseLikelihoodsNoPMD;
	mutable TBaseData _tmpBaseData;

public:
	TGenotypeLikelihoodCalculator();
	TGenotypeLikelihoodCalculator(TParameters & params, BAM::TReadGroups* ReadGroups, TLog* Logfile);
	~TGenotypeLikelihoodCalculator();

	void init(TParameters & params, BAM::TReadGroups* ReadGroups, TLog* Logfile);
	const TSequencingErrorModels& getSequencingErrorModels() const{ return _sequencingErrorModels; };

	bool hasPMD() const;
	bool recalibrationChangesQualities() const;

	double getErrorRate(const BAM::TBase & base) const;
	double getErrorWithPMD(const BAM::TBase & base) const;
	uint8_t getPhredInt(const BAM::TBase & base) const;
	uint8_t getPhredIntWithPMD(const BAM::TBase & base) const;
	void recalibrate(BAM::TBase & base) const;
	void recalibrateWithPMD(BAM::TBase & base) const;
	void recalibrate(std::vector<BAM::TBase> & bases) const;
	void recalibrateWithPMD(std::vector<BAM::TBase> & bases) const;

	void calculateGenotypeLikelihoods(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) const;
	double calculateLogPMDS(const BAM::TBase & base, const Base ref, const double pi) const; //TODO: move to PMDS class?
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
