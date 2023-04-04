/*
 * TgenotypeLikelihoods.h
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_
#define GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_

#include <stddef.h>
#include <algorithm>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "PMD/TModels.h"
#include "SequencingError/TModels.h"
#include "TSite.h"
#include "coretools/Types/probability.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods{

class TGenotypeLikelihoodCalculator{
protected:
	PMD::TModels _pmdModels;
	SequencingError::TModels _sequencingErrorModels;

public:
	TGenotypeLikelihoodCalculator(const BAM::TReadGroups* ReadGroups);

	const SequencingError::TModels& sequencingErrorModels() const { return _sequencingErrorModels; };
	SequencingError::TModels& sequencingErrorModels() { return _sequencingErrorModels; };
	const PMD::TModels& postMortemDamageModels() const { return _pmdModels; };
	PMD::TModels& postMortemDamageModels() { return _pmdModels; };

	bool hasPMD() const;
	bool recalibrationChangesQualities() const;

	coretools::Probability errorRate(const BAM::TSequencedBase & base) const;
	coretools::Probability errorWithPMD(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability phredInt(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability phredIntWithPMD(const BAM::TSequencedBase & base) const;
	void recalibrate(BAM::TSequencedBase & base) const;
	void recalibrateWithPMD(BAM::TSequencedBase & base) const;

	//are vector versions used??
	void recalibrate(std::vector<BAM::TSequencedBase> & bases) const;
	void recalibrateWithPMD(std::vector<BAM::TSequencedBase> & bases) const;
	double calculateLogPMDS(const BAM::TSequencedBase & base, const genometools::Base & ref, const coretools::Probability & pi) const; //TODO: move to PMDS class?

	//functions performed on sites
	TGenotypeLikelihoods calculateGenotypeLikelihoods(const TSite & site) const;
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
