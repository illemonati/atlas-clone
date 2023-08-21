/*
 * TgenotypeLikelihoods.h
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_
#define GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_

#include "oldPMD/TModels.h"
#include "SequencingError/TModels.h"
#include "TAlignment.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TSite.h"
#include "coretools/Types/probability.h"
#include "genometools/PhredProbabilityTypes.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods{

class TGenotypeLikelihoodCalculator{
private:
	oldPMD::TModels _pmdModels;
	SequencingError::TModels _sequencingErrorModels;

public:
	TGenotypeLikelihoodCalculator(const BAM::TReadGroups* ReadGroups);

	const SequencingError::TModels& sequencingErrorModels() const { return _sequencingErrorModels; };
	SequencingError::TModels& sequencingErrorModels() { return _sequencingErrorModels; };
	const oldPMD::TModels& postMortemDamageModels() const { return _pmdModels; };
	oldPMD::TModels& postMortemDamageModels() { return _pmdModels; };

	coretools::Probability errorWithPMD(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability phredIntWithPMD(const BAM::TSequencedBase & base) const;
	void recalibrateWithPMD(BAM::TSequencedBase & base) const;

	void recalibrateWithPMD(BAM::TAlignment& aln) const;
	double calculateLogPMDS(const BAM::TSequencedBase & base, const genometools::Base & ref, const coretools::Probability & pi) const; //TODO: move to PMDS class?

	//functions performed on sites
	TGenotypeLikelihoods calculateGenotypeLikelihoods(const TSite & site) const;
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
