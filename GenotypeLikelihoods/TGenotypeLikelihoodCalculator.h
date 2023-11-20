/*
 * TgenotypeLikelihoods.h
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_
#define GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_

#include "PMD/TModels.h"
#include "TReadGroupInfo.h"
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

class TErrorModels{
private:
	PMD::TModels _pmd;
	SequencingError::TModels _recal;

public:
	TErrorModels(const BAM::RGInfo::TReadGroupInfo& RGInfo);

	const SequencingError::TModels& sequencingErrorModels() const { return _recal; };
	SequencingError::TModels& sequencingErrorModels() { return _recal; };
	const PMD::TModels& postMortemDamageModels() const { return _pmd; };
	PMD::TModels& postMortemDamageModels() { return _pmd; };

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
