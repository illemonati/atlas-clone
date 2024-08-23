/*
 * TgenotypeLikelihoods.h
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_
#define GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_

#include "coretools/Types/probability.h"

#include "PMD/TModels.h"
#include "SequencingError/TModels.h"
#include "genometools/GenotypeContainers.h"
#include "TReadGroupInfo.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }
namespace BAM { class TAlignment; }

namespace GenotypeLikelihoods{
class TSite;

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
	coretools::PhredInt phredIntWithPMD(const BAM::TSequencedBase & base) const;
	void recalibrateWithPMD(BAM::TSequencedBase & base) const;

	void recalibrateWithPMD(BAM::TAlignment& aln) const;
	double calculateLogPMDS(const BAM::TSequencedBase & base, const genometools::Base & ref, const coretools::Probability & pi) const; //TODO: move to PMDS class?

	//functions performed on sites
	genometools::TGenotypeLikelihoods calculateGenotypeLikelihoods(const TSite & site) const;
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
