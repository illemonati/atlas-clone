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
#include "TReadGroupInfo.h"
#include "genometools/Genotypes/Containers.h"

namespace BAM {
class TReadGroups;
struct TSequencedData;
class TAlignment;
} // namespace BAM

namespace GenotypeLikelihoods {
class TSite;

class TErrorModels {
private:
	PMD::TModels _pmd;
	SequencingError::TModels _recal;

public:
	TErrorModels(BAM::RGInfo::TReadGroupInfo &RGInfo);

	const SequencingError::TModels &sequencingErrorModels() const { return _recal; };
	SequencingError::TModels &sequencingErrorModels() { return _recal; };
	const PMD::TModels &postMortemDamageModels() const { return _pmd; };
	PMD::TModels &postMortemDamageModels() { return _pmd; };

	coretools::Probability errorWithPMD(const BAM::TSequencedData &data) const;
	coretools::PhredInt phredIntWithPMD(const BAM::TSequencedData &data) const;
	void recalibrateWithPMD(BAM::TSequencedData &data) const;

	void recalibrateWithPMD(BAM::TAlignment &aln) const;
	double calculateLogPMDS(const BAM::TSequencedData &data, const genometools::Base &ref,
	                        const coretools::Probability &pi) const; // TODO: move to PMDS class?

	// functions performed on sites
	genometools::TGenotypeLikelihoods calculateGenotypeLikelihoods(const TSite &site) const;
	genometools::TBaseLikelihoods calculateBaseLikelihoods(const TSite &site) const;
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
