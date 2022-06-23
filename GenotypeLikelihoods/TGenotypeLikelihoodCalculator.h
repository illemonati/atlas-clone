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

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TSite.h"
#include "probability.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods{

	/*
class TGenotypeLikelihoodCalculator_simple{
private:
public:

	template <typename PMD, typename SEQERR>
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase & base, const PMD & PmdModels, const SEQERR & SequencingErrorModels) const{
		TBaseLikelihoods baseLikelihoodsNoPMD = SequencingErrorModels.getBaseLikelihoods(base);
		return PmdModels.getBaseLikelihoods(base, baseLikelihoodsNoPMD);
	};

	template<typename PMD, typename SEQERR>
	TGenotypeLikelihoods getGenotypeLikelihoods(const TSite &site, const PMD &PmdModels,
												 const SEQERR &SequencingErrorModels) const {
		if (site.empty()) { return TGenotypeLikelihoods{1.}; }
		std::vector<TBaseLikelihoods> baseLikelihoods;
		baseLikelihoods.reserve(site.depth());
		// calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
		for (const auto &d : site) {
			baseLikelihoods.push_back(getBaseLikelihoods(d, PmdModels, SequencingErrorModels));
		}

		// calculate genotype likelihoods
		return getGLH(baseLikelihoods, site.depth());
	};
};
	*/

class TGenotypeLikelihoodCalculator{
protected:
	bool _initialized;
	//TGenotypeLikelihoodCalculator_simple _calculator;

	//TGenotypeDistribution _genotypeDistribution;
	TPostMortemDamage _pmdModels;
	SequencingError::TModels _sequencingErrorModels;

public:
	TGenotypeLikelihoodCalculator();
	TGenotypeLikelihoodCalculator(const BAM::TReadGroups* ReadGroups);
	~TGenotypeLikelihoodCalculator() = default;

	void init(const BAM::TReadGroups* ReadGroups);
	const SequencingError::TModels& getSequencingErrorModels() const { return _sequencingErrorModels; };
	SequencingError::TModels& getSequencingErrorModelsMutable() { return _sequencingErrorModels; };
	const TPostMortemDamage& getPostMortemDamageModels() const { return _pmdModels; };
	TPostMortemDamage& getPostMortemDamageModelsMutable() { return _pmdModels; };

	bool hasPMD() const;
	bool recalibrationChangesQualities() const;

	coretools::Probability getErrorRate(const BAM::TSequencedBase & base) const;
	coretools::Probability getErrorWithPMD(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability getPhredIntWithPMD(const BAM::TSequencedBase & base) const;
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
