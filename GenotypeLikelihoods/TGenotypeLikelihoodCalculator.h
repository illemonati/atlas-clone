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

class TGenotypeLikelihoodCalculator_simple{
private:
	mutable std::vector<TBaseData> _baseLikelihoods;
	mutable TBaseData _baseLikelihoodsNoPMD;

public:

	template <typename PMD, typename SEQERR>
	void fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseData & BaseLikelihoods, const PMD & PmdModels, const SEQERR & SequencingErrorModels) const{
		SequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD);
		PmdModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD, BaseLikelihoods);
	};

	template <typename PMD, typename SEQERR>
	void fillGenotypeLikelihoods(const TSite &site, TGenotypeLikelihoods &genotypeLikelihoods, const PMD & PmdModels, const SEQERR & SequencingErrorModels) const{
		//ensure base likelihoods have proper size
		if(_baseLikelihoods.size() < site.depth()){
			_baseLikelihoods.resize(site.depth());
		}

		if(site.empty()){
			genotypeLikelihoods.reset();
		} else {
			//calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
			for(size_t i=0; i<site.depth(); ++i){
				fillBaseLikelihoods(site[i], _baseLikelihoods[i], PmdModels, SequencingErrorModels);
			}

			//calculate genotype likelihoods
			genotypeLikelihoods.fill(_baseLikelihoods, site.depth());
		}
	};
};

class TGenotypeLikelihoodCalculator{
protected:
	bool _initialized;
	TGenotypeLikelihoodCalculator_simple _calculator;

	TGenotypeDistribution _genotypeDistribution;
	TPostMortemDamage _pmdModels;
	TSequencingErrorModels _sequencingErrorModels;

public:
	TGenotypeLikelihoodCalculator();
	TGenotypeLikelihoodCalculator(TParameters & params, const BAM::TReadGroups* ReadGroups, TLog* Logfile);
	~TGenotypeLikelihoodCalculator() = default;

	void init(TParameters & params, const BAM::TReadGroups* ReadGroups, TLog* Logfile);
	const TSequencingErrorModels& getSequencingErrorModels() const { return _sequencingErrorModels; };
	TSequencingErrorModels& getSequencingErrorModelsMutable() { return _sequencingErrorModels; };
	const TPostMortemDamage& getPostMortemDamageModels() const { return _pmdModels; };
	TPostMortemDamage& getPostMortemDamageModelsMutable() { return _pmdModels; };

	bool hasPMD() const;
	bool recalibrationChangesQualities() const;

	double getErrorRate(const BAM::TSequencedBase & base) const;
	double getErrorWithPMD(const BAM::TSequencedBase & base) const;
	uint8_t getPhredInt(const BAM::TSequencedBase & base) const;
	uint8_t getPhredIntWithPMD(const BAM::TSequencedBase & base) const;
	void recalibrate(BAM::TSequencedBase & base) const;
	void recalibrateWithPMD(BAM::TSequencedBase & base) const;

	//are vector versions used??
	void recalibrate(std::vector<BAM::TSequencedBase> & bases) const;
	void recalibrateWithPMD(std::vector<BAM::TSequencedBase> & bases) const;
	double calculateLogPMDS(const BAM::TSequencedBase & base, const Base ref, const double pi) const; //TODO: move to PMDS class?

	//functions performed on sites
	void calculateGenotypeLikelihoods(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) const;
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPELIKELIHOODCALCULATOR_H_ */
