/*
 * TGenotypeModel.h
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_

#include "GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "probability.h"

namespace GenotypeLikelihoods {

//-------------------------------------------
// TGenotypeDistribution
// Base class.
//-------------------------------------------

class TGenotypeDistribution {
public:
	TGenotypeDistribution()                                                                            = default;
	virtual ~TGenotypeDistribution()                                                                   = default;
	virtual TGenotypeLikelihoods getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const = 0;
	virtual double weightedSum(const TGenotypeLikelihoods &likelihoods) const            = 0;
	virtual coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
														 genometools::Genotype genotype) const         = 0;
	virtual double normalize(TGenotypeLikelihoods &likelihoods) const                                  = 0;
	virtual std::string typeString() const noexcept                                                    = 0;
	virtual bool isInvariant() const noexcept                                                          = 0;
};

class THaploidDistribution final : public TGenotypeDistribution {
	static constexpr TGenotypeProbabilities _weights{{0.25, 0., 0., 0., 0.25, 0., 0., 0.25, 0., 0.25}};
public:
	static inline const std::string name = "haploid";
	THaploidDistribution()               = default;

	double weightedSum(const TGenotypeLikelihoods &likelihoods) const override;
	TGenotypeLikelihoods getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	virtual double normalize(TGenotypeLikelihoods &likelihoods) const override;
	std::string typeString() const noexcept override { return name; }
	bool isInvariant() const noexcept override {return true;}
};

class TDiploidDistribution final : public TGenotypeDistribution {
	constexpr static TGenotypeProbabilities _weights{}; // 0.1 everywhere
public:
	static inline const std::string name = "diploid";
	TDiploidDistribution()               = default;

	double weightedSum(const TGenotypeLikelihoods &likelihoods) const override;
	TGenotypeLikelihoods getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	virtual double normalize(TGenotypeLikelihoods &likelihoods) const override;
	std::string typeString() const noexcept override { return name; }
	bool isInvariant() const noexcept override {return false;}
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
