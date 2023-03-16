/*
 * TGenotypeModel.h
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_

#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "coretools/Types/probability.h"

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
	virtual coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
														 genometools::Genotype genotype) const         = 0;
	virtual double normalize_add(TGenotypeLikelihoods &likelihoods)                                    = 0;
	virtual void estimate()                                                                            = 0;
	virtual std::string typeString() const noexcept                                                    = 0;
	virtual std::string definition() const noexcept                                                    = 0;
	virtual bool isInvariant() const noexcept                                                          = 0;
};

class THaploidDistribution final : public TGenotypeDistribution {
	TBaseProbabilities _pi{};
	TBaseData _piSum{};
public:
	static inline const std::string name = "haploid";

	TGenotypeLikelihoods getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods) override;
	void estimate() override;
	std::string typeString() const noexcept override { return name; }
	std::string definition() const noexcept override;
	bool isInvariant() const noexcept override {return true;}
};

class TDiploidDistribution final : public TGenotypeDistribution {
	TGenotypeProbabilities _pi{};
	TGenotypeData _piSum{};
public:
	static inline const std::string name = "diploid";

	TGenotypeLikelihoods getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods) override;
	void estimate() override;
	std::string typeString() const noexcept override { return name; }
	std::string definition() const noexcept override;
	bool isInvariant() const noexcept override {return false;}
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
