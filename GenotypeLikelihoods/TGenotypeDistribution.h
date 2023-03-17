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
	virtual std::string_view typeString() const noexcept                                               = 0;
	virtual std::string definition() const noexcept                                                    = 0;
	virtual bool isInvariant() const noexcept                                                          = 0;
};

class THaploidDistribution final : public TGenotypeDistribution {
	TBaseProbabilities _pi{};
	TBaseData _piSum{};
public:
	static constexpr std::string_view name = "haploid";

	TGenotypeLikelihoods getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	std::string definition() const noexcept override;
	bool isInvariant() const noexcept override {return true;}
};

class TDiploidDistribution final : public TGenotypeDistribution {
	constexpr static double _het   = 6. / 4 * 0.001; // initial educated guess: theta = 0.001
	constexpr static auto _pi_init = TGenotypeProbabilities::normalize({1., _het, _het, _het, 1., _het, _het, 1., _het, 1.});
	TGenotypeProbabilities _pi     = _pi_init;
	TGenotypeData _piSum{};

public:
	static constexpr std::string_view name = "diploid";

	TGenotypeLikelihoods getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	std::string definition() const noexcept override;
	bool isInvariant() const noexcept override {return false;}
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
