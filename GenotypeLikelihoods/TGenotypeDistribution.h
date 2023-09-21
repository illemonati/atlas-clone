/*
 * TGenotypeModel.h
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_

#include "coretools/Containers/TMassFunction.h"
#include "coretools/Containers/TStrongArray.h"
#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "coretools/Types/probability.h"
#include "stattools/MLEInference/TNelderMead.h"

namespace GenotypeLikelihoods {

//-------------------------------------------
// TGenotypeDistribution
// Base class.
//-------------------------------------------

class TGenotypeDistribution {
public:
	TGenotypeDistribution()                                                                            = default;
	virtual ~TGenotypeDistribution()                                                                   = default;
	virtual TGenotypeLikelihoods P_dij(const TBaseLikelihoods &baseLikelihoods) const = 0;
	virtual coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
														 genometools::Genotype genotype) const         = 0;
	virtual double normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base ref)             = 0;
	virtual void estimate()                                                                            = 0;
	virtual std::string_view typeString() const noexcept                                               = 0;
	virtual void log() const                                                                           = 0;
	virtual bool isInvariant() const noexcept                                                          = 0;
};

class THaploidDistribution final : public TGenotypeDistribution {
	TBaseProbabilities _pi{};
	TBaseData _piSum{};
public:
	static constexpr std::string_view name = "haploid";

	TGenotypeLikelihoods P_dij(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	bool isInvariant() const noexcept override {return true;}
};

class TDiploidDistribution final : public TGenotypeDistribution {
	constexpr static double _het   = 6. / 4 * 0.001; // initial educated guess: theta = 0.001
	constexpr static auto _pi_init = TGenotypeProbabilities::normalize({1., _het, _het, _het, 1., _het, _het, 1., _het, 1.});
	TGenotypeProbabilities _pi     = _pi_init;
	TGenotypeData _piSum{};

public:
	static constexpr std::string_view name = "diploid";

	TGenotypeLikelihoods P_dij(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	bool isInvariant() const noexcept override {return false;}
};

class THKY85 final : public TGenotypeDistribution {
	double _mu      = 1;
	double _theta_r = 0.0001;
	double _theta_g = 0.0001;

	coretools::TStrongArray<TGenotypeProbabilities, genometools::Base> _pi;
	coretools::TStrongArray<TGenotypeData, genometools::Base> _likelihoodSum{};
	stattools::TNelderMead<3> _nelderMead;

public:
	static constexpr std::string_view name = "HKY85";
	THKY85();

	TGenotypeLikelihoods P_dij(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	bool isInvariant() const noexcept override {return false;}
};

class THKY85_mono final : public TGenotypeDistribution {
	double _mu    = 1;
	double _theta = 0.0001;

	coretools::TStrongArray<TBaseProbabilities, genometools::Base> _pi;
	coretools::TStrongArray<TBaseData, genometools::Base> _likelihoodSum{};
	stattools::TNelderMead<2> _nelderMead;

public:
	static constexpr std::string_view name = "HKY85_mono";
	THKY85_mono();

	TGenotypeLikelihoods P_dij(const TBaseLikelihoods &baseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
												 genometools::Genotype genotype) const override;
	double normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	bool isInvariant() const noexcept override {return true;}
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
