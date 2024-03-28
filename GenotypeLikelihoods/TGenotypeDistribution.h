/*
 * TGenotypeModel.h
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_

#include "GenotypeData.h"

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Types/probability.h"

#include "genometools/GenotypeTypes.h"

#include "genometools/VCF/TVcfWriter.h"
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
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const            = 0;
	virtual std::vector<double> pis() const                                                            = 0;
	virtual void reset()                                                                               = 0;
};

template<genometools::Ploidy P>
TGenotypeLikelihoods base2genotype(const TBaseLikelihoods &baseLikelihoods) {
	using genometools::Base;
	using coretools::average;
	constexpr auto P0 = coretools::P(0.);
	if constexpr (P == genometools::Ploidy::haploid) {
		return TGenotypeLikelihoods({baseLikelihoods[Base::A], P0, P0, P0, baseLikelihoods[Base::C], P0, P0,
									 baseLikelihoods[Base::G], P0, baseLikelihoods[Base::T]});
	} else {
		return TGenotypeLikelihoods(
			{baseLikelihoods[Base::A], average(baseLikelihoods[Base::A], baseLikelihoods[Base::C]),
			 average(baseLikelihoods[Base::A], baseLikelihoods[Base::G]),
			 average(baseLikelihoods[Base::A], baseLikelihoods[Base::T]), baseLikelihoods[Base::C],
			 average(baseLikelihoods[Base::C], baseLikelihoods[Base::G]),
			 average(baseLikelihoods[Base::C], baseLikelihoods[Base::T]), baseLikelihoods[Base::G],
			 average(baseLikelihoods[Base::G], baseLikelihoods[Base::T]), baseLikelihoods[Base::T]});
	}
}

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
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
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
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
};

class THKY85 final : public TGenotypeDistribution {
	constexpr static coretools::Probability _mu_init{0.5};
	constexpr static double _theta_r_init = 1e-2;
	constexpr static double _theta_g_init = 1e-4;

	coretools::Probability _mu = _mu_init;
	double _theta_r            = _theta_r_init;
	double _theta_g            = _theta_g_init;

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
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
};

class THKY85_mono final : public TGenotypeDistribution {
	constexpr static coretools::Probability _mu_init{0.5};
	constexpr static double _theta_init = 0.0001;

	coretools::Probability _mu = _mu_init;
	double _theta              = _theta_init;

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
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
