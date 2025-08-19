/*
 * TGenotypeModel.h
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_

#include "coretools/Types/probability.h"

#include "genometools/Genotypes/Containers.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Ploidy.h"

#include "stattools/MLEInference/TNelderMead.h"

namespace GenotypeLikelihoods {

class TGenotypeDistribution {
public:
	TGenotypeDistribution()          = default;
	virtual ~TGenotypeDistribution() = default;
	virtual genometools::TGenotypeLikelihoods P_dij(const genometools::TBaseLikelihoods &BaseLikelihoods) const = 0;
	virtual coretools::Probability getGenotypeLikelihood(const genometools::TBaseLikelihoods &BaseLikelihoods,
														 genometools::Genotype Genotype) const                  = 0;
	virtual double normalize_add(genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref)         = 0;
	virtual double add(const genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Fef)             = 0;
	virtual void estimate()                                                                                     = 0;
	virtual std::string_view typeString() const noexcept                                                        = 0;
	virtual void log() const                                                                                    = 0;
	virtual genometools::Ploidy ploidy() const noexcept                                                         = 0;
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const                     = 0;
	virtual std::vector<double> pis() const                                                                     = 0;
	virtual void reset()                                                                                        = 0;
};

template<genometools::Ploidy P>
genometools::TGenotypeLikelihoods base2genotype(const genometools::TBaseLikelihoods &BaseLikelihoods) {
	using genometools::Base;
	using coretools::average;
	constexpr auto P0 = coretools::P(0.);
	if constexpr (P == genometools::Ploidy::haploid) {
		return genometools::TGenotypeLikelihoods({BaseLikelihoods[Base::A], P0, P0, P0, BaseLikelihoods[Base::C], P0, P0,
									 BaseLikelihoods[Base::G], P0, BaseLikelihoods[Base::T]});
	} else {
		return genometools::TGenotypeLikelihoods(
			{BaseLikelihoods[Base::A], average(BaseLikelihoods[Base::A], BaseLikelihoods[Base::C]),
			 average(BaseLikelihoods[Base::A], BaseLikelihoods[Base::G]),
			 average(BaseLikelihoods[Base::A], BaseLikelihoods[Base::T]), BaseLikelihoods[Base::C],
			 average(BaseLikelihoods[Base::C], BaseLikelihoods[Base::G]),
			 average(BaseLikelihoods[Base::C], BaseLikelihoods[Base::T]), BaseLikelihoods[Base::G],
			 average(BaseLikelihoods[Base::G], BaseLikelihoods[Base::T]), BaseLikelihoods[Base::T]});
	}
}

class THaploidDistribution final : public TGenotypeDistribution {
	genometools::TBaseProbabilities _pi{};
	genometools::TBaseData _piSum{};
public:
	static constexpr std::string_view name = "haploid";

	genometools::TGenotypeLikelihoods P_dij(const genometools::TBaseLikelihoods &BaseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const genometools::TBaseLikelihoods &BaseLikelihoods,
												 genometools::Genotype Genotype) const override;
	double normalize_add(genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	double add(const genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	genometools::Ploidy ploidy() const noexcept override {return genometools::Ploidy::haploid;}
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
};

class TDiploidDistribution final : public TGenotypeDistribution {
	constexpr static double _het   = 6. / 4 * 0.001; // initial educated guess: theta = 0.001
	constexpr static auto _pi_init = genometools::TGenotypeProbabilities::normalize({1., _het, _het, _het, 1., _het, _het, 1., _het, 1.});
	genometools::TGenotypeProbabilities _pi     = _pi_init;
	genometools::TGenotypeData _piSum{};

public:
	static constexpr std::string_view name = "diploid";

	genometools::TGenotypeLikelihoods P_dij(const genometools::TBaseLikelihoods &BaseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const genometools::TBaseLikelihoods &BaseLikelihoods,
												 genometools::Genotype Genotype) const override;
	double normalize_add(genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	double add(const genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	genometools::Ploidy ploidy() const noexcept override {return genometools::Ploidy::diploid;}
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
};

class THKY85 final : public TGenotypeDistribution {
	constexpr static double _mu_init      = 0.5;
	constexpr static double _theta_r_init = 1e-3;
	constexpr static double _theta_g_init = 1e-4;

	double _mu      = _mu_init;
	double _theta_r = _theta_r_init;
	double _theta_g = _theta_g_init;

	coretools::TStrongArray<genometools::TGenotypeProbabilities, genometools::Base> _pi;
	coretools::TStrongArray<genometools::TGenotypeData, genometools::Base> _likelihoodSum{};
	stattools::TNelderMead<3> _nelderMead;

public:
	static constexpr std::string_view name = "HKY85";
	THKY85() : THKY85(_mu_init, _theta_r_init, _theta_g_init) {}
	THKY85(double Mu, double Theta_r, double Theta_g);

	genometools::TGenotypeLikelihoods P_dij(const genometools::TBaseLikelihoods &BaseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const genometools::TBaseLikelihoods &BaseLikelihoods,
												 genometools::Genotype Genotype) const override;
	double normalize_add(genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	double add(const genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	genometools::Ploidy ploidy() const noexcept override {return genometools::Ploidy::diploid;}
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
};

class THKY85_mono final : public TGenotypeDistribution {
	constexpr static coretools::Probability _mu_init{1./3};
	constexpr static double _theta_init = 0.0001;

	coretools::Probability _mu = _mu_init;
	double _theta              = _theta_init;

	coretools::TStrongArray<genometools::TBaseProbabilities, genometools::Base> _pi;
	coretools::TStrongArray<genometools::TBaseData, genometools::Base> _likelihoodSum{};
	stattools::TNelderMead<2> _nelderMead;

public:
	static constexpr std::string_view name = "HKY85_mono";
	THKY85_mono();

	genometools::TGenotypeLikelihoods P_dij(const genometools::TBaseLikelihoods &BaseLikelihoods) const override;
	coretools::Probability getGenotypeLikelihood(const genometools::TBaseLikelihoods &BaseLikelihoods,
												 genometools::Genotype Genotype) const override;
	double normalize_add(genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	double add(const genometools::TGenotypeLikelihoods &Likelihoods, genometools::Base Ref) override;
	void estimate() override;
	std::string_view typeString() const noexcept override { return name; }
	void log() const override;
	genometools::Ploidy ploidy() const noexcept override {return genometools::Ploidy::haploid;}
	virtual void addHeader(std::vector<std::string> &Header, std::string_view Prefix) const override;
	std::vector<double> pis() const override;
	virtual void reset() override;
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEDISTRIBUTION_H_ */
