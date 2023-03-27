
/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPMDTYPE_H_
#define TPMDTYPE_H_

#include "TGenotypeData.h"
#include "TPMDFunction.h"
#include "TSequencedBase.h"
#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods {

enum class ReadEnd : size_t { min = 0, forward5 = min, reverse5, forward3, reverse3, max };
class TPMDType {
	static constexpr size_t _N = coretools::index(genometools::Base::max) + 1;

public:
	using PMDTable      = std::vector<coretools::TStrongArray<
		coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, _N>, genometools::Base, _N>,
		GenotypeLikelihoods::ReadEnd>>;
	TPMDType()          = default;
	virtual ~TPMDType() = default;

	virtual bool hasDamage() const noexcept              = 0;
	virtual std::string functionString() const noexcept  = 0;
	virtual std::string_view typeString() const noexcept = 0;

	virtual void estimate(const PMDTable &PMDTable) = 0;

	virtual TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data,
											 const TBaseLikelihoods &baseLikelihoodsNoPMD) const = 0;
	virtual TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &data,
											const TBaseLikelihoods &baseLikelihoodsNoPMD) const  = 0;
	virtual TBaseProbabilities massFunction(genometools::Genotype g, const BAM::TSequencedBase &data,
											const TBaseLikelihoods &baseLikelihoodsNoPMD) const  = 0;

	virtual void simulate(BAM::TSequencedBase &data) const = 0;
};

//------------------------------------------------
// TPMDTypeNone
//------------------------------------------------
class TPMDTypeNone final : public TPMDType {
public:
	static constexpr TBaseMassFunctions massFunctions{
		{TBaseProbabilities::normalize({1., 0., 0., 0.}), TBaseProbabilities::normalize({0., 1., 0., 0.}),
		 TBaseProbabilities::normalize({0., 0., 1., 0.}), TBaseProbabilities::normalize({0., 0., 0., 1.})}};

	static constexpr std::string_view name = "none";
	TPMDTypeNone()                         = default;

	bool hasDamage() const noexcept override { return false; }
	std::string functionString() const noexcept override { return "none"; }
	std::string_view typeString() const noexcept override { return name; }

	void estimate(const PMDTable &) override {}

	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		// just copy
		return baseLikelihoodsNoPMD;
	}

	TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &, const TBaseLikelihoods &) const override {
		return massFunctions[b];
	}
	TBaseProbabilities massFunction(genometools::Genotype g, const BAM::TSequencedBase &, 
											const TBaseLikelihoods &baseLikelihoodsNoPMD) const  override {
		return massFunction(g, baseLikelihoodsNoPMD);
	}

	static TBaseProbabilities massFunction(genometools::Genotype g, const TBaseLikelihoods &baseLikelihoodsNoPMD);
	virtual void simulate(BAM::TSequencedBase &) const override {}
};


TPMDType *makeType(std::string_view pmdString);
}
#endif
