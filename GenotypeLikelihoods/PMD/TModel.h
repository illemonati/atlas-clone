
/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPMDTYPE_H_
#define TPMDTYPE_H_

#include "TGenotypeData.h"
#include "TFunction.h"
#include "TSequencedBase.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods::PMD {

enum class ReadEnd : size_t { min = 0, forward5 = min, reverse5, forward3, reverse3, max };
enum class Strand : size_t {min, Single=min, Double, max};
struct TModel {
	TModel()          = default;
	virtual ~TModel() = default;

	virtual bool hasDamage() const noexcept              = 0;
	virtual std::string functionString() const noexcept  = 0;
	virtual std::string_view typeString() const noexcept = 0;

	virtual void resize(size_t N)                                                                               = 0;
	virtual void estimate()                                                                                     = 0;
	virtual void add(genometools::Base from, BAM::TSequencedBase data)                                          = 0;
	virtual void writeTable(std::string_view name, std::array<coretools::TOutputFile, 2> &files) const noexcept = 0;

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
class TNoPMD final : public TModel {
public:
	static constexpr TBaseMassFunctions massFunctions{
		{TBaseProbabilities::normalize({1., 0., 0., 0.}), TBaseProbabilities::normalize({0., 1., 0., 0.}),
		 TBaseProbabilities::normalize({0., 0., 1., 0.}), TBaseProbabilities::normalize({0., 0., 0., 1.})}};

	static constexpr std::string_view name = "none";
	TNoPMD()                         = default;

	bool hasDamage() const noexcept override { return false; }
	std::string functionString() const noexcept override { return "none"; }
	std::string_view typeString() const noexcept override { return name; }

	virtual void resize(size_t) override {}
	void estimate() override {}
	void add(genometools::Base, BAM::TSequencedBase) override {}
	void writeTable(std::string_view, std::array<coretools::TOutputFile, 2>&) const noexcept override {}

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


TModel *makeType(std::string_view pmdString);
}
#endif
