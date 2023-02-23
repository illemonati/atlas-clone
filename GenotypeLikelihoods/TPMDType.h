
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

class TPMDType {
public:
	TPMDType()          = default;
	virtual ~TPMDType() = default;

	virtual bool hasDamage() const noexcept             = 0;
	virtual std::string functionString() const noexcept = 0;
	virtual std::string typeString() const noexcept     = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters)                   = 0;
	virtual void estimate(const PMDTable_RG &PMDTable, const TPMDEstimationParameters &EstimationParameters) = 0;

	virtual TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data,
												const TBaseLikelihoods &baseLikelihoodsNoPMD) const = 0;
	virtual TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &data, 
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

	static inline const std::string name = "none";
	TPMDTypeNone()                       = default;

	bool hasDamage() const noexcept override { return false; }
	std::string functionString() const noexcept override { return "none"; }
	std::string typeString() const noexcept override { return name; }

	void parseEstimationParameters(TPMDEstimationParameters &) override {}
	void estimate(const PMDTable_RG &, const TPMDEstimationParameters &) override {}

	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		// just copy
		return baseLikelihoodsNoPMD;
	}

	TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &, const TBaseLikelihoods &) const override {
		return massFunctions[b];
	}

	virtual void simulate(BAM::TSequencedBase &) const override {}
};

//------------------------------------------------------
// TPMDTypeDoubleStrand
//------------------------------------------------------
class TPMDTypeDoubleStrand final : public TPMDType {
private:
	std::unique_ptr<TPMDFunction> _pmdCT;
	std::unique_ptr<TPMDFunction> _pmdGA;

	coretools::Probability _probCT(const BAM::TSequencedBase &data) const noexcept {
		using coretools::Probability;
		if (data.distFrom3Prime < data.distFrom5Prime) { //from 3
			return !data.isReverseStrand() ? Probability{} : _pmdGA->prob(data.distFrom3Prime);
		} else { //from 5
			return !data.isReverseStrand() ? _pmdCT->prob(data.distFrom5Prime) : Probability{};
		}
	}

	coretools::Probability _probGA(const BAM::TSequencedBase &data) const noexcept {
		using coretools::Probability;
		if (data.distFrom3Prime < data.distFrom5Prime) { //from 3
			return !data.isReverseStrand() ? _pmdGA->prob(data.distFrom3Prime) : Probability{};
		} else { //from 5
			return !data.isReverseStrand() ? Probability{} : _pmdCT->prob(data.distFrom5Prime);
		}
	}

public:
	static inline const std::string name = "doubleStrand";
	TPMDTypeDoubleStrand(const std::vector<std::string> &Details);

	bool hasDamage() const noexcept override { return _pmdCT->hasDamage() || _pmdGA->hasDamage(); };
	std::string functionString() const noexcept override {
		return name + ":" + _pmdCT->string() + ":" + _pmdGA->string();
	}
	std::string typeString() const noexcept override { return name; }

	void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) override;
	void estimate(const PMDTable_RG &PMDTable, const TPMDEstimationParameters &EstimationParameters) override;

	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &data,
									   const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	virtual void simulate(BAM::TSequencedBase &data) const override;
};

//------------------------------------------------------
// TPMDTypeSingleStrand
//------------------------------------------------------
class TPMDTypeSingleStrand final : public TPMDType {
private:
	std::unique_ptr<TPMDFunction> _pmdCT3;
	std::unique_ptr<TPMDFunction> _pmdCT5;

	coretools::Probability _probCT(const BAM::TSequencedBase &data) const noexcept {
		if (data.isReverseStrand()) return coretools::Probability{};
		return data.distFrom3Prime < data.distFrom5Prime ? _pmdCT3->prob(data.distFrom3Prime)
														 : _pmdCT5->prob(data.distFrom5Prime);
	}

	coretools::Probability _probGA(const BAM::TSequencedBase &data) const noexcept {
		if (!data.isReverseStrand()) return coretools::Probability{};
		return data.distFrom3Prime < data.distFrom5Prime ? _pmdCT3->prob(data.distFrom3Prime)
														 : _pmdCT5->prob(data.distFrom5Prime);
	}

public:
	static inline const std::string name = "singleStrand";
	TPMDTypeSingleStrand(const std::vector<std::string> &Details);

	bool hasDamage() const noexcept override { return _pmdCT3->hasDamage() || _pmdCT5->hasDamage(); };
	std::string functionString() const noexcept override {
		return name + ":" + _pmdCT3->string() + ":" + _pmdCT5->string();
	}
	std::string typeString() const noexcept override { return name; }

	void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) override;
	void estimate(const PMDTable_RG &PMDTable, const TPMDEstimationParameters &EstimationParameters) override;

	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &data,
									   const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	virtual void simulate(BAM::TSequencedBase &data) const override;
};

TPMDType *makeType(std::string_view pmdString);
}
#endif
