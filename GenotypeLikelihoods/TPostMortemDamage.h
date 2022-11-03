/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPOSTMORTEMDAMAGE_H_
#define TPOSTMORTEMDAMAGE_H_

#include "TSequencedBase.h"
#include "coretools/Containers/TMassFunction.h"
#include "TReadGroupInfo.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"
#include <array>
#include <map>
#include <memory>
#include <stdint.h>
#include <string>
#include <vector>

#include <armadillo>

#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TPMDTables.h"
#include "TReadGroups.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM {
class TSequencedBase;
}

namespace GenotypeLikelihoods {

using TPMDEstimationParameters = std::map<std::string, double>;

//---------------------------------------------------------------
// TPMDFunction
//---------------------------------------------------------------
// pure abstract base class
class TPMDFunction {
public:
	TPMDFunction()          = default;
	virtual ~TPMDFunction() = default;

	virtual bool hasDamage() const noexcept = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) = 0;
	virtual void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
					   const TPMDEstimationParameters &EstimationParameters)               = 0;
	virtual std::string string() const noexcept                                            = 0;
	virtual coretools::Probability prob(uint16_t pos) const noexcept                                       = 0;
};

class TPMDFunctionNoPMD final : public TPMDFunction {
public:
	static inline const std::string name    = "none";
	static inline const std::string example = name;
	TPMDFunctionNoPMD(const std::string &string);

	bool hasDamage() const noexcept override { return false; }
	std::string string() const noexcept override { return name + "[]"; }

	void parseEstimationParameters(TPMDEstimationParameters &) override{};
	void learn(const TPMDTable &, const genometools::Base &, const genometools::Base &,
			   const TPMDEstimationParameters &) override{};

	coretools::Probability prob(uint16_t) const noexcept override { return 0.0; }
};

class TPMDFunctionExponential final : public TPMDFunction {
private:
	uint16_t _lastPosition;
	double _a, _b, _c;
	std::vector<coretools::Probability> _values;

	void _initialEstimatesOLS(const countVec &pmdCounts, const countVec &pmdSums, std::array<double, 3> &Parameters);
	void _fillFAndJacobian(arma::vec &F, arma::mat &J, const countVec &pmdCounts, const countVec &pmdSums,
						   const std::array<double, 3> &Parameters);
	void _estimateWithNewtonRaphson(const countVec &pmdCounts, const countVec &pmdSums,
									std::array<double, 3> &Parameters, uint32_t numNRIterations, double epsilon);
	double _calcLL(const countVec &pmdCounts, const countVec &pmdSums, const std::array<double, 3> &Parameters);
	void _fillPMDProbabilities();

public:
	static inline const std::string name    = "Exponential";
	static inline const std::string example = name + "[lastPosition,a,b,c]";
	static inline const std::string epsilon = "PMDExponentialEpsilon";
	static inline const std::string numNR   = "PMDExponentialNumNR";
	TPMDFunctionExponential(const std::string &string);

	bool hasDamage() const noexcept override { return _lastPosition > 0; }
	std::string string() const noexcept override {
		return name + "[" + coretools::str::toString(_lastPosition) + ',' +
			   coretools::str::concatenateString(std::vector{_a, _b, _c}, ",") + "]";
	}

	void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) override;
	void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
			   const TPMDEstimationParameters &EstimationParameters) override;

	coretools::Probability prob(uint16_t pos) const noexcept override;
};

class TPMDFunctionEmpiric final : public TPMDFunction {
private:
	std::vector<coretools::Probability> _values;

public:
	static inline const std::string name    = "Empiric";
	static inline const std::string example = name + "[p1,p2,...]";
	TPMDFunctionEmpiric(const std::string &string);

	bool hasDamage() const noexcept override { return _values.size() + _values.back().get() != 1.0; }
	std::string string() const noexcept override {
		return name + "[" + coretools::str::concatenateString(_values, ",") + "]";
	}

	void parseEstimationParameters(TPMDEstimationParameters &) override{};
	void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
			   const TPMDEstimationParameters &EstimationParameters) override;

	coretools::Probability prob(uint16_t pos) const noexcept override;
};

//------------------------------------------------
// TPMDType
// pure abstract base class
//------------------------------------------------
class TPMDType {
public:
	TPMDType()          = default;
	virtual ~TPMDType() = default;

	virtual bool hasDamage() const noexcept             = 0;
	virtual std::string functionString() const noexcept = 0;
	virtual std::string typeString() const noexcept     = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters)                   = 0;
	virtual void estimate(const PMDTable_RG &PMDTable, const TPMDEstimationParameters &EstimationParameters) = 0;

	virtual TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &data,
												const TBaseLikelihoods &baseLikelihoodsNoPMD) const = 0;
	virtual TBaseProbabilities getMassFunction(genometools::Base b, const BAM::TSequencedBase &data, 
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

	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		// just copy
		return baseLikelihoodsNoPMD;
	}

	TBaseProbabilities getMassFunction(genometools::Base b, const BAM::TSequencedBase &, const TBaseLikelihoods &) const override {
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
		return !data.isReverseStrand() ? _pmdCT->prob(data.distFrom5Prime) : _pmdGA->prob(data.distFrom3Prime);
	}

	coretools::Probability _probGA(const BAM::TSequencedBase &data) const noexcept {
		return !data.isReverseStrand() ? _pmdGA->prob(data.distFrom3Prime) : _pmdCT->prob(data.distFrom5Prime);
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

	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &data,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	TBaseProbabilities getMassFunction(genometools::Base b, const BAM::TSequencedBase &data,
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

	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &data,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	TBaseProbabilities getMassFunction(genometools::Base b, const BAM::TSequencedBase &data,
									   const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	virtual void simulate(BAM::TSequencedBase &data) const override;
};

//------------------------------------------------------
// TPostMortemDamage
//------------------------------------------------------
class TPostMortemDamage {
private:
	std::vector<std::unique_ptr<TPMDType>> _pmdObjects;
	bool _hasPMD = false;

	void _initializeFromString(const std::string &pmdString);
	std::vector<uint16_t> _initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename);
	void _setHasDamage();

public:
	TPostMortemDamage() = default;
	TPostMortemDamage(const std::string &pmdString, const BAM::TReadGroups &ReadGroups,
					  std::vector<uint16_t> &ReadGroupsWithoutPMD) {
		ReadGroupsWithoutPMD = initialize(pmdString, ReadGroups);
	}
	constexpr bool hasPMD() const noexcept { return _hasPMD; };
	const TPMDType &operator[](uint16_t ReadGroupIndex) const noexcept { return *_pmdObjects[ReadGroupIndex]; }
	TPMDType &operator[](uint16_t ReadGroupIndex) noexcept { return *_pmdObjects[ReadGroupIndex]; }

	std::vector<uint16_t> initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups);
	void initialize(BAM::RGInfo::TReadGroupInfo & RgInfo);
	void writeToFile(const BAM::TReadGroups &ReadGroups, const std::string filename) const;
	void writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
	                 const std::string filename) const;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &data,
	                                    const TBaseLikelihoods &baseLikelihoodsNoPMD) const;
	TBaseProbabilities getMassFunction(genometools::Base b, const BAM::TSequencedBase &data,
									   const TBaseLikelihoods &baseLikelihoodsNoPMD) const;

	std::string functionString() const noexcept {
		std::string r;
		for (auto &p: _pmdObjects) {
			r.append(p->functionString()).append(1, ';');
		}
		return r;
	}

};

}; // namespace GenotypeLikelihoods

#endif /* TPOSTMORTEMDAMAGE_H_ */
