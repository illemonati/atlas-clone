/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPOSTMORTEMDAMAGE_H_
#define TPOSTMORTEMDAMAGE_H_

#include "TMassFunction.h"
#include "TStrongArray.h"
#include "probability.h"
#include <array>
#include <map>
#include <memory>
#include <stdint.h>
#include <string>
#include <vector>

#include <armadillo>

#include "GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TPMDTables.h"
#include "TReadGroups.h"
#include "stringFunctions.h"

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
	virtual double prob(uint16_t pos) const noexcept                                       = 0;
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

	double prob(uint16_t) const noexcept override { return 0.0; }
};

class TPMDFunctionExponential final : public TPMDFunction {
private:
	uint16_t _lastPosition;
	double _a, _b, _c;
	std::vector<double> _probs;

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

	double prob(uint16_t pos) const noexcept override;
};

class TPMDFunctionEmpiric final : public TPMDFunction {
private:
	std::vector<double> _parameters;

public:
	static inline const std::string name    = "Empiric";
	static inline const std::string example = name + "[p1,p2,...]";
	TPMDFunctionEmpiric(const std::string &string);

	bool hasDamage() const noexcept override { return _parameters.size() + _parameters.back() != 1.0; }
	std::string string() const noexcept override {
		return name + "[" + coretools::str::concatenateString(_parameters, ",") + "]";
	}

	void parseEstimationParameters(TPMDEstimationParameters &) override{};
	void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
			   const TPMDEstimationParameters &EstimationParameters) override;

	double prob(uint16_t pos) const noexcept override;
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

	virtual TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base,
												const TBaseLikelihoods &baseLikelihoodsNoPMD) const = 0;
	virtual TBaseMassFunctions getMassFunctions(const BAM::TSequencedBase &base) const              = 0;

	virtual void simulate(BAM::TSequencedBase &base) const   = 0;
	virtual void simulate(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
						  const bool &IsReverseStrand) const = 0;
};

//------------------------------------------------
// TPMDTypeNone
//------------------------------------------------
class TPMDTypeNone final : public TPMDType {
	static constexpr TBaseMassFunctions _massFunctions{
		{TBaseProbabilities{{1., 0., 0., 0.}}, {{0., 1., 0., 0.}}, {{0., 0., 1., 0.}}, {{0., 0., 1., 0.}}}};

public:
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

	TBaseMassFunctions getMassFunctions(const BAM::TSequencedBase &) const override { return _massFunctions; }

	void simulate(BAM::TSequencedBase &) const override {}
	void simulate(genometools::Base &, uint16_t, uint16_t, const bool &) const override {}
};

//------------------------------------------------------
// TPMDTypeDoubleStrand
//------------------------------------------------------
class TPMDTypeDoubleStrand final : public TPMDType {
private:
	std::unique_ptr<TPMDFunction> _pmdCT;
	std::unique_ptr<TPMDFunction> _pmdGA;

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

	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	TBaseMassFunctions getMassFunctions(const BAM::TSequencedBase &base) const override;

	void simulate(BAM::TSequencedBase &base) const override;
	void simulate(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
				  const bool &IsReverseStrand) const override;
};

//------------------------------------------------------
// TPMDTypeSingleStrand
//------------------------------------------------------
class TPMDTypeSingleStrand final : public TPMDType {
private:
	std::unique_ptr<TPMDFunction> _pmdCT3;
	std::unique_ptr<TPMDFunction> _pmdCT5;

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

	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base,
										const TBaseLikelihoods &baseLikelihoodsNoPMD) const override;

	TBaseMassFunctions getMassFunctions(const BAM::TSequencedBase &base) const override;

	void simulate(BAM::TSequencedBase &base) const override;
	void simulate(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
				  const bool &IsReverseStrand) const override;
};

//------------------------------------------------------
// TPostMortemDamage
//------------------------------------------------------
class TPostMortemDamage {
private:
	std::vector<std::unique_ptr<TPMDType>> _pmdObjects;
	bool _hasPMD = false;

	void _initializeFromString(const std::string &pmdString);
	void _initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename,
							 std::vector<uint16_t> &ReadGroupsWithoutPMD);
	void _setHasDamage();

public:
	TPostMortemDamage() = default;
	TPostMortemDamage(const std::string &pmdString, const BAM::TReadGroups &ReadGroups,
					  std::vector<uint16_t> &ReadGroupsWithoutPMD) {
		initialize(pmdString, ReadGroups, ReadGroupsWithoutPMD);
	}
	constexpr bool hasPMD() const noexcept { return _hasPMD; };
	const TPMDType &operator[](uint16_t ReadGroupIndex) const noexcept { return *_pmdObjects[ReadGroupIndex]; }
	TPMDType &operator[](uint16_t ReadGroupIndex) noexcept { return *_pmdObjects[ReadGroupIndex]; }

	void initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups,
					std::vector<uint16_t> &ReadGroupsWithoutPMD);
	void writeToFile(const BAM::TReadGroups &ReadGroups, const std::string filename) const;
	void writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
	                 const std::string filename) const;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base,
	                                    const TBaseLikelihoods &baseLikelihoodsNoPMD) const;
};

}; // namespace GenotypeLikelihoods

#endif /* TPOSTMORTEMDAMAGE_H_ */
