/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPOSTMORTEMDAMAGE_H_
#define TPOSTMORTEMDAMAGE_H_

#include "GenotypeTypes.h"
#include "TPMDTables.h"
#include "TSite.h"
//#include "auxiliaryTools.h"

#include <algorithm>
#include <math.h>
#include <memory>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

namespace GenotypeLikelihoods {

using coretools::TLog;
using coretools::TParameters;
using coretools::TRandomGenerator;

// Existing Functions
extern const std::string PMDFunctionName_none;
extern const std::string PMDFunctionName_empiric;
extern const std::string PMDFunctionName_exponential;
// extern const std:.string PMDFunctionName_Skoglund;

// Existing types
extern const std::string PMDTypeName_none;
extern const std::string PMDTypeName_singleStrand;
extern const std::string PMDTypeName_doubleStrand;

// Estimation Parameters
extern const std::string PMDEstimationExponential_epsilon;
extern const std::string PMDEstimationExponential_numNR;

using TPMDEstimationParameters = std::map<std::string, double>;

//---------------------------------------------------------------
// TPMDFunction
//---------------------------------------------------------------
// pure abstract base class
class TPMDFunction {
public:
	TPMDFunction()          = default;
	virtual ~TPMDFunction() = default;

	virtual bool hasDamage() const noexcept      = 0;
	virtual std::string name() const noexcept    = 0;
	virtual std::string example() const noexcept = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters, TParameters &Params,
					       TLog *Logfile)                    = 0;
	virtual void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
			   const TPMDEstimationParameters &EstimationParameters) = 0;
	virtual std::string string() const noexcept = 0;
	virtual double prob(uint16_t pos) const noexcept = 0;
};

class TPMDFunctionNoPMD : public TPMDFunction {
private: 
	std::vector<double> _parameters;
public:
	TPMDFunctionNoPMD(const std::string &string);
	~TPMDFunctionNoPMD() = default;

	bool hasDamage() const noexcept override { return false; };
	std::string name() const noexcept override { return PMDFunctionName_none; };
	std::string example() const noexcept override { return PMDFunctionName_none; };
	std::string string() const noexcept override { return name() + "[" + coretools::concatenateString(_parameters, ",") + "]"; };


	void parseEstimationParameters(TPMDEstimationParameters &, TParameters &, TLog *) override{};
	void learn(const TPMDTable &, const genometools::Base &, const genometools::Base &,
		   const TPMDEstimationParameters &) override{};

	double prob(uint16_t) const noexcept override { return 0.0; };
};

/*
class TPMDFunctionSkoglund:public TPMDFunction{
public:
    TPMDFunctionSkoglund(const std::string & string);
    ~TPMDFunctionSkoglund() = default;

    bool hasDamage() const override { return true; };
    std::string name() const override { return PMDFunctionName_Skoglund; };
    std::string example(){ return std::string(PMDFunctionName_Skoglund) + "[p,c]"; };

    void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog*
Logfile); void learn(const TPMDTable & table, const Base & from, const Base & to);

    double prob(uint16_t pos) const override;
};
*/

class TPMDFunctionExponential : public TPMDFunction {
private:
	std::vector<double> _parameters;
	uint16_t _lastPosition;
	std::vector<double> _probs;

	void _initialEstimatesOLS(const countVec &pmdCounts, const countVec &pmdSums, std::vector<double> &Parameters);
	void _fillFAndJacobian(arma::vec &F, arma::mat &J, const countVec &pmdCounts, const countVec &pmdSums,
			       const std::vector<double> &Parameters);
	void _estimateWithNewtonRaphson(const countVec &pmdCounts, const countVec &pmdSums, std::vector<double> &Parameters,
					uint32_t numNRIterations, double epsilon);
	double _calcLL(const countVec &pmdCounts, const countVec &pmdSums, const std::vector<double> &Parameters);
	void _fillPMDProbabilities();

public:
	TPMDFunctionExponential(const std::string &string);
	~TPMDFunctionExponential() = default;

	bool hasDamage() const noexcept override { return true; };
	std::string name() const noexcept override { return PMDFunctionName_exponential; };
	std::string example() const noexcept override { return std::string(PMDFunctionName_exponential) + "[a,b,c]"; };
	std::string string() const noexcept override { return name() + "[" + coretools::concatenateString(_parameters, ",") + "]"; };

	void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters, TParameters &Params,
				       TLog *Logfile) override;
	void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
		   const TPMDEstimationParameters &EstimationParameters) override;

	double prob(uint16_t pos) const noexcept override;
};

class TPMDFunctionEmpiric : public TPMDFunction {
private:
	std::vector<double> _parameters;
public:
	TPMDFunctionEmpiric(const std::string &string);
	~TPMDFunctionEmpiric(){};

	bool hasDamage() const noexcept override { return true; };
	std::string name() const noexcept override { return PMDFunctionName_empiric; };
	std::string example() const noexcept override { return std::string(PMDFunctionName_empiric) + "[p1,p2,...]"; };
	std::string string() const noexcept override { return name() + "[" + coretools::concatenateString(_parameters, ",") + "]"; };

	void parseEstimationParameters(TPMDEstimationParameters &, TParameters &, TLog *) override{};
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
	virtual std::string type() const noexcept           = 0;
	virtual std::string functionString() const noexcept = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters, TParameters &Params,
					       TLog *Logfile) = 0;
	virtual void estimate(const TPMDTableReadGroup &PMDTable, const TPMDEstimationParameters &EstimationParameters) = 0;

	virtual void fillBaseLikelihoods(const BAM::TSequencedBase &base, const TBaseProbabilities &baseLikelihoodsNoPMD,
					 TBaseProbabilities &baseLikelihoods) const = 0;

	virtual void simulatePMD(BAM::TSequencedBase &base, TRandomGenerator &RandomGenerator) const   = 0;
	virtual void simulatePMD(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
				 const bool &IsReverseStrand, TRandomGenerator &RandomGenerator) const = 0;
};

//------------------------------------------------
// TPMDTypeNone
//------------------------------------------------
class TPMDTypeNone : public TPMDType {
public:
	TPMDTypeNone()  = default;
	~TPMDTypeNone() = default;

	bool hasDamage() const noexcept override { return false; };
	std::string type() const noexcept override { return PMDTypeName_none; };
	std::string functionString() const noexcept override { return "none"; };

	void parseEstimationParameters(TPMDEstimationParameters &, TParameters &, TLog *) override {}
	void estimate(const TPMDTableReadGroup &, const TPMDEstimationParameters &) override {}

	void fillBaseLikelihoods(const BAM::TSequencedBase &, const TBaseProbabilities &baseLikelihoodsNoPMD,
				 TBaseProbabilities &baseLikelihoods) const override {
		// just copy
		baseLikelihoods = baseLikelihoodsNoPMD;
	};

	void simulatePMD(BAM::TSequencedBase &, TRandomGenerator &) const override{};
	void simulatePMD(genometools::Base &, uint16_t, uint16_t, const bool &, TRandomGenerator &) const override{};
};

//------------------------------------------------------
// TPMDTypeDoubleStrand
//------------------------------------------------------
class TPMDTypeDoubleStrand : public TPMDType {
private:
	std::unique_ptr<TPMDFunction> _pmdCT;
	std::unique_ptr<TPMDFunction> _pmdGA;

public:
	TPMDTypeDoubleStrand(const std::vector<std::string> &Details);
	~TPMDTypeDoubleStrand() = default;

	bool hasDamage() const noexcept override { return _pmdCT->hasDamage() || _pmdGA->hasDamage(); };
	std::string type() const noexcept override { return PMDTypeName_doubleStrand; };
	std::string functionString() const noexcept override;

	void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters, TParameters &Params,
				       TLog *Logfile) override;
	void estimate(const TPMDTableReadGroup &PMDTable, const TPMDEstimationParameters &EstimationParameters) override;

	void fillBaseLikelihoods(const BAM::TSequencedBase &base, const TBaseProbabilities &baseLikelihoodsNoPMD,
				 TBaseProbabilities &baseLikelihoods) const override;

	void simulatePMD(BAM::TSequencedBase &base, TRandomGenerator &RandomGenerator) const override;
	void simulatePMD(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
			 const bool &IsReverseStrand, TRandomGenerator &RandomGenerator) const override;
};

//------------------------------------------------------
// TPMDTypeSingleStrand
//------------------------------------------------------
class TPMDTypeSingleStrand : public TPMDType {
private:
	std::unique_ptr<TPMDFunction> _pmdCT3;
	std::unique_ptr<TPMDFunction> _pmdCT5;

public:
	TPMDTypeSingleStrand(const std::vector<std::string> &Details);
	~TPMDTypeSingleStrand() = default;

	bool hasDamage() const noexcept override { return _pmdCT3->hasDamage() || _pmdCT5->hasDamage(); };
	std::string type() const noexcept override { return PMDTypeName_singleStrand; };
	std::string functionString() const noexcept override;

	void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters, TParameters &Params,
				       TLog *Logfile) override;

	void estimate(const TPMDTableReadGroup &PMDTable, const TPMDEstimationParameters &EstimationParameters) override;

	void fillBaseLikelihoods(const BAM::TSequencedBase &base, const TBaseProbabilities &baseLikelihoodsNoPMD,
				 TBaseProbabilities &baseLikelihoods) const override;

	void simulatePMD(BAM::TSequencedBase &base, TRandomGenerator &RandomGenerator) const override;
	void simulatePMD(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
			 const bool &IsReverseStrand, TRandomGenerator &RandomGenerator) const override;
};

//------------------------------------------------------
// TPostMortemDamage
//------------------------------------------------------
class TPostMortemDamage {
private:
	std::vector<std::shared_ptr<TPMDType>> _pmdObjects;
	bool _hasPMD = false;

	void _createPMDType(const std::string &pmdString, std::shared_ptr<TPMDType> &ptr);
	void _initializeFromString(const std::string &pmdString, TLog *logfile);
	void _initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename, TLog *logfile,
				 std::vector<uint16_t> &ReadGroupsWithoutPMD);
	void _setHasDamage();

public:
	TPostMortemDamage() = default;
	TPostMortemDamage(const std::string &pmdString, const BAM::TReadGroups &ReadGroups, TLog *Logfile,
			  std::vector<uint16_t> &ReadGroupsWithoutPMD) {
		initialize(pmdString, ReadGroups, Logfile, ReadGroupsWithoutPMD);
	}

	constexpr bool hasPMD() const noexcept { return _hasPMD; };
	const TPMDType &operator[](uint16_t ReadGroupIndex) const noexcept { return *_pmdObjects[ReadGroupIndex]; };
	TPMDType &operator[](uint16_t ReadGroupIndex) noexcept { return *_pmdObjects[ReadGroupIndex]; };

	void initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups, TLog *Logfile,
	                std::vector<uint16_t> &ReadGroupsWithoutPMD);
	void writeToFile(const BAM::TReadGroups &ReadGroups, const std::string filename) const ;
	void writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
	                 const std::string filename) const;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, const TBaseProbabilities &baseLikelihoodsNoPMD,
	                         TBaseProbabilities &baseLikelihoods) const;
};

}; // namespace GenotypeLikelihoods

#endif /* TPOSTMORTEMDAMAGE_H_ */
