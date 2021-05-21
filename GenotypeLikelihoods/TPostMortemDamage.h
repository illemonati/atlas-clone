	/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPOSTMORTEMDAMAGE_H_
#define TPOSTMORTEMDAMAGE_H_


#include "TPMDTables.h"
#include "TSite.h"
//#include "auxiliaryTools.h"

#include <math.h>
#include <algorithm>
#include <memory>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


namespace GenotypeLikelihoods{

//Existing Functions
extern const std::string PMDFunctionName_none;
extern const std::string PMDFunctionName_empiric;
extern const std::string PMDFunctionName_exponential;
//extern const std:.string PMDFunctionName_Skoglund;

//Existing types
extern const std::string PMDTypeName_none;
extern const std::string PMDTypeName_singleStrand;
extern const std::string PMDTypeName_doubleStrand;

//Estimation Parameters
extern const std::string PMDEstimationExponential_epsilon;
extern const std::string PMDEstimationExponential_numNR;

//-------------------------------------
// TPMDEstimationParameters
//-------------------------------------
class TPMDEstimationParameters{
private:
	std::map<std::string, double> _parameters;

public:
	TPMDEstimationParameters();
	~TPMDEstimationParameters() = default;

	bool exists(const std::string & Parameter) const{
		return _parameters.find(Parameter) == _parameters.end();
	};

	bool add(const std::string & Parameter, const double & Value){
		return(_parameters.emplace(Parameter, Value).second);
	};

	double operator[](const std::string & Parameter) const{
		auto it =  _parameters.find(Parameter);
		if(it == _parameters.end()){
			throw std::runtime_error("TPMDEstimationParameters::double operator[](const std::string & Parameter) const: Parameter '" + Parameter + "' not set!");
		}
		return it->second;
	}
};

//---------------------------------------------------------------
//TPMDFunction
//---------------------------------------------------------------
//pure abstract base class
class TPMDFunction{
protected:
	std::vector<double> _parameters;
	void _parseParameters(const std::string & s);

public:
	TPMDFunction() = default;
	virtual ~TPMDFunction() = default;

	virtual bool hasDamage() const = 0;
	virtual std::string name() const = 0;
	virtual std::string example() = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile) = 0;
	virtual void learn(const TPMDTable & Table, const BAM::Base & from, const BAM::Base & to, const TPMDEstimationParameters & EstimationParameters) = 0;
	std::string string(){ return name() + "[" + concatenateString(_parameters, ",") + "]"; };

	virtual double prob(const uint16_t & pos) const = 0;
};

class TPMDFunctionNoPMD: public TPMDFunction{
public:
	TPMDFunctionNoPMD(const std::string & string);
	~TPMDFunctionNoPMD() = default;

	bool hasDamage() const override { return false; };
	std::string name() const override { return PMDFunctionName_none; };
	std::string example(){ return PMDFunctionName_none; };

	void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile){};
	void learn(const TPMDTable & Table, const BAM::Base & from, const BAM::Base & to, const TPMDEstimationParameters & EstimationParameters){};

	double prob(const uint16_t & pos) const override { return 0.0; };
};

/*
class TPMDFunctionSkoglund:public TPMDFunction{
public:
	TPMDFunctionSkoglund(const std::string & string);
	~TPMDFunctionSkoglund() = default;

	bool hasDamage() const override { return true; };
	std::string name() const override { return PMDFunctionName_Skoglund; };
	std::string example(){ return std::string(PMDFunctionName_Skoglund) + "[p,c]"; };

	void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile);
	void learn(const TPMDTable & table, const Base & from, const Base & to);

	double prob(const uint16_t & pos) const override;
};
*/

class TPMDFunctionExponential:public TPMDFunction{
private:
	uint16_t _lastPosition;
	std::vector<double> _probs;

	void _initialEstimatesOLS(const countVec & pmdCounts, const countVec& pmdSums, std::vector<double> & Parameters);
	void _fillFAndJacobian(arma::vec & F, arma::mat & J, const countVec & pmdCounts, const countVec& pmdSums, const std::vector<double> & Parameters);
	void _estimateWithNewtonRaphson(const countVec & pmdCounts, const countVec& pmdSums, std::vector<double> & Parameters, const uint32_t & numNRIterations, const double & epsilon);
	double _calcLL(const countVec & pmdCounts, const countVec& pmdSums, const std::vector<double> & Parameters);
	void _fillPMDProbabilities();

public:
	TPMDFunctionExponential(const std::string & string);
	~TPMDFunctionExponential() = default;

	bool hasDamage() const override { return true; };
	std::string name() const override { return PMDFunctionName_exponential; };
	std::string example(){ return std::string(PMDFunctionName_exponential) + "[a,b,c]"; };

	void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile);
	void learn(const TPMDTable & Table, const BAM::Base & from, const BAM::Base & to, const TPMDEstimationParameters & EstimationParameters);

	double prob(const uint16_t & pos) const override;
};

class TPMDFunctionEmpiric:public TPMDFunction{
public:
	TPMDFunctionEmpiric(const std::string & string);
	~TPMDFunctionEmpiric(){};

	bool hasDamage() const override { return true; };
	std::string name() const override { return PMDFunctionName_empiric; };
	std::string example(){ return std::string(PMDFunctionName_empiric) + "[p1,p2,...]"; };

	void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile){};
	void learn(const TPMDTable & Table, const BAM::Base & from, const BAM::Base & to, const TPMDEstimationParameters & EstimationParameters);

	double prob(const uint16_t & pos) const override;
};

//------------------------------------------------
// TPMDType
// pure abstract base class
//------------------------------------------------
class TPMDType{
protected:
	void _initializeFunction(const std::string & pmdString, std::unique_ptr<TPMDFunction> & ptr);

public:
	TPMDType() = default;
	virtual ~TPMDType() = default;

	virtual bool hasDamage() const = 0;
	virtual std::string type() const = 0;
	virtual std::string functionString() const = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile){};
	virtual void estimate(const TPMDTableReadGroup & PMDTable, const TPMDEstimationParameters & EstimationParameters){};

	virtual void fillBaseLikelihoods(const BAM::TSequencedBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const = 0;

	virtual void simulatePMD(BAM::TSequencedBase & base, TRandomGenerator & RandomGenerator) const = 0;
	virtual void simulatePMD(BAM::Base & base, const uint16_t & DistFrom5Prime, const uint16_t & DistFrom3Prime, const bool & IsReverseStrand, TRandomGenerator & RandomGenerator) const = 0;
};

//------------------------------------------------
// TPMDTypeNone
//------------------------------------------------
class TPMDTypeNone: public TPMDType{
public:
	TPMDTypeNone() = default;
	~TPMDTypeNone() = default;

	bool hasDamage() const override { return false; };
	std::string type() const override { return PMDTypeName_none; };
	std::string functionString() const override { return "none"; };

	void fillBaseLikelihoods(const BAM::TSequencedBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const override {
		//just copy
		baseLikelihoods = baseLikelihoodsNoPMD;
	};

	void simulatePMD(BAM::TSequencedBase & base, TRandomGenerator & RandomGenerator) const override {};
	void simulatePMD(BAM::Base & base, const uint16_t & DistFrom5Prime, const uint16_t & DistFrom3Prime, const bool & IsReverseStrand, TRandomGenerator & RandomGenerator) const override {};
};

//------------------------------------------------------
//TPMDTypeDoubleStrand
//------------------------------------------------------
class TPMDTypeDoubleStrand:public TPMDType{
private:
	std::unique_ptr<TPMDFunction> _pmdCT;
	std::unique_ptr<TPMDFunction> _pmdGA;

public:
	TPMDTypeDoubleStrand(const std::vector<std::string> & Details);
	~TPMDTypeDoubleStrand() = default;

	bool hasDamage() const override { return _pmdCT->hasDamage() || _pmdGA->hasDamage(); };
	std::string type() const override { return PMDTypeName_doubleStrand; };
	std::string functionString() const override;

	void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile) override;
	void estimate(const TPMDTableReadGroup & PMDTable, const TPMDEstimationParameters & EstimationParameters) override;

	void fillBaseLikelihoods(const BAM::TSequencedBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const override;

	void simulatePMD(BAM::TSequencedBase & base, TRandomGenerator & RandomGenerator) const override;
	void simulatePMD(BAM::Base & base, const uint16_t & DistFrom5Prime, const uint16_t & DistFrom3Prime, const bool & IsReverseStrand, TRandomGenerator & RandomGenerator) const override;
};

//------------------------------------------------------
//TPostMortemDamage
//------------------------------------------------------
class TPostMortemDamage{
private:
	std::vector< std::shared_ptr<TPMDType> > _pmdObjects;
	bool _hasPMD;

	void _createPMDType(const std::string & pmdString, std::shared_ptr<TPMDType> & ptr);
	void _initializeFromString(const std::string & pmdString, TLog* logfile);
	void _initializeFromFile(const BAM::TReadGroups & ReadGroups, const std::string & filename, TLog* logfile, std::vector<uint16_t> & ReadGroupsWithoutPMD);
	void _setHasDamage();

public:
	TPostMortemDamage();

	bool hasPMD() const{ return _hasPMD; };
	const TPMDType& operator[](const uint16_t & ReadGroupIndex) const { return *_pmdObjects[ReadGroupIndex]; };
	TPMDType& operator[](const uint16_t & ReadGroupIndex){ return *_pmdObjects[ReadGroupIndex]; };

	void initialize(const std::string & pmdString, const BAM::TReadGroups & ReadGroups, TLog* Logfile, std::vector<uint16_t> & ReadGroupsWithoutPMD);
	void writeToFile(const BAM::TReadGroups & ReadGroups, const std::string filename);
	void writeToFile(const BAM::TReadGroups & ReadGroups, const BAM::TReadGroupMap & ReadGroupMap, const std::string filename);
	void fillBaseLikelihoods(const BAM::TSequencedBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const;
};

}; //end namespace


#endif /* TPOSTMORTEMDAMAGE_H_ */
