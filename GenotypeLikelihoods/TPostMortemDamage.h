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
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


namespace GenotypeLikelihoods{

//Existing Functions
#define PMDFunctionName_none "none"
#define PMDFunctionName_empiric "Empiric"
#define PMDFunctionName_exponential "Exponential"
//#define PMDFunctionName_Skoglund "Skoglund"

//Existing types
#define PMDTypeName_none "none"
#define PMDTypeName_singleStrand "singleStrand"
#define PMDTypeName_doubleStrand "doubleStrand"

//Estimation Parameters
#define PMDEstimationExponential_epsilon "PMDExponentialEpsilon"
#define PMDEstimationExponential_numNR "PMDExponentialNumNR"

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
	virtual void learn(const TPMDTable & Table, const Base & from, const Base & to, const TPMDEstimationParameters & EstimationParameters) = 0;
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
	void learn(const TPMDTable & Table, const Base & from, const Base & to, const TPMDEstimationParameters & EstimationParameters){};

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
	void learn(const TPMDTable & Table, const Base & from, const Base & to, const TPMDEstimationParameters & EstimationParameters);

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
	void learn(const TPMDTable & Table, const Base & from, const Base & to, const TPMDEstimationParameters & EstimationParameters);

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
	virtual void estimate(const TPMDTables & PMDTables, const TPMDEstimationParameters & EstimationParameters){};

	virtual void fillBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const = 0;
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

	void fillBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const override {
		//just copy
		baseLikelihoods = baseLikelihoodsNoPMD;
	};
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

	void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile);
	void estimate(const TPMDTableReadGroup & PMDTable, const TPMDEstimationParameters & EstimationParameters);

	void fillBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const override;
};

//------------------------------------------------------
//TPostMortemDamage
//------------------------------------------------------
class TPostMortemDamage{
private:
	std::vector< std::unique_ptr<TPMDType> > _pmdObjects;
	bool _hasPMD;

	void _createPMDType(const std::string & pmdString, std::unique_ptr<TPMDType> & ptr);
	void _initializeFromString(const std::string & pmdString, TLog* logfile);
	void _initializeFromFile(const BAM::TReadGroups & ReadGroups, const std::string & filename, TLog* logfile);
	void _setHasDamage();

public:
	TPostMortemDamage();
	bool hasPMD() const{ return _hasPMD; };

	void initialize(const std::string & pmdString, const BAM::TReadGroups & ReadGroups, TLog* Logfile);
	void initialize(TParameters & params, const BAM::TReadGroups & ReadGroups, TLog* logfile);

	void parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile);
	void estimate(const TPMDTables & PMDTables, const BAM::TReadGroups & ReadGroups, TLog* logfile, const TPMDEstimationParameters & EstimationParameters);

	void writeToFile(const BAM::TReadGroups & ReadGroups, const std::string filename);

	void calculateBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const;
};

}; //end namespace


#endif /* TPOSTMORTEMDAMAGE_H_ */
