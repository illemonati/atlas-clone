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

#define PMDFunctionName_none "none"
#define PMDFunctionName_empiric "Empiric"
#define PMDFunctionName_exponential "Exponential"
#define PMDFunctionName_Skoglund "Skoglund"

#define PMDTypeName_none "none"
#define PMDTypeName_singleStrand "singleStrand"
#define PMDTypeName_doubleStrand "doubleStrand"


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

	virtual void learn(const TPMDTable & table, const Base & from, const Base & to) = 0;
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

	void learn(const TPMDTable & table, const Base & from, const Base & to){};

	double prob(const uint16_t & pos) const override { return 0.0; };
};

class TPMDFunctionSkoglund:public TPMDFunction{
public:
	TPMDFunctionSkoglund(const std::string & string);
	~TPMDFunctionSkoglund() = default;

	bool hasDamage() const override { return true; };
	std::string name() const override { return PMDFunctionName_Skoglund; };
	std::string example(){ return std::string(PMDFunctionName_Skoglund) + "[p,c]"; };

	void learn(const TPMDTable & table, const Base & from, const Base & to);

	double prob(const uint16_t & pos) const override;
};

class TPMDFunctionExponential:public TPMDFunction{
public:
	TPMDFunctionExponential(const std::string & string);
	~TPMDFunctionExponential() = default;

	bool hasDamage() const override { return true; };
	std::string name() const override { return PMDFunctionName_exponential; };
	std::string example(){ return std::string(PMDFunctionName_exponential) + "[a,b,c]"; };

	void learn(const TPMDTable & table, const Base & from, const Base & to);

	double prob(const uint16_t & pos) const override;
};

class TPMDFunctionEmpiric:public TPMDFunction{
public:
	TPMDFunctionEmpiric(const std::string & string);
	~TPMDFunctionEmpiric(){};

	bool hasDamage() const override { return true; };
	std::string name() const override { return PMDFunctionName_empiric; };
	std::string example(){ return std::string(PMDFunctionName_empiric) + "[p1,p2,...]"; };

	void learn(const TPMDTable & table, const Base & from, const Base & to);

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
	std::unique_ptr<TPMDFunction> pmdCT;
	std::unique_ptr<TPMDFunction> pmdGA;

public:
	TPMDTypeDoubleStrand(const std::string & string);
	~TPMDTypeDoubleStrand() = default;

	bool hasDamage() const override { return pmdCT->hasDamage() || pmdGA->hasDamage(); };
	std::string type() const override { return PMDTypeName_doubleStrand; };
	std::string functionString() const override;

	void fillBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const override;
};

//------------------------------------------------------
//TPostMortemDamage
//------------------------------------------------------
class TPostMortemDamage{
private:
	std::vector< std::unique_ptr<TPMDType> > _pmdObjects;
	bool _hasPMD;

	void _createPMDType(const std::string & type, const std::string & functions, std::unique_ptr<TPMDType> & ptr);
	void _initializeFromFile(BAM::TReadGroups & ReadGroups, const std::string filename, TLog* logfile);


public:
	TPostMortemDamage();
	bool hasPMD() const{ return _hasPMD; };

	void initialize(TParameters & params, BAM::TReadGroups & ReadGroups, TLog* logfile);

	void writeToFile(const BAM::TReadGroups & ReadGroups, const std::string filename);

	void calculateBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const;
};

}; //end namespace


#endif /* TPOSTMORTEMDAMAGE_H_ */
