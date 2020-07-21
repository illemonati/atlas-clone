/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPOSTMORTEMDAMAGE_H_
#define TPOSTMORTEMDAMAGE_H_

#include <math.h>
#include "TReadGroups.h"
#include "TGenotypeMap.h"
#include "TSite.h"
#include "auxiliaryTools.h"
#include <algorithm>
#include "TBase.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


namespace GenotypeLikelihoods{

enum PMDType {pmdCT=0, pmdGA, pmdGT, pmdCA};


//---------------------------------------------------------------
//TPMDTable
//---------------------------------------------------------------
class TPMDTable{
private:
	long*** counts; //they are [read group][reference][read]
	long** sums;
	bool sumsCalculated;
	int maxLength;
	TGenotypeMap genoMap;

	void calculateSums();
	void deleteSums();
	void fillFAndJacobian(arma::vec & F, arma::mat & J, Base & from, Base & to, double* oldParams);
	void fillF(arma::vec & F, Base & from, Base & to, double* oldParams);
	double calcLL(Base & from, Base & to, double* oldParams);

public:
	TPMDTable(int MaxLength);
	~TPMDTable();
	void empty();
	void add(const int & pos, const Base & ref, const Base & read);
	void writeTable(std::ofstream & out, std::string prefix);
	void writeTableWithCounts(std::ofstream & out, std::string prefix);
	std::string getPMDString(const Base first, const Base second);
	std::string fitExponentialModel(Base from, Base to, int & numNRIterations, double & eps, std::string readGroupName, int maxReadLength, TLog* logfile);
};

class TPMDTables{
private:
	BAM::TReadGroupMap* readGroupMap;
	BAM::TReadGroups* readGroups;
	int maxReadLength;
	int origNumReadGroups;
	int numReadGroups;
	TPMDTable** forward;
	TPMDTable** reverse;
	bool _initialized;

public:
	TPMDTables();
	TPMDTables(BAM::TReadGroups* ReadGroups, int maxLengthForInference, int MaxReadLength, BAM::TReadGroupMap* ReadGroupMapObject);
	~TPMDTables();

	void initialize(BAM::TReadGroups* ReadGroups, int maxLengthForInference, int MaxReadLength, BAM::TReadGroupMap* ReadGroupMapObject);

	void addFromFivePrime(const uint16_t readGroup, const uint16_t pos, const Base & ref, const Base & read);
	void addFromThreePrime(const uint16_t readGroup, const uint16_t pos, const Base & ref, const Base & read);
	void writePMDFile(std::string filename);
	void writeTable(std::string filename);
	void writeTableWithCounts(std::string filename);
	void fitExponentialModel(int numNRIterations, double eps, std::string & filename, TLog* logfile);
};


//---------------------------------------------------------------
//TPMDFunction
//---------------------------------------------------------------
//Note: Base class is to be used when there is no PMD!
class TPMDFunction{
protected:
	std::string functionName;

	virtual void setName(){ functionName = "noPMD"; };
public:
	TPMDFunction(){ setName(); };
	TPMDFunction(TPMDFunction & other){
		setName();
	};
	virtual ~TPMDFunction(){};
	virtual void getCopy(TPMDFunction* & pointer){
		pointer = new TPMDFunction();
	};
	virtual double getProb(const uint16_t & pos) const{
		return 0.0;
	};
	virtual std::string getString() const{ return "P(pmd|pos) = 0.0"; };
	virtual std::string getFunctionName() const{ return functionName; };
	virtual bool hasDamage() const{ return false; };
};

class TPMDSkoglund:public TPMDFunction{
private:
	double lambda, c;

protected:
	void setName(){ functionName = "Skoglund"; };

public:
	TPMDSkoglund(double & Lambda, double & C);
	TPMDSkoglund(TPMDSkoglund & other){
		setName();
		lambda = other.lambda;
		c = other.lambda;
	};
	~TPMDSkoglund(){};
	void getCopy(TPMDFunction* & pointer){
		pointer = new TPMDSkoglund(lambda, c);
	};
	double getProb(const uint16_t & pos) const;
	std::string getString() const;
	bool hasDamage() const{ return true; };
};

class TPMDExponential:public TPMDFunction{
private:
	double a,b,c;

protected:
	void setName(){ functionName = "Exponential"; };

public:
	TPMDExponential(double & A, double & B, double & C);
	TPMDExponential(TPMDExponential & other){
		setName();
		a = other.a;
		b = other.b;
		c = other.c;
	};
	~TPMDExponential(){};
	void getCopy(TPMDFunction* & pointer){
		pointer = new TPMDExponential(a, b, c);
	};
	double getProb(const uint16_t & pos) const;
	std::string getString() const;
	bool hasDamage() const{ return true; };
};

class TPMDEmpiric:public TPMDFunction{
private:
	int length;
	std::vector<double> probs;
	double last;

protected:
	void setName(){ functionName = "Empiric"; };

public:
	TPMDEmpiric(std::string & values, std::string & example);
	TPMDEmpiric(std::vector<double> Probs);
	~TPMDEmpiric(){};
	void getCopy(TPMDFunction* & pointer){
		pointer = new TPMDEmpiric(probs);
	};
	double getProb(const uint16_t & pos) const;
	std::string getString() const;
	bool hasDamage() const{ return true; };
};

//------------------------------------------------------
//TPMDDoubleStrand
//------------------------------------------------------
class TPMDDoubleStrand{
private:
	TPMDFunction* myFunctions[2];
	bool functionsInitialized[2];

public:
	TPMDDoubleStrand();
	TPMDDoubleStrand(TParameters & params, TLog* logfile);
	TPMDDoubleStrand(const TPMDDoubleStrand & other);
	TPMDDoubleStrand(TPMDDoubleStrand && other);

	~TPMDDoubleStrand();

	void initialize(TParameters & params, TLog* logfile);
	void initialize(const TPMDDoubleStrand & other);
	void initializeFunction(std::string pmdString, PMDType type);

	//for getProb: distance is zero based!!! TODO: remove these functions
	double getProb(const uint16_t pos, PMDType type) const{ return myFunctions[type]->getProb(pos); };
	double getProbFivePrime(const uint16_t pos) const{ return myFunctions[pmdCT]->getProb(pos); };
	double getProbThreePrime(const uint16_t pos) const{ return myFunctions[pmdGA]->getProb(pos); };

	std::string getFunctionString(PMDType type) const{ return myFunctions[type]->getString(); };
	bool functionInitialized(PMDType type) const{
		return functionsInitialized[type];
	};
	bool hasDamage() const{ return myFunctions[pmdCT]->hasDamage() || myFunctions[pmdGA]->hasDamage(); };
	bool hasDamageCT() const{ return myFunctions[pmdCT]->hasDamage(); };
	bool hasDamageGA() const{ return myFunctions[pmdGA]->hasDamage(); };

	void fillBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const;

//	double getProbPMD(int readGroup, Base & ref, Base & read, double & pmdCT, double & pmdGA, double & errorRate);
//	double getProbNoPMD(int readGroup, Base & ref, Base & read, double & pmdCT, double & pmdGA, double & errorRate);
};

//------------------------------------------------------
//TPostMortemDamage
//------------------------------------------------------
class TPostMortemDamage{
private:
	std::vector<TPMDDoubleStrand> _pmdObjects;
	bool _hasPMD;

	PMDType getEnumPMDType(std::string pmdType);
	void initializeFromFile(BAM::TReadGroups & ReadGroups, const std::string filename, TLog* logfile);

public:
	TPostMortemDamage();
	bool hasPMD() const{ return _hasPMD; };
	void initialize(TParameters & params, BAM::TReadGroups & ReadGroups, TLog* logfile);
	void calculateBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const;
};

}; //end namespace


#endif /* TPOSTMORTEMDAMAGE_H_ */
