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
#include <algorithm>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


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
//	std::string getPMDStringCT();
//	std::string getPMDStringGA();
	std::string fitExponentialModel(Base from, Base to, int & numNRIterations, double & eps, std::string readGroupName, int maxReadLength, TLog* logfile);
};

class TPMDTables{
public:
	TReadGroupMap& readGroupMapObject;
	TReadGroups& readGroups;
	int maxReadLength;
	int origNumReadGroups;
	int numReadGroups;
	TPMDTable** forward;
	TPMDTable** reverse;

	TPMDTables(TReadGroups& ReadGroups, int maxLengthForInference, int MaxReadLength, TReadGroupMap & ReadGroupMapObject);
	~TPMDTables();
//	void initializeReadGroupMap(BamTools::SamHeader* bamHeader, TParameters & params, TLog* logfile);
	void addFromFivePrime(const int readGroup, const int pos, const Base & ref, const Base & read);
	void addFromThreePrime(const int readGroup, const int pos, const Base & ref, const Base & read);
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
	virtual double getProb(int & pos){
		return 0.0;
	};
	virtual std::string getString(){ return "P(pmd|pos) = 0.0"; };
	virtual std::string getFunctionName(){ return functionName; };
	virtual bool hasDamage(){ return false; };
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
	double getProb(int & pos);
	std::string getString();
	bool hasDamage(){ return true; };
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
	double getProb(int & pos);
	std::string getString();
	bool hasDamage(){ return true; };
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
	double getProb(int & pos);
	std::string getString();
	bool hasDamage(){ return true; };
};

//------------------------------------------------------
//TPMD
//------------------------------------------------------
class TPMD{
private:
	TPMDFunction* myFunctions[4];
	bool functionsInitialized[4];

public:
	TPMD(){
		for(int pmdType=0; pmdType<4; ++pmdType){
			myFunctions[pmdType] = NULL;
			functionsInitialized[pmdType] = false;
		}
	};

	TPMD(TParameters & params, TLog* logfile){
		TPMD();
		initialize(params, logfile);
	};

	TPMD(TPMD & other){initialize(other);};
	~TPMD(){
		if(functionsInitialized[pmdCT]) delete myFunctions[pmdCT];
		if(functionsInitialized[pmdGA]) delete myFunctions[pmdGA];
	};
	void initialize(TParameters & params, TLog* logfile);
	void initialize(TPMD & other);
	void initializeFunction(std::string pmdString, PMDType type);
	//for getProb: distance is zero based!!!
	double getProb(int pos, PMDType type){ return myFunctions[type]->getProb(pos); };
	double getProbFivePrime(int pos){ return myFunctions[pmdCT]->getProb(pos); };
	double getProbThreePrime(int pos){ return myFunctions[pmdGA]->getProb(pos); };
	std::string getFunctionString(PMDType type){ return myFunctions[type]->getString(); };
	bool functionInitialized(PMDType type){
		return functionsInitialized[type];
	};
	bool hasDamage(){ return myFunctions[pmdCT]->hasDamage() & myFunctions[pmdGA]->hasDamage(); };
	bool hasDamageCT(){ return myFunctions[pmdCT]->hasDamage(); };
	bool hasDamageGA(){ return myFunctions[pmdGA]->hasDamage(); };


//	double getProbPMD(int readGroup, Base & ref, Base & read, double & pmdCT, double & pmdGA, double & errorRate);
//	double getProbNoPMD(int readGroup, Base & ref, Base & read, double & pmdCT, double & pmdGA, double & errorRate);
};



#endif /* TPOSTMORTEMDAMAGE_H_ */
