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

enum PMDType {pmdCT=0, pmdGA};
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
	void add(int & pos, Base & ref, Base & read);
	void writeTable(std::ofstream & out, std::string prefix);
	void writeTableWithCounts(std::ofstream & out, std::string prefix);
	std::string getPMDStringCT();
	std::string getPMDStringGA();
	std::string fitExponentialModel(Base from, Base to, int & numNRIterations, double & eps);
};

class TPMDTables{
public:
	TReadGroups* readGroups;
	TPMDTable** forward;
	TPMDTable** reverse;

	TPMDTables(TReadGroups* ReadGroups, int maxLength);
	~TPMDTables();
	void addForward(int readGroup, int pos, Base & ref, Base & read);
	void addReverse(int readGroup, int pos, Base & ref, Base & read);
	void writePMDFile(std::string filename);
	void writeTable(std::string filename);
	void writeTableWithCounts(std::string filename);
	void fitExponentialModel(int numNRIterations, double eps, std::string & filename);
};


//---------------------------------------------------------------
//TPMDFunction
//---------------------------------------------------------------
//Note: Base class is to be used when there is no PMD!
class TPMDFunction{
public:
	TPMDFunction(){};
	virtual ~TPMDFunction(){};
	virtual double getProb(int & pos){
		return 0.0;
	};
	virtual std::string getString(){ return "P(pmd|pos) = 0.0"; };
};

class TPMDSkoglund:public TPMDFunction{
private:
	double lambda, c;

public:
	TPMDSkoglund(double & Lambda, double & C);
	~TPMDSkoglund(){};
	double getProb(int & pos);
	std::string getString();
};

class TPMDExponential:public TPMDFunction{
private:
	double a,b,c;

public:
	TPMDExponential(double & A, double & B, double & C);
	double getProb(int & pos);
	std::string getString();
};

class TPMDEmpiric:public TPMDFunction{
private:
	int length;
	double* probs;
	double last;

public:
	TPMDEmpiric(std::string & values, std::string & example);
	double getProb(int & pos);
	std::string getString();
};

//------------------------------------------------------
//TPMD
//------------------------------------------------------
class TPMD{
private:
	TPMDFunction* myFunctions[2];
	bool functionsInitialized[2];

public:
	TPMD(){
		myFunctions[pmdCT] = NULL;
		myFunctions[pmdGA] = NULL;
		functionsInitialized[pmdCT] = false;
		functionsInitialized[pmdGA] = false;
	};

	~TPMD(){
		if(functionsInitialized[pmdCT]) delete myFunctions[pmdCT];
		if(functionsInitialized[pmdGA]) delete myFunctions[pmdGA];
	};
	void initializeFunction(std::string & pmdString, PMDType type);
	//for getProb: distance is zero based!!!
	double getProb(int pos, PMDType type){ return myFunctions[type]->getProb(pos); };
	double getProbCT(int pos){ return myFunctions[pmdCT]->getProb(pos); };
	double getProbGA(int pos){ return myFunctions[pmdGA]->getProb(pos); };
	std::string getFunctionString(PMDType type){ return myFunctions[type]->getString(); };
	bool functionInitialized(PMDType type){
		return functionsInitialized[type];
	};
};



#endif /* TPOSTMORTEMDAMAGE_H_ */
