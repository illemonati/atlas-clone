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

public:
	TPMDTable(int MaxLength);
	~TPMDTable();
	void empty();
	void add(int & pos, Base & ref, Base & read);
	void writeTable(std::ofstream & out, std::string prefix);
	void writeTableWithCounts(std::ofstream & out, std::string prefix);
	std::string getPMDStringCT();
	std::string getPMDStringGA();
	long*** getPointerToCounts();
	long** getPointerToSums();
	int getMaxLength(){return maxLength;};
};

//---------------------------------------------------------------
//TPMDFitModel
//---------------------------------------------------------------
class TPMDFitModel{
protected:
	TLog* logfile;
	int numNRIterations;
	double eps;
	Base from, to;
	int numParams;

	//data
	std::string readGroupName;
	long*** counts;
	long** sums;
	int maxLength;
	int lastPositionToConsiderPlusOne;

	void findLastPositionToConsider();
	virtual void fitWithOLS(double* params){ throw "Not implemente din base class TPMDFitModel!"; };
	virtual void fillFAndJacobian(arma::vec & F, arma::mat & J, double* oldParams){ throw "Not implemente din base class TPMDFitModel!"; };
	virtual double calcLL(double* oldParams){ throw "Not implemente din base class TPMDFitModel!"; };
	virtual void runNewtonRaphson(double* params){ throw "Not implemente din base class TPMDFitModel!"; };
	virtual std::string assembleModelString(double* params){ throw "Not implemente din base class TPMDFitModel!"; };

public:
	TPMDFitModel(int NumNRIterations, double Eps, TLog* Logfile);
	virtual ~TPMDFitModel(){
		clear();
	};
	void clear();
	void setData(TPMDTable* table, std::string ReadGroupName);
	std::string fit(Base From, Base To);
};

class TPMDFitExponentialModel:public TPMDFitModel{
protected:
	void fitWithOLS(double* params);
	void fillFAndJacobian(arma::vec & F, arma::mat & J, double* oldParams);
	double calcLL(double* oldParams);
	void runNewtonRaphson(double* params);
	std::string assembleModelString(double* params);

public:
	TPMDFitExponentialModel(int NumNRIterations, double Eps, TLog* Logfile);
};

//---------------------------------------------------------------
//TPMDTables
//---------------------------------------------------------------
class TPMDTables{
private:
	void fitModel(TPMDFitModel & model, std::string & filename);

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
};

class TPMDSkoglund:public TPMDFunction{
private:
	double lambda, c;

protected:
	virtual void setName(){ functionName = "Skoglund"; };

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

	TPMD(TParameters & params, TLog* logfile){
		myFunctions[pmdCT] = NULL;
		myFunctions[pmdGA] = NULL;
		functionsInitialized[pmdCT] = false;
		functionsInitialized[pmdGA] = false;

		initialize(params, logfile);
	};
	TPMD(TPMD & other){initialize(other);};
	~TPMD(){
		if(functionsInitialized[pmdCT]) delete myFunctions[pmdCT];
		if(functionsInitialized[pmdGA]) delete myFunctions[pmdGA];
	};
	void initialize(TParameters & params, TLog* logfile);
	void initialize(TPMD & other);
	void initializeFunction(std::string & pmdString, PMDType type, TLog* logfile);
	//for getProb: distance is zero based!!!
	double getProb(int pos, PMDType type){ return myFunctions[type]->getProb(pos); };
	double getProbCT(int pos){ return myFunctions[pmdCT]->getProb(pos); };
	double getProbGA(int pos){ return myFunctions[pmdGA]->getProb(pos); };
	std::string getFunctionString(PMDType type){ return myFunctions[type]->getString(); };
	bool functionInitialized(PMDType type){
		return functionsInitialized[type];
	};
//	double getProbPMD(int readGroup, Base & ref, Base & read, double & pmdCT, double & pmdGA, double & errorRate);
//	double getProbNoPMD(int readGroup, Base & ref, Base & read, double & pmdCT, double & pmdGA, double & errorRate);
};



#endif /* TPOSTMORTEMDAMAGE_H_ */
