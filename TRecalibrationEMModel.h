/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include "TRecalibrationEMAuxiliaryTools.h"
#include "TFile.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

#define qualfuncPosFuncContext_name "qualFuncPosFuncContext"
#define qualfuncPosFunc_name "qualFuncPosFunc"
#define qualfuncPosSpecificContext_name "qualFuncPosSpecificContext"
#define noRecal_name "none"


//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel_Base{
protected:
	std::string _name;
	unsigned int _numParameters;
	int _myShift;
	TRecalibrationEMQualityPositionMap _qualPosMap;
	long _numSitesAdded;

	double* _betas; //betas of the model
	double* _oldBetas; //use during estimation
	bool _initialized;

	void _allocateBetaMemory();
	void _freeBetaMemory();
	double _calcEpsilon(double & eta);

public:

	TRecalibrationEMModel_Base();
	TRecalibrationEMModel_Base(int Shift);
	virtual ~TRecalibrationEMModel_Base(){ _freeBetaMemory(); };

	std::string name(){ return _name; };
	int numParameters(){ return _numParameters; };
	void proposeNewParameters(double & lambda, arma::mat & JxF);
	void rejectProposedParameters();
	virtual double calcEpsilon(const TRecalibrationEMReadData & data){ throw "double calcEpsilon(TRecalibrationEMReadData & data) not defined for TRecalibrationEMModel_Base!"; };
	virtual void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){ throw "void addToFandJacobian(...) not defined for TRecalibrationEMModel_Base!"; };
	virtual void writeParametersToFile(TOutputFilePlain & out){ throw "void writeParametersToFile(TOutputFilePlain) not defined for TRecalibrationEMModel_Base!"; };
	virtual double getErrorRate(TBase & base){ throw "double getErrorRate(TBase & base) not defined for TRecalibrationEMModel_Base!"; };
};

class TRecalibrationEMModel_noRecal:public TRecalibrationEMModel_Base{
public:
	TRecalibrationEMModel_noRecal(int Shift);
	~TRecalibrationEMModel_noRecal(){};
	virtual double getErrorRate(TBase & base);
	void writeParametersToFile(TOutputFilePlain & out);
};

class TRecalibrationEMModel_qualFuncPosFuncContext:public TRecalibrationEMModel_Base{
protected:
	void _initialize(int NumReadGroups);

public:
	TRecalibrationEMModel_qualFuncPosFuncContext(int Shift);
	TRecalibrationEMModel_qualFuncPosFuncContext(std::vector<std::string> & vec, int Shift);
	virtual ~TRecalibrationEMModel_qualFuncPosFuncContext(){};

	virtual double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	virtual void writeParametersToFile(TOutputFilePlain & out);
	virtual double getErrorRate(TBase & base);
};

class TRecalibrationEMModel_qualFuncPosFunc:public TRecalibrationEMModel_qualFuncPosFuncContext{
public:
	TRecalibrationEMModel_qualFuncPosFunc(int Shift);
	TRecalibrationEMModel_qualFuncPosFunc(std::vector<std::string> & vec, int Shift);
	~TRecalibrationEMModel_qualFuncPosFunc(){};

	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	void writeParametersToFile(TOutputFilePlain & out);
	double getErrorRate(TBase & base);
};

class TRecalibrationEMModel_qualFuncPosSpecificContext:public TRecalibrationEMModel_Base{
private:
	int maxPos;
	int numParamsWithoutPositions;

public:
	TRecalibrationEMModel_qualFuncPosSpecificContext(int Shift, int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecificContext(std::vector<std::string> & vec, int Shift);
	~TRecalibrationEMModel_qualFuncPosSpecificContext(){};

	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	void writeParametersToFile(TOutputFilePlain & out);
	double getErrorRate(TBase & base);
};


//--------------------------------------------------------------------
// TRecalibrationEMModels
// Object containing a vector of recal models
//--------------------------------------------------------------------
class TRecalibrationEMModels{
private:
	std::vector<TRecalibrationEMModel_Base*> models;
	int totNumParameters;
	TRecalibrationEMReadGroupIndex readGroupIndex;
	TLog* logfile;

	//Newton Raphson Parameters
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;

	void _addModel(std::string & modelTag, std::vector<std::string> & values, bool verbose);
	void _writeParameters(TOutputFilePlain & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate);

public:
	TRecalibrationEMModels(int numReadGroups, TLog* Logfile);
	~TRecalibrationEMModels();

	void addSingleModelForAllReadGroups(std::string modelTag, std::vector<std::string> & values, bool verbose);
	void addModel(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<std::string> & values, bool verbose);
	void addModel(int readGroupId, bool isSecondMate, std::string modelTag, int maxPos);

	inline TRecalibrationEMModel_Base* operator[](int index){ return models[index]; };
	int numModels(){ return models.size(); };
	inline double calcEpsilon(const TRecalibrationEMReadData & data){
		return models[ readGroupIndex.index(data) ]->calcEpsilon(data);
	};

	bool hasReadGroupsWithoutModel();
	void addNoRecalModelIfMissing();
	void reportReadGroupsNotUsed(TReadGroups & readGroups);
	void setEMParamsToZero();
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weight, const double & weightJacobian);
	bool solveJxF(const int numSites);
	void proposeNewParameters(double lambda);
	void rejectProposedParameters();
	double getSteepestGradient();
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();
	void writeHeader(TOutputFilePlain & out);
	void writeParameters(TOutputFilePlain & out, TReadGroups & readGroups);
	inline double getErrorRate(TBase & base){
		return models[ readGroupIndex.index(base.readGroup, base.isSecondMate) ]->getErrorRate(base);
	};
};


#endif /* TRECALIBRATIONEMMODEL_H_ */
