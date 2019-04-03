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


#define noRecal_name "none"
#define qualFuncPosFunc_name "qualFuncPosFunc"
#define qualFuncPosFuncContext_name "qualFuncPosFuncContext"
#define qualFuncPosSpecific_name "qualFuncPosSpecific"
#define qualFuncPosSpecificContext_name "qualFuncPosSpecificContext"
#define qualFuncPosSpecificContextNew_name "qualFuncPosSpecificContextNew"

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

	void _parseParameterString(std::vector<std::string> & vec, std::vector<double>* values);
	void parseQualityParameters(double* betaPointer, std::vector<double> & values); //default function assuming quadratic model
	void parsePositionParameters(double* betaPointer, std::vector<double> & values); //default function assuming quadratic model
	void _allocateBetaMemory();
	void _freeBetaMemory();
	double _calcEpsilon(double & eta);

public:

	TRecalibrationEMModel_Base();
	TRecalibrationEMModel_Base(int Shift);
	virtual ~TRecalibrationEMModel_Base(){ _freeBetaMemory(); };

	std::string name(){ return _name; };
	int numParameters(){ return _numParameters; };
	int shift(){ return _myShift; };
	void setShift(int shift){ _myShift = shift;};
	virtual bool checkParameterRange(int maxPos){ return true; }; //check if parameters are in correct range
	virtual void proposeNewParameters(double & lambda, arma::mat & JxF);
	void rejectProposedParameters();
	virtual double calcEpsilon(const TRecalibrationEMReadData & data){ throw "double calcEpsilon(TRecalibrationEMReadData & data) not defined for TRecalibrationEMModel_Base!"; };
	virtual void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){ throw "void addToFandJacobian(...) not defined for TRecalibrationEMModel_Base!"; };
	void writeParametersToFile(TOutputFilePlain & out);
	std::string getModelString();
	virtual std::string getQualityString(); //default function assuming quadratic model
	virtual std::string getPositionString(); //default function assuming quadratic model
	virtual std::string getContextString(); //default function assuming context specific intercepts
	virtual double getErrorRate(TBase & base){ throw "double getErrorRate(TBase & base) not defined for TRecalibrationEMModel_Base!"; };
	virtual void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual){ throw "void fillTransformationTableForSimulation(...) not defined for TRecalibrationEMModel_Base!"; };
};

class TRecalibrationEMModel_noRecal:public TRecalibrationEMModel_Base{
public:
	TRecalibrationEMModel_noRecal(int Shift);
	~TRecalibrationEMModel_noRecal(){};
	double getErrorRate(TBase & base);
	std::string getQualityString(){return "-";};
	std::string getPositionString(){return "-";};
	std::string getContextString(){return "-";};
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosFunc:public TRecalibrationEMModel_Base{
public:
	TRecalibrationEMModel_qualFuncPosFunc(int Shift);
	TRecalibrationEMModel_qualFuncPosFunc(std::vector<std::string> & vec, int Shift);
	~TRecalibrationEMModel_qualFuncPosFunc(){};

	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	std::string getContextString();
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosFuncContext:public TRecalibrationEMModel_Base{
public:
	TRecalibrationEMModel_qualFuncPosFuncContext(int Shift);
	TRecalibrationEMModel_qualFuncPosFuncContext(std::vector<std::string> & vec, int Shift);
	~TRecalibrationEMModel_qualFuncPosFuncContext(){};

	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosSpecific:public TRecalibrationEMModel_Base{
private:
	int _maxPosPlusOne;
	int _numParamsWithoutPositions;

public:
	TRecalibrationEMModel_qualFuncPosSpecific(int Shift, int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecific(std::vector<std::string> & vec, int Shift);
	~TRecalibrationEMModel_qualFuncPosSpecific(){};

	bool checkParameterRange(int maxPos){ if(_maxPosPlusOne == maxPos + 1) return true; else return false; }; //check if parameters are in correct range
	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	std::string getPositionString();
	std::string getContextString();
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosSpecificContext:public TRecalibrationEMModel_Base{
private:
	int _maxPosPlusOne;
	int _numParamsWithoutPositions;

public:
	TRecalibrationEMModel_qualFuncPosSpecificContext(int Shift, int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecificContext(std::vector<std::string> & vec, int Shift);
	~TRecalibrationEMModel_qualFuncPosSpecificContext(){};

	bool checkParameterRange(int maxPos){ if(_maxPosPlusOne == maxPos + 1) return true; else return false; }; //check if parameters are in correct range
	void proposeNewParameters(double & lambda, arma::mat & JxF);
	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	std::string getPositionString();
	std::string getContextString();
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosSpecificContextNew:public TRecalibrationEMModel_Base{
private:
	int _maxPosMinusOne;
	int _maxPosPlusOne;
	int _numParamsWithoutPositions;

public:
	TRecalibrationEMModel_qualFuncPosSpecificContextNew(int Shift, int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecificContextNew(std::vector<std::string> & vec, int Shift);
	~TRecalibrationEMModel_qualFuncPosSpecificContextNew(){};

	bool checkParameterRange(int maxPos){ if(_maxPosMinusOne == maxPos + 1) return true; else return false; }; //check if parameters are in correct range
	void proposeNewParameters(double & lambda, arma::mat & JxF);
	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	std::string getPositionString();
	std::string getContextString();
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

//--------------------------------------------------------------------
// Global function to create models
//--------------------------------------------------------------------
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, std::vector<std::string> & values, int shift, bool verbose, TLog* logfile);
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, int maxPos, int shift, bool verbose, TLog* logfile);

//--------------------------------------------------------------------
// TRecalibrationEMModels
// Object containing a vector of recal models
//--------------------------------------------------------------------
class TRecalibrationEMModels{
private:
	std::vector<TRecalibrationEMModel_Base*> models;
	unsigned int totNumParameters;
	TRecalibrationEMReadGroupIndex readGroupIndex;
	TLog* logfile;

	//Newton Raphson Parameters
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;

	void _addModel(std::string & modelTag, std::vector<std::string> & values, bool verbose);
	void _createModelsFromString(std::string & string, TReadGroups & readGroups);
	void _createModelsFromFile(std::string filename, TReadGroups & readGroups, TReadGroupMap & readGroupMap);

	void _writeParameters(TOutputFilePlain & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate);

public:
	TRecalibrationEMModels(int numReadGroups, TLog* Logfile);
	~TRecalibrationEMModels();

	void addSameModelForAllReadGroups(std::string modelTag, std::vector<std::string> & values, bool verbose);
	void addModel(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<std::string> & values, bool verbose);
	void addModel(int readGroupId, bool isSecondMate, std::string modelTag, int maxPos);
	void addModelIfItDoesNotExist(int readGroupId, bool isSecondMate, std::string modelTag, int maxPos);
	void removeModel(int readGroupId, bool isSecondMate);
	void createModels(std::string string, TReadGroups & readGroups, TReadGroupMap & readGroupMap);

	inline TRecalibrationEMModel_Base* operator[](int index){ return models[index]; };
	int numModels(){ return models.size(); };
	bool modelExists(int readGroupId, bool isSecondMate){ return readGroupIndex.inUse(readGroupId, isSecondMate); };
	inline double calcEpsilon(const TRecalibrationEMReadData & data){
		return models[ readGroupIndex.index(data) ]->calcEpsilon(data);
	};

	bool hasReadGroupsWithoutModel();
	void addNoRecalModelIfMissing();
	void reportReadGroupsNotUsed(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void reportReadGroupsNotUsed(TReadGroups & readGroups);
	void reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups);
	void warningForMissingReadGroups(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void warningForMissingReadGroups(TReadGroups & readGroups);

	void setEMParamsToZero();
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weight, const double & weightJacobian);
	bool solveJxF(const int numSites);
	void proposeNewParameters(double lambda);
	void scaleParameters();
	void rejectProposedParameters();
	double getSteepestGradient();
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();
	void writeHeader(TOutputFilePlain & out);
	void writeParameters(TOutputFilePlain & out, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	inline double getErrorRate(TBase & base){
		return models[ readGroupIndex.index(base.readGroup, base.isSecondMate) ]->getErrorRate(base);
	};
};


#endif /* TRECALIBRATIONEMMODEL_H_ */
