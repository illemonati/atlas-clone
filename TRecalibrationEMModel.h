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
#define qualSpecificPosSpecific_name "qualSpecificPosSpecific"
#define qualSpecificPosSpecificContext_name "qualSpecificPosSpecificContext"


//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel_Base{
protected:
	std::string _name;
	unsigned int _numParameters;
	TRecalibrationEMQualityPositionMap _qualPosMap;

	double* _betas; //betas of the model
	double* _oldBetas; //use during estimation
	bool _initialized;

	//Newton Raphson Parameters
	double _Q, _oldQ;
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	unsigned int _numSitesAdded;
	bool _NRconverged;
	bool _NRStepAccepted;

	void _parseParameterString(std::vector<std::string> & vec, std::vector<double>* values);
	void parseQualityParameters(double* betaPointer, std::vector<double> & values); //default function assuming quadratic model
	void parsePositionParameters(double* betaPointer, std::vector<double> & values); //default function assuming quadratic model
	void _allocateBetaMemory();
	void _freeBetaMemory();
	double _calcEpsilon(double & eta);

public:

	TRecalibrationEMModel_Base();
	virtual ~TRecalibrationEMModel_Base(){ _freeBetaMemory(); };

	std::string name(){ return _name; };
	int numParameters(){ return _numParameters; };
	void setEMParamsToZero();
	virtual void checkParameterRange(std::vector<int> & Qualities, int maxPos){}; //check if parameters are in correct range

	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	double curQ(){ return _Q; };
	bool solveJxF();
	virtual void proposeNewParameters(double & lambda);
	bool acceptProposedParametersBasedOnQ();
	void rejectProposedParameters();

	double getSteepestGradient();
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();

	virtual double calcEpsilon(const TRecalibrationEMReadData & data){ throw "double calcEpsilon(TRecalibrationEMReadData & data) not defined for TRecalibrationEMModel_Base!"; };
	virtual void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){ throw "void addToFandJacobian(...) not defined for TRecalibrationEMModel_Base!"; };
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
	TRecalibrationEMModel_noRecal();
	~TRecalibrationEMModel_noRecal(){};
	double getErrorRate(TBase & base);
	std::string getQualityString(){return "-";};
	std::string getPositionString(){return "-";};
	std::string getContextString(){return "-";};
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosFunc:public TRecalibrationEMModel_Base{
public:
	TRecalibrationEMModel_qualFuncPosFunc();
	TRecalibrationEMModel_qualFuncPosFunc(std::vector<std::string> & vec);
	~TRecalibrationEMModel_qualFuncPosFunc(){};

	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	std::string getContextString();
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosFuncContext:public TRecalibrationEMModel_Base{
public:
	TRecalibrationEMModel_qualFuncPosFuncContext();
	TRecalibrationEMModel_qualFuncPosFuncContext(std::vector<std::string> & vec);
	~TRecalibrationEMModel_qualFuncPosFuncContext(){};

	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualFuncPosSpecific:public TRecalibrationEMModel_Base{
private:
	int _maxPosPlusOne;
	int _numParamsWithoutPositions;

public:
	TRecalibrationEMModel_qualFuncPosSpecific(int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecific(std::vector<std::string> & vec);
	~TRecalibrationEMModel_qualFuncPosSpecific(){};

	void checkParameterRange(std::vector<int> & Qualities, int maxPos);
	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
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
	TRecalibrationEMModel_qualFuncPosSpecificContext(int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecificContext(std::vector<std::string> & vec);
	~TRecalibrationEMModel_qualFuncPosSpecificContext(){};

	void checkParameterRange(std::vector<int> & Qualities, int maxPos);
	void proposeNewParameters(double & lambda);
	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
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
	TRecalibrationEMModel_qualFuncPosSpecificContextNew(int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecificContextNew(std::vector<std::string> & vec);
	~TRecalibrationEMModel_qualFuncPosSpecificContextNew(){};

	void checkParameterRange(std::vector<int> & Qualities, int maxPos);
	void proposeNewParameters(double & lambda);
	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	std::string getPositionString();
	std::string getContextString();
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

class TRecalibrationEMModel_qualSpecficPosSpecific:public TRecalibrationEMModel_Base{
private:
	int _maxPosPlusOne;
	int _numQualities;
	int _maxQualPlusOne;
	int* _qualityIndex;
	TQualityMap qualMap;

public:
	TRecalibrationEMModel_qualSpecficPosSpecific(std::vector<int> & Qualities, int MaxQual, int MaxPos);
	TRecalibrationEMModel_qualSpecficPosSpecific(std::vector<std::string> & vec);
	~TRecalibrationEMModel_qualSpecficPosSpecific(){};

	void checkParameterRange(std::vector<int> & Qualities, int maxPos);
	double calcEpsilon(const TRecalibrationEMReadData & data);
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	std::string getQualityString();
	std::string getPositionString();
	std::string getContextString();
	double getErrorRate(TBase & base);
	void fillTransformationTableForSimulation(int*** transformedQuality, int MaxPos, int MaxQual);
};

//--------------------------------------------------------------------
// Global function to create models
//--------------------------------------------------------------------
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, std::vector<std::string> & values, bool verbose, TLog* logfile);
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, int maxPos, bool verbose, TLog* logfile);

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
	void addModelIfItDoesNotExist(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<int> & Qualities, int maxPos);
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
	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	double curQ();
	bool solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void rejectProposedParameters();

	double getSteepestGradient();
	void writeHeader(TOutputFilePlain & out);
	void writeParameters(TOutputFilePlain & out, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	inline double getErrorRate(TBase & base){
		return models[ readGroupIndex.index(base.readGroup, base.isSecondMate) ]->getErrorRate(base);
	};
};


#endif /* TRECALIBRATIONEMMODEL_H_ */
