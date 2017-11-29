/*
 * TSimulatorQualityTransformation.h
 *
 *  Created on: Nov 28, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATORQUALITYTRANSFORMATION_H_
#define TSIMULATORQUALITYTRANSFORMATION_H_

#include "TBase.h"
#include "TLog.h"
#include "TRandomGenerator.h"
#include "TGenotypeMap.h"
#include "TParameters.h"

//-------------------------------
//TSimulatorQualityDist
//-------------------------------
//Class of a fixed quality
class TSimulatorQualityDist{
protected:
	int _max;

	//tmp variables
	int tmpInt;

public:
	TSimulatorQualityDist();
	TSimulatorQualityDist(std::string & s);
	virtual ~TSimulatorQualityDist(){};
	int max(){ return _max; };
	virtual int sample(){ return _max; };
	virtual void sample(int* qualities, int & len);
	virtual void printDetails(TLog* logfile);
};

//Class of a normal distribution
class TSimulatorQualityDistNormal:public TSimulatorQualityDist{
private:
	double mean, sd;
	int _max, _maxPlusOne;
	int _min;

	//densities
	int size;
	double* densities;
	double* cumulDensities;
	bool densitiesInitialized;

	TRandomGenerator* randomGenerator;

public:
	TSimulatorQualityDistNormal(std::string & s, TRandomGenerator* RandomGenerator);
	~TSimulatorQualityDistNormal(){
		if(densitiesInitialized){
			delete[] densities;
			delete[] cumulDensities;
		}
	};
	void parseFunctionString(std::string & s);
	void fillDensities();
	int sample();
	void sample(int* qualities, int & len);
	void printDetails(TLog* logfile);
};

//-----------------------------------------------
//TSimulatorQualityTransformation
//-----------------------------------------------
class TSimulatorQualityTransformation{
protected:
	TSimulatorQualityDist* qualityDist;
	TRandomGenerator* randomGenerator;
	TQualityMap qualityMap;

	//tmp vars
	int p;

public:
	TSimulatorQualityTransformation(TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorQualityTransformation(){};
	virtual void simulateQualitiesAndErrors(Base* bases, int* qualities, int & len);
	virtual void printDetails(TLog* logfile);
};


//-------------------------------------
//TSimulatorQualityTransformationRecal
//-------------------------------------
class TSimulatorQualityTransformationRecal:public TSimulatorQualityTransformation{
private:
	std::vector<double> betas;
	int maxReadLengthPlusOne;
	int maxQualPlusOne;
	int numContextPlusOne;
	int*** transformedQuality; //index are [qual][pos][context]
	Base previousBase;
	TGenotypeMap genoMap;

	//private functions
	void fillTransformationTable(int maxReadLength);
	void clearTransformationTable();
	void simulateQualitiesAndErrors(Base* bases, int* qualities, int & len);
	inline int transformQual(const int & qual, const int & pos, const int & context);

public:
	TSimulatorQualityTransformationRecal(const std::string & s, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	~TSimulatorQualityTransformationRecal(){
		clearTransformationTable();
	};
	void printDetails(TLog* logfile);
};

/*
//------------------------------------
//TSimulatorQualityTransformationBQSR
//------------------------------------
class TSimulatorQualityTransformationBQSR:public TSimulatorQualityTransformation{
private:
	int maxPos;

	//quality params
	int phi1;
	double phi2;
	bool trueQualToFakeQualTableInitialized = false;
	double* trueQualToFakeQual;
	double kappa, lambda;

	//position params
	double revIntercept = 1.0;
	double intercept = 1.0;
	double m = 0.0;

	//optimization algorithm params
	double* w;
	bool weightsInitialized;
	std::vector< std::vector<double> > QBetaQBetaP;
	bool betaQBetaPInitialized = false;

	//quality functions
	void parseBQSRQualInput(TParameters & params);
	double returnFakeError(int & trueQual);
	void setFakeQualityDistribution();
	void initializeTrueQualToFakeQualTable();
	int sampleFakeQuality();

	//position functions
	void calculateSlopeIntercept();
	double returnBetaPp(int & pos);

	//optimization algorithm functions
	void fillWeights(double & kappa_cur, double & lambda_cur);
	void fillQBetaQBetaP();
	double returnCurMean();
	double returnCurSD(double & kappa);
	double returnDelta(double & curMean, double & curSD);

public:
	TSimulatorQualityTransformationBQSR(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorQualityTransformationBQSR(){
		if(trueQualToFakeQualTableInitialized)
			delete[] trueQualToFakeQual;
		if(weightsInitialized)
			delete[] w;
	};
};
*/


#endif /* TSIMULATORQUALITYTRANSFORMATION_H_ */
