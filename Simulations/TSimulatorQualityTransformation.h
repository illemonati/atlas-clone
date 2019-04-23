/*
 * TSimulatorQualityTransformation.h
 *
 *  Created on: Nov 28, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATORQUALITYTRANSFORMATION_H_
#define TSIMULATORQUALITYTRANSFORMATION_H_

#include "../TBase.h"
#include "../TLog.h"
#include "../TRandomGenerator.h"
#include "../TQualityMap.h"
#include "../TParameters.h"
#include "TSimulatorReadLength.h"
#include "TRecalibrationEMModel.h"

//-------------------------------
//TSimulatorQualityDist
//-------------------------------
//Class of a fixed quality
class TSimulatorQualityDist{
protected:
	int _min;
	int _max;
	int _maxPlusOne;
	double _mean, _sd;


	//tmp variables
	int tmpInt;

public:
	TSimulatorQualityDist();
	TSimulatorQualityDist(std::string & s);
	virtual ~TSimulatorQualityDist(){};
	int min(){ return _min; };
	int max(){ return _max; };
	double mean(){return _mean; };
	double sd(){return _sd; };
	virtual int sample(){ return _max; };
	virtual void sample(int* qualities, const int & len);
	virtual void printDetails(TLog* logfile);
};

//Class of a binned distribution
class TSimulatorQualityDistBinned:public TSimulatorQualityDist{
private:
	TRandomGenerator* randomGenerator;
	std::vector<int> qualBins;
	int numQualBins;

public:
	TSimulatorQualityDistBinned(std::string & s, TRandomGenerator* RandomGenerator);
	void sample(int* qualities, const int & len);
	void printDetails(TLog* logfile);
};


//Class of a normal distribution
class TSimulatorQualityDistNormal:public TSimulatorQualityDist{
public:
//private:
	//densities
	int size;
	double* densities;
	double* cumulDensities;
	bool densitiesInitialized;

	TRandomGenerator* randomGenerator;

//public:
	TSimulatorQualityDistNormal(std::string & s, TRandomGenerator* RandomGenerator);
	TSimulatorQualityDistNormal(double & mean, double & sd, int min, int max, TRandomGenerator* RandomGenerator);;
	~TSimulatorQualityDistNormal(){
		if(densitiesInitialized){
			delete[] densities;
			delete[] cumulDensities;
		}
	};
	void parseFunctionString(std::string & s);
	void fillDensities();
	int sample();
	void sample(int* qualities, const int & len);
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
	virtual void simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len, const bool isReverseStrand);
	virtual void printDetails(TLog* logfile);
};


//-------------------------------------
//TSimulatorQualityTransformationRecal
//-------------------------------------
class TSimulatorQualityTransformationRecal:public TSimulatorQualityTransformation{
private:
	TRecalibrationEMModel_Base* model;
	int maxReadLengthPlusOne;
	int maxQualPlusOne;
	int numContext;
	int*** transformedQuality; //index are [qual][pos][context]
	TGenotypeMap genoMap;

	//private functions
	void fillTransformationTable(const std::string & modelTag, std::vector<std::string> & values, int maxReadLength);
	void clearTransformationTable();
	void transformQualities(Base* bases, int* qualities, const int & len, bool & isReverseStrand);

public:
	TSimulatorQualityTransformationRecal(std::string string, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	TSimulatorQualityTransformationRecal(const std::string & modelTag, std::vector<std::string> & values, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	~TSimulatorQualityTransformationRecal(){
		clearTransformationTable();
	};

	int getTransformedQuality(int qual, int pos, int context){
		return(transformedQuality[qual][pos][context]);
	};

	void simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len, const bool isReverseStrand);
	void printDetails(TLog* logfile);
};


//------------------------------------
//TSimulatorQualityTransformationBQSR
//------------------------------------
class TSimulatorQualityTransformationBQSR:public TSimulatorQualityTransformation{
private:
	//quality params
	TSimulatorReadLength* readLengthDist;
	int phi1;
	double phi2;
	int maxReadLength;
	int minPhredInt, maxPhredInt, maxPhredIntPlusOne;
	double meanPhred, sdPhred;
	double trueQual;
	TSimulatorQualityDistNormal* fakePhredDist;
	double kappa, lambda;

	//position params
	double revIntercept;
	double intercept;
	double m;

	//optimization algorithm params
	double* w;
	bool weightsInitialized;
	std::vector< std::vector<double> > errorBetaQBetaP;
	std::vector< std::vector<double> > QBetaQBetaP;

	//quality functions
	void parseBQSRQualInput(TParameters & params);
	double returnTrueError(const int & trueQual);
	void setFakePhredDistribution(TLog* logfile);
	int sampleFakePhredInt();

	//position functions
	void calculateSlopeIntercept();
	double returnBetaPp(const int & pos);

	//optimization algorithm functions
	void fillWeights(double & kappa_cur, double & lambda_cur);
	void fillQBetaQBetaP();
	double returnCurMean();
	double returnCurSD(double & kappa);
	double returnDelta(double & curMean, double & curSD);

public:
	TSimulatorQualityTransformationBQSR(const std::string & s, TSimulatorReadLength* ReadLengthDist, TLog* logfile, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorQualityTransformationBQSR(){
		delete fakePhredDist;
		if(weightsInitialized)
			delete[] w;
	};
	void simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len, const bool isReverseStrand);
};


//---------------------------------------------------------
//TSimulatorQualityTransformParameters
//---------------------------------------------------------
struct TSimulatorQualityTransformParameters{
	std::string type; //either none, recal or BQSR
	std::string parameters_firstMate;
	std::string parameters_secondMate;

	TSimulatorQualityTransformParameters(){};
	TSimulatorQualityTransformParameters(std::string Type, std::string Parameters_firstMate, std::string Parameters_secondMate){
		type = Type;
		parameters_firstMate = Parameters_firstMate;
		parameters_secondMate = Parameters_secondMate;
	};

};


#endif /* TSIMULATORQUALITYTRANSFORMATION_H_ */
