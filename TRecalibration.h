/*
 * TRecalibration.h
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#ifndef TRecalIBRATION_H_
#define TRecalIBRATION_H_

#include "bamtools/api/BamReader.h"
#include "TSite.h"
#include <omp.h>
#include "TReadGroups.h"
#include "bamtools/api/SamHeader.h"
#include "QualityTables.h"
#include "TQualityMap.h"

//---------------------------------------------------------------
//TQualityIndex
//---------------------------------------------------------------
class TQualityIndex{
	//Note: quality as stored in bases ranges from 33 to max!
public:
	int minPhredInt, maxPhredInt, numQ, last, first;
	int* index;

	TQualityIndex(int MinQ, int MaxQ){
		//here minQ and maxQ are not in ascii. by default minQ = 0 and maxQ = 100
		minPhredInt = MinQ;
		maxPhredInt = MaxQ;
		numQ = maxPhredInt - minPhredInt + 1;
		last = numQ - 1;
		first = 0;

		//fill index
		index = new int[numQ + 33];
		for(int i=0; i < maxPhredInt + 1; ++i){
			if(i < minPhredInt) index[i] = 0;
			else index[i] = i - minPhredInt;
		}
	};

	~TQualityIndex(){
		delete[] index;
	};

	int& getIndex(const int & quality){
		if(quality < 33) throw "Quality is negative!";
		if(quality > maxPhredInt + 33) return last;
		return index[quality - 33];
	};

	int getPhredIntFromIndex(const int & index){
		if(index < 0) throw "Quality index is negative!";
		if(index > numQ) return maxPhredInt;
		return minPhredInt + index;
	};
};

//---------------------------------------------------------------
//TRecalibration: default = no recalibration
//---------------------------------------------------------------
class TRecalibration{
protected:
	TReadGroupMap& readGroupMapObject;
	TQualityMap qualityMap;

public:
	TRecalibration(TReadGroupMap& ReadGroupMap);

	virtual ~TRecalibration(){};

	virtual bool recalibrationChangesQualities(){
		return false;
	};

/*	char getQualityAsChar(const TBase & base, int & minOutQuality, int & maxOutQuality){
		int qual = getphredInt(base) + 33;
		if(qual > maxOutQuality) qual = maxOutQuality;
		if(qual < minOutQuality) qual = minOutQuality;
		return qual;
	};*/

	void calcEmissionProbabilities(TSite & site);
	virtual double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	double getErrorRateFromBase(const TBase & base);
	virtual int getQuality(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	virtual int getQualityFromBase(const TBase & base);
//	virtual int getphredInt(const TBase & base){
//		return base.phredInt;
//	};
//	virtual int getQuality(const TBase & base){
//		return base.quality;
//	};

	virtual bool requiresEstimation(){ return false;};
	int findReadGroupIndex(std::string & name, BamTools::SamReadGroupDictionary & readGroups);
};

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
class TRecalibrationEMQualityPositionMap{
public:
	int maxQual;
	int maxPos;
	double* eta;
	double* etaSquared;
	double* position;
	double* positionSquared;
	bool initialized;

	TRecalibrationEMQualityPositionMap(){
		initialized = false;
		initialize(500, 255); //TODO: think about default!
	};

	~TRecalibrationEMQualityPositionMap(){
		clear();
	};

	void clear(){
		if(initialized){
			delete[] eta;
			delete[] etaSquared;
			delete[] position;
			delete[] positionSquared;
			initialized = false;
		}
	};

	void initialize(int MaxQual, int MaxPos){
		clear();
		maxQual = MaxQual;
		maxPos = MaxPos;
		eta = new double[maxQual+1];
		etaSquared = new double[maxQual+1];
		position = new double[maxPos+1];
		positionSquared = new double[maxPos+1];
		initialized = true;

		//fill qualities. Use TQualityMap for conversion
		TQualityMap qualiMap;
		for(int q=0; q<=maxQual; q++){
			double eps = qualiMap.qualityToError(q);
			if(eps < 0.0000000001) eps = 0.0000000001;
			else if(eps > 0.9999999999) eps = 0.9999999999;

			eta[q] = log(eps / (1.0 - eps));
			etaSquared[q] = eta[q] * eta[q];
		}

		//fill positions
		for(int p = 0; p<=maxPos; p++){
			position[p] = p;
			positionSquared[p] = p * p;
		}
	};
};

struct TRecalibrationEMReadData{
	uint8_t quality;
	uint8_t position;
	float D[4];
	uint8_t context;
	uint8_t readGroup;

	void setD(Base base, double PMD_CT, double PMD_GA){
		switch(base){
			case A: D[0] = 0.0; //geno = AA
					D[1] = 1.0; //geno = CC
					D[2] = 1.0 - PMD_GA; //geno = GG
					D[3] = 1.0; //geno = TT
					break;
			case C: D[0] = 1.0; //geno = AA
					D[1] = PMD_CT; //geno = CC
					D[2] = 1.0; //geno = GG
					D[3] = 1.0; //geno = TT
					break;
			case G: D[0] = 1.0; //geno = AA
					D[1] = 1.0; //geno = CC
					D[2] = PMD_GA; //geno = GG
					D[3] = 1.0; //geno = TT
					break;
			case T: D[0] = 1.0; //geno = AA
					D[1] = 1.0 - PMD_CT; //geno = CC
					D[2] = 1.0; //geno = GG
					D[3] = 0.0; //geno = TT
					break;
			case N:
					D[0] = 0.0;
					D[1] = 0.0;
					D[2] = 0.0;
					D[3] = 0.0;
					break;
		}
	};
};

//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel{
protected:
	int numReadGroups;
	int totNumParams;
	int* readGroupShifts;
	TRecalibrationEMQualityPositionMap qualPosMap;

	double** betas; //betas of the model
	double** oldBetas; //use during estimation
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	bool initialized;
	bool EMParamsInitialized;

	void _initialize(int NumReadGroups);
	double _calcEpsilon(double & eta);

public:
	int numParams;
	long numSitesAdded;

	TRecalibrationEMModel();
	TRecalibrationEMModel(int NumReadGroups);
	virtual ~TRecalibrationEMModel(){
		delete[] readGroupShifts;
		for(int r=0; r<numReadGroups; ++r){
			delete[] betas[r];
			delete[] oldBetas[r];
		}
		delete[] betas;
		delete[] oldBetas;
	};

	bool setParams(std::vector<std::string> & vec, int & rg);
	void initializeEMParams();
	void setEMParamsToZero();
	virtual double calcEpsilon(TRecalibrationEMReadData & data);
	virtual void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data);
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	void rejectProposedParameters();
	double getSteepestGradient();
	virtual void writeHeader(std::ofstream & out);
	void writeParametersToFile(std::ofstream & out, const uint8_t & readGroup);
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();
	virtual double getErrorRate(int rg, double originalErrorRate, const uint8_t & posInRead, const uint8_t & context);
};

class TRecalibrationEMModelNoContext:public TRecalibrationEMModel{
public:
	TRecalibrationEMModelNoContext(int NumReadGroups);

	double calcEpsilon(TRecalibrationEMReadData & data);
	void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data);
	void writeHeader(std::ofstream & out);
	double getErrorRate(int rg, double originalErrorRate, const int & posInRead, const uint8_t & context);
};

class TRecalibrationEMModelPositionSpecific:public TRecalibrationEMModel{
private:
	int maxPos;
public:
	TRecalibrationEMModelPositionSpecific(int NumReadGroups, int maxPos);

	double calcEpsilon(TRecalibrationEMReadData & data);
	void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data);
	void writeHeader(std::ofstream & out);
	double getErrorRate(int rg, double originalErrorRate, const int & posInRead, const uint8_t & context);
};

//--------------------------------------------------------------------
// TRecalibrationEM
//--------------------------------------------------------------------
class TRecalibrationEMSite{
public:
	TRecalibrationEMReadData* data;
	float P_g_given_d_oldBeta[4];
	unsigned int numReads;

	TRecalibrationEMSite();
	~TRecalibrationEMSite();
	TRecalibrationEMSite(TSite & site, int* readGroupMap, TQualityMap & qualiMap);
	double dePhred(double quality){
		double tmp = pow(10.0, quality / -10.0);
		if(tmp < 0.0000000001) return 0.0000000001;
		if(tmp > 0.9999999999) return 0.9999999999;
		return tmp;
	};
	void calcEpsilon(TRecalibrationEMModel* & model, float* & epsilon);
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel* & model, float* & freqs, float* & epsilon);
	double calcLL(TRecalibrationEMModel* & model, float* & freqs, float* & epsilon);
	double calcQ(TRecalibrationEMModel* & model, float* & epsilon);
	void addToJacobianAndF(TRecalibrationEMModel* & model, float* & epsilon);
};

class TRecalibrationEMWindow{
public:
	std::vector<TRecalibrationEMSite*> sites;
	float* freqs; //base frequencies
	int* readGroupMap;

	TRecalibrationEMWindow(TBaseFrequencies* baseFreqs, int* ReadGroupMap);
	virtual ~TRecalibrationEMWindow(){
		delete[] freqs;
		for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
			delete *site;
		}
		sites.clear();
	};
	unsigned int getMaxDepth();
	virtual void addSite(TSite & site, TQualityMap & qualiMap);
	long numSites();
	long numSitesDepthTwoOrMore();
	long cumulativeDepth();
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	double calcLL(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	double calcQ(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	void addToJacobianAndF(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	void setEuqalBaseFrequencies();
};

class TRecalibrationEM:public TRecalibration{
public:
	TLog* logfile;
	BamTools::SamHeader* bamHeader;
	std::string* readGroupNames;
	TRecalibrationEMModel* model;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;

	//variables for EM
	bool estimatetionRequired;
	bool equalBaseFrequencies;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	unsigned int maxDepth; //sites with higher depth will be ignored
	float* tmpEpsilon;
	bool tmpEpsilonInitialized;

	TRecalibrationEM(BamTools::SamHeader* BamHeader, std::string &name, TParameters & params, TLog* Logfile, TReadGroupMap& ReadGroupMap);
	~TRecalibrationEM(){
		delete[] readGroupNames;
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
		delete model;
		if(tmpEpsilonInitialized)
			delete[] tmpEpsilon;
	};
	bool recalibrationChangesQualities(){ return true; };
	bool requiresEstimation(){ return estimatetionRequired;};
	void addNewWindow(TBaseFrequencies* freqs);
	void addSite(TSite & site);
	long numSites();
	long numSitesDepthTwoOrMore();
	long cumulativeDepth();
	void prepareWindowsforEM();
	void runNewtonRaphson(int & maxNewtonraphsonIteratios, double & maxFThreshold, TLog* logfile, bool & writeTmpTables, std::string debugFilename);
	void runEM(std::string outputName, bool & writeTmpTables);
	void writeCurrentEstimates(std::string filename, double & LL);
	void writeHeader(std::ofstream & out);
	void writeParams(std::ofstream & out, double & LL);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
	double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	int getQuality(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
};


#endif /* TRecalIBRATION_H_ */
