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
//	virtual double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	virtual double getErrorRate(TBase & base);
//	double getErrorRateFromBase(const TBase & base);
	virtual int getQuality(TBase & base);
//	virtual int getQualityFromBase(const TBase & base, TQualityMap & qualMap);
//	virtual int getphredInt(const TBase & base){
//		return base.phredInt;
//	};
//	virtual int getQuality(const TBase & base){
//		return base.quality;
//	};

	virtual bool requiresEstimation(){ return false;};
	int findReadGroupIndex(std::string & name, BamTools::SamReadGroupDictionary & readGroups);
};

//--------------------------------------------------------------------
// TRecalibrationEMQualityPositionMap
// Look-up tables for position and quality. Only indexes will be stored.
//--------------------------------------------------------------------
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


//--------------------------------------------------------------------
// TRecalibrationEMReadData
// Per site data storage
//--------------------------------------------------------------------
struct TRecalibrationEMReadData{
	uint8_t quality;
	uint8_t position;
	float D[4];
	uint8_t context;
	uint16_t readGroup;
	bool isSecond;

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
// TRecalibrationEMDataTable
// Object to store for which qualities and positions data is available.
//--------------------------------------------------------------------
class TRecalibrationEMDataTable{
public:
	int numReadGroups;
	int maxQual;
	bool*** qualities; //qualities[readGroup][first/second][quality]
	unsigned int** maxPos; //maxPos[readGroup][first/second]

	TRecalibrationEMDataTable(int NumReadGroups, int MaxQual){
		numReadGroups = NumReadGroups;
		maxQual = MaxQual;

		qualities = new bool**[numReadGroups];
		maxPos = new unsigned int*[numReadGroups];
		for(int rg = 0; rg<numReadGroups; ++rg){
			qualities[rg] = new bool*[2];
			qualities[rg][0] = new bool[maxQual];
			qualities[rg][1] = new bool[maxQual];

			maxPos[rg] = new unsigned int[2];
		}

		clear();
	};

	~TRecalibrationEMDataTable(){
		for(int rg = 0; rg<numReadGroups; ++rg){
			delete[] qualities[rg][0];
			delete[] qualities[rg][1];
			delete[] qualities[rg];

			delete[] maxPos[rg];
		}

		delete[] qualities;
		delete[] maxPos;
	};

	void clear(){
		for(int rg = 0; rg<numReadGroups; ++rg){
			for(int i=0; i<maxQual; ++i){
				qualities[rg][0][i] = 0;
				qualities[rg][1][i] = 0;
			}
			maxPos[rg][0] = 0;
			maxPos[rg][1] = 0;
		}
	};

	void add(TRecalibrationEMReadData & data){
		++qualities[data.readGroup][(int) data.isSecond][data.quality];
		if(maxPos[data.readGroup][data.isSecond] < data.position)
				maxPos[data.readGroup][data.isSecond] = data.position;
	};

};

//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel{
protected:
	int numReadGroups;
	uint8_t defaultNumParams;
	uint8_t* numParamsPerReadGroup;
	//int numParameterSets; //is equal to numReadGroups * 2 if all read groups have parameters for first and second mates
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
	void _allocateBetaMemory();
	double _calcEpsilon(double & eta);

public:

	long numSitesAdded;

	TRecalibrationEMModel();
	virtual ~TRecalibrationEMModel(){
		delete[] readGroupShifts;
		for(int r=0; r<numReadGroups; ++r){
			delete[] betas[r];
			delete[] oldBetas[r];
		}
		delete[] betas;
		delete[] oldBetas;
	};

	virtual void initializeWithValues(std::vector<std::string> & vec, int NumReadGroups);
	virtual void initializeWithHeader(std::vector<std::string> & vec, int NumReadGroups);
	virtual void initializeWithDataTable(TRecalibrationEMDataTable & dataTable);
	void initializeEMParams();
	bool setParams(std::vector<std::string> & vec, int & rg);
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
	virtual double getErrorRate(TBase & base);
};

class TRecalibrationEMModelNoContext:public TRecalibrationEMModel{
public:
	TRecalibrationEMModelNoContext();

	double calcEpsilon(TRecalibrationEMReadData & data);
	void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data);
	void writeHeader(std::ofstream & out);
	double getErrorRate(TBase & base);
};

class TRecalibrationEMModelPositionSpecific:public TRecalibrationEMModel{
private:
	unsigned int* maxPos;

	void _initializeFromVectorSize(std::vector<std::string> & vec, int NumReadGroups);

public:
	TRecalibrationEMModelPositionSpecific();

	void initializeWithValues(std::vector<std::string> & vec, int NumReadGroups);
	void initializeWithHeader(std::vector<std::string> & vec, int NumReadGroups);
	void initializeWithDataTable(TRecalibrationEMDataTable & dataTable);
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
	TRecalibrationEMSite(TSite & site, TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap);

	void addToDataTable(TRecalibrationEMDataTable & dataTable);
/*
	double dePhred(double quality){
		double tmp = pow(10.0, quality / -10.0);
		if(tmp < 0.0000000001) return 0.0000000001;
		if(tmp > 0.9999999999) return 0.9999999999;
		return tmp;
	};

	*/
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
	TReadGroupMap* readGroupMapObject;

	TRecalibrationEMWindow(TBaseFrequencies* baseFreqs, TReadGroupMap* ReadGroupMap);
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
	void addToDataTable(TRecalibrationEMDataTable & dataTable);
	long cumulativeDepth();
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	double calcLL(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	double calcQ(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	void addToJacobianAndF(TRecalibrationEMModel* & model, float* & tmpEpsilon);
	void setEuqalBaseFrequencies();
};

class TRecalibrationEM:public TRecalibration{
private:
	TLog* logfile;
	BamTools::SamHeader* bamHeader;
	std::string* readGroupNames;
	TRecalibrationEMModel* model;
	bool modelInitialized;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;

	//variables for EM
	bool equalBaseFrequencies;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	unsigned int maxDepth; //sites with higher depth will be ignored
	float* tmpEpsilon;
	bool tmpEpsilonInitialized;

	void _initializeModel(std::string modelTag);

public:
	TRecalibrationEM(BamTools::SamHeader* BamHeader, TLog* Logfile, TReadGroupMap& ReadGroupMap);
	~TRecalibrationEM(){
		if(modelInitialized)
			delete model;
		delete[] readGroupNames;
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
		delete model;
		if(tmpEpsilonInitialized)
			delete[] tmpEpsilon;
	};

	void initialize(std::string string);
	void initializeRecalibrationParametersFromString(std::string & string);
	void initializeRecalibrationParametersFromFile(std::string filename);
	void initializeForParameterEstimation(TParameters & args);
	bool recalibrationChangesQualities(){ return true; };
	void addNewWindow(TBaseFrequencies* freqs);
	void addSite(TSite & site, TQualityMap & qualiMap);
	long numSites();
	long numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTable & dataTable);
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
	double getErrorRate(TBase & base);
	int getQuality(TBase & base);
};


#endif /* TRecalIBRATION_H_ */
