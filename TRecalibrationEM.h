/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEM_H_
#define TRECALIBRATIONEM_H_

#include "TRecalibration.h"
#include "TRecalibrationEMModel.h"

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
	unsigned int** countsPerReadGroup;
	unsigned int totalCounts;

	TRecalibrationEMDataTable(int NumReadGroups, int MaxQual);
	~TRecalibrationEMDataTable();

	void clear();
	void add(TRecalibrationEMReadData & data);
	void assembleCountsPerReadGroup();
};

//--------------------------------------------------------------------
// TRecalibrationEMReadGroupIndex
// Object to map read group ID and mate to an internal index used to store recal parameters
//--------------------------------------------------------------------
class TRecalibrationEMReadGroupIndex{
private:
	bool initialized;
	int** readGroupIndex;
	bool** readGroupInUse;

	void _init(int NumReadGroups);
	void _free();


public:
	int numReadGroups;
	int numIndexes;

	TRecalibrationEMReadGroupIndex();
	~TRecalibrationEMReadGroupIndex();

	void initialize(int NumReadGroups);
	void initialize(TRecalibrationEMDataTable & dataTable);

	void setAllAsUsed();
	void setAllToSingleIndex();
	int setAsUsed(int readGroup, bool isSecondMate);

	bool inUse(const int & readGroup, const bool & isSecondMate){
		return readGroupInUse[readGroup][isSecondMate];
	};

	int index(const int & readGroup, const bool & isSecondMate){
		return readGroupIndex[readGroup][isSecondMate];
	};

	bool nextNotInUse(std::pair<int, bool> & pair);
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimationParameters
// Object to store all parameters related to the EM algorithm
//--------------------------------------------------------------------
class TRecalibrationEMEstimationParameters{
public:
	bool equalBaseFrequencies;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	unsigned int minRequiredObservationsPerReadGroup;

	TRecalibrationEMEstimationParameters(TParameters & args, TLog* logfile);
};


//--------------------------------------------------------------------
// TRecalibrationEMSite
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
	void calcEpsilon(TRecalibrationEMModel_Base* & model, float* & epsilon);
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel_Base* & model, float* & freqs, float* & epsilon);
	double calcLL(TRecalibrationEMModel_Base* & model, float* & freqs, float* & epsilon);
	double calcQ(TRecalibrationEMModel_Base* & model, float* & epsilon);
	void addToJacobianAndF(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMModel_Base* & model, float* & epsilon);
};

//--------------------------------------------------------------------
// TRecalibrationEMWindow
//--------------------------------------------------------------------
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
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel_Base* & model, float* & tmpEpsilon);
	double calcLL(TRecalibrationEMModel_Base* & model, float* & tmpEpsilon);
	double calcQ(TRecalibrationEMModel_Base* & model, float* & tmpEpsilon);
	void addToJacobianAndF(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMModel_Base* & model, float* & tmpEpsilon);
	void setEuqalBaseFrequencies();
};

//--------------------------------------------------------------------
// TRecalibrationEM
//--------------------------------------------------------------------
class TRecalibrationEM:public TRecalibration{
private:
	TLog* logfile;
	BamTools::SamHeader* bamHeader;
	std::string* readGroupNames;
	std::vector<TRecalibrationEMModel_Base*> models;
	bool modelsInitialized;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;
	TRecalibrationEMReadGroupIndex readGroupIndex;

	//variables for estimation
	bool estimationParametersInitialized;
	TRecalibrationEMEstimationParameters* EMParameters;
	std::string modelTagForEstimation;
	float* tmpEpsilon;
	bool tmpEpsilonInitialized;
	unsigned int maxDepth; //sites with higher depth will be ignored

	void _addModel(std::string modelTag, std::vector<std::string> & values, bool verbose);
	void _addModelForEstimation(std::string modelTag, int maxPos);
	void _printReadGroupNameToLogfile(std::pair<int, bool> & missingReadGroupInfo);

public:
	TRecalibrationEM(BamTools::SamHeader* BamHeader, TLog* Logfile, TReadGroupMap& ReadGroupMap);
	~TRecalibrationEM(){
		if(modelsInitialized){
			for(int i=0; i<readGroupIndex.numIndexes; ++i)
				delete models[i];
			delete[] models;
		}
		delete[] readGroupNames;
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
		delete models;
		if(tmpEpsilonInitialized)
			delete[] tmpEpsilon;
	};

	void initialize(std::string string);
	void initializeRecalibrationParametersFromString(std::string & string);
	void initializeRecalibrationParametersFromFile(std::string filename);
	void performEstimation(std::string outputName, bool & writeTmpTables);
	void runEM(std::string outputName, bool & writeTmpTables);
	void initializeForParameterEstimation(TParameters & args);
	bool recalibrationChangesQualities(){ return true; };
	void addNewWindow(TBaseFrequencies* freqs);
	void addSite(TSite & site, TQualityMap & qualiMap);
	long numSites();
	long numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTable & dataTable);
	long cumulativeDepth();
	void prepareWindowsforEM();
	void runNewtonRaphson(int & maxNewtonraphsonIteratios, double & maxFThreshold, TLog* logfile);
	void writeCurrentEstimates(std::string filename, double & LL);
	void writeHeader(std::ofstream & out);
	void writeParams(std::ofstream & out, double & LL);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
	double getErrorRate(TBase & base);
	int getQuality(TBase & base);
};




#endif /* TRECALIBRATIONEM_H_ */
