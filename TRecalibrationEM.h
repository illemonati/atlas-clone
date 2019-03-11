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
	unsigned int minRequiredObservations;

	TRecalibrationEMEstimationParameters(TParameters & args, TLog* logfile);
};

//--------------------------------------------------------------------
// TRecalibrationEMSite
// Object to store all data at one site
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
	inline void calcEpsilon(TRecalibrationEMModels & models, float* & epsilon){
		for(unsigned int k=0; k<numReads; ++k)
			epsilon[k] = models.calcEpsilon(data[k]);
	};
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, float* & freqs, float* & epsilon);
	double calcLL(TRecalibrationEMModels & models, float* & freqs, float* & epsilon);
	double calcQ(TRecalibrationEMModels & models, float* & epsilon);
	void addToJacobianAndF(TRecalibrationEMModels & models, float* & epsilon);
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
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, float* & tmpEpsilon);
	double calcLL(TRecalibrationEMModels & models, float* & tmpEpsilon);
	double calcQ(TRecalibrationEMModels & models, float* & tmpEpsilon);
	void addToJacobianAndF(TRecalibrationEMModels & models, float* & tmpEpsilon);
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
	TRecalibrationEMModels* models;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;

	//variables for estimation
	bool estimationParametersInitialized;
	TRecalibrationEMEstimationParameters* EMParameters;
	std::string modelTagForEstimation;
	float* tmpEpsilon;
	bool tmpEpsilonInitialized;
	unsigned int maxDepth; //sites with higher depth will be ignored

	void runEM(int numSitesWithData, std::string outputName, bool & writeTmpTables);
	void runNewtonRaphson(int numSitesWithData);
	void prepareWindowsforEM();


public:
	TRecalibrationEM(BamTools::SamHeader* BamHeader, TLog* Logfile, TReadGroupMap& ReadGroupMap);
	~TRecalibrationEM(){
		delete models;
		delete[] readGroupNames;
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
		if(tmpEpsilonInitialized)
			delete[] tmpEpsilon;

		if(estimationParametersInitialized)
			delete EMParameters;
	};

	bool recalibrationChangesQualities(){ return true; };

	void initialize(std::string string);
	void initializeRecalibrationParametersFromString(std::string & string);
	void initializeRecalibrationParametersFromFile(std::string filename);
	void initializeForParameterEstimation(TParameters & args);

	void performEstimation(std::string outputName, bool & writeTmpTables);

	void addNewWindow(TBaseFrequencies* freqs);
	void addSite(TSite & site, TQualityMap & qualiMap);
	long numSites();
	long numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTable & dataTable);
	long cumulativeDepth();

	void writeCurrentEstimates(std::string filename);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);

	inline double getErrorRate(TBase & base){
		return models->getErrorRate(base);
	};
	inline int getQuality(TBase & base){
		double q = models->getErrorRate(base);
		return qualityMap.errorToQuality(q);
	};

};




#endif /* TRECALIBRATIONEM_H_ */
