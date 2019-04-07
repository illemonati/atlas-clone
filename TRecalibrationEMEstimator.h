/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMESTIMATOR_H_
#define TRECALIBRATIONEMESTIMATOR_H_

#include "TRecalibrationEMModel.h"
#include "TSite.h"


//--------------------------------------------------------------------
// TRecalibrationEMSite
// Object to store all data at one site
//--------------------------------------------------------------------
class TRecalibrationEMSite{
public:
	TRecalibrationEMReadData* data;
	double P_g_given_d_oldBeta[4];
	unsigned int numReads;

	TRecalibrationEMSite();
	~TRecalibrationEMSite();
	TRecalibrationEMSite(TSite & site, TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap);

	void addToDataTable(TRecalibrationEMDataTable & dataTable);
	inline void calcEpsilon(TRecalibrationEMModels & models, double* & epsilon){
		for(unsigned int k=0; k<numReads; ++k)
			epsilon[k] = models.calcEpsilon(data[k]);
	};
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, double* & freqs, double* & epsilon);
	double calcLL(TRecalibrationEMModels & models, double* & freqs, double* & epsilon);
	double calcQ(TRecalibrationEMModels & models, double* & epsilon);
	void addToQ(TRecalibrationEMModels & models);
	void addToJacobianAndF(TRecalibrationEMModels & models, double* & epsilon);
};

//--------------------------------------------------------------------
// TRecalibrationEMWindow
//--------------------------------------------------------------------
class TRecalibrationEMWindow{
public:
	std::vector<TRecalibrationEMSite*> sites;
	double* freqs; //base frequencies
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
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, double* & tmpEpsilon);
	double calcLL(TRecalibrationEMModels & models, double* & tmpEpsilon);
	double calcQ(TRecalibrationEMModels & models, double* & tmpEpsilon);
	void addToQ(TRecalibrationEMModels & models);
	void addToJacobianAndF(TRecalibrationEMModels & models, double* & tmpEpsilon);
	void setEuqalBaseFrequencies();
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator{
private:
	TLog* logfile;
	TReadGroups* _readGroups;
	TReadGroupMap* _readGroupMap;
	TRecalibrationEMModels* models;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;

	//variables for estimation
	bool equalBaseFrequencies;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	unsigned int minRequiredObservations;
	std::string modelTagForEstimation;
	double* tmpEpsilon;
	bool tmpEpsilonInitialized;
	unsigned int maxDepth; //sites with higher depth will be ignored

	void _runEM(int numSitesWithData, std::string outputName, bool & writeTmpTables);
	void _runNewtonRaphson(int numSitesWithData);
	void _initializeTmpEpsilon();


public:
	TRecalibrationEMEstimator(TParameters & args, TReadGroups* ReadGroups, TLog* Logfile, TReadGroupMap* ReadGroupMap);
	~TRecalibrationEMEstimator(){
		delete models;
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
		if(tmpEpsilonInitialized)
			delete[] tmpEpsilon;
	};

	//functions to add data
	void addNewWindow(TBaseFrequencies* freqs);
	void addSite(TSite & site, TQualityMap & qualiMap);
	long numSites();
	long numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTable & dataTable);
	long cumulativeDepth();

	//function to estimate
	void initializeFromString(const std::string string);
	void performEstimation(std::string outputName, bool & writeTmpTables);
	void writeCurrentEstimates(std::string filename);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
};




#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
