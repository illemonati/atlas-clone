/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMESTIMATOR_H_
#define TRECALIBRATIONEMESTIMATOR_H_

#include "../GenotypeLikelihoods/auxiliaryTools.h"
#include "../GenotypeLikelihoods/TSequencingErrorModel.h"
#include "../TSite.h"

namespace GenotypeLikelihoods{

//--------------------------------------------------------------------
// TRecalibrationEMSite
// Object to store all data at one site
//--------------------------------------------------------------------
class TRecalibrationEMSite{
private:
	TRecalibrationEMReadData* data;
	Base trueBase;
	double P_g_given_d_oldBeta[4];

	inline void calcEpsilon(TSequencingErrorModels & models, double* & epsilon){
		for(unsigned int k=0; k<numReads; ++k)
			epsilon[k] = models.calcEpsilon(data[k]);
	};

	inline double calcB(const double & D){
		return 1.33333333333333333333 * D - 1.0;
	};

	void _save(TSite & site, TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap);
	void _addToJacobianAndF(TSequencingErrorModels & models, double* & epsilon);
	void _addToJacobianAndFKnownGenotype(TSequencingErrorModels & models, double* & epsilon);

public:
	unsigned int numReads;

	TRecalibrationEMSite();
	~TRecalibrationEMSite();
	TRecalibrationEMSite(TSite & site, TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap);
	TRecalibrationEMSite(TSite & site, TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap, const Base TrueBase);

	void addToDataTable(TRecalibrationEMDataTables & dataTable);
	double fill_P_g_given_d_beta_AND_calcLL(TSequencingErrorModels & models, double* & freqs, double* & epsilon);
	double calcLL(TSequencingErrorModels & models, double* & freqs, double* & epsilon);
	void addToQ(TSequencingErrorModels & models);
	void addToJacobianAndF(TSequencingErrorModels & models, double* & epsilon);
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
	void addSite(TSite & site, TQualityMap & qualiMap);
	void addSite(TSite & site, TQualityMap & qualiMap, const Base TrueBase);
	long numSites();
	long numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTables & dataTable);
	long cumulativeDepth();
	double fill_P_g_given_d_beta_AND_calcLL(TSequencingErrorModels & models, double* & tmpEpsilon);
	double calcLL(TSequencingErrorModels & models, double* & tmpEpsilon);
	//double calcQ(TRecalibrationEMModels & models, double* & tmpEpsilon);
	void addToQ(TSequencingErrorModels & models);
	void addToJacobianAndF(TSequencingErrorModels & models, double* & tmpEpsilon);
	void setEuqalBaseFrequencies();
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator{
protected:
	TLog* logfile;
	TReadGroups* _readGroups;
	TReadGroupMap* _readGroupMap;
	TSequencingErrorModels* models;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;

	//variables for estimation
	bool equalBaseFrequencies;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	unsigned int minRequiredObservations;
	std::string recalFile; //file name in case a file with model is provided
	TSequencingErrorCovariateDefinition covariateDefitionForEstimation;
	double* tmpEpsilon;
	bool tmpEpsilonInitialized;
	unsigned int maxDepth; //sites with higher depth will be ignored

	void _initializeModels();
	void _runEM(std::string outputName, bool & writeTmpTables);
	void _runNewtonRaphson();
	void _initializTmpVariablesForEstimation();

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
	void addSite(TSite & site, TQualityMap & qualiMap, const Base TrueBase);
	long numSites();
	long numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTables & dataTable);
	long cumulativeDepth();

	//function to estimate
	void initializeFromFile(const std::string string);
	void performEstimation(std::string outputName, bool & writeTmpTables);
	void performEstimationKnownGenotypes(std::string outputName, bool & writeTmpTables);
	void writeCurrentEstimates(std::string filename);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
};

}; //end namespace

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
