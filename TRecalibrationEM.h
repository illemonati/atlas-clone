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
/*
	double dePhred(double quality){
		double tmp = pow(10.0, quality / -10.0);
		if(tmp < 0.0000000001) return 0.0000000001;
		if(tmp > 0.9999999999) return 0.9999999999;
		return tmp;
	};

	*/
	void calcEpsilon(TRecalibrationEMModel_Base* & model, float* & epsilon);
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel_Base* & model, float* & freqs, float* & epsilon);
	double calcLL(TRecalibrationEMModel_Base* & model, float* & freqs, float* & epsilon);
	double calcQ(TRecalibrationEMModel_Base* & model, float* & epsilon);
	void addToJacobianAndF(TRecalibrationEMModel_Base* & model, float* & epsilon);
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
	void addToJacobianAndF(TRecalibrationEMModel_Base* & model, float* & tmpEpsilon);
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
	TRecalibrationEMModel_Base* model;
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
	void runNewtonRaphson(int & maxNewtonraphsonIteratios, double & maxFThreshold, TLog* logfile);
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




#endif /* TRECALIBRATIONEM_H_ */
