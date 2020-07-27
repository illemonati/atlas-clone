/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMESTIMATOR_H_
#define TRECALIBRATIONEMESTIMATOR_H_

#include "auxiliaryTools.h"
#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TSite.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"

namespace GenotypeLikelihoods{


/*
//--------------------------------------------------------------------
// TRecalibrationEMSite
// Object to store all data at one site
//--------------------------------------------------------------------
class TRecalibrationEMSite{
private:
	std::vector<BAM::TBase> data;
	//TRecalibrationEMReadData* data;
	Base trueBase;
	double P_bbar_given_d_oldTheta[4];

	inline void calcEpsilon(TSequencingErrorModels & models, double* epsilon){
		for(auto& b : data){
			epsilon[k] = models.getErrorRate(b);
		}
	};

	inline double calcB(const double & D){
		return 1.33333333333333333333 * D - 1.0;
	};

	void _save(TSite & site, BAM::TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap);
	void _addToJacobianAndF(TSequencingErrorModels & models, double* epsilon);
	void _addToJacobianAndFKnownGenotype(TSequencingErrorModels & models, double* epsilon);

public:
	//unsigned int numReads;

	TRecalibrationEMSite();
	TRecalibrationEMSite(TSite & site, BAM::TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap);
	TRecalibrationEMSite(TSite & site, BAM::TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap, const Base TrueBase);

	void addToDataTable(TRecalibrationEMDataTables & dataTable);
	double fill_P_g_given_d_beta_AND_calcLL(TSequencingErrorModels & models, double* freqs, double* epsilon);
	double calcLL(TSequencingErrorModels & models, double* freqs, double* epsilon);
	void addToQ(TSequencingErrorModels & models);
	void addToJacobianAndF(TSequencingErrorModels & models, double* epsilon);
};

*/

//--------------------------------------------------------------------
// TRecalibrationEMWindow
//--------------------------------------------------------------------
class TRecalibrationEMWindow{
private:
	std::vector<TSiteStorage> sites;
	TGenotypeDistribution* genoDist;


	//tmp variables



public:

	double freqs[4]; //base frequencies
	BAM::TReadGroupMap* readGroupMapObject;

	TRecalibrationEMWindow(const TBaseData & baseFreqs, BAM::TReadGroupMap* ReadGroupMap);

	uint32_t getMaxDepth();
	void addSite(TSite & site);
	size_t size();
	size_t numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTables & dataTables);
	uint64_t cumulativeDepth();


	void setGenotypeDistribution(TGenotypeDistribution* GenoDist){ genoDist = GenoDist; };
	const TBaseData& expectedBaseFrequencies(){ return genoDist->baseFrequencies(); };


	//------------ OLD -----------------
	double fill_P_g_given_d_beta_AND_calcLL(TSequencingErrorModels & models, double* tmpEpsilon);
	double calcLL(TSequencingErrorModels & models, double* tmpEpsilon);
	//double calcQ(TRecalibrationEMModels & models, double* & tmpEpsilon);
	void addToQ(TSequencingErrorModels & models);
	void addToJacobianAndF(TSequencingErrorModels & models, double* tmpEpsilon);
	void setEuqalBaseFrequencies();
	//------------ END OLD -----------------

	std::vector<TSiteStorage>::iterator begin(){ return sites.begin(); };
	std::vector<TSiteStorage>::iterator end(){ return sites.end(); };
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator:public TGenotypeLikelihoodCalculator{
protected:
	std::vector<TRecalibrationEMWindow> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;
	TRecalibrationEMDataTables dataTables;

	//variables for estimation
	bool equalBaseFrequencies;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	unsigned int minRequiredObservations;
	std::string recalFile; //file name in case a file with model is provided
	TSequencingErrorCovariateDefinition covariateDefitionForEstimation;
	unsigned int maxDepth; //sites with higher depth will be ignored

	void _initializeModels();
	void _runEM(std::string outputName, bool & writeTmpTables);
	void _runNewtonRaphson();

	//functions to estimate theta_epsilon (sequencing error rates)
	void _calculate_EMWeights_epsilon(std::vector<TBaseData> & EMWeights);
	double _calculate_Q_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d);
	void _calculate_J_F_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d);

	void _updateEM_theta_epsilon();


public:
	TRecalibrationEMEstimator(TParameters & args, BAM::TReadGroups* ReadGroups, TLog* Logfile, BAM::TReadGroupMap* ReadGroupMap);
	~TRecalibrationEMEstimator(){
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
	};

	//functions to add data
	void addNewWindow(const TBaseData & freqs);
	void addSite(TSite & site);
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
