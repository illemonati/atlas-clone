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

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator:public TGenotypeLikelihoodCalculator{
private:
	std::vector<TSite> _sites;
	TGenotypeDistribution* _genoDist;
	TRecalibrationEMDataTables _dataTables;

	//variables for estimation
	bool equalBaseFrequencies;
	int _numEMIterations;
	double _minDeltaLL;
	int _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	unsigned int _minRequiredObservations;
	std::string _recalFile; //file name in case a file with model is provided
	TSequencingErrorCovariateDefinition _covariateDefitionForEstimation;

	void _initializeModels();
	void _runEM(std::string outputName, bool & writeTmpTables);
	void _fillRelevantBaseFrequencies(TBaseData & baseFreq, const Genotype genotype);

	//functions to estimate theta_epsilon (sequencing error rates)
	void _calculate_EMWeights_epsilon(std::vector<TBaseData> & EMWeights);
	double _calculate_Q_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d);
	void _calculate_J_F_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d);
	void _updateEM_theta_epsilon();

	//function to calculate LL (currently uses haploid model)
	double _calculateLL_fullModel();

public:
	TRecalibrationEMEstimator(TParameters & args, BAM::TReadGroups* ReadGroups, TLog* Logfile, BAM::TReadGroupMap* ReadGroupMap);

	//functions to add data
	void addSite(const TSite & site);
	size_t numSites();
	size_t numSitesDepthTwoOrMore();
	void addToDataTable(TRecalibrationEMDataTables & dataTable);
	size_t cumulativeDepth();

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
