/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMESTIMATOR_H_
#define TRECALIBRATIONEMESTIMATOR_H_

#include "RecalEstimatorTools.h"
#include "auxiliaryTools.h"
#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TSite.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"

namespace GenotypeLikelihoods{

namespace RecalEstimator{

//--------------------------------------------------------------------------------------
// TModelIndex
// Object to map read group ID and mate to an internal index used to store recal models
//--------------------------------------------------------------------------------------
class TModelIndex{
private:
	std::vector< std::array< std::shared_ptr<TSequencingErrorModelRecal>, 2> > _index;

public:
	TModelIndex(const BAM::TReadGroups& ReadGroups){
		_index.resize(ReadGroups.size());
	};
	~TModelIndex() = default;

	void set(uint16_t ReadGroupId, const bool & IsSecondMate, std::shared_ptr<TSequencingErrorModelRecal> & Model, const BAM::TReadGroupMap & ReadGroupMap){
		//also set for all models pooled
		for(auto& r : ReadGroupMap.readGroupsPooledWith(ReadGroupId)){
			_index[r][IsSecondMate] = Model;
		}
	};

	bool inUse(const BAM::TSequencedBase & base) const{
		return _index[base.readGroupID][(int) base.isSecondMate()].get();
	};

	TSequencingErrorModelRecal& operator()(const BAM::TSequencedBase & base) const{
		return *_index[base.readGroupID][base.isSecondMate()];
	};
};

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
// store pointers to models so that they can be estimated
//--------------------------------------------------------------------------

class TSequencingErrorModelVectorForEstimation{
private:
	//vector of pointers to models that require estimation
	std::vector< std::shared_ptr<TSequencingErrorModelRecal> > _models;
	TModelIndex _modelIndex;

public:
	TSequencingErrorModelVectorForEstimation(TSequencingErrorModels & SequencingErrorModels,
										     const RecalEstimatorTools::TRecalDataTables & DataTables,
											 const BAM::TReadGroups & ReadGroups,
											 const BAM::TReadGroupMap & ReadGroupMap,
											 uint32_t MinRequiredObservations,
											 coretools::TLog* Logfile);
	~TSequencingErrorModelVectorForEstimation() = default;

	size_t size() const { return _models.size(); };
	void fillBaseLikelihoods(const BAM::TSequencedBase & base,  TBaseLikelihoods & baseLikelihoods) const;

	//functions to estimate rho
	void prepareRhoEstimationFromEMWeights();
	void addBaseForRhoEstimation(BAM::TSequencedBase & base, const TBaseLikelihoods & EMWeights);
	void estimateRho();

	//functions to estimate beta
	void setNewtonRaphsonParamsToZero();
	void addToFandJacobian(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d);
	void addToQ(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d);
	void setQToZero();
	double curQ();
	bool solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();

	void print();
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator{
private:
	coretools::TLog* _logfile;
	std::vector<TSite> _sites;
	std::unique_ptr<TGenotypeDistribution> _genoDist;
	const BAM::TReadGroupMap* _readGroupMap;
	const BAM::TReadGroups* _readGroups;

	TGenotypeLikelihoodCalculator_simple _genotypeLikelihoodCalculator;
	std::unique_ptr<TSequencingErrorModelVectorForEstimation> _modelsToEstimate;
	RecalEstimatorTools::TRecalDataTables _dataTables;

	//variables for estimation
	bool equalBaseFrequencies;
	int _numEMIterations;
	double _minDeltaLL;
	int _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	unsigned int _minRequiredObservations;
	bool _writeTmpTables;
	std::string _recalFile; //file name in case a file with model is provided

	size_t _numSitesDepthTwoOrMore();
	void _initializeModels(TSequencingErrorModels & SequencingErrorModels);
	void _runEM(std::string outputName, const TPostMortemDamage & PmdModels);
	void _fillRelevantBaseFrequencies(TBaseProbabilities & baseFreq, const genometools::Genotype & genotype);

	//functions to estimate theta_epsilon (sequencing error rates)
	void _calculate_EMWeights_epsilon(std::vector<TBaseLikelihoods> & EMWeights, const TPostMortemDamage & PmdModels);
	double _calculate_Q_beta(const std::vector<TBaseLikelihoods> & EM_weights_bbar_given_d);
	void _calculate_J_F_beta(const std::vector<TBaseLikelihoods> & EM_weights_bbar_given_d);
	void _updateEM_theta_epsilon(const TPostMortemDamage & PmdModels);

	//function to calculate LL (currently uses haploid model)
	double _calculateLL_fullModel(const TPostMortemDamage & PmdModels);

public:
	TRecalibrationEMEstimator(coretools::TParameters & args, coretools::TLog* Logfile, const BAM::TReadGroups* ReadGroups, const BAM::TReadGroupMap* ReadGroupMap);

	//functions to add data
	void addSite(const TSite & site);

	//function to estimate
	void performEstimation(std::string outputName, TSequencingErrorModels & SequencingErrorModels, const TPostMortemDamage & PmdModels);

	void writeCurrentEstimates(std::string filename);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
};

}; //end namespace RecalEstimator

}; //end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
