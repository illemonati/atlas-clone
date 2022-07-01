/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMESTIMATOR_H_
#define TRECALIBRATIONEMESTIMATOR_H_

#include <array>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TReadGroups.h"
#include "TSite.h"

namespace BAM {
class TSequencedBase;
}
namespace GenotypeLikelihoods {
class TPostMortemDamage;
}
namespace GenotypeLikelihoods {
namespace SequencingError {
class TModelRecal;
}
} // namespace GenotypeLikelihoods
namespace GenotypeLikelihoods {
namespace SequencingError {
class TModels;
}
} // namespace GenotypeLikelihoods

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
// store pointers to models so that they can be estimated
//--------------------------------------------------------------------------

class TModelVectorForEstimation {
private:
	// vector of pointers to models that require estimation
	std::vector<TModelRecal *> _models;                    // non-owning
	std::vector<std::array<TModelRecal *, 2>> _modelIndex; // non-owning
public:
	TModelVectorForEstimation() = default;
	void reset(TModels &SequencingErrorModels, const RecalEstimatorTools::TRecalDataTables &DataTables,
			   const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
			   uint32_t MinRequiredObservations);

	size_t size() const { return _models.size(); };
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const;

	// functions to estimate rho
	void prepareRhoEstimationFromEMWeights();
	void addBaseForRhoEstimation(BAM::TSequencedBase &base, const TBaseLikelihoods &EMWeights);
	void estimateRho();

	// functions to estimate beta
	void setNewtonRaphsonParamsToZero();
	void addToFandJacobian(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d);
	void addToQ(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d);
	void setQToZero();
	double curQ();
	void solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();

	void writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const;

	std::string getModelsDefinition();
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator {
private:
	std::vector<TSite> _sites;
	std::unique_ptr<TGenotypeDistribution> _genoDist;
	std::vector<TGenotypeLikelihoods> _P_g_I_d;
	const BAM::TReadGroupMap *_readGroupMap;
	const BAM::TReadGroups *_readGroups;

	//TGenotypeLikelihoodCalculator_simple _genotypeLikelihoodCalculator;
	TModelVectorForEstimation _modelsToEstimate;
	RecalEstimatorTools::TRecalDataTables _dataTables;

	// variables for estimation
	bool equalBaseFrequencies;
	int _numEMIterations;
	double _minDeltaLL;
	int _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	unsigned int _minRequiredObservations;
	bool _writeTmpTables;
	std::string _recalFile; // file name in case a file with model is provided

	size_t _numSitesDepthTwoOrMore();
	void _initializeModels(TModels &SequencingErrorModels);
	void _runEM(const std::string &outputName, const TPostMortemDamage &PmdModels);

	// functions to estimate theta_epsilon (sequencing error rates)
	void _updateRho(const TPostMortemDamage &PmdModels);
	double _calculate_Q_beta();
	void _calculate_J_F_beta();
	void _updateEM_theta_epsilon(const TPostMortemDamage &PmdModels);
	double _calculate_LL_updateWeights(const TPostMortemDamage &PmdModels);

public:
	TRecalibrationEMEstimator(const BAM::TReadGroups *ReadGroups, const BAM::TReadGroupMap *ReadGroupMap);

	// functions to add data
	void addSite(const TSite &site);

	// function to estimate
	void performEstimation(const std::string &outputName, TModels &SequencingErrorModels,
						   const TPostMortemDamage &PmdModels);

	void writeCurrentEstimates(const std::string &filename);
	double calcLL();
	void calcLikelihoodSurface(const std::string &filename, int numMarginalGridPoints);
	void calcQSurface(const std::string &filename, int numMarginalGridPoints);
};

}; // namespace SequencingError

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
