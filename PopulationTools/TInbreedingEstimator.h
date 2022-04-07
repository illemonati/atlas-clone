//
// Created by caduffm on 3/8/22.
//

#ifndef ATLAS_TINBREEDINGESTIMATOR_H
#define ATLAS_TINBREEDINGESTIMATOR_H

#include <limits>

#include "TPopulationLikelihoods.h"

#include "TLog.h"
#include "TNamesPositions.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"

#include "TDAGBuilder.h"
#include "TMCMC.h"
#include "TParameterTyped.h"
#include "TPriorBeta.h"
#include "TPriorBernouilli.h"
#include "TPriorUniform.h"

namespace PopulationTools {

typedef coretools::Probability TypeF;                // F in [0, 1]
typedef coretools::ZeroOpenOneClosed<double> TypeP;  // p in (0, 1]
typedef coretools::WeakType<double> TypeLogGamma;    // log_gamma in (-inf, inf) -> gamma = exp(log_gamma) is > 0
typedef coretools::ZeroOneOpen<double> TypePi;       // pi in (0, 1)
typedef coretools::WeakType<bool> TypeFModel;        // FModel = 0,1
typedef coretools::WeakType<bool> TypePModel;        // FModel = 0,1
typedef genometools::HighPrecisionPhredIntProbability PhredType;
typedef genometools::TSampleLikelihoods<PhredType> TypeGTL;

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

class TInbreedingEstimatorPrior
    : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeGTL, 2> {
private:
	using typename TBase<stattools::TObservationBase, TypeGTL, 2>::Storage;

	// parameters
	stattools::TParameterTyped<TypeF, 1> *_F;
	stattools::TParameterTyped<TypeP, 1> *_p;
	stattools::TParameterTyped<TypeFModel, 1> *_FModel;
	stattools::TParameterTyped<TypeFModel, 1> *_pModel;

	// dimensions
	size_t _numLoci    = 0;
	size_t _numSamples = 0;

	// proposal probability for model switch F
	coretools::Probability _q_FModel_To_HWE        = 0;
	coretools::LogProbability _log_q_FModel_To_HWE = 0;

	// propose new values of F after model switch
	coretools::StrictlyPositive<double> _lambdaNewF = coretools::StrictlyPositive<double>::min();

	// proposal probability for model switch p
	coretools::Probability _q_PModel_To_NullModel        = 0;
	coretools::LogProbability _log_q_PModel_To_NullModel = 0;

	// propose new values of p after model switch
	coretools::StrictlyPositive<double> _lambdaNewP = coretools::StrictlyPositive<double>::min();

	// initial estimates for p
	const std::vector<double> &_initialEstimatesP;

	// calculate LL
	[[nodiscard]] coretools::LogProbability
	_calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
	                               const genometools::THardyWeinbergGenotypeProbabilities &Probs) const;
	[[nodiscard]] coretools::LogProbability _calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
	                                                                       double P) const;
	[[nodiscard]] coretools::LogProbability _calculateLLSumOverIndividuals(const Storage &Data, size_t Locus, double P,
	                                                                       double F) const;

	// DAG
	void _readCommandLineArguments();

	// initial values
	void _setInitialF();
	void _setInitialP();

	// common update function
	double _calculateLLRatio_UpdateP(const Storage &Data, size_t Locus);
	[[nodiscard]] double _calculateProbabilityOfProposingThisFOrP(double Value, coretools::StrictlyPositive<double> Lambda) const;

	// update functions for F
	void _updateFAndFModel(const Storage &Data);
	void _updateRegularF(const Storage &Data);
	void _updateFToHWE(const Storage &Data);
	void _updateHWEToF(const Storage &Data);
	[[nodiscard]] double _getRandomNewF() const;

	// update functions for p
	void _updateP(const Storage &Data, size_t Locus);
	void _updatePToNull(const Storage &Data, size_t Locus);
	void _updateRegularP(const Storage &Data, size_t Locus);
	void _updateNullToP(const Storage &Data, size_t Locus);
	[[nodiscard]] double _getRandomNewP() const;

public:
	TInbreedingEstimatorPrior(stattools::TParameterTyped<TypeF, 1> *F, stattools::TParameterTyped<TypeP, 1> *P,
	                          stattools::TParameterTyped<TypeFModel, 1> *FModel, stattools::TParameterTyped<TypeFModel, 1> *PModel,
	                          const std::vector<double> &InitialEstimatesP);
	void initializeStorageOfPriorParameters() override;

	[[nodiscard]] std::string name() const override;

	// estimate initial values
	void estimateInitialPriorParameters() override;

	// full log densities
	[[nodiscard]] double getSumLogPriorDensity(const Storage &Data) const override;

	// update all hyperprior parameters
	void updateParams() override;

	// simulate values
	void simulateUnderPrior() override;
};

//------------------------------------------
// TInbreedingEstimatorModel
//------------------------------------------

class TInbreedingEstimatorModel {
private:
	// parameters
	stattools::TParameterTyped<TypeF, 1> _F;
	stattools::TParameterTyped<TypeFModel, 1> _FModel;
	stattools::TParameterTyped<TypePi, 1> _pi;
	stattools::TParameterTyped<TypePModel, 1> _pModel;
	stattools::TParameterTyped<TypeLogGamma, 1> _log_gamma;
	stattools::TParameterTyped<TypeP, 1> _p;
	stattools::TObservationTyped<TypeGTL, 2> _observation;

public:
	TInbreedingEstimatorModel(const std::string &Filename, std::shared_ptr<stattools::TDAGBuilder> &DAGBuilder,
	                          const genometools::TPopulationLikelihoods<stattools::TValueFixed<TypeGTL>> &Likelihoods);
};

//------------------------------------------
// TInbreedingEstimator
//------------------------------------------

class TInbreedingEstimator {
private:
	// data
	genometools::TPopulationLikelihoods<stattools::TValueFixed<TypeGTL>> _likelihoods;

	// read input data
	void _readData();

public:
	void runEstimation(coretools::TParameters &Parameters);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimateInbreeding : public coretools::TTask {
public:
	TTask_estimateInbreeding() {
		_explanation = "Estimating the inbreeding coefficient";
		_citations.emplace("Burger et al. (2020) Current Biology");
	};

	void run() override {
		using namespace coretools::instances;
		TInbreedingEstimator inbreedingEstimator;
		inbreedingEstimator.runEstimation(parameters());
	};
};

} // end namespace PopulationTools

#endif // ATLAS_TINBREEDINGESTIMATOR_H
