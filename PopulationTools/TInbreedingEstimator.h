//
// Created by caduffm on 3/8/22.
//

#ifndef ATLAS_TINBREEDINGESTIMATOR_H
#define ATLAS_TINBREEDINGESTIMATOR_H

#include <string>
#include <vector>

#include "coretools/Types/probability.h"

#include "genometools/PhredProbabilityTypes.h"
#include "genometools/THardyWeinbergGenotypeProbabilities.h"
#include "genometools/TSampleLikelihoods.h"
#include "genometools/VCF/TPopulationLikelihoods.h"

#include "coretools/Types/TStringHash.h"
#include "stattools/Deterministics/TLinkFunctions.h"
#include "stattools/ParametersObservations/TObservation.h"
#include "stattools/ParametersObservations/TParameter.h"
#include "stattools/Priors/TPriorBernouilli.h"
#include "stattools/Priors/TPriorBeta.h"
#include "stattools/Priors/TPriorUniform.h"

namespace stattools {
class TObservationBase;
}

namespace PopulationTools {

class TInbreedingEstimatorPrior; // forward declaration

//------------------------------------------
// Types and parameter specifications
//------------------------------------------

namespace inbr {

using TypeF        = coretools::Probability;      // F in [0, 1]
using TypeP        = coretools::Probability;      // p in [0, 1]
using TypeLogGamma = coretools::Unbounded;        // log_gamma in (-inf, inf) -> gamma = exp(log_gamma) is > 0
using TypeGamma    = coretools::StrictlyPositive; // gamma in (0, inf]
using TypePi       = coretools::ZeroOneOpen;      // pi in (0, 1)
using TypeFModel   = coretools::Boolean;          // FModel = 0,1
using TypePModel   = coretools::Boolean;          // PModel = 0,1
using PhredType    = genometools::HighPrecisionPhredIntProbability;
using TypeGTL      = genometools::TSampleLikelihoods<PhredType>;

constexpr static size_t NumDimParams = 1;
constexpr static size_t NumDimGTL    = 2;

using BoxOnFModel = stattools::prior::TUniformFixed<stattools::TParameterBase, TypeFModel, NumDimParams>;
using SpecFModel  = stattools::ParamSpec<TypeFModel, stattools::Hash<coretools::toHash("fModel")>, BoxOnFModel>;

using BoxOnF = stattools::prior::TUniformFixed<stattools::TParameterBase, TypeF, NumDimParams>;
using SpecF  = stattools::ParamSpec<TypeF, stattools::Hash<coretools::toHash("f")>, BoxOnF,
                                    stattools::RJMCMC<SpecFModel, coretools::probdist::TBetaDistr>>;

using BoxOnPi = stattools::prior::TUniformFixed<stattools::TParameterBase, TypePi, NumDimParams>;
using SpecPi  = stattools::ParamSpec<TypePi, stattools::Hash<coretools::toHash("pi")>, BoxOnPi>;

using BoxOnPModel = stattools::prior::TBernouilliInferred<stattools::TParameterBase, TypePModel, NumDimParams, SpecPi>;
using SpecPModel  = stattools::ParamSpec<TypePModel, stattools::Hash<coretools::toHash("pModel")>, BoxOnPModel>;

using BoxOnLogGamma = stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogGamma, NumDimParams>;
using SpecLogGamma  = stattools::ParamSpec<TypeLogGamma, stattools::Hash<coretools::toHash("logGamma")>, BoxOnLogGamma>;

using BoxOnGamma = stattools::det::TExp<stattools::TParameterBase, TypeGamma, NumDimParams, SpecLogGamma>;
using SpecGamma  = stattools::ParamSpec<TypeGamma, stattools::Hash<coretools::toHash("gamma")>, BoxOnGamma>;

using BoxOnP = stattools::prior::TBetaSymmetricZeroMixtureInferred<stattools::TParameterBase, TypeP, NumDimParams,
                                                                   SpecGamma, SpecPModel, TInbreedingEstimatorPrior>;
using SpecP =
    stattools::ParamSpec<TypeP, stattools::Hash<coretools::toHash("p")>, BoxOnP,
                         stattools::RJMCMC<SpecPModel, coretools::probdist::TBetaDistr>, stattools::Parallelize<stattools::MarkovOrder::different>>;

using BoxOnObs = TInbreedingEstimatorPrior;
using SpecObs  = stattools::TObservation<TypeGTL, NumDimGTL, BoxOnObs>;

} // namespace inbr

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

class TInbreedingEstimatorPrior
    : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, inbr::TypeGTL, 2> {
private:
	using Base = TStochasticBase<stattools::TObservationBase, inbr::TypeGTL, 2>;

	using typename Base::Storage;
	using typename Base::UpdatedStorage;
	using BoxType = TInbreedingEstimatorPrior;

	using TypeParamF      = stattools::TParameter<inbr::SpecF, BoxType>;
	using TypeParamP      = stattools::TParameter<inbr::SpecP, BoxType>;
	using TypeParamFModel = stattools::TParameter<inbr::SpecFModel, BoxType>;
	using TypeParamPModel = stattools::TParameter<inbr::SpecPModel, BoxType>;

	// parameters
	TypeParamF *_F;
	TypeParamP *_p;
	TypeParamFModel *_FModel;
	TypeParamPModel *_pModel;

	// dimensions
	size_t _numLoci    = 0;
	size_t _numSamples = 0;

	// openMP
	size_t _numThreads = 1;

	// initial estimates for p
	const std::vector<double> &_initialEstimatesP;

	// calculate LL
	[[nodiscard]] coretools::LogProbability
	_calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
	                               const genometools::THardyWeinbergGenotypeProbabilities &Probs) const;
	[[nodiscard]] coretools::LogProbability _calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
	                                                                       coretools::Probability P) const;
	[[nodiscard]] coretools::LogProbability _calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
	                                                                       coretools::Probability P,
	                                                                       coretools::Probability F) const;

	// DAG
	void _readCommandLineArguments();

	// initial values
	void _setInitialF();
	void _setInitialP();

	// common update function
	[[nodiscard]] double _calculateLLRatio_UpdateP(size_t Locus) const;

	// update functions for F
	[[nodiscard]] double _calculateLLRatio_FToHWE() const;
	[[nodiscard]] double _calculateLLRatio_HWEToF() const;

	// simulate values
	void _simulateUnderPrior(Storage *) override;

public:
	TInbreedingEstimatorPrior(TypeParamF *F, TypeParamP *P, TypeParamFModel *FModel, TypeParamPModel *PModel,
	                          const std::vector<double> &InitialEstimatesP);
	void initialize() override;

	[[nodiscard]] std::string name() const override;

	// estimate initial values
	void guessInitialValues() override;

	// updates
	[[nodiscard]] double calculateLLRatio(TypeParamF *, size_t) const;
	[[nodiscard]] double calculateLLRatio(TypeParamP *, size_t) const;
	[[nodiscard]] double calculateLLRatio(TypeParamFModel *, size_t) const;
	[[nodiscard]] double calculateLLRatio(TypeParamPModel *, size_t) const;
	void updateTempVals(TypeParamF *, size_t, bool);
	void updateTempVals(TypeParamP *, size_t, bool);
	void updateTempVals(TypeParamFModel *, size_t, bool);
	void updateTempVals(TypeParamPModel *, size_t, bool);

	// full log densities
	[[nodiscard]] double getSumLogPriorDensity(const Storage &Data) const override;
};

//------------------------------------------
// TInbreedingEstimatorModel
//------------------------------------------

class TInbreedingEstimatorModel {
private:
	// parameters and boxes
	inbr::BoxOnFModel _boxOnFModel;
	stattools::TParameter<inbr::SpecFModel, inbr::BoxOnObs> _FModel;

	inbr::BoxOnF _boxOnF;
	stattools::TParameter<inbr::SpecF, inbr::BoxOnObs> _F;

	inbr::BoxOnPi _boxOnPi;
	stattools::TParameter<inbr::SpecPi, inbr::BoxOnPModel> _pi;

	inbr::BoxOnPModel _boxOnPModel;
	stattools::TParameter<inbr::SpecPModel, inbr::BoxOnObs> _pModel;

	inbr::BoxOnLogGamma _boxOnLogGamma;
	stattools::TParameter<inbr::SpecLogGamma, inbr::BoxOnGamma> _logGamma;

	inbr::BoxOnGamma _boxOnGamma;
	stattools::TParameter<inbr::SpecGamma, inbr::BoxOnP> _gamma;

	inbr::BoxOnP _boxOnP;
	stattools::TParameter<inbr::SpecP, inbr::BoxOnObs> _p;

	inbr::BoxOnObs _boxOnObs;
	stattools::TObservation<inbr::TypeGTL, inbr::NumDimGTL, inbr::BoxOnObs> _observation;

public:
	TInbreedingEstimatorModel(
	    const std::string &Filename,
	    const genometools::TPopulationLikelihoods<stattools::TValueFixed<inbr::TypeGTL>> &Likelihoods);
};

//------------------------------------------
// TInbreedingEstimator
//------------------------------------------

class TInbreedingEstimator {
private:
	// data
	genometools::TPopulationLikelihoods<stattools::TValueFixed<inbr::TypeGTL>> _likelihoods;

	// read input data
	void _readData();

public:
	void run();
};

} // end namespace PopulationTools

#endif // ATLAS_TINBREEDINGESTIMATOR_H
