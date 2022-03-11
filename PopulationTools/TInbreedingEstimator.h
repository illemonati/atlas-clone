//
// Created by caduffm on 3/8/22.
//

#ifndef ATLAS_TINBREEDINGESTIMATOR_H
#define ATLAS_TINBREEDINGESTIMATOR_H

#include "TParameters.h"
#include "TLog.h"
#include "TPopulationLikelihoods.h"
#include "TRandomGenerator.h"
#include <limits>
#include "TDAGBuilder.h"
#include "TPriorFixed.h"
#include "TPriorInferred.h"
#include "TMCMC.h"
#include "TTask.h"
#include "TNamesPositions.h"

namespace PopulationTools {

typedef coretools::Probability TypeF; // F in [0, 1]
typedef coretools::ZeroOpenOneClosed<double> TypeP; // p in (0, 1]
typedef coretools::WeakType<double> TypeC; // c in (-inf, inf) -> gamma = exp(c) is > 0
typedef coretools::WeakType<bool> TypeZ; // z = 0,1
typedef genometools::HighPrecisionPhredIntProbability PhredType;
typedef genometools::TSampleLikelihoods<PhredType> TypeGTL;

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

class InbreedingPrior : public stattools::MCMCPrior {
public:
    inline static const std::string inbreeding = "inbreeding";
};

class TInbreedingEstimatorPrior : public stattools::prior::TBaseNonIID<stattools::TValueFixed, TypeGTL, 2> {
private:
    // parameters
    std::shared_ptr<stattools::TParameterTyped<TypeF, 1>> _F;
    std::shared_ptr<stattools::TParameterTyped<TypeP, 1>> _p;
    std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>> _z;

    // dimensions
    size_t _numLoci;
    size_t _numSamples;

    // proposal probability for model switch F
    coretools::Probability _q_FModel_To_HWE;
    coretools::LogProbability _log_q_FModel_To_HWE;

    // propose new values after model switch
    coretools::StrictlyPositive<double> _lambdaNewF;

    // initial estimates for p
    std::vector<double>& _initialEstimatesP;

    // calculate LL
    coretools::LogProbability _calculateLLSumOverIndividuals(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t Locus, const genometools::THardyWeinbergGenotypeProbabilities & Probs) const;
    coretools::LogProbability _calculateLLSumOverIndividuals(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t Locus, double P) const;
    coretools::LogProbability _calculateLLSumOverIndividuals(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t Locus, double P, double F) const;

    // DAG
    void _registerPriorParameters() override;
    void _readCommandLineArguments();

    // initial values
    void _setInitialF();
    void _setInitialP();

    // log-densities
    double _getPriorDensity_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t index) const override;
    double _getLogPriorDensity_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t index) const override;
    double _getPriorDensityOld_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t index) const override;
    double _getLogPriorDensityOld_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t index) const override;
    double _getLogPriorRatio_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t index) const override;
    double _getExpectedValueFromPriorParameters(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & data, size_t index) const override;

    // update functions for F
    void _updateFAndZ(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data);
    void _updateRegularF(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data);
    void _updateFToHWE(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data);
    void _updateHWEToF(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data);
    double _calculateProbabilityOfProposingThisF(double F) const;
    double _getRandomNewF() const;

    // update functions for p
    void _updateP(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t Locus);

public:
    TInbreedingEstimatorPrior(std::shared_ptr<stattools::TParameterTyped<TypeF, 1>>  F,
                              std::shared_ptr<stattools::TParameterTyped<TypeP, 1>>  P,
                              std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>>  Z,
                              const std::vector<double> & InitialEstimatesP);
    ~TInbreedingEstimatorPrior() override = default;
    void initializeStorageOfPriorParameters() override;

    // estimate initial values
    void estimateInitialPriorParameters() override;

    // full log densities
    double getSumLogPriorDensity(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data) const override;

    // update all hyperprior parameters
    void updateParams() override;

    // simulate values
    void simulateUnderPrior() override;
};

//------------------------------------------
// TInbreedingEstimator
//------------------------------------------

class TInbreedingEstimator {
private:
    std::shared_ptr<stattools::TDAGBuilder> _dagBuilder;
    genometools::TPopulationLikelihoods<PhredType> _likelihoods;

    // filenames
    std::string _prefix;
    std::string _filename;

    // define DAG
    void _defineDAG();
    auto _defineFAndZ();
    auto _definePAndC();
    void _defineObservation(const std::shared_ptr<stattools::TParameterTyped<TypeF, 1>> & F, const std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>> & Z, const std::shared_ptr<stattools::TParameterTyped<TypeP, 1>> & P);

    // read input data
    void _readData();

public:
    TInbreedingEstimator(coretools::TParameters &Parameters, coretools::TLog *Logfile,
                                  coretools::TRandomGenerator *RandomGenerator);

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
        TInbreedingEstimator inbreedingEstimator(parameters(), &logfile(), &randomGenerator());
        inbreedingEstimator.runEstimation(parameters());
    };
};

}; // end namespace PopulationTools

#endif //ATLAS_TINBREEDINGESTIMATOR_H
