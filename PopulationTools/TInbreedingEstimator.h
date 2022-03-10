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
typedef genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> TypeGTL;

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

class InbreedingPrior : public stattools::MCMCPrior {
public:
    inline static const std::string inbreeding = "inbreeding";
};

class THardyWeinbergWithInbreedingGenotypeProbabilities : public genometools::THardyWeinbergGenotypeProbabilities {
public:
    THardyWeinbergWithInbreedingGenotypeProbabilities(const coretools::Probability p = 0.5, const coretools::Probability F = 0.0) noexcept{ set(p, F); }

    constexpr void set(const coretools::Probability p, const coretools::Probability F) noexcept {
        using BG = genometools::BiallelicGenotype;
        _probabilities[index(BG::haploidFirst)]  = p.complement(); // (1-p)
        _probabilities[index(BG::haploidSecond)] = p;              // p
        _probabilities[index(BG::homoFirst)]     = F.complement() * _probabilities[index(BG::haploidFirst)] * _probabilities[index(BG::haploidFirst)]
                                                        + F*_probabilities[index(BG::haploidFirst)]; // (1-F)*(1-p)^2 + F(1-p)
        _probabilities[index(BG::het)]           = F.complement() * 2.0 * _probabilities[index(BG::haploidFirst)] * _probabilities[index(BG::haploidSecond)]; // (1-F)2p(1-p)
        _probabilities[index(BG::homoSecond)]    = 1.0 - _probabilities[index(BG::homoFirst)] - _probabilities[index(BG::het)];
    }
};

class TInbreedingEstimatorPrior : public stattools::prior::TBaseNonIID<stattools::TValueFixed, TypeGTL, 2> {
private:
    // parameters
    std::shared_ptr<stattools::TParameterTyped<TypeF, 1>> _F;
    std::shared_ptr<stattools::TParameterTyped<TypeP, 1>> _p;
    std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>> _z;


    // proposal probability for model switch F
    double _q_FModel_To_HWE;
    double _log_q_FModel_To_HWE;

    // propose new values after model switch
    double _lambdaNewF;

    // observation
    TPopulationLikelihoods _likelihoods; // TODO: how to formalize those to observation class?

    // calculate LL
    static coretools::LogProbability _calculateLLSumOverIndividuals(TSampleLikelihoods* data, size_t N, const genometools::THardyWeinbergGenotypeProbabilities & Probs);
    static coretools::LogProbability _calculateLLSumOverIndividuals(TSampleLikelihoods* data, size_t N, double P);
    static coretools::LogProbability _calculateLLSumOverIndividuals(TSampleLikelihoods* data, size_t N, double P, double F);

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
    void _updateFAndZ();
    void _updateRegularF();
    void _updateFToHWE();
    void _updateHWEToF();
    double _calculateProbabilityOfProposingThisF(double F) const;
    double _getRandomNewF() const;

    // update functions for p
    void _updateP(TSampleLikelihoods* data, size_t Locus);

public:
    TInbreedingEstimatorPrior(const std::shared_ptr<stattools::TParameterTyped<TypeF, 1>> & F,
                              const std::shared_ptr<stattools::TParameterTyped<TypeP, 1>> & P,
                              const std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>> & Z);
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

    TPopulationLikelihoods _likelihoods;

    // filenames
    std::string _prefix;
    std::string _filename;

    // define DAG
    void _defineDAG();
    void _defineF();
    void _defineP();

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
