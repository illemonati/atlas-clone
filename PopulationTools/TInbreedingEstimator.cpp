//
// Created by caduffm on 3/8/22.
//

#include "TInbreedingEstimator.h"

#include <utility>

namespace PopulationTools {

// TODO:
//  - if F is in zero-model: don't add to meanVar; and don't adjust proposal kernel based on this

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

TInbreedingEstimatorPrior::TInbreedingEstimatorPrior(std::shared_ptr<stattools::TParameterTyped<TypeF, 1>> F,
                                                     std::shared_ptr<stattools::TParameterTyped<TypeP, 1>> P,
                                                     std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>> Z,
                                                     const std::vector<double> & InitialEstimatesP) :
                                                     _F(std::move(F)), _p(std::move(P)), _z(std::move(Z)), _initialEstimatesP(InitialEstimatesP){
    // set parameters
    _q_FModel_To_HWE = 0.0;
    _log_q_FModel_To_HWE = 0.0;
    _lambdaNewF = 0.0;
    _numLoci = 0;
    _numSamples = 0;

    _registerPriorParameters();
}

void TInbreedingEstimatorPrior::_registerPriorParameters() {
    _priorParameters.push_back(_F);
    _priorParameters.push_back(_p);
    _priorParameters.push_back(_z);
}

void TInbreedingEstimatorPrior::initializeStorageOfPriorParameters() {
    // get shared ptr to data
    auto data = this->_getSharedPtrParam(this->_parametersBelow[0]);
    const auto & lociNames = data->getDimensionName(0);

    // set dimensions
    _numLoci = data->dimensions()[0];
    _numSamples = data->dimensions()[1];

    // F and z: one value
    _F->initializeStorageSingleElementBasedOnPrior();
    _z->initializeStorageSingleElementBasedOnPrior();

    // p: linear of length L
    _p->initializeStorageBasedOnPrior({_numLoci}, {lociNames});

    // read command line parameters
    _readCommandLineArguments();
}

void TInbreedingEstimatorPrior::_readCommandLineArguments(){
    using namespace coretools::instances;
    // parameters for model switch of F
    _q_FModel_To_HWE = parameters().getParameterWithDefault("probMovingToModelNoF", 0.1);
    _log_q_FModel_To_HWE = log(_q_FModel_To_HWE);
    logfile().list("Will propose move to model without F with probability ", _q_FModel_To_HWE, ". (parameter 'probMovingToModelNoF')");

    // Read lambda for proposing new F
    _lambdaNewF = parameters().getParameterWithDefault("lambdaF", 100.0);
    logfile().list("Lambda of exponential distribution used for the proposal of new F after move to model with F is set to ", coretools::str::toString(_lambdaNewF), ". (parameter 'lambdaF')");
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t Locus, const genometools::THardyWeinbergGenotypeProbabilities & Probs) const {
    // sum over all individuals of log sum_g P(d|g)P(g|p,F)
    coretools::LogProbability sum = 0.0;
    for (size_t i = 0; i < _numSamples; i++){
        const auto linearIndex = Data->getIndex({Locus, i});
        sum += Data->value(linearIndex).HWESum<coretools::LogProbability>(Probs);
    }
    return sum;
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &Data, size_t Locus, double P) const {
    return _calculateLLSumOverIndividuals(Data, Locus, genometools::THardyWeinbergGenotypeProbabilities(P));
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &Data, size_t Locus, double P, double F) const {
    return _calculateLLSumOverIndividuals(Data, Locus, genometools::THardyWeinbergWithInbreedingGenotypeProbabilities(P, F));
}

void TInbreedingEstimatorPrior::updateParams() {
    // get shared ptr to data
    auto data = this->_getSharedPtrParam(this->_parametersBelow[0]);

    _updateFAndZ(data);

    for (size_t l = 0; l < _numLoci; l++) {
        _updateP(data, l);
    }
}

void TInbreedingEstimatorPrior::_updateFAndZ(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data){
    using namespace coretools::instances;
    if (_F->isUpdated()) {
        if (_z->value()) { // we're in F-model
            if (_z->isUpdated() && randomGenerator().getRand() < _q_FModel_To_HWE){
                _updateFToHWE(Data); // propose model switch
            } else {
                _updateRegularF(Data); // update F
            }
        } else { // we're in HWE-model
            if (_z->isUpdated()) {
                _updateHWEToF(Data); // propose model switch
            }
        }
    }
}

double TInbreedingEstimatorPrior::_calculateProbabilityOfProposingThisF(double F) const{
    // When switching from HWE to F-Model: propose values from truncated exponential distribution with a certain lambdaF
    // -> in Hastings ratio, for the HWE-Model (without F), we need to know the probability of picking exactly that F
    // corresponds to log f(F) in Hastings ratios
    return coretools::TExponentialDistr::logDensity(F, _lambdaNewF);
}

double TInbreedingEstimatorPrior::_getRandomNewF() const{
    using namespace coretools::instances;
    return randomGenerator().getExponentialRandomTruncated(_lambdaNewF, 0.0, 1.0);
}

// TODO: Limit analysis to a single population?

void TInbreedingEstimatorPrior::_updateFToHWE(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data) {
    // propose model switch: from F-Model to HWE-Model
    _F->set(0.0);
    _z->set(false);

    // calculate RJ-MCMC term of Hastings ratio
    double log_f_F = _calculateProbabilityOfProposingThisF(_F->oldValue());
    double logH = log_f_F - _log_q_FModel_To_HWE - _F->getLogPriorDensityOld();

    // calculate LL Hastings-Ratio
    for (size_t l = 0; l < _numLoci; l++){
        logH +=  _calculateLLSumOverIndividuals(Data, l, _p->value(l)) // HWE model
                 - _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->oldValue()); // F-Model
    }

    if (_F->acceptOrReject(logH, false)){
        // accepted model switch -> we're now in null model
        // -> stop adding to posterior mean/var
        _F->doTrackMeanVar(false);
    } else {
        // rejected model switch -> we're still in F-model
        _z->set(true);
    }
}

void TInbreedingEstimatorPrior::_updateRegularF(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data) {
    // update F without model switch: F -> F'
    if (_F->update()) {
        double logH = _F->getLogPriorRatio();
        for (size_t l = 0; l < _numLoci; l++){
            logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->value())
                    - _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->oldValue());
        }
        _F->acceptOrReject(logH);
    }
}

void TInbreedingEstimatorPrior::_updateHWEToF(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data){
    // propose model switch: from HWE-Model to F-Model
    _F->set(_getRandomNewF());
    _z->set(true);

    // calculate RJ-MCMC term of Hastings ratio
    double log_f_F = _calculateProbabilityOfProposingThisF(_F->value()); // "compensating" the effect of no F in HWE model
    double logH = _log_q_FModel_To_HWE + _F->getLogPriorDensity() - log_f_F;

    // calculate LL Hastings-Ratio
    for (size_t l = 0; l < _numLoci; l++){
        logH +=  _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->value()) // F-Model
                 - _calculateLLSumOverIndividuals(Data, l, _p->value(l)); // HWE-Model
    }

    if (_F->acceptOrReject(logH, false)){
        // accepted model switch -> we're now in F-model
        // -> start adding to posterior mean/var
        _F->doTrackMeanVar(true); // F is in non-zero model -> add to posterior mean/var
    } else {
        // rejected model switch -> we're still in null-model
        _z->set(false);
    }
}

void TInbreedingEstimatorPrior::_updateP(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t Locus){
    if (_p->updateSpecificIndex(Locus)){
        double logH = _p->getLogPriorRatio();
        if (_z->value()){ // F-Model
            logH += _calculateLLSumOverIndividuals(Data, Locus, _p->value(Locus), _F->value())
                    - _calculateLLSumOverIndividuals(Data, Locus, _p->oldValue(Locus), _F->value());
        } else { // HWE-Model
            logH += _calculateLLSumOverIndividuals(Data, Locus, _p->value(Locus))
                    - _calculateLLSumOverIndividuals(Data, Locus, _p->oldValue(Locus));
        }

        _p->acceptOrReject(Locus, logH);
    }
}

void TInbreedingEstimatorPrior::_setInitialF(){
    // set initial F
    if (!_F->hasFixedInitialValue()){
        if (_z->hasFixedInitialValue() && _z->value() == 0){
            // user wants to start in zero-model
            _F->set(0.0);
        } else {
            _F->set(_getRandomNewF());
        }
    }
    // set initial z
    _F->value() == 0.0 ? _z->set(false) : _z->set(true);
}

void TInbreedingEstimatorPrior::_setInitialP(){
    using namespace coretools::instances;
    logfile().list("Initializing allele frequencies to values estimated from sample genotype likelihoods");
    logfile().warning("This task is not implemented for multiple populations! Considering all samples to be from one population.");

    // now read values into _p
    if (!_p->hasFixedInitialValue()){
        for (size_t l = 0; l < _numLoci; l++){
            _p->set(l, _initialEstimatesP[l]);
        }
    }
    // clear initialEstimatesP
    _initialEstimatesP.clear();
}

void TInbreedingEstimatorPrior::estimateInitialPriorParameters() {
    _setInitialF();
    _setInitialP();
}

double TInbreedingEstimatorPrior::getSumLogPriorDensity(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data) const {
    double sum = 0.;
    for (size_t l = 0; l < _numLoci; l++){
        if (_z->value()){
            sum += _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->value()); // F-Model
        } else {
            sum += _calculateLLSumOverIndividuals(Data, l, _p->value(l)); // HWE-Model
        }
    }
    return sum;
}

void TInbreedingEstimatorPrior::simulateUnderPrior() {
    // TODO: ? In here or in Atlas simulator?
}

double TInbreedingEstimatorPrior::_getLogPriorDensity_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> & Data, size_t Index) const {
    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": not implemented.");
}
double TInbreedingEstimatorPrior::_getPriorDensity_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &, size_t) const {
    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": not implemented.");
}
double TInbreedingEstimatorPrior::_getPriorDensityOld_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &, size_t) const {
    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": not implemented.");
}
double TInbreedingEstimatorPrior::_getLogPriorDensityOld_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &, size_t) const {
    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": not implemented.");
}
double TInbreedingEstimatorPrior::_getLogPriorRatio_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &, size_t) const {
    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": not implemented.");
}
double TInbreedingEstimatorPrior::_getExpectedValueFromPriorParameters(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &, size_t) const{
    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": not implemented.");
}

//------------------------------------------
// TInbreedingEstimator
//------------------------------------------

TInbreedingEstimator::TInbreedingEstimator(coretools::TParameters &, coretools::TLog *, coretools::TRandomGenerator *) {}

auto TInbreedingEstimator::_defineFAndZ(){
    auto priorOnF = std::make_shared<stattools::prior::TUniformFixed<stattools::TValueUpdated, TypeF, 1>>();
    stattools::TParameterDefinition defF;
    defF.setAllFiles(_filename);
    auto F = std::make_shared<stattools::TParameterTyped<TypeF, 1>>("F", priorOnF, defF);
    _dagBuilder->addToDAG(F);

    auto priorOnZ = std::make_shared<stattools::prior::TUniformFixed<stattools::TValueUpdated, TypeZ, 1>>();
    stattools::TParameterDefinition defZ;
    defZ.setAllFiles(_filename);
    auto z = std::make_shared<stattools::TParameterTyped<TypeZ, 1>>("z", priorOnZ, defZ);
    _dagBuilder->addToDAG( z);

    return std::make_tuple(F, z);
}

auto TInbreedingEstimator::_definePAndC(){
    // c
    auto priorOnC = std::make_shared<stattools::prior::TUniformFixed<stattools::TValueUpdated, TypeC, 1>>();
    stattools::TParameterDefinition defC;
    defC.setAllFiles(_filename);
    auto c = std::make_shared<stattools::TParameterTyped<TypeC, 1>>("c", priorOnC, defC);
    _dagBuilder->addToDAG(c);

    // p
    auto priorOnP = std::make_shared<stattools::prior::TUniformFixed<stattools::TValueUpdated, TypeP, 1>>(); // std::make_shared<stattools::prior::TSymmetricBetaInferred<TypeP, 1, TypeC>>(c);
    stattools::TParameterDefinition defP;
    defP.setAllFiles(_filename);
    defP.setJumpSizeForAll(false);
    auto p = std::make_shared<stattools::TParameterTyped<TypeP, 1>>("p", priorOnP, defP);
    _dagBuilder->addToDAG(p);

    return p;
}

void TInbreedingEstimator::_defineObservation(const std::shared_ptr<stattools::TParameterTyped<TypeF, 1>> & F,
                                              const std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>> & Z,
                                              const std::shared_ptr<stattools::TParameterTyped<TypeP, 1>> & P){
    auto inbreedingPrior = std::make_shared<TInbreedingEstimatorPrior>(F, P, Z, _likelihoods.alleleFrequencies());
    stattools::TObservationDefinition def;
    auto observation = std::make_shared<stattools::TObservationTyped<TypeGTL, 2>>("genotypeLikelihoods", inbreedingPrior, def);
    observation->setStorage(_likelihoods.getStorage());
    _dagBuilder->addToDAG(observation);
}

void TInbreedingEstimator::_readData(){
    //read data
    using namespace coretools::instances;
    _likelihoods.doSaveAlleleFrequencies();
    if(parameters().parameterExists("trueAlleleFreq")){
        _likelihoods.doSaveTrueAlleleFrequencies();
        logfile().list("Will save true allele frequencies for population likelihoods.");
    }
    _likelihoods.readData(parameters(), &logfile());
}

void TInbreedingEstimator::_defineDAG(){
    _readData();

    // define F
    auto [F, z] = _defineFAndZ();

    // define C and P
    auto p = _definePAndC();

    // define observation
    _defineObservation(F, z, p);
}

void TInbreedingEstimator::runEstimation(coretools::TParameters &Parameters) {
    // read file names
    std::string vcfFileName = Parameters.getParameter("vcf");
    vcfFileName = coretools::str::extractBeforeLast(vcfFileName, ".vcf");
    _prefix = Parameters.getParameterWithDefault<std::string>("prefix", "");
    _filename = Parameters.getParameterWithDefault<std::string>("out", vcfFileName);

    // build DAG
    _dagBuilder = std::make_shared<stattools::TDAGBuilder>();
    _defineDAG();
    _dagBuilder->buildDAG();

    // run MCMC
    stattools::TMCMC mcmc;
    mcmc.runMCMC(_prefix, _dagBuilder);
}


} // end namespace PopulationTools