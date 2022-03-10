//
// Created by caduffm on 3/8/22.
//

#include "TInbreedingEstimator.h"

namespace PopulationTools {

// TODO:
//  - if F is in zero-model: don't add to meanVar; and don't adjust proposal kernel based on this

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

TInbreedingEstimatorPrior::TInbreedingEstimatorPrior(const std::shared_ptr<stattools::TParameterTyped<TypeF, 1>> &F,
                                                     const std::shared_ptr<stattools::TParameterTyped<TypeP, 1>> &P,
                                                     const std::shared_ptr<stattools::TParameterTyped<TypeZ, 1>> &Z) {
    // set parameters
    _F = F;
    _p = P;
    _z = Z;
    _q_FModel_To_HWE = 0.0;
    _log_q_FModel_To_HWE = 0.0;
    _lambdaNewF = 0.0;

    _registerPriorParameters();
}

void TInbreedingEstimatorPrior::_registerPriorParameters() {
    _priorParameters.push_back(_F);
    _priorParameters.push_back(_p);
    _priorParameters.push_back(_z);
}

void TInbreedingEstimatorPrior::initializeStorageOfPriorParameters() {
    // assemble dimension names from _likelihoods
    auto lociNames = std::make_shared<genometools::TNamesPositions>(_likelihoods.positions());

    // F and z: one value
    _F->initializeStorageSingleElementBasedOnPrior();
    _z->initializeStorageSingleElementBasedOnPrior();

    // p: linear of length L
    _p->initializeStorageBasedOnPrior({_likelihoods.getNumLoci()}, {lociNames});

    // read command line parameters
    _readCommandLineArguments();
}

void TInbreedingEstimatorPrior::_readCommandLineArguments(){
    using namespace coretools::instances;
    // parameters for model switch of F
    _q_FModel_To_HWE = parameters().getParameterWithDefault("probMovingToModelNoF", 0.1);
    if (_q_FModel_To_HWE < 0.0 || _q_FModel_To_HWE > 1.0){
        throw "Invalid argument 'probMovingToModelNoF': Must be in interval [0,1] (not " + coretools::str::toString(_q_FModel_To_HWE) + ")!";
    }
    _log_q_FModel_To_HWE = log(_q_FModel_To_HWE);
    logfile().list("Will propose move to model without F with probability ", _q_FModel_To_HWE, ". (parameter 'probMovingToModelNoF')");

    // Read lambda for proposing new F
    auto lambda = parameters().getParameterWithDefault<std::string>("lambdaF", "100.0");
    logfile().list("Lambda of exponential distribution used for the proposal of new F after move to model with F is set to ", lambda, ". (parameter 'lambdaF')");
}

double TInbreedingEstimatorPrior::_getLogPriorDensity_vec(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &, size_t) const {
    // TODO: Harmonize with observation!
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(TSampleLikelihoods* data, size_t N, const genometools::THardyWeinbergGenotypeProbabilities & Probs){
    // sum over all individuals of log sum_g P(d|g)P(g|p,F)
    coretools::LogProbability sum = 0.0;
    for (size_t i = 0; i < N; i++){
        sum += data[i].HWESum<coretools::LogProbability>(Probs);
    }
    return sum;
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(TSampleLikelihoods* data, size_t N, double P){
    return _calculateLLSumOverIndividuals(data, N, genometools::THardyWeinbergGenotypeProbabilities(P));
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(TSampleLikelihoods* data, size_t N, double P, double F){
    return _calculateLLSumOverIndividuals(data, N, THardyWeinbergWithInbreedingGenotypeProbabilities(P, F));
}

void TInbreedingEstimatorPrior::updateParams() {
    _updateFAndZ();

    size_t l = 0;
    for (_likelihoods.begin(); !_likelihoods.end(); _likelihoods.next(), ++l) {
        _updateP(_likelihoods.curData(), l);
    }
}

void TInbreedingEstimatorPrior::_updateFAndZ(){
    using namespace coretools::instances;
    if (_F->isUpdated()) {
        if (_z->value()) { // we're in F-model
            if (_z->isUpdated() && randomGenerator().getRand() < _q_FModel_To_HWE){
                _updateFToHWE(); // propose model switch
            } else {
                _updateRegularF(); // update F
            }
        } else { // we're in HWE-model
            if (_z->isUpdated()) {
                _updateHWEToF(); // propose model switch
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

void TInbreedingEstimatorPrior::_updateFToHWE() {
    // propose model switch: from F-Model to HWE-Model
    _F->set(0.0);
    _z->set(false);

    // calculate RJ-MCMC term of Hastings ratio
    double log_f_F = _calculateProbabilityOfProposingThisF(_F->oldValue()); // "compensating" the effect of no F in HWE model
    double logH = log_f_F - _log_q_FModel_To_HWE + _F->getLogPriorDensity(); // TODO: If _F didn't have a uniform prior, how would this go into Hastings ratio in model switch?

    // calculate LL Hastings-Ratio
    size_t l = 0;
    for (_likelihoods.begin(); !_likelihoods.end(); _likelihoods.next(), ++l){
        logH +=  _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l)) // HWE model
                 - _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l), _F->oldValue()); // F-Model
    }

    if (!_F->acceptOrReject(logH)){
        // rejected -> re-set model
        _z->set(true);
    }
}

void TInbreedingEstimatorPrior::_updateRegularF() {
    // update F without model switch: F -> F'
    if (_F->update()) {
        double logH = _F->getLogPriorRatio();
        size_t l = 0;
        for (_likelihoods.begin(); !_likelihoods.end(); _likelihoods.next(), ++l) {
            logH += _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l), _F->value())
                    - _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l), _F->oldValue());
        }
        _F->acceptOrReject(logH);
    }
}

void TInbreedingEstimatorPrior::_updateHWEToF(){
    // propose model switch: from HWE-Model to F-Model
    _F->set(_getRandomNewF());
    _z->set(true);

    // calculate RJ-MCMC term of Hastings ratio
    double log_f_F = _calculateProbabilityOfProposingThisF(_F->value()); // "compensating" the effect of no F in HWE model
    double logH = _log_q_FModel_To_HWE - log_f_F + _F->getLogPriorRatio(); // TODO: If _F didn't have a uniform prior, how would this go into Hastings ratio in model switch?

    // calculate LL Hastings-Ratio
    size_t l = 0;
    for (_likelihoods.begin(); !_likelihoods.end(); _likelihoods.next(), ++l){
        logH +=  _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l), _F->value()) // F-Model
                 - _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l)); // HWE-Model
    }

    if (!_F->acceptOrReject(logH)){
        // rejected -> re-set modelF
        _z->set(false);
    }
}

void TInbreedingEstimatorPrior::_updateP(TSampleLikelihoods* data, size_t Locus){
    if (_p->updateSpecificIndex(Locus)){
        double logH = _p->getLogPriorRatio();
        if (_z->value()){ // F-Model
            logH += _calculateLLSumOverIndividuals(data, _likelihoods.curSampleSize(), _p->value(Locus), _F->value()) - _calculateLLSumOverIndividuals(data, _likelihoods.curSampleSize(), _p->oldValue(Locus), _F->value());
        } else { // HWE-Model
            logH += _calculateLLSumOverIndividuals(data, _likelihoods.curSampleSize(), _p->value(Locus)) - _calculateLLSumOverIndividuals(data, _likelihoods.curSampleSize(), _p->oldValue(Locus));
        }

        _p->acceptOrReject(Locus, logH);
    }
}

void TInbreedingEstimatorPrior::_setInitialF(){
    if (_F->hasFixedInitialValue()){
        _F->value() == 0.0 ? _z->set(false) : _z->set(true);
    } else {
        // F: decide in which model to start
        using namespace coretools::instances;
        if (parameters().parameterExists("startInZeroModel")){
            _z->set(false);
            _F->set(0.0);
        } else {
            _z->set(true);
            _F->set(_getRandomNewF());
        }
    }
}

void TInbreedingEstimatorPrior::_setInitialP(){
    using namespace coretools::instances;
    std::vector<double> freqs;
    if(parameters().parameterExists("trueAlleleFreq")){
        freqs = _likelihoods.donateTrueAlleleFrequencies();
        if (freqs.size() != _likelihoods.getNumLoci()){
            throw coretools::str::toString("Failed to initialize p with true allele frequencies (parameter 'trueAlleleFreq'): Did not receive one allele frequency per locus! Number of loci from VCF: ", _likelihoods.getNumLoci(), " vs number of allele frequencies ", freqs.size(), ")");
        }
        logfile().list("Initializing allele frequencies to true values read from trueAlleleFreq file");
        logfile().warning("This task is not implemented for multiple populations! Considering all samples to be from one population.");
    } else {
        freqs = _likelihoods.donateAlleleFrequencies();
        logfile().list("Initializing allele frequencies to values estimated from sample genotype likelihoods");
        logfile().warning("This task is not implemented for multiple populations! Considering all samples to be from one population.");
    }

    // now read values into _p
    if (!_p->hasFixedInitialValue()){
        for (size_t l = 0; l < _likelihoods.getNumLoci(); l++){
            _p->set(l, freqs[l]);
        }
    }
}

void TInbreedingEstimatorPrior::estimateInitialPriorParameters() {
    _setInitialF();
    _setInitialP();
}

double TInbreedingEstimatorPrior::getSumLogPriorDensity(const std::shared_ptr<const stattools::TParameterObservationTypedBase<stattools::TValueFixed, TypeGTL, 2>> &) const {
    double sum = 0.;
    size_t l = 0;
    for (_likelihoods.begin(); !_likelihoods.end(); _likelihoods.next(), ++l){
        if (_z->value()){
            sum += _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l), _F->value()); // F-Model
        } else {
            sum += _calculateLLSumOverIndividuals(_likelihoods.curData(), _likelihoods.curSampleSize(), _p->value(l)); // HWE-Model
        }
    }
    return sum;
}

void TInbreedingEstimatorPrior::simulateUnderPrior() {
    // TODO: ? In here or in Atlas simulator?
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

void TInbreedingEstimator::_defineF(){
    auto priorOnF = std::make_shared<stattools::prior::TUniformFixed<stattools::TValueUpdated, TypeF, 1>>();
    stattools::TParameterDefinition defF;
    defF.setAllFiles(_filename);
    auto F = std::make_shared<stattools::TParameterTyped<TypeF, 1>>("F", priorOnF, defF);
    _dagBuilder->addToDAG(F);
}

void TInbreedingEstimator::_defineP(){
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
}

void TInbreedingEstimator::_readData(){
    //read data
    using namespace coretools::instances;
    _likelihoods.doSaveAlleleFrequencies();
    if(parameters().parameterExists("trueAlleleFreq")){
        _likelihoods.doSaveTrueAlleleFrequencies();
        logfile().list("Will save true allele frequencies in population_likelihoods object");
    }
    _likelihoods.readData(parameters(), &logfile());
}

void TInbreedingEstimator::_defineDAG(){
    _readData();
    _defineF();
    _defineP();

    // TODO: Create TInbreedingPrior to build the complete DAG!
}

void TInbreedingEstimator::runEstimation(coretools::TParameters &Parameters) {
    // read file names
    std::string vcfFileName = _likelihoods.getVCFName();
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