//
// Created by caduffm on 3/8/22.
//

#include "TInbreedingEstimator.h"

#include <array>
#include <exception>
#include <math.h>
#include <stdexcept>

#include "TDAGBuilder.h"
#include "TMCMC.h"
#include "TParameterObservationTypedBase.h"
#include "TPriorBernouilli.h"
#include "TPriorBeta.h"
#include "TPriorUniform.h"
#include "TRandomGenerator.h"
#include "mathFunctions.h"
#include "stringFunctions.h"

namespace stattools {
class TParameterBase;
}
namespace stattools {
class TParameterObservationBase;
}

namespace PopulationTools {

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

TInbreedingEstimatorPrior::TInbreedingEstimatorPrior(stattools::TParameterTyped<TypeF, 1> *F,
                                                     stattools::TParameterTyped<TypeP, 1> *P,
                                                     stattools::TParameterTyped<TypeFModel, 1> *FModel,
                                                     stattools::TParameterTyped<TypePModel, 1> *PModel,
                                                     const std::vector<double> &InitialEstimatesP)
    : _F(F), _p(P), _FModel(FModel), _pModel(PModel), _initialEstimatesP(InitialEstimatesP) {
	_priorParameters.push_back(_F);
	_priorParameters.push_back(_p);
	_priorParameters.push_back(_FModel);
	_priorParameters.push_back(_pModel);

	this->_flagPriorParameters();
}

std::string TInbreedingEstimatorPrior::name() const { return "inbreeding"; }

void TInbreedingEstimatorPrior::initializeStorageOfPriorParameters() {
	// get shared ptr to data
	auto data             = this->_storageBelow[0];
	const auto &lociNames = data->getDimensionName(0);

	// set dimensions
	_numLoci    = data->dimensions()[0];
	_numSamples = data->dimensions()[1];

	// F and FModel: one value
	_F->initializeStorageSingleElementBasedOnPrior();
	_FModel->initializeStorageSingleElementBasedOnPrior();

	// p and pModel: linear of length L
	_p->initializeStorageBasedOnPrior({_numLoci}, {lociNames});
	_pModel->initializeStorageBasedOnPrior({_numLoci}, {lociNames});

	// read command line parameters
	_readCommandLineArguments();
}

void TInbreedingEstimatorPrior::_readCommandLineArguments() {
	using namespace coretools::instances;
	// parameters for model switch of F
	_q_FModel_To_HWE     = parameters().getParameterWithDefault("probMovingToModelNoF", 0.1);
	_log_q_FModel_To_HWE = log(_q_FModel_To_HWE);
	logfile().list("Will propose move to model without F with probability ", _q_FModel_To_HWE,
	               ". (parameter 'probMovingToModelNoF')");

	// Read lambda for proposing new F
	_lambdaNewF = parameters().getParameterWithDefault("lambdaF", 100.0);
	logfile().list("Lambda of exponential distribution used for the proposal of "
	               "new F after move to model with F is set to ",
	               coretools::str::toString(_lambdaNewF), ". (parameter 'lambdaF')");

	// parameters for model switch of p
	_q_PModel_To_NullModel     = parameters().getParameterWithDefault("probMovingToModelP0", 0.1);
	_log_q_PModel_To_NullModel = log(_q_PModel_To_NullModel);
	logfile().list("Will propose move to monomorphic model with probability ", _q_PModel_To_NullModel,
	               ". (parameter 'probMovingToModelP0')");

	// Read lambda for proposing new p
	_lambdaNewP = parameters().getParameterWithDefault("lambdaP", 100.0);
	logfile().list("Lambda of exponential distribution used for the proposal of "
	               "new p after move to polymorphic model is set to ",
	               coretools::str::toString(_lambdaNewP), ". (parameter 'lambdaP')");
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(
    const Storage &Data, size_t Locus, const genometools::THardyWeinbergGenotypeProbabilities &Probs) const {
	// sum over all individuals of log sum_g P(d|g)P(g|p,F)
	coretools::LogProbability sum = 0.0;
	for (size_t i = 0; i < _numSamples; i++) {
		const auto linearIndex = Data.getIndex({Locus, i});
		sum += Data[linearIndex].HWESum<coretools::LogProbability>(Probs);
	}
	return sum;
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
                                                                                    coretools::Probability P) const {
	return _calculateLLSumOverIndividuals(Data, Locus, genometools::THardyWeinbergGenotypeProbabilities(P));
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
                                                                                    coretools::Probability P,
                                                                                    coretools::Probability F) const {
	return _calculateLLSumOverIndividuals(Data, Locus,
	                                      genometools::THardyWeinbergWithInbreedingGenotypeProbabilities(P, F));
}

void TInbreedingEstimatorPrior::updateParams() {
	// get ptr to data
	auto data = *this->_storageBelow[0];

	_updateFAndFModel(data);

	for (size_t l = 0; l < _numLoci; l++) { _updateP(data, l); }
}

void TInbreedingEstimatorPrior::_updateFAndFModel(const Storage &Data) {
	using namespace coretools::instances;
	if (_F->isUpdated()) {
		if (_FModel->value()) { // we're in F-model
			if (_FModel->isUpdated() && randomGenerator().getRand() < _q_FModel_To_HWE) {
				_updateFToHWE(Data); // propose model switch
			} else {
				_updateRegularF(Data); // update F
			}
		} else { // we're in HWE-model
			if (_FModel->isUpdated()) {
				_updateHWEToF(Data); // propose model switch
			}
		}
	}
}

double
TInbreedingEstimatorPrior::_calculateProbabilityOfProposingThisFOrP(double Value,
                                                                    coretools::StrictlyPositive<double> Lambda) const {
	// When switching from One to Null-Model: propose values from truncated
	// exponential distribution with a certain lambda
	// -> in Hastings ratio, for the Null-Model, we need to know the
	// probability of picking exactly that Value corresponds to log f(Value) in Hastings
	// ratios
	// applies to model switches of both F and p
	return coretools::TExponentialDistr::logDensity(Value, Lambda);
}

double TInbreedingEstimatorPrior::_getRandomNewF() const {
	using namespace coretools::instances;
	return randomGenerator().getExponentialRandomTruncated(_lambdaNewF, 0.0, 1.0);
}

double TInbreedingEstimatorPrior::_getRandomNewP() const {
	using namespace coretools::instances;
	return randomGenerator().getExponentialRandomTruncated(_lambdaNewP, 0.0, 1.0);
}

void TInbreedingEstimatorPrior::_updateFToHWE(const Storage &Data) {
	// propose model switch: from F-Model to HWE-Model
	_F->set(0.0);
	_FModel->set(false);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_F = _calculateProbabilityOfProposingThisFOrP(_F->oldValue(), _lambdaNewF);
	double logH    = log_f_F - _log_q_FModel_To_HWE - _F->getLogPriorDensityOld();

	// calculate LL Hastings-Ratio
	for (size_t l = 0; l < _numLoci; l++) {
		if (_pModel->value(l)) {                                          // only if polyploid
			logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l)) // HWE model
			        - _calculateLLSumOverIndividuals(Data, l, _p->value(l),
			                                         _F->oldValue()); // F-Model
		}
	}

	if (_F->acceptOrReject_DontCountRejection(logH)) {
		// accepted model switch -> we're now in null model
		// -> stop adding to posterior mean/var
		_F->doTrackMeanVar(false);
	} else {
		// rejected model switch -> we're still in F-model
		_FModel->set(true);
	}
}

void TInbreedingEstimatorPrior::_updateRegularF(const Storage &Data) {
	// update F without model switch: F -> F'
	if (_F->update()) {
		double logH = _F->getLogPriorRatio();
		for (size_t l = 0; l < _numLoci; l++) {
			if (_pModel->value(l)) { // only if polyploid
				logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->value()) -
				        _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->oldValue());
			}
		}
		_F->acceptOrReject(logH);
	}
}

void TInbreedingEstimatorPrior::_updateHWEToF(const Storage &Data) {
	// propose model switch: from HWE-Model to F-Model
	_F->set(_getRandomNewF());
	_FModel->set(true);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_F = _calculateProbabilityOfProposingThisFOrP(
	    _F->value(), _lambdaNewF); // "compensating" the effect of no F in HWE model
	double logH = _log_q_FModel_To_HWE + _F->getLogPriorDensity() - log_f_F;

	// calculate LL Hastings-Ratio
	for (size_t l = 0; l < _numLoci; l++) {
		if (_pModel->value(l)) { // only if polyploid
			logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l),
			                                       _F->value())              // F-Model
			        - _calculateLLSumOverIndividuals(Data, l, _p->value(l)); // HWE-Model
		}
	}

	if (_F->acceptOrReject_DontCountRejection(logH)) {
		// accepted model switch -> we're now in F-model
		// -> start adding to posterior mean/var
		_F->doTrackMeanVar(true); // F is in non-zero model -> add to posterior mean/var
	} else {
		// rejected model switch -> we're still in null-model
		_FModel->set(false);
	}
}

void TInbreedingEstimatorPrior::_updateP(const Storage &Data, size_t Locus) {
	using namespace coretools::instances;
	if (_p->isUpdated()) {
		if (_pModel->value(Locus)) { // we're in the polyploid model
			if (_pModel->isUpdated() && randomGenerator().getRand() < _q_PModel_To_NullModel) {
				_updatePToNull(Data, Locus); // propose model switch
			} else {
				_updateRegularP(Data, Locus); // update p
			}
		} else { // we're in null-model
			if (_pModel->isUpdated()) {
				_updateNullToP(Data, Locus); // propose model switch
			}
		}
	}
}

double TInbreedingEstimatorPrior::_calculateLLRatio_UpdateP(const Storage &Data, size_t Locus) {
	if (_FModel->value()) { // F-Model
		return _calculateLLSumOverIndividuals(Data, Locus, _p->value(Locus), _F->value()) -
		       _calculateLLSumOverIndividuals(Data, Locus, _p->oldValue(Locus), _F->value());
	} else { // HWE-Model
		return _calculateLLSumOverIndividuals(Data, Locus, _p->value(Locus)) -
		       _calculateLLSumOverIndividuals(Data, Locus, _p->oldValue(Locus));
	}
}

void TInbreedingEstimatorPrior::_updatePToNull(const Storage &Data, size_t Locus) {
	// propose model switch: from p-Model to null-Model
	_p->set(Locus, 0.0);
	_pModel->set(Locus, false);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_p = _calculateProbabilityOfProposingThisFOrP(_p->oldValue(Locus), _lambdaNewP);
	double logH    = log_f_p - _log_q_PModel_To_NullModel - _p->getLogPriorDensityOld(Locus);
	logH += _pModel->getLogPriorRatio(Locus);
	logH += _calculateLLRatio_UpdateP(Data, Locus);

	if (_p->acceptOrReject_DontCountRejection(Locus, logH)) {
		// accepted model switch -> we're now in null model
		// -> stop adding to posterior mean/var
		_p->doTrackMeanVar(false); // TODO: should be specific for one locus
	} else {
		// rejected model switch -> we're still in P-model
		_pModel->set(Locus, true);
	}
}

void TInbreedingEstimatorPrior::_updateRegularP(const Storage &Data, size_t Locus) {
	if (_p->updateSpecificIndex(Locus)) {
		const double logH = _calculateLLRatio_UpdateP(Data, Locus) + _p->getLogPriorRatio(Locus);
		_p->acceptOrReject(Locus, logH);
	}
}

void TInbreedingEstimatorPrior::_updateNullToP(const Storage &Data, size_t Locus) {
	// propose model switch: from null-Model to p-Model
	_p->set(Locus, _getRandomNewP());
	_pModel->set(Locus, true);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_p = _calculateProbabilityOfProposingThisFOrP(
	    _p->value(Locus), _lambdaNewP); // "compensating" the effect of no p in Null-model
	double logH = _log_q_PModel_To_NullModel + _p->getLogPriorDensity(Locus) - log_f_p;
	logH += _pModel->getLogPriorRatio(Locus);
	logH += _calculateLLRatio_UpdateP(Data, Locus);

	if (_p->acceptOrReject_DontCountRejection(Locus, logH)) {
		// accepted model switch -> we're now in P-model
		// -> start adding to posterior mean/var
		_p->doTrackMeanVar(true); // p is in non-zero model -> add to posterior mean/var
	} else {
		// rejected model switch -> we're still in null-model
		_pModel->set(Locus, false);
	}
}

// todo: posterior mean of p: average over models; also for F

void TInbreedingEstimatorPrior::_setInitialF() {
	// set initial F
	if (!_F->hasFixedInitialValue()) {
		if (_FModel->hasFixedInitialValue() && _FModel->value() == 0) {
			// user wants to start in zero-model
			_F->set(0.0);
		} else {
			_F->set(_getRandomNewF());
		}
	}
	// set initial z
	_F->value() == 0.0 ? _FModel->set(false) : _FModel->set(true);
}

void TInbreedingEstimatorPrior::_setInitialPAndPModel() {
	using namespace coretools::instances;
	logfile().list("Initializing allele frequencies to values estimated from "
	               "sample genotype likelihoods");

	// now read values into _p and _pModel
	if (!_p->hasFixedInitialValue()) {
		for (size_t l = 0; l < _numLoci; l++) {
			if (_pModel->hasFixedInitialValue()) {
				if (_pModel->value(l) == 0) {
					_p->set(l, 0.0);
				} else {
					_p->set(l, std::max(_initialEstimatesP[l], TypeP::min().get())); // start in 1-model: prevent p = 0
				}
			} else {
				if (_initialEstimatesP[l] == 0.0) {
					_p->set(l, 0.0);
					_pModel->set(l, false);
				} else {
					_p->set(l, _initialEstimatesP[l]);
					_pModel->set(l, true);
				}
			}
		}
	}
}

void TInbreedingEstimatorPrior::estimateInitialPriorParameters() {
	_setInitialF();
	_setInitialPAndPModel();
}

double TInbreedingEstimatorPrior::getSumLogPriorDensity(const Storage &Data) const {
	double sum = 0.;
	for (size_t l = 0; l < _numLoci; l++) {
		if (_FModel->value()) {
			sum += _calculateLLSumOverIndividuals(Data, l, _p->value(l),
			                                      _F->value()); // F-Model
		} else {
			sum += _calculateLLSumOverIndividuals(Data, l, _p->value(l)); // HWE-Model
		}
	}
	return sum;
}

void TInbreedingEstimatorPrior::_simulateUnderPrior(Storage *) {
	throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": Not implemented. Use Atlas simulator instead.");
}

//------------------------------------------
// TInbreedingEstimatorModel
//------------------------------------------

TInbreedingEstimatorModel::TInbreedingEstimatorModel(
    const std::string &Filename, std::shared_ptr<stattools::TDAGBuilder> &DAGBuilder,
    const genometools::TPopulationLikelihoods<stattools::TValueFixed<TypeGTL>> &Likelihoods)
    : _F("F", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeF, 1>>(), {Filename}),
      _FModel("FModel", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeFModel, 1>>(),
              {Filename}),
      _pi("pi", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypePi, 1>>(), {Filename}),
      _pModel("pModel",
              std::make_shared<stattools::prior::TBernouilliInferred<stattools::TParameterBase, TypePModel, 1, TypePi>>(
                  &_pi),
              {Filename}),
      _log_gamma("log_gamma",
                 std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogGamma, 1>>(),
                 {Filename}),
      _p("p",
         std::make_shared<stattools::prior::TBetaSymmetricZeroMixtureInferred<stattools::TParameterBase, TypeP, 1,
                                                                              TypeLogGamma, TypePModel, true>>(
             &_log_gamma, &_pModel),
         {Filename}),
      _observation(
          "genotypeLikelihoods",
          std::make_shared<TInbreedingEstimatorPrior>(&_F, &_p, &_FModel, &_pModel, Likelihoods.alleleFrequencies()),
          {}) {

	_p.getDefinition().setJumpSizeForAll(false);

	DAGBuilder->addToDAG(_F);
	DAGBuilder->addToDAG(_FModel);
	DAGBuilder->addToDAG(_pi);
	DAGBuilder->addToDAG(_pModel);
	DAGBuilder->addToDAG(_log_gamma);
	DAGBuilder->addToDAG(_p);
	DAGBuilder->addToDAG(_observation);

	// set storage
	_observation.storage() = Likelihoods.getStorage();
}

//------------------------------------------
// TInbreedingEstimator
//------------------------------------------

void TInbreedingEstimator::_readData() {
	// read data
	using namespace coretools::instances;
	_likelihoods.doSaveAlleleFrequencies();
	if (parameters().parameterExists("trueAlleleFreq")) {
		_likelihoods.doSaveTrueAlleleFrequencies();
		logfile().list("Will save true allele frequencies for population likelihoods.");
	}
	_likelihoods.readData(parameters(), &logfile());
}

void TInbreedingEstimator::runEstimation(coretools::TParameters &Parameters) {
	// read file names
	std::string vcfFileName = Parameters.getParameter("vcf");
	vcfFileName             = coretools::str::extractBeforeLast(vcfFileName, ".vcf");
	auto prefix             = Parameters.getParameterWithDefault<std::string>("prefix", "");
	auto filename           = Parameters.getParameterWithDefault<std::string>("out", vcfFileName);

	// read data
	_readData();

	// build DAG
	auto dagBuilder = std::make_shared<stattools::TDAGBuilder>();
	TInbreedingEstimatorModel model(filename, dagBuilder, _likelihoods);
	dagBuilder->buildDAG();

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(prefix, dagBuilder);
}

} // end namespace PopulationTools