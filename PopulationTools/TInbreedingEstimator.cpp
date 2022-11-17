//
// Created by caduffm on 3/8/22.
//

#include "TInbreedingEstimator.h"

#include <array>
#include <cmath>
#include <exception>
#include <stdexcept>

#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/MCMC/TMCMC.h"
#include "stattools/ParametersObservations/TParameterObservationTypedBase.h"
#include "stattools/Priors/TPriorBernouilli.h"
#include "stattools/Priors/TPriorBeta.h"
#include "stattools/Priors/TPriorUniform.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Strings/stringFunctions.h"

namespace stattools {
class TParameterBase;
}
namespace stattools {
class TParameterObservationBase;
}

namespace PopulationTools {

using namespace coretools::instances;

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

TInbreedingEstimatorPrior::TInbreedingEstimatorPrior(stattools::TParameterTyped<TypeF, 1> *F,
                                                     stattools::TParameterTyped<TypeP, 1> *P,
                                                     stattools::TParameterTyped<TypeFModel, 1> *FModel,
                                                     stattools::TParameterTyped<TypePModel, 1> *PModel,
                                                     const std::vector<double> &InitialEstimatesP)
    : _F(F), _p(P), _FModel(FModel), _pModel(PModel), _initialEstimatesP(InitialEstimatesP) {
	this->addPriorParameter({_F, _p, _FModel, _pModel});

	// read command line parameters
	_readCommandLineArguments();
}

std::string TInbreedingEstimatorPrior::name() const { return "inbreeding"; }

void TInbreedingEstimatorPrior::initializeInferred() {
	assert(this->_storageBelow.size() == 1);
	auto data             = this->_storageBelow[0];
	const auto &lociNames = data->getDimensionName(0);

	// set dimensions
	_numLoci    = data->dimensions()[0];
	_numSamples = data->dimensions()[1];

	// F and FModel: one value
	_F->initStorage();
	_FModel->initStorage();

	// p and pModel: linear of length L
	_p->initStorage({_numLoci}, {lociNames});
	_pModel->initStorage({_numLoci}, {lociNames});
}

void TInbreedingEstimatorPrior::_readCommandLineArguments() {
	// parameters for model switch of F
	_q_FModel_To_HWE     = parameters().getParameterWithDefault("probMovingToModelNoF", 0.1);
	_log_q_FModel_To_HWE = log(_q_FModel_To_HWE);
	logfile().list("Will propose move to model without F with probability ", _q_FModel_To_HWE,
	               ". (use 'probMovingToModelNoF' to change)");

	// Read lambda for proposing new F
	_lambdaNewF = parameters().getParameterWithDefault("lambdaF", 100.0);
	logfile().list("Setting lambda of exponential distribution used for the proposal of "
	               "new F after move to model with F to ",
	               _lambdaNewF, ". (use 'lambdaF' to change)");

	// parameters for model switch of p
	_q_PModel_To_NullModel     = parameters().getParameterWithDefault("probMovingToModelP0", 0.1);
	_log_q_PModel_To_NullModel = log(_q_PModel_To_NullModel);
	logfile().list("Will propose move to monomorphic model with probability ", _q_PModel_To_NullModel,
	               ". (use 'probMovingToModelP0' to change)");

	// Read lambda for proposing new p
	_lambdaNewP = parameters().getParameterWithDefault("lambdaP", 100.0);
	logfile().list("Lambda of exponential distribution used for the proposal of "
	               "new p after move to polymorphic model is set to ",
	               _lambdaNewP, ". (use 'lambdaP' to change)");
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
	// get reference to data
	const auto &data = *this->_storageBelow[0];

	_updateF(data);
	for (size_t l = 0; l < _numLoci; l++) { _updateP(data, l); }
}

void TInbreedingEstimatorPrior::_updateF(const Storage &Data) {
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

double calcLogProposalProb(double Value, coretools::StrictlyPositive<double> Lambda) {
	// When switching from One to Null-Model: propose values from truncated
	// exponential distribution with a certain lambda
	// -> in Hastings ratio, for the Null-Model, we need to know the
	// probability of picking exactly that Value corresponds to log f(Value) in Hastings
	// ratios
	// applies to model switches of both F and p
	return coretools::probdist::TExponentialDistr::logDensity(Value, Lambda);
}

double getRandomExpTrunc(double Lambda) {
	// When switching from null to One-Model: propose new value from truncated exponential distribution
	// used for both F and p
	auto val = coretools::instances::randomGenerator().getExponentialRandomTruncated(Lambda, 0.0, 1.0);
	while (val == 0.0 || val == 1.0) { // make sure not to propose exactly 0 or 1
		val = coretools::instances::randomGenerator().getExponentialRandomTruncated(Lambda, 0.0, 1.0);
	}
	return val;
}

void TInbreedingEstimatorPrior::_updateFToHWE(const Storage &Data) {
	// propose model switch: from F-Model to HWE-Model
	_F->set(0.0);
	_FModel->set(false);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_F = calcLogProposalProb(_F->oldValue(), _lambdaNewF);
	double logH    = log_f_F - _log_q_FModel_To_HWE - _F->getLogPriorDensityOld();

	// calculate LL Hastings-Ratio
	for (size_t l = 0; l < _numLoci; l++) {
		if (_pModel->value(l)) {                                                             // only if polyploid
			logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l))                    // HWE model
			        - _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->oldValue()); // F-Model
		}
	}

	if (!_F->acceptOrReject_DontCountRejection(logH)) {
		// rejected model switch -> we're still in F-model
		_FModel->set(true);
	}
	_FModel->addToMeanVar();
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
	_F->set(getRandomExpTrunc(_lambdaNewF));
	_FModel->set(true);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_F = calcLogProposalProb(_F->value(), _lambdaNewF); // "compensating" the effect of no F in HWE model
	double logH    = _log_q_FModel_To_HWE + _F->getLogPriorDensity() - log_f_F;

	// calculate LL Hastings-Ratio
	for (size_t l = 0; l < _numLoci; l++) {
		if (_pModel->value(l)) {                                                       // only if polyploid
			logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->value()) // F-Model
			        - _calculateLLSumOverIndividuals(Data, l, _p->value(l));           // HWE-Model
		}
	}

	if (!_F->acceptOrReject_DontCountRejection(logH)) {
		// rejected model switch -> we're still in null-model
		_FModel->set(false);
	}
	_FModel->addToMeanVar();
}

void TInbreedingEstimatorPrior::_updateP(const Storage &Data, size_t Locus) {
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
	double log_f_p = calcLogProposalProb(_p->oldValue(Locus), _lambdaNewP);
	double logH    = log_f_p - _log_q_PModel_To_NullModel - _p->getLogPriorDensityOld(Locus);
	logH += _pModel->getLogPriorRatio(Locus);
	logH += _calculateLLRatio_UpdateP(Data, Locus);

	if (!_p->acceptOrReject_DontCountRejection(Locus, logH)) {
		// rejected model switch -> we're still in P-model
		_pModel->set(Locus, true);
	}
	_pModel->addToMeanVar();
}

void TInbreedingEstimatorPrior::_updateRegularP(const Storage &Data, size_t Locus) {
	if (_p->updateSpecificIndex(Locus)) {
		double logH = std::numeric_limits<double>::lowest();
		if (_p->value(Locus) > 0.0 && _p->value(Locus) < 1.0) {
			// prevent jumping to zero-model in here: always reject if p==0 or p==1 (invalid for Beta-distribution)
			logH = _calculateLLRatio_UpdateP(Data, Locus) + _p->getLogPriorRatio(Locus);
		}
		_p->acceptOrReject(Locus, logH);
	}
}

void TInbreedingEstimatorPrior::_updateNullToP(const Storage &Data, size_t Locus) {
	// propose model switch: from null-Model to p-Model
	_p->set(Locus, getRandomExpTrunc(_lambdaNewP));
	_pModel->set(Locus, true);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_p = calcLogProposalProb(_p->value(Locus), _lambdaNewP);
	double logH    = _log_q_PModel_To_NullModel + _p->getLogPriorDensity(Locus) - log_f_p;
	logH += _pModel->getLogPriorRatio(Locus);
	logH += _calculateLLRatio_UpdateP(Data, Locus);

	if (!_p->acceptOrReject_DontCountRejection(Locus, logH)) {
		// rejected model switch -> we're still in null-model
		_pModel->set(Locus, false);
	}
	_pModel->addToMeanVar();
}

void TInbreedingEstimatorPrior::_setInitialF() {
	// set initial F
	if (!_F->hasFixedInitialValue()) {
		if (_FModel->hasFixedInitialValue() && _FModel->value() == 0) {
			// user wants to start in zero-model
			_F->set(0.0);
		} else {
			_F->set(getRandomExpTrunc(_lambdaNewF));
		}
	}
	// set initial z
	_F->value() == 0.0 ? _FModel->set(false) : _FModel->set(true);
}

void TInbreedingEstimatorPrior::_setInitialP() {
	logfile().list("Initializing allele frequencies to values estimated from "
	               "sample genotype likelihoods.");

	// now read values into _p and _pModel
	if (!_p->hasFixedInitialValue()) {
		for (size_t l = 0; l < _numLoci; l++) {
			if (_pModel->hasFixedInitialValue()) {
				if (_pModel->value(l) == 0) { // user wants to start in 0-Model
					_p->set(l, 0.0);
				} else { // user wants to start in 1-model: prevent p = 0 and p = 1
					auto val = std::max(_initialEstimatesP[l], TypeP::min().get());
					val      = std::min(val, TypeP::max().get());
					_p->set(l, val);
				}
			} else {
				_p->set(l, _initialEstimatesP[l]);
				if (_initialEstimatesP[l] == 0.0 || _initialEstimatesP[l] == 1.0) {
					_pModel->set(l, false);
				} else {
					_pModel->set(l, true);
				}
			}
		}
	} else { // _p has fixed initial value -> set pModel
		if (!_pModel->hasFixedInitialValue()) {
			for (size_t l = 0; l < _numLoci; l++) {
				(_p->value(l) == 0.0 || _p->value(l) == 1.0) ? _pModel->set(l, false) : _pModel->set(l, true);
			}
		}
	}
}

void TInbreedingEstimatorPrior::estimateInitialPriorParameters() {
	_setInitialF();
	_setInitialP();
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
    const std::string &Filename, stattools::TDAGBuilder &DAGBuilder,
    const genometools::TPopulationLikelihoods<stattools::TValueFixed<TypeGTL>> &Likelihoods)
    : _F("F", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeF, 1>>(),
         {Filename + "_F"}),
      _FModel("withInbreeding",
              std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeFModel, 1>>(),
              {Filename + "_F"}),
      _pi("pi", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypePi, 1>>(),
          {Filename + "_p"}),
      _pModel("isPolymorph",
              std::make_shared<stattools::prior::TBernouilliInferred<stattools::TParameterBase, TypePModel, 1, TypePi>>(
                  &_pi),
              {Filename + "_p"}),
      _log_gamma("log_gamma",
                 std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogGamma, 1>>(),
                 {Filename + "_p"}),
      _p("p",
         std::make_shared<stattools::prior::TBetaSymmetricZeroMixtureInferred<stattools::TParameterBase, TypeP, 1,
                                                                              TypeLogGamma, TypePModel, true>>(
             &_log_gamma, &_pModel),
         {Filename + "_p"}),
      _observation(
          "genotypeLikelihoods",
          std::make_shared<TInbreedingEstimatorPrior>(&_F, &_p, &_FModel, &_pModel, Likelihoods.alleleFrequencies()),
          Likelihoods.getStorage(), {}) {

	_p.getDefinition().setJumpSizeForAll(false);
	_p.getDefinition().setEqualNumberOfUpdates(false);

	DAGBuilder.addToDAG({&_F, &_FModel, &_pi, &_pModel, &_log_gamma, &_p});
	DAGBuilder.addToDAG(&_observation);
}

//------------------------------------------
// TInbreedingEstimator
//------------------------------------------

void TInbreedingEstimator::_readData() {
	// read data
	_likelihoods.doSaveAlleleFrequencies();
	if (parameters().parameterExists("trueAlleleFreq")) {
		_likelihoods.doSaveTrueAlleleFrequencies();
		logfile().list("Will save true allele frequencies for population likelihoods.");
	}
	_likelihoods.readData();
}

void TInbreedingEstimator::run() {
	// read file names
	std::string vcfFileName = parameters().getParameter("vcf");
	vcfFileName             = coretools::str::extractBeforeLast(vcfFileName, ".vcf");
	auto prefix             = parameters().getParameterWithDefault<std::string>("prefix", "");
	auto filename           = parameters().getParameterWithDefault<std::string>("out", vcfFileName);

	// read data
	_readData();

	// build DAG
	stattools::TDAGBuilder dagBuilder;
	TInbreedingEstimatorModel model(filename, dagBuilder, _likelihoods);
	dagBuilder.buildDAG();

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(prefix, &dagBuilder);
}

} // end namespace PopulationTools
