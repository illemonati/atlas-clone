//
// Created by caduffm on 3/8/22.
//

#include "TInbreedingEstimator.h"

#include <utility>

namespace PopulationTools {

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

TInbreedingEstimatorPrior::TInbreedingEstimatorPrior(stattools::TParameterTyped<TypeF, 1> *F,
                                                     stattools::TParameterTyped<TypeP, 1> *P,
                                                     stattools::TParameterTyped<TypeZ, 1> *Z,
                                                     const std::vector<double> &InitialEstimatesP)
    : _F(F), _p(P), _z(Z), _initialEstimatesP(InitialEstimatesP) {
	_priorParameters.push_back(_F);
	_priorParameters.push_back(_p);
	_priorParameters.push_back(_z);

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

	// F and z: one value
	_F->initializeStorageSingleElementBasedOnPrior();
	_z->initializeStorageSingleElementBasedOnPrior();

	// p: linear of length L
	_p->initializeStorageBasedOnPrior({_numLoci}, {lociNames});

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
                                                                                    double P) const {
	return _calculateLLSumOverIndividuals(Data, Locus, genometools::THardyWeinbergGenotypeProbabilities(P));
}

coretools::LogProbability TInbreedingEstimatorPrior::_calculateLLSumOverIndividuals(const Storage &Data, size_t Locus,
                                                                                    double P, double F) const {
	return _calculateLLSumOverIndividuals(Data, Locus,
	                                      genometools::THardyWeinbergWithInbreedingGenotypeProbabilities(P, F));
}

void TInbreedingEstimatorPrior::updateParams() {
	// get ptr to data
	auto data = *this->_storageBelow[0];

	_updateFAndZ(data);

	for (size_t l = 0; l < _numLoci; l++) { _updateP(data, l); }
}

void TInbreedingEstimatorPrior::_updateFAndZ(const Storage &Data) {
	using namespace coretools::instances;
	if (_F->isUpdated()) {
		if (_z->value()) { // we're in F-model
			if (_z->isUpdated() && randomGenerator().getRand() < _q_FModel_To_HWE) {
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

double TInbreedingEstimatorPrior::_calculateProbabilityOfProposingThisF(double F) const {
	// When switching from HWE to F-Model: propose values from truncated
	// exponential distribution with a certain lambdaF
	// -> in Hastings ratio, for the HWE-Model (without F), we need to know the
	// probability of picking exactly that F corresponds to log f(F) in Hastings
	// ratios
	return coretools::TExponentialDistr::logDensity(F, _lambdaNewF);
}

double TInbreedingEstimatorPrior::_getRandomNewF() const {
	using namespace coretools::instances;
	return randomGenerator().getExponentialRandomTruncated(_lambdaNewF, 0.0, 1.0);
}

void TInbreedingEstimatorPrior::_updateFToHWE(const Storage &Data) {
	// propose model switch: from F-Model to HWE-Model
	_F->set(0.0);
	_z->set(false);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_F = _calculateProbabilityOfProposingThisF(_F->oldValue());
	double logH    = log_f_F - _log_q_FModel_To_HWE - _F->getLogPriorDensityOld();

	// calculate LL Hastings-Ratio
	for (size_t l = 0; l < _numLoci; l++) {
		logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l)) // HWE model
		        - _calculateLLSumOverIndividuals(Data, l, _p->value(l),
		                                         _F->oldValue()); // F-Model
	}

	if (_F->acceptOrReject_DontCountRejection(logH)) {
		// accepted model switch -> we're now in null model
		// -> stop adding to posterior mean/var
		_F->doTrackMeanVar(false);
	} else {
		// rejected model switch -> we're still in F-model
		_z->set(true);
	}
}

void TInbreedingEstimatorPrior::_updateRegularF(const Storage &Data) {
	// update F without model switch: F -> F'
	if (_F->update()) {
		double logH = _F->getLogPriorRatio();
		for (size_t l = 0; l < _numLoci; l++) {
			logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->value()) -
			        _calculateLLSumOverIndividuals(Data, l, _p->value(l), _F->oldValue());
		}
		_F->acceptOrReject(logH);
	}
}

void TInbreedingEstimatorPrior::_updateHWEToF(const Storage &Data) {
	// propose model switch: from HWE-Model to F-Model
	_F->set(_getRandomNewF());
	_z->set(true);

	// calculate RJ-MCMC term of Hastings ratio
	double log_f_F =
	    _calculateProbabilityOfProposingThisF(_F->value()); // "compensating" the effect of no F in HWE model
	double logH = _log_q_FModel_To_HWE + _F->getLogPriorDensity() - log_f_F;

	// calculate LL Hastings-Ratio
	for (size_t l = 0; l < _numLoci; l++) {
		logH += _calculateLLSumOverIndividuals(Data, l, _p->value(l),
		                                       _F->value())              // F-Model
		        - _calculateLLSumOverIndividuals(Data, l, _p->value(l)); // HWE-Model
	}

	if (_F->acceptOrReject_DontCountRejection(logH)) {
		// accepted model switch -> we're now in F-model
		// -> start adding to posterior mean/var
		_F->doTrackMeanVar(true); // F is in non-zero model -> add to posterior mean/var
	} else {
		// rejected model switch -> we're still in null-model
		_z->set(false);
	}
}

void TInbreedingEstimatorPrior::_updateP(const Storage &Data, size_t Locus) {
	if (_p->updateSpecificIndex(Locus)) {
		double logH = _p->getLogPriorRatio();
		if (_z->value()) { // F-Model
			logH += _calculateLLSumOverIndividuals(Data, Locus, _p->value(Locus), _F->value()) -
			        _calculateLLSumOverIndividuals(Data, Locus, _p->oldValue(Locus), _F->value());
		} else { // HWE-Model
			logH += _calculateLLSumOverIndividuals(Data, Locus, _p->value(Locus)) -
			        _calculateLLSumOverIndividuals(Data, Locus, _p->oldValue(Locus));
		}

		_p->acceptOrReject(Locus, logH);
	}
}

void TInbreedingEstimatorPrior::_setInitialF() {
	// set initial F
	if (!_F->hasFixedInitialValue()) {
		if (_z->hasFixedInitialValue() && _z->value() == 0) {
			// user wants to start in zero-model
			_F->set(0.0);
		} else {
			_F->set(_getRandomNewF());
		}
	}
	// set initial z
	_F->value() == 0.0 ? _z->set(false) : _z->set(true);
}

void TInbreedingEstimatorPrior::_setInitialP() {
	using namespace coretools::instances;
	logfile().list("Initializing allele frequencies to values estimated from "
	               "sample genotype likelihoods");

	// now read values into _p
	if (!_p->hasFixedInitialValue()) {
		for (size_t l = 0; l < _numLoci; l++) { _p->set(l, std::max(_initialEstimatesP[l], TypeP::min().get())); }
	}
}

void TInbreedingEstimatorPrior::estimateInitialPriorParameters() {
	_setInitialF();
	_setInitialP();
}

double TInbreedingEstimatorPrior::getSumLogPriorDensity(const Storage &Data) const {
	double sum = 0.;
	for (size_t l = 0; l < _numLoci; l++) {
		if (_z->value()) {
			sum += _calculateLLSumOverIndividuals(Data, l, _p->value(l),
			                                      _F->value()); // F-Model
		} else {
			sum += _calculateLLSumOverIndividuals(Data, l, _p->value(l)); // HWE-Model
		}
	}
	return sum;
}

void TInbreedingEstimatorPrior::simulateUnderPrior() {
	throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": Not implemented. Use Atlas simulator instead.");
}

//------------------------------------------
// TInbreedingEstimatorModel
//------------------------------------------

TInbreedingEstimatorModel::TInbreedingEstimatorModel(
    const std::string &Filename, std::shared_ptr<stattools::TDAGBuilder> &DAGBuilder,
    const genometools::TPopulationLikelihoods<stattools::TValueFixed<TypeGTL>> &Likelihoods)
    : _F("F", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeF, 1>>(), {Filename}),
      _z("z", std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeZ, 1>>(), {Filename}),
      _log_alpha("log_alpha",
                 std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogAlpha, 1>>(),
                 {Filename}),
      _log_beta("log_beta",
                std::make_shared<stattools::prior::TUniformFixed<stattools::TParameterBase, TypeLogAlpha, 1>>(),
                {Filename}),
      _p("p",
         std::make_shared<stattools::prior::TBetaInferred<stattools::TParameterBase, TypeP, 1, TypeLogAlpha,
                                                          TypeLogBeta, true, true>>(&_log_alpha, &_log_beta),
         {Filename}),
      _observation("genotypeLikelihoods",
                   std::make_shared<TInbreedingEstimatorPrior>(&_F, &_p, &_z, Likelihoods.alleleFrequencies()), {}) {

	_p.getDefinition().setJumpSizeForAll(false);

	DAGBuilder->addToDAG(_F);
	DAGBuilder->addToDAG(_z);
	DAGBuilder->addToDAG(_log_alpha);
	DAGBuilder->addToDAG(_log_beta);
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