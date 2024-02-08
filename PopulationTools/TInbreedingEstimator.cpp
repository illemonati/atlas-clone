//
// Created by caduffm on 3/8/22.
//

#include "TInbreedingEstimator.h"

#include <array>
#include <cmath>
#include <stdexcept>

#include "coretools/Math/mathFunctions.h"
#include "stattools/DAG/TDAGBuilder.h"
#include "stattools/MCMC/TMCMC.h"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace PopulationTools {

using namespace coretools::instances;

//------------------------------------------
// TInbreedingEstimatorPrior
//------------------------------------------

TInbreedingEstimatorPrior::TInbreedingEstimatorPrior(TypeParamF *F, TypeParamP *P, TypeParamFModel *FModel,
                                                     TypeParamPModel *PModel,
                                                     const std::vector<double> &InitialEstimatesP)
    : _F(F), _p(P), _FModel(FModel), _pModel(PModel), _initialEstimatesP(InitialEstimatesP) {
	this->addPriorParameter({_F, _p, _FModel, _pModel});

	// read command line parameters
	_readCommandLineArguments();
}

std::string TInbreedingEstimatorPrior::name() const { return "inbreeding"; }

void TInbreedingEstimatorPrior::initialize() {
	assert(this->_storageBelow.size() == 1);
	auto data             = this->_storageBelow[0];
	const auto &lociNames = data->getDimensionName(0);

	// set dimensions
	_numLoci    = data->dimensions()[0];
	_numSamples = data->dimensions()[1];

	// F and FModel: one value
	_F->initStorage(this);
	_FModel->initStorage(this);

	// p and pModel: linear of length L
	_p->initStorage(this, {_numLoci}, {lociNames});
	_pModel->initStorage(this, {_numLoci}, {lociNames});
}

void TInbreedingEstimatorPrior::_readCommandLineArguments() {
	// Read number of threads
	_numThreads = coretools::getNumThreads();
	coretools::instances::logfile().list("Running with ", _numThreads, +" threads (argument 'numThreads')");
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

double TInbreedingEstimatorPrior::calculateLLRatio(TypeParamFModel *, size_t) const {
	// RJ-MCMC update of FModel
	if (_FModel->value() == 0) { // 1 -> 0
		return _calculateLLRatio_FToHWE();
	} else {
		return _calculateLLRatio_HWEToF();
	}
}

void TInbreedingEstimatorPrior::updateTempVals(TypeParamFModel *, size_t, bool) { /* empty: no temporary values */
}

double TInbreedingEstimatorPrior::_calculateLLRatio_FToHWE() const {
	// RJ-MCMC update: F-model (1) to HWE-model (0)
	const auto &data = *this->_storageBelow[0];

	double LLRatio = 0.0;
#pragma omp parallel for num_threads(_numThreads) reduction(+ : LLRatio)
	for (size_t l = 0; l < _numLoci; l++) {
		if (_pModel->value(l)) {                                                                // only if polyploid
			LLRatio += _calculateLLSumOverIndividuals(data, l, _p->value(l))                    // HWE model
			           - _calculateLLSumOverIndividuals(data, l, _p->value(l), _F->oldValue()); // F-Model
		}
	}
	return LLRatio;
}

double TInbreedingEstimatorPrior::_calculateLLRatio_HWEToF() const {
	// RJ-MCMC update: HWE-model (0) to F-model (1)
	const auto &data = *this->_storageBelow[0];

	double LLRatio = 0.0;
#pragma omp parallel for num_threads(_numThreads) reduction(+ : LLRatio)
	for (size_t l = 0; l < _numLoci; l++) {
		if (_pModel->value(l)) {                                                          // only if polyploid
			LLRatio += _calculateLLSumOverIndividuals(data, l, _p->value(l), _F->value()) // F-Model
			           - _calculateLLSumOverIndividuals(data, l, _p->value(l));           // HWE-Model
		}
	}

	return LLRatio;
}

double TInbreedingEstimatorPrior::calculateLLRatio(TypeParamF *, size_t) const {
	// regular Metropolis-Hastings update of F
	const auto &data = *this->_storageBelow[0];

	double LL = 0.0;
#pragma omp parallel for num_threads(_numThreads) reduction(+ : LL)
	for (size_t l = 0; l < _numLoci; l++) {
		if (_pModel->value(l)) { // only if polyploid
			LL += _calculateLLSumOverIndividuals(data, l, _p->value(l), _F->value()) -
			      _calculateLLSumOverIndividuals(data, l, _p->value(l), _F->oldValue());
		}
	}
	return LL;
}

void TInbreedingEstimatorPrior::updateTempVals(TypeParamF *, size_t, bool) { /* empty: no temporary values */
}

double TInbreedingEstimatorPrior::calculateLLRatio(TypeParamPModel *, size_t Locus) const {
	// RJ-MCMC update of pModel
	return _calculateLLRatio_UpdateP(Locus);
}

void TInbreedingEstimatorPrior::updateTempVals(TypeParamPModel *, size_t, bool) { /* empty: no temporary values */
}

double TInbreedingEstimatorPrior::_calculateLLRatio_UpdateP(size_t Locus) const {
	const auto &data = *this->_storageBelow[0];

	if (_FModel->value()) { // F-Model
		return _calculateLLSumOverIndividuals(data, Locus, _p->value(Locus), _F->value()) -
		       _calculateLLSumOverIndividuals(data, Locus, _p->oldValue(Locus), _F->value());
	} else { // HWE-Model
		return _calculateLLSumOverIndividuals(data, Locus, _p->value(Locus)) -
		       _calculateLLSumOverIndividuals(data, Locus, _p->oldValue(Locus));
	}
}

double TInbreedingEstimatorPrior::calculateLLRatio(TypeParamP *, size_t Locus) const {
	// prevent jumping to zero-model in here: always reject if p==0 or p==1 (invalid for Beta-distribution)
	if (_p->value(Locus) > 0.0 && _p->value(Locus) < 1.0) { return _calculateLLRatio_UpdateP(Locus); }
	return std::numeric_limits<double>::lowest();
}

void TInbreedingEstimatorPrior::updateTempVals(TypeParamP *, size_t, bool) {}

void TInbreedingEstimatorPrior::_setInitialF() {
	// set initial F
	if (!_F->hasFixedInitialValue()) {
		if (_FModel->hasFixedInitialValue() && _FModel->value() == 0) {
			// user wants to start in zero-model
			_F->set(0.0);
		} else {
			_F->set(_F->proposeNewValueRJMCMC());
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
					auto val = std::max(_initialEstimatesP[l], inbr::TypeP::min().get());
					val      = std::min(val, inbr::TypeP::max().get());
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

void TInbreedingEstimatorPrior::guessInitialValues() {
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

void TInbreedingEstimatorPrior::_simulateUnderPrior(Storage *) { DEVERROR("Use Atlas simulator instead."); }

//------------------------------------------
// TInbreedingEstimatorModel
//------------------------------------------

TInbreedingEstimatorModel::TInbreedingEstimatorModel(
    const std::string &Filename,
    const genometools::TPopulationLikelihoods<stattools::TValueFixed<inbr::TypeGTL>> &Likelihoods)
    : _boxOnFModel(), _FModel("withInbreeding", &_boxOnFModel, {Filename + "_F"}), _boxOnF(),
      _F("F", &_boxOnF, {Filename + "_F"}, &_FModel), _boxOnPi(), _pi("pi", &_boxOnPi, {Filename + "_p"}),
      _boxOnPModel(&_pi), _pModel("isPolymorph", &_boxOnPModel, {Filename + "_p"}), _boxOnLogGamma(),
      _logGamma("log_gamma", &_boxOnLogGamma, {Filename + "_p"}), _boxOnGamma(&_logGamma),
      _gamma("gamma", &_boxOnGamma, {Filename + "_p"}), _boxOnP(&_gamma, &_pModel),
      _p("p", &_boxOnP, {Filename + "_p"}, &_pModel, {1}),
      _boxOnObs(&_F, &_p, &_FModel, &_pModel, Likelihoods.alleleFrequencies()),
      _observation("genotypeLikelihoods", &_boxOnObs, Likelihoods.getStorage(), {}) {

	_p.getDefinition().setJumpSizeForAll(false);
	_p.getDefinition().setRJProposalDistrParams("1,100"); // parameters for Beta distribution
	_F.getDefinition().setRJProposalDistrParams("1,100"); // parameters for Beta distribution
}

//------------------------------------------
// TInbreedingEstimator
//------------------------------------------

void TInbreedingEstimator::_readData() {
	// read data
	_likelihoods.doSaveAlleleFrequencies();
	if (parameters().exists("trueAlleleFreq")) {
		_likelihoods.doSaveTrueAlleleFrequencies();
		logfile().list("Will save true allele frequencies for population likelihoods.");
	}
	_likelihoods.readData();
}

void TInbreedingEstimator::run() {
	// read file names
	std::string vcfFileName = parameters().get("vcf");
	vcfFileName             = coretools::str::extractBeforeLast(vcfFileName, ".vcf");
	auto prefix             = parameters().get<std::string>("prefix", "");
	auto filename           = parameters().get<std::string>("out", vcfFileName);

	// read data
	_readData();

	// create model
	TInbreedingEstimatorModel model(filename, _likelihoods);

	// run MCMC
	stattools::TMCMC mcmc;
	mcmc.runMCMC(prefix);
}

} // end namespace PopulationTools
