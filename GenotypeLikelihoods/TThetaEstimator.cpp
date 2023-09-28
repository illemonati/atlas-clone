/*
 * TThetaEstimator.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#include "TThetaEstimator.h"

#include <armadillo>
#include <exception>
#include <memory>
#include <stdlib.h>

#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSite.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/weakTypes.h"

namespace GenotypeLikelihoods {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using coretools::str::toString;

//---------------------------------------------------------------
// TThetaEstimator_base
//---------------------------------------------------------------

TGenotypeProbabilities getPGenotype(double expTheta, const TBaseProbabilities &baseFrequencies) {
	using namespace genometools;
	// assumes that base frequencies are set!
	TGenotypeLikelihoods lGeno;
	for (Base b = Base::min; b < Base::max; ++b) {
		// homozygous genotypes
		Genotype hom = genotype(b, b);
		lGeno[hom]   = baseFrequencies[b] * (expTheta + baseFrequencies[b].get() * (1.0 - expTheta));
		// heterozygous genotypes
		for (Base c = coretools::next(b); c < Base::max; ++c) {
			Genotype het = genotype(b, c);
			lGeno[het]   = 2.0 * baseFrequencies[b].get() * baseFrequencies[c].get() * (1.0 - expTheta);
		}
	}
	return TGenotypeProbabilities::normalize(lGeno);
};

GenotypeLikelihoods::TGenotypeProbabilities getPGenotype(const Theta &thisTheta) {
	return getPGenotype(thisTheta.expMinusTheta, thisTheta.baseFreq);
};

TThetaEstimator_base::TThetaEstimator_base()
	: _useTmpFile(coretools::instances::parameters().parameterExists("useTmpFile")),
	  _minSitesWithData(coretools::instances::parameters().getParameterWithDefault<int>("minSitesWithData", 1000)),
	  _extraVerbose(coretools::instances::parameters().parameterExists("extraVerbose")) {
	_initDataStorage();
};

TThetaEstimator_base::TThetaEstimator_base(const TThetaEstimator_base &other) {
	_data            = nullptr;

	// data
	_useTmpFile = other._useTmpFile;
	_initDataStorage();

	// minimum window size
	_minSitesWithData = other._minSitesWithData;

	// extra verbosity
	_extraVerbose   = other._extraVerbose;
	_equalBaseFreqs = other._equalBaseFreqs;
};

void TThetaEstimator_base::_initDataStorage() {
	_data.reset();

	// file or vector?
	if (_useTmpFile) {
		_tmpFileName =
			"temporaryDataForThetaEstimation_" + randomGenerator().getRandomAlphaNumericString(10) + ".tmp.gz";
		logfile().list("Will write temporary data to file '" + _tmpFileName + "'.");
		_data = std::make_unique<TThetaEstimatorDataFile>(_tmpFileName);
	} else {
		_data = std::make_unique<TThetaEstimatorDataVector>();
	}
};

void TThetaEstimator_base::_readParametersRegardingInitialSearch() {
	using coretools::str::toString;
	logfile().startIndent("Parameters of the initial theta search:");
	_initialTheta = parameters().getParameterWithDefault("initTheta", 0.01);
	logfile().list("Will start with an initial theta of " + toString(_initialTheta) + ".");
	_initThetaNumSearchIterations = parameters().getParameterWithDefault("initThetaNumSearchIterations", 10);
	if (_initThetaNumSearchIterations > 0) {
		logfile().list("Will run " + toString(_initThetaNumSearchIterations) +
					   " iterations of a crude search for an initial theta.");
		_initThetaSearchFactor = parameters().getParameterWithDefault("initThetaSearchFactor", 100);
		logfile().list("The initial search factor will be " + toString(_initThetaSearchFactor) + ".");
	} else {
		_initThetaSearchFactor = 0;
	}
	logfile().endIndent();
};

void TThetaEstimator_base::_findGoodStartingTheta(TThetaEstimatorData *thisData, Theta &thisTheta,
												 const std::string &tag) {
	logfile().listFlush("Estimating initial parameters" + tag + " ...");

	// set base frequencies to initial base frequencies
	thisTheta.baseFreq = _equalBaseFreqs ? TBaseProbabilities{} : thisData->baseFrequencies();

	// variables
	double initTheta = _initialTheta;
	double oldTheta  = initTheta;
	double expTheta  = exp(-initTheta);

	// calc initial LL
	const auto pGenotype   = getPGenotype(expTheta, thisTheta.baseFreq);
	thisTheta.LL = thisData->calcLogLikelihood(pGenotype);

	// run iterations
	double oldLL  = thisTheta.LL;
	double factor = _initThetaSearchFactor;
	int numUpdates;
	for (int i = 0; i < _initThetaNumSearchIterations; ++i) {
		// first test increase in theta
		numUpdates = -1;
		do {
			++numUpdates;
			oldLL    = thisTheta.LL;
			oldTheta = initTheta;
			initTheta *= factor;
			expTheta     = exp(-initTheta);
			thisTheta.LL = thisData->calcLogLikelihood(getPGenotype(expTheta, thisTheta.baseFreq));
		} while (oldLL < thisTheta.LL);
		if (numUpdates == 0) {
			// then test decrease in theta
			initTheta    = oldTheta;
			thisTheta.LL = oldLL;
			// maybe smaller?
			do {
				oldLL    = thisTheta.LL;
				oldTheta = initTheta;
				initTheta /= factor;
				expTheta     = exp(-initTheta);
				thisTheta.LL = thisData->calcLogLikelihood(getPGenotype(expTheta, thisTheta.baseFreq));
			} while (oldLL < thisTheta.LL);
		}
		factor       = sqrt(factor);
		initTheta    = oldTheta;
		thisTheta.LL = oldLL;
	}
	// return previous
	thisTheta.set(oldTheta);
	thisTheta.LL = oldLL;

	// check if values make sense. If theta < 1/(10*windowsize), set it to 1/(10*windowsize)
	if (thisTheta.theta < 0.1 / (double)thisData->size()) {
		thisTheta.set(0.1 / (double)thisData->size());
	} else if (thisTheta.theta > 1.0) {
		thisTheta.set(1.0);
	}

	// conclude
	logfile().done();
	logfile().conclude("Initial base frequencies: " + thisTheta.string());
	logfile().conclude("Initial theta = ", thisTheta.theta);
}

//-------------------------------------------------------
// TThetaEstimator
//-------------------------------------------------------
TThetaEstimator::TThetaEstimator() : TThetaEstimator_base() {
	_initialTheta                 = 0.0;
	_initThetaSearchFactor        = -1;
	_initThetaNumSearchIterations = -1;

	// EM
	_equalBaseFreqs = parameters().parameterExists("equalBaseFreqs");
	if (_equalBaseFreqs) logfile().list("Will assume equal base frequencies. (parameter 'equalBaseFreqs')");
	else logfile().list("Will estimate base frequencies. (use 'equalBaseFreqs' to assume equal base frequencies)");

	logfile().startIndent("Parameters of EM algorithm to infer theta:");
	_numIterations = parameters().getParameterWithDefault<int>("iterations", 100);
	logfile().list("Will run up to " + toString(_numIterations) + " iterations.");
	if (_equalBaseFreqs) {
		_numThetaOnlyUpdates = _numIterations;
	} else {
		_numThetaOnlyUpdates = parameters().getParameterWithDefault<int>("iterationsThetaOnly", 10);
		logfile().list("In each iteration, theta will be updated " + toString(_numThetaOnlyUpdates) + " times.");
	}
	_maxEpsilon = parameters().getParameterWithDefault("maxEps", 0.000001);
	logfile().list("Will run EM until deltaLL < " + toString(_maxEpsilon) + ".");

	//NR
	_NewtonRaphsonNumIterations = parameters().getParameterWithDefault<int>("NRiterations", 10);
	logfile().list("Will run Newton-Raphson algorithm up to " + toString(_NewtonRaphsonNumIterations) + " times.");

	_NewtonRaphsonMaxF = parameters().getParameterWithDefault("maxF", 0.00001);
	logfile().list("Will run Newton-Raphson algorithm until max(F) < " + toString(_NewtonRaphsonMaxF) + ".");
	logfile().endIndent();

	// params regarding initial search
	_readParametersRegardingInitialSearch();
}

void TThetaEstimator::clear() {
	_data->clear();
	_estimationSuccessful = false;
};

void TThetaEstimator::add(const GenotypeLikelihoods::TSite &site,
						  const GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods) {
	_data->add(site, genotypeLikelihoods);
};

void TThetaEstimator::add(const TWindow &window, const TGenotypeLikelihoodCalculator &glCalculator) {
	for (std::vector<TSite>::const_iterator it = window.cbegin(); it != window.cend(); ++it) {
		add(*it, glCalculator.calculateGenotypeLikelihoods(*it));
	}
};

bool TThetaEstimator::_NRAllParams(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) {
	using namespace genometools;
	using coretools::Probability;
	using coretools::index;
	// calculate substitution probabilities

	// Calculate all genotype probabilities for all sites
	const GenotypeLikelihoods::TGenotypeData P_G = _data->P_G(pGenotype); // see paper

	// prepare storage
	arma::mat::fixed<6,6> Jacobian;
	arma::vec::fixed<6> F;

	double rho = _theta.expMinusTheta / (1.0 - _theta.expMinusTheta);
	double mu  = _data->sizeWithData();

	TBaseProbabilities &baseFreq = _theta.baseFreq; // for cleaner code
	for (int n = 0; n < _NewtonRaphsonNumIterations; ++n) {
		// i) calculate F (Note: index is zero based!)
		F(4) = _data->sizeWithData();
		F(5) = 0.0;
		for (Base k = Base::min; k < Base::max; ++k) {
			const auto hom = genotype(k, k);
			double tmpSum  = 0.0;
			for (Base l = Base::min; l < Base::max; ++l) {
				if (l != k) {
					const auto het = genotype(k, l);
					tmpSum += P_G[het];
				}
			}
			F(coretools::index(k)) =
				P_G[hom] * (1.0 + baseFreq[k].get() / (rho + baseFreq[k].get())) + tmpSum - mu * baseFreq[k].get();
			F(4) -= P_G[hom] * (rho + 1.0) / (rho + baseFreq[k].get());
			F(5) += baseFreq[k].get();
		}
		F(5) = F(5) - 1.0;

		// ii) fill Jacobian (Note: index is zero based!)
		Jacobian.zeros();
		TBaseData tmp;
		for (Base k = Base::min; k < Base::max; ++k) {
			const auto hom = genotype(Base(k), Base(k));
			tmp[k]         = P_G[hom] / ((baseFreq[k].get() + rho) * (baseFreq[k].get() + rho));
		}

		for (Base k = Base::min; k < Base::max; ++k) {
			const auto i   = coretools::index(k);
			Jacobian(i, i) = tmp[k] * rho - mu;
			Jacobian(i, 4) = -tmp[k];
			Jacobian(5, i) = 1.0;
			Jacobian(4, i) = tmp[k] * (rho + 1.0);
			Jacobian(i, 5) = -baseFreq[k].get();
			Jacobian(4, 4) += tmp[k] * (1.0 - baseFreq[k].get());
		}

		// iii) now estimate new parameters
		mu = _data->sizeWithData();

		arma::vec::fixed<6> JxF;
		if (solve(JxF, Jacobian, F)) {
			// baseFreq = TBaseProbabilities(baseFreq, JxF, std::minus<>());
			/*for(Base k = Base::min; k < Base::max; ++k){
				baseFreq[k] = baseFreq[k].get() - JxF(index(k));
				}*/
			baseFreq.for_each_index(
				[&JxF](Probability p, Base i) { return static_cast<Probability>(p.get() - JxF(coretools::index(i))); });
			rho -= JxF(4);
			mu -= JxF(5);

			// check if we break
			double maxF = 0.0;
			for (int i = 0; i < 6; ++i) {
				if (F(i) > maxF) maxF = F(i);
			}

			if (maxF < _NewtonRaphsonMaxF) {
				_theta.set(-log(rho / (1.0 + rho)));
				return true;
			}
		} else {
			return false;
		}
	}
	return true;
};

void TThetaEstimator::_NROnlyTheta(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) {
	using namespace genometools;

	// Calculate all genotype probabilities for all sites
	const GenotypeLikelihoods::TGenotypeData P_G = _data->P_G(pGenotype); // see paper

	double rho = _theta.expMinusTheta / (1.0 - _theta.expMinusTheta);

	for (int n = 0; n < _NewtonRaphsonNumIterations; ++n) {

		// i) calculate F() (Note: index is zero based!)
		// ii) fill Jacobian (Note: index is zero based!)
		double F        = _data->sizeWithData();
		double Jacobian = 0.0;
		for (Base k = Base::min; k < Base::max; ++k) {
			const auto hom = genotype(k, k);

			const double tmp = (_theta.baseFreq[k].get() + rho);
			F -= P_G[hom] * (rho + 1.0) / tmp;

			const double tmpSum = P_G[hom] / (tmp * tmp);
			Jacobian += tmpSum * (1.0 - _theta.baseFreq[k].get());
		}

		// iii) now estimate new parameters
		rho = rho - F / Jacobian;

		// check if we break
		if (F < _NewtonRaphsonMaxF) {
			_theta.set(-log(rho / (1.0 + rho)));
			break;
		}
	}
};

void TThetaEstimator::_runEMForTheta() {
	// increase initialTheta if we fail to calculate inverse of Jacobian.
	//  this may be the case if initialTheta is smaller than true theta and likelihood is very flat
	double startingTheta = _initialTheta;
	_theta.LL             = -9e100;
	while (!_equalBaseFreqs && !_NRAllParams(getPGenotype(_theta))) {
		// solve did not work -> start with higher theta!
		startingTheta *= 2.0;
		_theta.set(startingTheta);
		if (startingTheta > 1.0) UERROR("Failed to estimate Theta, issues calculating inverse of Jacobian!");
	}

	// start EM loop
	for (int iter = 0; iter < _numIterations; ++iter) {
		// update only theta: most difficult parameter and it is much faster to update only this one alone.
		int i           = 0;
		double oldTheta = 0.0;
		double oldLL    = _theta.LL;
		do {
			oldTheta = _theta.theta;
			_NROnlyTheta(getPGenotype(_theta));
			++i;
		} while (i < _numThetaOnlyUpdates && _theta.theta != oldTheta);

		// update all params
		const auto pGenotype = getPGenotype(_theta);
		if (!_equalBaseFreqs) _NRAllParams(pGenotype);

		// e) do we break EM? Check LL
		_theta.LL = _data->calcLogLikelihood(pGenotype);
		if ((_theta.LL - oldLL) < _maxEpsilon) break;

		// maybe theta = 0?
		if (_theta.theta < 0.1 / (double)_data->size()) {
			//(theta is somewhere between 1/numLoci and 0,
			oldTheta = _theta.theta;

			// test with theta = 0.0
			_theta.set(0.0);
			_theta.LL   = _data->calcLogLikelihood(getPGenotype(_theta));

			// theta is between zero and a very small number -> don't care about exact value
			if (_theta.LL < oldLL) {
				_theta.set(oldTheta);
				_theta.LL = oldLL;
			}
			break;
		}

		if (_extraVerbose) logfile().list(toString(iter) + ") current theta = " + toString(_theta.theta));
	}
	if (_extraVerbose) logfile().list("EM converged, current theta = " + toString(_theta.theta));
}

void TThetaEstimator::_estimateConfidenceInterval() {
	using namespace genometools;
	// we estimate an approximate confidence interval for theta using the Fisher information
	// This function assumes that EM has already been run!

	// calculate P(g|theta, pi)
	const auto pGenotype = getPGenotype(_theta);

	// calclate d/dtheta P(g|theta, pi)
	TGenotypeData deriv_pGenotype;

	for (Base k = Base::min; k < Base::max; ++k) {
		// homozygous genotype
		const auto hom = genotype(k, k);
		deriv_pGenotype[hom] =
			(_theta.baseFreq[k].get() * _theta.baseFreq[k].get() - _theta.baseFreq[k].get()) * _theta.expMinusTheta;
		// heterozygous genotypes
		for (Base l = coretools::next(k); l < Base::max; ++l) {
			const auto het       = genotype(k, l);
			deriv_pGenotype[het] = 2.0 * _theta.baseFreq[k].get() * _theta.baseFreq[l].get() * _theta.expMinusTheta;
		}
	}

	double FisherInfo = _data->fisherInfo(pGenotype, deriv_pGenotype);

	// estimate confidence interval
	// TODO: Fisher Info can be negative -> SQRT will be nan!
	_theta.thetaConfidence = 1.96 / sqrt(FisherInfo);
}

void TThetaEstimator::_calcExpectedHet(){
	using namespace genometools;
	//calculating epxected heterozygosity under the Felsenstein model
	double hom = 0.0;
	for (Base k = Base::min; k < Base::max; ++k) {
		hom += _theta.baseFreq[k] * (_theta.expMinusTheta + _theta.baseFreq[k] * (1.0 - _theta.expMinusTheta));
	}
	_expectedHet = 1.0 - hom;
}

//------------------------------------------------------------
// Functions to run estimation-
//------------------------------------------------------------
bool TThetaEstimator::estimateTheta() {
	_estimationSuccessful = false;
	if (_data->sizeWithData() < _minSitesWithData) {
		logfile().warning("Can not estimate theta: only " + toString(_data->sizeWithData()) +
						  " sites with data in this region (minSitesWithData = " + toString(_minSitesWithData) + ")!");
		return false;
	}

	// estimate starting parameters
	_findGoodStartingTheta(_data.get(), _theta, "");

	// Run EM
	if (_extraVerbose) {
		logfile().startIndent("Running EM to find ML estimate:");
		_runEMForTheta();
		logfile().endIndent();
	} else {
		logfile().listFlush("Running EM to find ML estimate ...");
		_runEMForTheta();
		logfile().done();
	}
	logfile().conclude("theta was estimated at ", _theta.theta);
	_calcExpectedHet();
	logfile().conclude("implies expected heterozygosity of ", _expectedHet);

	// confidence intervals
	logfile().listFlush("Estimating approximate confidence intervals from Fisher-Information ...");
	_estimateConfidenceInterval();
	logfile().done();
	logfile().conclude("95% confidence intervals are theta +- " + toString(_theta.thetaConfidence));
	_estimationSuccessful = true;
	return true;
}

void TThetaEstimator::setTheta(const double Theta) { _theta.set(Theta); };

void TThetaEstimator::setBaseFreq(const GenotypeLikelihoods::TBaseProbabilities &BaseFreq) {
	_theta.baseFreq = BaseFreq;
};

void TThetaEstimator::addToHeader(std::vector<std::string> &header, const std::string &prefix) {
	_data->addToHeader(header, prefix);
	header.push_back(prefix + "pi(A)");
	header.push_back(prefix + "pi(C)");
	header.push_back(prefix + "pi(G)");
	header.push_back(prefix + "pi(T)");
	header.push_back(prefix + "theta_MLE");
	header.push_back(prefix + "theta_C95_l");
	header.push_back(prefix + "theta_C95_u");
	header.push_back(prefix + "LL");
	header.push_back(prefix + "expHet_MLE");
}

void TThetaEstimator::writeEstimateFrequenciesAndTheta(coretools::TOutputFile &out) {
	if (_estimationSuccessful) {
		// base frequencies
		for (genometools::Base k = genometools::Base::min; k < genometools::Base::max; ++k) {
			out << _theta.baseFreq[k];
		}

		// theta estimates
		out << _theta.theta;
		out << _theta.theta - _theta.thetaConfidence;
		out << _theta.theta + _theta.thetaConfidence;
		out << _theta.LL;
		out << _expectedHet;
	} else {
		// frequencies
		out << "-"
			<< "-"
			<< "-"
			<< "-";

		// theta estimates
		out << "-"
			<< "-"
			<< "-"
			<< "-"
			<< "-";
	}
};

void TThetaEstimator::writeResultsToFile(coretools::TOutputFile &out) {
	// number of sites
	_data->writeSite(out);

	// estimated params
	writeEstimateFrequenciesAndTheta(out);
};

void TThetaEstimator::calcLikelihoodSurface(coretools::TOutputFile &out, size_t steps) {
	// write header
	out.writeHeader({"log10(theta)", "theta", "LL"});
	out.precision(12);

	// calculate likelihood surface
	double minLogTheta = -5.0;
	double maxLogTheta = -1.0;
	double stepSize    = (maxLogTheta - minLogTheta) / ((double)steps - 1.0);

	for (size_t i = 0; i < steps; ++i) {
		// calc theta and expMinusTheta
		_theta.setLog(minLogTheta + stepSize * i);

		// calculate	substitution probabilities and Likelihood
		_theta.LL   = _data->calcLogLikelihood(getPGenotype(_theta));

		// write results
		out.writeln(_theta.logTheta, _theta.theta, _theta.LL);
	}
};

void TThetaEstimator::bootstrapTheta() {
	logfile().listFlush("Bootstrapping sites ...");

	_data->bootstrap();
	logfile().done();

	// estimate theta
	estimateTheta();

	// clean up
	_data->clearBootstrap();
};

//---------------------------------------------------------------
// TThetaEstimatorRatio
//---------------------------------------------------------------
TThetaEstimatorRatio::TThetaEstimatorRatio() : TThetaEstimator_base() {
	_clearCounters();

	// data2
	if (_useTmpFile) {
		_data2 = std::make_unique<TThetaEstimatorDataFile>(_tmpFileName + "2.tmp.gz");
	} else
		_data2 = std::make_unique<TThetaEstimatorDataVector>();

	// MCMC params
	logfile().startIndent("Parameters of MCMC algorithm:");
	_burnin = parameters().getParameterWithDefault<int>("burnin", 1000);
	logfile().list("Will run a burnin of " + toString(_burnin) + " iterations.");
	_numIterations = parameters().getParameterWithDefault<int>("iterations", 10000);
	logfile().list("Will run MCMC for " + toString(_numIterations) + " iterations.");
	_thinning = parameters().getParameterWithDefault<int>("thinning", 1);
	if (_thinning < 1 || _thinning > _numIterations) UERROR("Thinning must be > 1 and < number iterations!");
	if (_thinning > 1) {
		if (_thinning == 2)
			logfile().list("Will print every second iterations to the output file (thinning = 2)");
		else if (_thinning == 3)
			logfile().list("Will print every third iterations to the output file (thinning = 3)");
		else
			logfile().list("Will print every " + toString(_thinning) +
						   "th iterations to the output file (thinning = " + toString(_thinning) + ")");
	}

	// normal prior on ratio phi = log(theta_1 / theta_2)
	_phiPriorMean          = parameters().getParameterWithDefault("phiPriorMean", 0.0);
	_phiPriorVar           = parameters().getParameterWithDefault("phiPriorVar", 1.0);
	_phiPriorOneOverTwoVar = 1.0 / 2.0 / _phiPriorVar;
	logfile().list("Will assume a normal prior on phi ~ N(" + toString(_phiPriorMean) + ", " + toString(_phiPriorVar) +
				   ").");

	// proposal kernel
	_sdProposalKernelTheta1    = parameters().getParameterWithDefault("sdProposalTheta", 0.1);
	_sdProposalKernelTheta2    = _sdProposalKernelTheta1;
	_sdProposalKernelBaseFreq1 = parameters().getParameterWithDefault("sdProposalFreq", 0.01);
	_sdProposalKernelBaseFreq2 = _sdProposalKernelBaseFreq1;
	logfile().list("Will use initial proposal kernel standard deviations of " + toString(_sdProposalKernelTheta1) +
				   " and " + toString(_sdProposalKernelBaseFreq1) + " for thetas and base frequencies, respectively.");
	logfile().endIndent();

	// params regarding initial search
	_readParametersRegardingInitialSearch();
};

void TThetaEstimatorRatio::_clearCounters() {
	_numAcceptedTheta1    = 0;
	_numAcceptedTheta2    = 0;
	_numAcceptedBaseFreq1 = 0;
	_numAcceptedBaseFreq2 = 0;
};

void TThetaEstimatorRatio::_concludeAcceptanceRate(int numAccepted, int length, const std::string &name) {
	double acceptanceRate = (double)numAccepted / (double)length;
	logfile().conclude("Acceptance rate " + name + " = " + toString(acceptanceRate));
};

void TThetaEstimatorRatio::_concludeAcceptanceRateUpdateProposal(int numAccepted, int length, double &sd,
																const std::string &name) {
	double acceptanceRate = (double)numAccepted / (double)length;
	sd *= acceptanceRate * 3.0;
	logfile().conclude("Acceptance rate " + name + " = " + toString(acceptanceRate) + " (updated proposal sd to " +
					   toString(sd) + ")");
};

void TThetaEstimatorRatio::_concludeAcceptanceRates(int length) {
	_concludeAcceptanceRate(_numAcceptedTheta1, length, "theta 1");
	_concludeAcceptanceRate(_numAcceptedTheta2, length, "theta 2");
	_concludeAcceptanceRate(_numAcceptedBaseFreq1, length, "base frequencies 1");
	_concludeAcceptanceRate(_numAcceptedBaseFreq2, length, "base frequencies 1");
};

void TThetaEstimatorRatio::_concludeAcceptanceRatesUpdateProposal(int length) {
	_concludeAcceptanceRateUpdateProposal(_numAcceptedTheta1, length, _sdProposalKernelTheta1, "theta 1");
	_concludeAcceptanceRateUpdateProposal(_numAcceptedTheta2, length, _sdProposalKernelTheta2, "theta 2");
	_concludeAcceptanceRateUpdateProposal(_numAcceptedBaseFreq1, length, _sdProposalKernelBaseFreq1, "base frequencies 1");
	_concludeAcceptanceRateUpdateProposal(_numAcceptedBaseFreq2, length, _sdProposalKernelBaseFreq2, "base frequencies 1");
};

void TThetaEstimatorRatio::estimateRatio(std::string ouputName) {
	logfile().startIndent("Running MCMC to estimate phi = log(theta1 / theta2):");

	// check if there is sufficient data
	logfile().list(toString(_data->sizeWithData()) + " sites with data available for region 1.");
	if (_data->sizeWithData() < _minSitesWithData) UERROR("Not enough sites for region 1!");
	logfile().list(toString(_data2->sizeWithData()) + " sites with data available for region 2.");
	if (_data2->sizeWithData() < _minSitesWithData) UERROR("Not enough sites for region 2!");

	// get good starting values
	_findGoodStartingTheta(_data.get(), _theta, " region 1");
	_findGoodStartingTheta(_data2.get(), _theta2, " region 2");

	// first run burnin
	_clearCounters();
	if (_burnin > 0) {
		int oldProg                = 0;
		std::string progressString = "Running burnin of length " + toString(_burnin) + " ...";
		logfile().listFlush(progressString + "(0%)");
		for (int i = 0; i < _burnin; ++i) {
			_oneMCMCIteration();

			// print progress
			int prog = (double)i / (double)_burnin * 100.0;
			if (prog > oldProg) {
				oldProg = prog;
				logfile().overListFlush(progressString, "(", oldProg, "%)");
			}
		}
		logfile().overList(progressString + "done!   ");
		_concludeAcceptanceRatesUpdateProposal(_burnin);
	}

	// open MCMC output file
	ouputName += "_thetaRatioMCMC.txt.gz";
	logfile().list("Will write MCMC chain to file '" + ouputName + "'.");
	gz::ogzstream out(ouputName.c_str());
	if (!out) UERROR("Failed to open file '", ouputName, "' for writing!");

	// write header
	out << "log_theta_1\tlog_theta_2\tlog_phi\n";

	// now run chain with sampling
	_clearCounters();
	int oldProg                = 0;
	std::string progressString = "Running MCMC chain of length " + toString(_numIterations) + " ...";
	logfile().listFlush(progressString + "(0%)");
	for (int i = 0; i < _numIterations; ++i) {
		_oneMCMCIteration();

		// print to file
		if (i % _thinning == 0) {
			out << _theta.logTheta << "\t" << _theta2.logTheta << "\t" << _theta.logTheta - _theta2.logTheta << "\n";
		}

		// print progress
		int prog = (double)i / (double)_numIterations * 100.0;
		if (prog > oldProg) {
			oldProg = prog;
			logfile().overListFlush(progressString, "(", oldProg, "%)");
		}
	}
	logfile().overList(progressString + " done!   ");
	_concludeAcceptanceRates(_numIterations);

	// clean up
	logfile().endIndent();
	out.close();
};

bool TThetaEstimatorRatio::_updateTheta(TThetaEstimatorData *thisData, Theta &thisTheta, double otherLogThetaMean,
									   double thisSdProposalKernel) {
	// propose
	double newLogTheta = randomGenerator().getNormalRandom(thisTheta.logTheta, thisSdProposalKernel);
	double newExpTheta = exp(-exp(newLogTheta)); // we update log(theta) but need exp(-theta)

	// calc LL
	const auto _pGenotype   = getPGenotype(newExpTheta, thisTheta.baseFreq);
	double newLL = thisData->calcLogLikelihood(_pGenotype);

	// calc hastings ratio with prior
	// we use a uniform prior on log(theta)
	// and normal on log(phi) = log(theta) - log(theta2)
	double logH               = newLL - thisTheta.LL;
	double newLogPhiMinusMean = newLogTheta - otherLogThetaMean;
	double oldLogPhiMinusMean = thisTheta.logTheta - otherLogThetaMean;
	logH += _phiPriorOneOverTwoVar * (newLogPhiMinusMean * newLogPhiMinusMean - oldLogPhiMinusMean * oldLogPhiMinusMean);

	// accept or reject
	if (log(randomGenerator().getRand()) < logH) {
		thisTheta.setLog(newLogTheta, newLL);
		return true;
	} else
		return false;
};

bool TThetaEstimatorRatio::_updateBaseFrequencies(TThetaEstimatorData *thisData, Theta &thisTheta,
												 double thisSdProposalKernel) {
	using namespace genometools;
	// propose: select one frequency at random and shift this one
	// make sure frequencies are not outside [0,1]
	const auto thisBase = Base(randomGenerator().sample(4));
	TBaseLikelihoods tmpBaseLik;
	double tmp = thisTheta.baseFreq[thisBase].get() + randomGenerator().getNormalRandom(0.0, thisSdProposalKernel);
	if (tmp > 1.0) {
		tmpBaseLik[thisBase] = 2.0 - tmp;
	} else if (tmp < 0.0) {
		tmpBaseLik[thisBase] = -tmp;
	} else {
		tmpBaseLik[thisBase] = tmp;
	}

	// now scale all others so the sum will be 1.0
	double scale = (double)tmpBaseLik[thisBase].complement() / (double)thisTheta.baseFreq[thisBase].complement();
	for (Base k = Base::min; k < Base::max; ++k) {
		if (k != thisBase) { tmpBaseLik[k] = thisTheta.baseFreq[thisBase].get() * scale; }
	}

	const auto tmpBaseFreq = TBaseProbabilities::normalize(tmpBaseLik);

	// calc LL & hastings ratio (use uniform prior, i.e. all combinations are equally likely)
	const auto _pGenotype   = getPGenotype(thisTheta.expMinusTheta, tmpBaseFreq);
	double newLL = thisData->calcLogLikelihood(_pGenotype);
	double logH  = newLL - thisTheta.LL;

	// accept or reject
	if (log(randomGenerator().getRand()) < logH) {
		thisTheta.baseFreq = tmpBaseFreq;
		thisTheta.LL       = newLL;
		return true;
	} else
		return false;
};

void TThetaEstimatorRatio::_oneMCMCIteration() {
	// update theta
	_numAcceptedTheta1 += _updateTheta(_data.get(), _theta, _theta2.logTheta + _phiPriorMean, _sdProposalKernelTheta1);
	_numAcceptedTheta2 += _updateTheta(_data2.get(), _theta2, _theta.logTheta - _phiPriorMean, _sdProposalKernelTheta2);

	// update base frequencies
	_numAcceptedBaseFreq1 += _updateBaseFrequencies(_data.get(), _theta, _sdProposalKernelBaseFreq1);
	_numAcceptedBaseFreq2 += _updateBaseFrequencies(_data2.get(), _theta2, _sdProposalKernelBaseFreq2);
};

}; // namespace GenotypeLikelihoods
