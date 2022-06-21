/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMEstimator.h"
#include <algorithm>
#include <armadillo>
#include <exception>
#include <functional>
#include <iostream>
#include <math.h>
#include <stdexcept>

#include "GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TGenotypeData.h"
#include "TLog.h"
#include "TParameters.h"
#include "TPostMortemDamage.h"
#include "TSequencedBase.h"
#include "TSequencingErrorModel.h"
#include "TSequencingErrorModels.h"
#include "probability.h"
#include "stringFunctions.h"

#include "devtools.h"

namespace GenotypeLikelihoods {

namespace SequencingError {

using coretools::instances::logfile;
using coretools::instances::parameters;

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
//--------------------------------------------------------------------------

TModelVectorForEstimation::TModelVectorForEstimation(TModels &SequencingErrorModels,
													 const RecalEstimatorTools::TRecalDataTables &DataTables,
													 const BAM::TReadGroups &ReadGroups,
													 const BAM::TReadGroupMap &ReadGroupMap,
													 uint32_t MinRequiredObservations)
	: _modelIndex(ReadGroups.size()) {
	using MS = RecalEstimatorTools::ModelStatusTypes;
	// Copy models that are 1) in use after pooling and 2) have data.
	// Note: data table is already pooled!

	// prepare storage for reporting
	RecalEstimatorTools::TModelStati modelStati;

	// copy models
	for (auto &r : ReadGroupMap.readGroupsInUse()) {
		// add to model stati
		modelStati.add(r);

		// loop over mates
		for (uint8_t mate = 0; mate < 2; ++mate) {
			const RecalEstimatorTools::TRecalDataTable &table = DataTables[r][mate];

			// check if there is sufficient data for _first mate
			if (table.size() > 0) {
				// check if model is estimatable
				if (SequencingErrorModels(r, mate).estimatable()) {
					// copy model and update index
					auto model = SequencingErrorModels.getRecal(r, mate);
					model->checkOrInit(table);
					_models.push_back(model);
					modelStati[r][MS::copied].set(mate);
					for (auto &rr : ReadGroupMap.readGroupsPooledWith(r)) { _modelIndex[rr][mate] = model; }

					// check if there is limited data
					if (table.size() < MinRequiredObservations) { modelStati[r][MS::littleData].set(mate); }
				} else {
					modelStati[r][MS::dataButNoRecal].set(mate);
				}
			} else if (SequencingErrorModels(r, mate).estimatable()) {
				modelStati[r][MS::noData].set(mate);
			}
		}
	}

	// report models that will be estimated
	modelStati.report(MS::copied, "Read groups for which models will be estimated:", ReadGroups);
	modelStati.report(MS::noData, "Read groups excluded because they have no data:", ReadGroups);
	modelStati.report(MS::dataButNoRecal, "Read groups with data but no recal model:", ReadGroups);
	if (modelStati.num(MS::copied) == 0) { throw "No recal models need estimation!"; }
	modelStati.report(MS::littleData, "Read groups with very little data (consider pooling):", ReadGroups);
};

TBaseLikelihoods TModelVectorForEstimation::getBaseLikelihoods(const BAM::TSequencedBase &base) const {
	return _modelIndex[base.readGroupID][base.isSecondMate()]->getBaseLikelihoods(base);
};

// functions to estimate rho
//-------------------------------------------------------------------
// functions to estimate rho
void TModelVectorForEstimation::prepareRhoEstimationFromEMWeights() {
	for (auto &model : _models) { model->prepareRhoEstimationFromEMWeights(); }
};

void TModelVectorForEstimation::addBaseForRhoEstimation(BAM::TSequencedBase &base, const TBaseLikelihoods &EMWeights) {
	_modelIndex[base.readGroupID][base.isSecondMate()]->addBaseForRhoEstimation(base, EMWeights);
};

void TModelVectorForEstimation::estimateRho() {
	for (auto &model : _models) { model->estimateRho(); }
};

// functions to estimate beta
//-------------------------------------------------------------------

void TModelVectorForEstimation::addToFandJacobian(const BAM::TSequencedBase &base,
												  const TBaseLikelihoods &EM_weights_bbar_given_d) {
	_modelIndex[base.readGroupID][base.isSecondMate()]->addToFandJacobian(base, EM_weights_bbar_given_d);
};

void TModelVectorForEstimation::addToQ(const BAM::TSequencedBase &base,
									   const TBaseLikelihoods &EM_weights_bbar_given_d) {
	_modelIndex[base.readGroupID][base.isSecondMate()]->addToQ(base, EM_weights_bbar_given_d);
};

void TModelVectorForEstimation::setNewtonRaphsonParamsToZero() {
	for (auto &model : _models) { model->setNewtonRaphsonParamsToZero(); }
};

void TModelVectorForEstimation::setQToZero() {
	for (auto &model : _models) { model->setQToZero(); }
};

double TModelVectorForEstimation::curQ() {
	double Q = 0.0;
	for (auto &model : _models) { Q += model->curQ(); }
	return Q;
};

void TModelVectorForEstimation::solveJxF() {
	for (auto &model : _models) {
		if (!model->solveJxF()) {
			UERROR("Issue solving JxF! This may be due to a lack of data. Consider adding more sites. Jacobian: ");
		}
	}
};

void TModelVectorForEstimation::proposeNewParameters(double lambda) {
	for (auto &model : _models) { model->proposeNewParameters(lambda); }
};

unsigned int TModelVectorForEstimation::acceptProposedParametersBasedOnQ() {
	unsigned int numAccepted = 0;
	for (auto &model : _models) { numAccepted += (unsigned int)model->acceptProposedParametersBasedOnQ(); }
	return numAccepted;
};

void TModelVectorForEstimation::adjustParametersPostEstimation() {
	for (auto &model : _models) { model->adjustParametersPostEstimation(); }
};

double TModelVectorForEstimation::getSteepestGradient() {
	double maxF = 0.0;
	for (auto &model : _models) {
		double tmp = model->getSteepestGradient();
		if (fabs(tmp) > maxF) maxF = fabs(tmp);
	}
	return maxF;
};

void TModelVectorForEstimation::writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const {
	TModelNoRecal _noRecal;
	// open file and write header
	coretools::TOutputFile out(Filename);
	out.writeHeader({"readGroup", "mate", "covariates", "rho"});

	// add models
	for (uint16_t r = 0; r < ReadGroups.size(); ++r) {
		for (uint8_t mate = 0; mate < 2; ++mate) {
			out << ReadGroups.getName(r) << std::array{"first", "second"}[mate];
			if (_modelIndex[r][mate])
				out << _modelIndex[r][mate]->getCovariateDefinition() << _modelIndex[r][mate]->getRhoDefinition();
			else 
				out << _noRecal.getCovariateDefinition() << _noRecal.getRhoDefinition();
			out << std::endl;
		}
	}
}

std::string TModelVectorForEstimation::getModelsDefinition() {
	std::string s;
	for (auto &model : _models) {
		s += model->getCovariateDefinition() + '\t' + model->getRhoDefinition() + '\n';
	}
	return s;
};

//---------------------------------------------------------------
// TRecalibrationEMEstimator
//---------------------------------------------------------------
TRecalibrationEMEstimator::TRecalibrationEMEstimator(const BAM::TReadGroups *ReadGroups,
													 const BAM::TReadGroupMap *ReadGroupMap) {
	// read groups
	_readGroups   = ReadGroups;
	_readGroupMap = ReadGroupMap;

	// genotype distribution: currently only allow for haploid
	_genoDist = std::make_unique<TGenotypeDistribution_haploid>();

	// estimation parameters
	logfile().startIndent("Settings regarding the EM algorithm:");

	_minRequiredObservations = 10000; // constant for reporting
	_numEMIterations         = parameters().getParameterWithDefault<int>("iterations", 200);
	logfile().list("Will perform at max ", _numEMIterations, " EM iterations.");
	_minDeltaLL = parameters().getParameterWithDefault<double>("minDeltaLL", 0.000001);
	logfile().list("Will stop EM when deltaLL < ", _minDeltaLL, ".");
	_NewtonRaphsonNumIterations = parameters().getParameterWithDefault<int>("NRiterations", 20);
	logfile().list("Will conduct at max ", _NewtonRaphsonNumIterations, " Newton-Raphson iterations");
	_NewtonRaphsonMaxF = parameters().getParameterWithDefault<double>("maxF", 0.0001);
	logfile().list("Will stop Newton-Raphson when F < ", _NewtonRaphsonMaxF, ".");

	// base frequency model
	equalBaseFrequencies = true;
	if (parameters().parameterExists("estimateBaseFrequencies")) {
		equalBaseFrequencies = false;
		logfile().list("Will estimate the base frequencies. (parameter ''estimateBaseFrequencies)");

		// TODO: implement estimation of genotype distribution!
		throw std::runtime_error("Estimation of genotype distribution not yet implemented!");
	} else if (equalBaseFrequencies) {
		logfile().list("Will assume equal base frequencies {0.25, 0.25, 0.25, 0.25}. (use 'estimateBaseFrequencies' to "
					   "estimate them)");
	}

	// write tmp tables?
	if (parameters().parameterExists("writeTmpTables")) {
		_writeTmpTables = true;
		logfile().list(
			"Will write intermediate estimates of EM and Newton-Raphson to file. (parameter 'writeTmpTables')");
	} else {
		_writeTmpTables = false;
		logfile().list("Will not write intermediate estimates. (use 'writeTmpTables' to get this debug information)");
	}

	logfile().endIndent();
};

//-----------------------------------------------
// Function to initialize models used for estimation
//-----------------------------------------------
size_t TRecalibrationEMEstimator::_numSitesDepthTwoOrMore() {
	size_t _numSites = 0;
	for (auto &s : _sites) {
		if (s.depth() > 1) ++_numSites;
	}
	return _numSites;
};

void TRecalibrationEMEstimator::_initializeModels(SequencingError::TModels &SequencingErrorModels) {
	using coretools::str::toString;
	// count data available for recal
	logfile().listFlush("Counting data available for recal ...");
	// Note: data tables pool read groups!
	_dataTables.initialize(_readGroups, _readGroupMap);
	_dataTables.add(_sites);
	logfile().done();

	logfile().conclude("Number of sites with data: " + toString(_sites.size()));
	size_t numSitesDepthTwoOrMore = _numSitesDepthTwoOrMore();
	logfile().conclude("Number of sites with depth > 1: " + toString(numSitesDepthTwoOrMore));
	logfile().conclude("Number of bases: " + toString(_dataTables.size()));
	if (numSitesDepthTwoOrMore < 100) throw "Less than 100 sites with depth >= 2 available - aborting estimation!";

	// identify models with data that can be estimated
	logfile().startIndent("Identifying models to estimate:");
	_modelsToEstimate = std::make_unique<TModelVectorForEstimation>(SequencingErrorModels, _dataTables, *_readGroups,
																	*_readGroupMap, _minRequiredObservations);
	logfile().endIndent();
};

//----------------------------
// Functions to add data
//----------------------------
void TRecalibrationEMEstimator::addSite(const TSite &site) {
	if (!site.empty()) { _sites.emplace_back(site); }
};

//----------------------------
// Functions for estimation
//----------------------------
void TRecalibrationEMEstimator::performEstimation(const std::string &outputName,
												  SequencingError::TModels &SequencingErrorModels,
												  const TPostMortemDamage &PmdModels) {
	// initialize models
	_initializeModels(SequencingErrorModels);

	// run EM
	_runEM(outputName, PmdModels);

	// writing final estimates
	std::string filename = outputName + "_recalibrationEM.txt";
	logfile().listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename);
	logfile().done();
};

std::vector<TBaseLikelihoods> TRecalibrationEMEstimator::_calculate_EMWeights_epsilon(const TPostMortemDamage &PmdModels) {
	using genometools::Base;
	// make sure EM-weight storage is of appropriate size
	std::vector<TBaseLikelihoods> EMWeights;
	EMWeights.reserve(_dataTables.size());
	// loop over all bases and calculate EM-weights
	for (auto &s : _sites) {
		// get relevant base frequencies P(b): from known genotype or distribution if genotype is unknown
		const TBaseLikelihoods baseFreq{getBaseFrequences(s.genotype)};

		// calculate weights per base
		for (auto &b : s) {
			// calculate P(bbar) = \sum_b P(bbar|b)P(b)
			const auto PMD = PmdModels.getBaseLikelihoods(b, baseFreq);

			// calculate P(d|bbar)
			EMWeights.push_back(_modelsToEstimate->getBaseLikelihoods(b));
			auto &EMi = EMWeights.back();

			// calculate P(d|bbar) \propto P(d|bbar)P(bbar)
			EMi *= PMD;
			normalize(EMi);
		}
	}
	return EMWeights;
};

std::vector<TBaseLikelihoods> TRecalibrationEMEstimator::_updateRho_getWeights(const TPostMortemDamage &PmdModels) {
	std::vector<TBaseLikelihoods> weights;
	weights.reserve(_dataTables.size());
	for (auto &s_i : _sites) {
		for (auto &d_ij : s_i) {
			const auto P_d_I_bbar = _modelsToEstimate->getBaseLikelihoods(d_ij);
			const auto P_bbar_I_b = PmdModels.getBaseLikelihoods(d_ij, P_d_I_bbar);
			const auto P_b_I_g    = getBaseFrequences(s_i.genotype);
		}
	}
	return weights;
}

double TRecalibrationEMEstimator::_calculate_Q_beta(const std::vector<TBaseLikelihoods> &EM_weights_bbar_given_d) {
	_modelsToEstimate->setQToZero();

	// loop over all bases and add them to Q
	size_t index = 0;
	for (auto &s : _sites) {
		for (auto &b : s) {
			_modelsToEstimate->addToQ(b, EM_weights_bbar_given_d[index]);
			++index;
		}
	}

	// return total Q
	return _modelsToEstimate->curQ();
};

void TRecalibrationEMEstimator::_calculate_J_F_beta(const std::vector<TBaseLikelihoods> &EM_weights_bbar_given_d) {
	logfile().listFlush("Calculating Jacobian and gradient ...");
	_modelsToEstimate->setNewtonRaphsonParamsToZero();

	size_t index = 0;
	for (auto &s : _sites) {
		for (auto &b : s) {
			_modelsToEstimate->addToFandJacobian(b, EM_weights_bbar_given_d[index]);
			++index;
		}
	}

	// solve J^-1 x F
	_modelsToEstimate->solveJxF();
	logfile().done();
};

void TRecalibrationEMEstimator::_updateEM_theta_epsilon(const TPostMortemDamage &PmdModels) {
	using coretools::str::toString;
	logfile().startIndent("Updating sequencing error models (theta_epsilon):");

	// 1) calculate EM weights
	//-------------------------
	logfile().listFlushDots("Calculating EM weights");
	const auto EM_weights_bbar_given_d = _calculate_EMWeights_epsilon(PmdModels);
	logfile().done();
	OUT(EM_weights_bbar_given_d.size());

	// 2) update rho
	//-------------------------
	/*
	logfile().listFlushDots("Updating rho");
	_sequencingErrorModels.prepareRhoEstimationFromEMWeights();
	size_t index = 0;
	for(auto& s : _sites){
		for(auto& b : s){
		_sequencingErrorModels.addBaseForRhoEstimation(b, EM_weights_bbar_given_d[index]);
		++index;
		}
	}
	_sequencingErrorModels.estimateRho();
	logfile().done();
	*/

	// 3) Calculate Q_beta at current location
	//-------------------------
	logfile().listFlushDots("Calculating Q_beta at current location");
	double curQ = _calculate_Q_beta(EM_weights_bbar_given_d);
	logfile().done();
	logfile().conclude("Q_beta = ", curQ);

	// 4) Use Newton-Raphson to optimize
	//-------------------------
	logfile().startIndent("Optimizing Q_beta using a Newton-Raphson algorithm:");

	for (int i = 0; i < _NewtonRaphsonNumIterations; ++i) {
		logfile().startIndent("Running Newton-Raphson iteration " + toString(i + 1) + ":");

		// a) fill Jacobin and F
		_calculate_J_F_beta(EM_weights_bbar_given_d);

		// b) update parameters using backtracking
		double lambda           = 1.0;
		int log2_lambda         = 0;
		size_t numUpdatedModels = 0;
		size_t numUpdatedModels_old;

		while (numUpdatedModels < _modelsToEstimate->size() && lambda > 1.0E-20) {
			// propose move
			logfile().listFlush("Proposing move with log2(lambda) = " + toString(log2_lambda) + " ... ");
			_modelsToEstimate->proposeNewParameters(lambda);

			OUT(_modelsToEstimate->getModelsDefinition());

			// calculate Q at new location
			double Q = _calculate_Q_beta(EM_weights_bbar_given_d);

			OUT(curQ);
			OUT(Q);

			// check if we accept or backtrack
			numUpdatedModels_old = numUpdatedModels;
			numUpdatedModels     = _modelsToEstimate->acceptProposedParametersBasedOnQ();
			logfile().write(toString(numUpdatedModels) + "/" + toString(_modelsToEstimate->size()) +
							" models converged.");

			if (numUpdatedModels > numUpdatedModels_old) {
				logfile().conclude("Q_beta was increased from " + toString(curQ) + " to " + toString(Q));
			}

			// backtrack
			lambda = lambda / 2.0; // backtrack;
			--log2_lambda;
		}

		if (numUpdatedModels < _modelsToEstimate->size()) {
			logfile().conclude("Some models did not improve even with log2(lambda) = " + toString(log2_lambda) +
							   ", aborting Newton-Raphson.");
		}

		// c) adjust parameters post estimation
		_modelsToEstimate->adjustParametersPostEstimation();

		// d) get largest gradient (F) to check if we break NR optimization
		double maxF = _modelsToEstimate->getSteepestGradient();
		logfile().conclude("max(F) = " + toString(maxF));
		logfile().endIndent();
		if (maxF < _NewtonRaphsonMaxF || numUpdatedModels == 0) break;

		logfile().endIndent();
	}
	logfile().endIndent();
	logfile().endIndent();
};

double TRecalibrationEMEstimator::_calculateLL_fullModel(const TPostMortemDamage &PmdModels) {
	double LL                              = 0.0;
	const TGenotypeProbabilities &genoFreq = _genoDist->genotypeFrequencies();
	for (auto &s : _sites) {
		// calculate genotype likelihoods
		const TGenotypeLikelihoods genotypeLikelihoods = _genotypeLikelihoodCalculator.getGenotypeLikelihoods(s, PmdModels, *_modelsToEstimate);

		// weight by genotype prior
		if (s.genotype == genometools::Genotype::NN) {
			LL += log(weightedSum(genotypeLikelihoods, genoFreq));
		} else {
			LL += log(genotypeLikelihoods[s.genotype].get());
		}
	}

	return LL;
};

void TRecalibrationEMEstimator::_runEM(const std::string &outputName, const TPostMortemDamage &PmdModels) {
	using coretools::str::toString;
	// run EM
	logfile().startNumbering("Running EM algorithm:");

	std::ofstream out;
	std::string filename;

	// calculate initial LL
	double oldLL = _calculateLL_fullModel(PmdModels);
	logfile().conclude("Initial log Likelihood = " + toString(oldLL));

	// running iterations
	for (int iter = 0; iter < _numEMIterations; ++iter) {
		logfile().number("EM Iteration:");
		logfile().addIndent();

		// update theta_epsilon (sequencing errors)
		_updateEM_theta_epsilon(PmdModels);

		// calculate LL
		double LL = _calculateLL_fullModel(PmdModels);
		logfile().conclude("Current Log Likelihood = " + toString(LL));

		// check if we break based on LL
		double deltaLL = LL - oldLL;
		logfile().conclude("delta LL = " + toString(deltaLL));
		if (iter > 0 && deltaLL < _minDeltaLL) {
			logfile().conclude("EM has converged (delta LL < " + toString(_minDeltaLL) + ")");
			break;
		}
		oldLL = LL;

		// write current estimates to file
		if (_writeTmpTables) {
			filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
			logfile().listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename);
			logfile().done();
		}

		// end loop
		logfile().endIndent();
		if (iter == _numEMIterations - 1) logfile().warning("EM has not converged after maximum number of iterations!");
	}

	// finalize
	logfile().endNumbering();
};

//----------------------------
// Other functions
//----------------------------
void TRecalibrationEMEstimator::writeCurrentEstimates(const std::string & Filename) {
	// open file and write header
	_modelsToEstimate->writeRecalFile(*_readGroups, Filename);
};

double TRecalibrationEMEstimator::calcLL() {
	throw std::runtime_error("double TRecalibrationEMEstimator::calcLL() not yet implemented!");
};

void TRecalibrationEMEstimator::calcLikelihoodSurface(const std::string &, int) {
	throw std::runtime_error("void TRecalibrationEMEstimator::calcLikelihoodSurface(std::string filename, int "
							 "numMarginalGridPoints) not yet implemented!");
};

void TRecalibrationEMEstimator::calcQSurface(const std::string &, int) {
	throw std::runtime_error("TRecalibrationEMEstimator::calcQSurface(std::string filename, int numMarginalGridPoints) "
							 "not yet implemented!");
};

}; // namespace SequencingError

}; // end namespace GenotypeLikelihoods
