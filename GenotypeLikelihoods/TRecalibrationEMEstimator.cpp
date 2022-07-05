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
#include <numeric>
#include <stdexcept>

#include "GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TLog.h"
#include "TParameters.h"
#include "TPostMortemDamage.h"
#include "TSequencedBase.h"
#include "TSequencingErrorModel.h"
#include "TSequencingErrorModels.h"
#include "algorithms.h"
#include "probability.h"
#include "stringFunctions.h"

#include "devtools.h"

namespace GenotypeLikelihoods {

namespace SequencingError {

using coretools::instances::logfile;
using coretools::instances::parameters;
using genometools::Genotype;
using genometools::Base;

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
//--------------------------------------------------------------------------

void TModelVectorForEstimation::reset(TModels &SequencingErrorModels,
									  const RecalEstimatorTools::TRecalDataTables &DataTables,
									  const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
									  uint32_t MinRequiredObservations) {
	using MS = RecalEstimatorTools::ModelStatusTypes;
	// Copy models that are 1) in use after pooling and 2) have data.
	// Note: data table is already pooled!

	// prepare storage for reporting
	_modelIndex.clear();
	_models.clear();
	_modelIndex.resize(ReadGroups.size());
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

TBaseLikelihoods TModelVectorForEstimation::getBaseLikelihoods(const BAM::TSequencedBase &data) const {
	return _modelIndex[data.readGroupID][data.isSecondMate()]->getBaseLikelihoods(data);
};

//-------------------------------------------------------------------
// functions to estimate rho
void TModelVectorForEstimation::prepareRhoEstimationFromEMWeights() {
	for (auto &model : _models) { model->prepareRhoEstimationFromEMWeights(); }
};

void TModelVectorForEstimation::addBaseForRhoEstimation(BAM::TSequencedBase &data, const TBaseLikelihoods &EMWeights) {
	_modelIndex[data.readGroupID][data.isSecondMate()]->addBaseForRhoEstimation(data, EMWeights);
};

void TModelVectorForEstimation::estimateRho() {
	for (auto &model : _models) { model->estimateRho(); }
};

// functions to estimate beta
//-------------------------------------------------------------------

void TModelVectorForEstimation::addToQFJ(const BAM::TSequencedBase &data, coretools::Probability p_g_I_d, coretools::Probability p_bbar_I_gd) {
		_modelIndex[data.readGroupID][data.isSecondMate()]->addToQFJ(data, p_g_I_d, p_bbar_I_gd);
};

void TModelVectorForEstimation::setQJF_0() {
	for (auto &model : _models) {
		model->setQFJ_0();
	}
};

double TModelVectorForEstimation::curQ() {
	return std::accumulate(_models.begin(), _models.end(), 0.0, [](auto tot, const auto &val) { return tot + val->curQ(); });
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
	if (_models.empty()) return std::string{};

	return std::accumulate(_models.begin() + 1, _models.end(), _models.front()->getCovariateDefinition(),
						   [](auto tot, auto &model) { return tot + '\n' + model->getCovariateDefinition(); });
};

std::string TModelVectorForEstimation::getRhoDefinition() {
	if (_models.empty()) return std::string{};

	return std::accumulate(_models.begin() + 1, _models.end(), _models.front()->getRhoDefinition(),
						   [](auto tot, auto &model) { return tot + '\n' + model->getRhoDefinition(); });
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
	_genoDist = std::make_unique<THaploidDistribution>();
	logfile().list("Will use a ", _genoDist->typeString(), " genotype distribution.");

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
	_modelsToEstimate.reset(SequencingErrorModels, _dataTables, *_readGroups, *_readGroupMap, _minRequiredObservations);
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
	const std::string filename = outputName + "_recalibrationEM.txt";
	logfile().listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename);
	logfile().done();
};

void TRecalibrationEMEstimator::_estimateRho_updatePij(const TPostMortemDamage &PmdModels) {
	using genometools::genotype;
	_Pijs.clear();
	for (size_t i = 0; i < _sites.size(); ++i) {
		for (auto &d_ij : _sites[i]) {
			_Pijs.emplace_back(0.);
			auto& Pij = _Pijs.back();
			const auto L_eps = _modelsToEstimate.getBaseLikelihoods(d_ij);
			for (auto a = Base::min; a < Base::max; ++a) {
				const auto g_aa = genotype(a, a);
				const TBaseLikelihoods lk_a{PmdModels.getMassFunction(a, d_ij, L_eps)};

				Pij[g_aa] = lk_a[d_ij.base];

				_modelsToEstimate.addBaseForRhoEstimation(d_ij, _Pis[i][g_aa] * lk_a);
				if (!_genoDist->isInvariant()) {
					for (auto b = genometools::next(a); b < Base::max; ++b) {
						const auto g_ab = genotype(a, b);
						constexpr auto p05 = coretools::Probability(0.5);
						const TBaseLikelihoods lk_ab =
							p05 * lk_a + p05 * TBaseLikelihoods(PmdModels.getMassFunction(b, d_ij, L_eps));

						Pij[genotype(a, b)] = lk_ab[d_ij.base];

						_modelsToEstimate.addBaseForRhoEstimation(d_ij, _Pis[i][g_ab] * lk_ab);
					}
				}
			}
		}
	}
	_modelsToEstimate.estimateRho();
}

double TRecalibrationEMEstimator::_calculate_Q_updateJF() {
	_modelsToEstimate.setQJF_0();

	size_t ij = 0;
	for (size_t i = 0; i < _sites.size(); ++i) {
		const auto& Pi = _Pis[i];
		for (auto &d_ij : _sites[i]) {
			const auto &Pij = _Pijs[ij++];
			for (auto a = Base::min; a < Base::max; ++a) {
				const auto g_aa = genometools::genotype(a, a);
				_modelsToEstimate.addToQFJ(d_ij, Pi[g_aa], Pij[g_aa]);
				if (!_genoDist->isInvariant()) {
					for (auto b = genometools::next(a); b < Base::max; ++b) {
						const auto g_ab = genometools::genotype(a, b);
						_modelsToEstimate.addToQFJ(d_ij, Pi[g_ab], Pij[g_ab]);
					}
				}
			}
		}
	}

	// return total Q
	return _modelsToEstimate.curQ();
};

void TRecalibrationEMEstimator::_updateEpsilon(const TPostMortemDamage &PmdModels) {
	using coretools::str::toString;
	logfile().startIndent("Updating sequencing error models (theta_epsilon):");

	// 1) Weights are calculated during LL calculation

	// 2) update rho
	//-------------------------
	logfile().listFlushDots("Updating rho");
	_estimateRho_updatePij(PmdModels);
	logfile().done();
	logfile().conclude("rho = ", _modelsToEstimate.getRhoDefinition());

	// 3) Calculate Q_beta at current location
	//-------------------------
	logfile().listFlushDots("Calculating Q_beta at current location");
	const double curQ = _calculate_Q_updateJF();
	logfile().done();
	logfile().conclude("Q_beta = ", curQ);

	// 4) Use Newton-Raphson to optimize
	//-------------------------
	logfile().startIndent("Optimizing Q_beta using a Newton-Raphson algorithm:");

	const auto nTot = _modelsToEstimate.size();

	for (int i = 0; i < _NewtonRaphsonNumIterations; ++i) {
		logfile().startIndent("Running Newton-Raphson iteration " + toString(i + 1) + ":");

		// a) fill Jacobin and F
		logfile().listFlush("Calculating Jacobian and gradient ...");
		_modelsToEstimate.solveJxF();
		logfile().done();

		// b) update parameters using backtracking
		double lambda   = 1.0;
		int log2_lambda = 0;
		size_t nUpdated = 0;

		while (nUpdated < nTot && lambda > 1.0E-20) {
			// propose move
			logfile().listFlush("Proposing move with log2(lambda) = " + toString(log2_lambda) + " ... ");
			_modelsToEstimate.proposeNewParameters(lambda);
			const auto pModels = _modelsToEstimate.getModelsDefinition();

			// calculate Q at new location
			const double Q = _calculate_Q_updateJF();

			// check if we accept or backtrack
			const auto numUpdatedModels_old = nUpdated;
			nUpdated                = _modelsToEstimate.acceptProposedParametersBasedOnQ();
			logfile().write(toString(nUpdated) + "/" + toString(nTot) +
			                " models converged.");

			logfile().conclude("Model = ", pModels);
			logfile().conclude("Q_beta = ", Q);

			if (nUpdated > numUpdatedModels_old) {
				logfile().conclude("Q_beta was increased from " + toString(curQ) + " to " + toString(Q));
			}

			// backtrack
			lambda = lambda / 2.0; // backtrack;
			--log2_lambda;
		}

		if (nUpdated < nTot) {
			logfile().conclude("Some models did not improve even with log2(lambda) = " + toString(log2_lambda) +
							   ", aborting Newton-Raphson.");
		}

		// c) adjust parameters post estimation
		_modelsToEstimate.adjustParametersPostEstimation();

		// d) get largest gradient (F) to check if we break NR optimization
		double maxF = _modelsToEstimate.getSteepestGradient();
		logfile().conclude("max(F) = " + toString(maxF));
		logfile().endIndent();
		if (maxF < _NewtonRaphsonMaxF || nUpdated == 0) break;

		logfile().endIndent();
	}
	logfile().endIndent();
	logfile().endIndent();
};

double TRecalibrationEMEstimator::_calculate_LL_updatePi(const TPostMortemDamage &PmdModels) {
	_Pis.clear();
	double LL = 0.0;
	for (auto &s_i : _sites) {
		if (s_i.genotype == Genotype::NN) {
			_Pis.emplace_back(1.); // Start at 1,1,1,1,1,1,1,1
			auto &L = _Pis.back();
			for (auto &d_ij : s_i) {
				const auto L_eps = _modelsToEstimate.getBaseLikelihoods(d_ij);
				const auto L_D   = PmdModels.getBaseLikelihoods(d_ij, L_eps);
				L *= _genoDist->getGenotypeLikelihoods(L_D);
			}
			LL += log(_genoDist->normalize(L));
			// unknow genotype
		} else {
			// know genotype, probability of correct genotype is 1
			_Pis.emplace_back(0.); 
			_Pis.back()[s_i.genotype] = 1; 
			double L = 1.;
			for (auto &d_ij : s_i) {
				const auto L_eps = _modelsToEstimate.getBaseLikelihoods(d_ij);
				const auto L_D   = PmdModels.getBaseLikelihoods(d_ij, L_eps);
				L *= _genoDist->getGenotypeLikelihood(L_D, s_i.genotype);
			}
			LL += log(L);
		}
	}

	return LL;
};

void TRecalibrationEMEstimator::_runEM(const std::string &outputName, const TPostMortemDamage &PmdModels) {
	using coretools::str::toString;
	// run EM
	logfile().startNumbering("Running EM algorithm:");

	// calculate initial LL
	double oldLL = _calculate_LL_updatePi(PmdModels);
	logfile().conclude("Initial log Likelihood = " + toString(oldLL));

	// running iterations
	for (int i = 0; i < _numEMIterations; ++i) {
		logfile().number("EM Iteration:");
		logfile().addIndent();

		// update theta_epsilon (sequencing errors)
		_updateEpsilon(PmdModels);

		// calculate LL
		const double LL = _calculate_LL_updatePi(PmdModels);
		logfile().conclude("Current Log Likelihood = " + toString(LL));

		// check if we break based on LL
		double deltaLL = LL - oldLL;
		logfile().conclude("delta LL = " + toString(deltaLL));
		if (i > 0 && deltaLL < _minDeltaLL) {
			logfile().conclude("EM has converged (delta LL < " + toString(_minDeltaLL) + ")");
			break;
		}
		oldLL = LL;

		// write current estimates to file
		if (_writeTmpTables) {
			std::string filename = outputName + "_recalibrationEM_Loop" + toString(i) + ".txt";
			logfile().listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename);
			logfile().done();
		}

		// end loop
		logfile().endIndent();
		if (i == _numEMIterations - 1) logfile().warning("EM has not converged after maximum number of iterations!");
	}

	// finalize
	logfile().endNumbering();
};

//----------------------------
// Other functions
//----------------------------
void TRecalibrationEMEstimator::writeCurrentEstimates(const std::string & Filename) {
	// open file and write header
	_modelsToEstimate.writeRecalFile(*_readGroups, Filename);
};

}; // namespace SequencingError

}; // end namespace GenotypeLikelihoods
