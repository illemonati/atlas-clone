/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "SequencingError/TRecalEstimator.h"
#include <algorithm>
#include <armadillo>
#include <exception>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <stdexcept>

#include "coretools/Containers/TStrongArray.h"
#include "genometools/GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "PMD/TModels.h"
#include "TSequencedBase.h"
#include "SequencingError/TModel.h"
#include "SequencingError/TModels.h"
#include "coretools/algorithms.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"

#include "coretools/devtools.h"

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
					TModelRecal *model = reinterpret_cast<TModelRecal*>(&SequencingErrorModels(r, mate));
					//auto model = SequencingErrorModels.getRecal(r, mate);
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
	if (modelStati.num(MS::copied) == 0) { UERROR("No recal models need estimation!"); }
	modelStati.report(MS::littleData, "Read groups with very little data (consider pooling):", ReadGroups);
};

TBaseLikelihoods TModelVectorForEstimation::getBaseLikelihoods(const BAM::TSequencedBase &data) const {
	return model(data)->baseLikelihoods(data);
};

//-------------------------------------------------------------------
// functions to estimate rho

void TModelVectorForEstimation::addToRho(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) {
	model(data)->addToRho(data, P_g_I_d, P_bbar_I_d);
};

void TModelVectorForEstimation::estimateRho() {
	for (auto &model : _models) { model->estimateRho(); }
};

// functions to estimate beta
//-------------------------------------------------------------------

double TModelVectorForEstimation::Q() const {
	return std::accumulate(_models.begin(), _models.end(), 0.0, [](auto tot, const auto &val) { return tot + val->Q(); });
};

void TModelVectorForEstimation::solveJxF() {
	for (auto &model : _models) {
		model->solveJxF();
	}
};

void TModelVectorForEstimation::propose(double lambda) {
	for (auto &model : _models) { model->propose(lambda); }
};

unsigned int TModelVectorForEstimation::acceptOrReject() {
	unsigned int numAccepted = 0;
	for (auto &model : _models) { numAccepted += (unsigned int)model->acceptOrReject(); }
	return numAccepted;
};

void TModelVectorForEstimation::adjust() {
	for (auto &model : _models) { model->adjust(); }
};

double TModelVectorForEstimation::maxF() const {
	double maxF = 0.0;
	for (auto &model : _models) {
		maxF = std::max(maxF, model->maxF());
	}
	return maxF;
};

void TModelVectorForEstimation::writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const {
	TModelNoRecal _noRecal;
	// open file and write header
	coretools::TOutputFile out(Filename, {"readGroup", "mate", "covariates", "rho"});

	// add models
	for (uint16_t r = 0; r < ReadGroups.size(); ++r) {
		for (uint8_t mate = 0; mate < 2; ++mate) {
			out.write(ReadGroups.getName(r), std::array{"first", "second"}[mate]);
			if (_modelIndex[r][mate])
				out.writeln(_modelIndex[r][mate]->epsilonDefinition(), _modelIndex[r][mate]->rhoDefinition());
			else 
				out.writeln(_noRecal.epsilonDefinition(), _noRecal.rhoDefinition());
		}
	}
}

std::string TModelVectorForEstimation::getModelsDefinition() {
	if (_models.empty()) return std::string{};

	return std::accumulate(_models.begin() + 1, _models.end(), _models.front()->epsilonDefinition(),
						   [](auto tot, auto &model) { return tot + '\n' + model->epsilonDefinition(); });
};

std::string TModelVectorForEstimation::getRhoDefinition() {
	if (_models.empty()) return std::string{};

	return std::accumulate(_models.begin() + 1, _models.end(), _models.front()->rhoDefinition(),
						   [](auto tot, auto &model) { return tot + '\n' + model->rhoDefinition(); });
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
	const auto dist = parameters().getParameterWithDefault("genoDist", "haploid");
	if (dist == "haploid") {
		_genoDist = std::make_unique<THaploidDistribution>();
	} else if (dist == "diploid") {
		_genoDist = std::make_unique<TDiploidDistribution>();
	} else {
		UERROR("Genotype distribution ", dist, " does not exist. Use 'haploid' or 'diploid'!");
	}
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
		DEVERROR("Estimation of genotype distribution not yet implemented!");
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

	logfile().conclude("Number of sites with data: ", _sites.size());
	size_t numSitesDepthTwoOrMore = _numSitesDepthTwoOrMore();
	logfile().conclude("Number of sites with depth > 1: ", numSitesDepthTwoOrMore);
	logfile().conclude("Number of bases: ", _dataTables.size());
	if (numSitesDepthTwoOrMore < 100) UERROR("Less than 100 sites with depth >= 2 available - aborting estimation!");

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
												  const PMD::TModels &PmdModels) {
	// initialize models
	_initializeModels(SequencingErrorModels);

	// run EM
	_runEM(outputName, PmdModels);

	// writing final estimates
	const std::string filename = outputName + "_recal.txt";
	logfile().list("Writing final estimates to file '", filename, "'.");
	_modelsToEstimate.writeRecalFile(*_readGroups, filename);
	BAM::RGInfo::TReadGroupInfo r(*_readGroups);
	SequencingErrorModels.addToRGInfo(r);
	r.write(outputName + "_recal.json");
};

void TRecalibrationEMEstimator::_estimateRho_updatePbbar(const PMD::TModels &PmdModels) {
	using genometools::genotype;
	coretools::TStrongArray<size_t, Base> cnt{0};
	_P_bbar_I_gds.clear();
	for (size_t i = 0; i < _sites.size(); ++i) {
		for (const auto &d_ij : _sites[i]) {
			++cnt[d_ij.base];
			_P_bbar_I_gds.emplace_back(0.);
			auto &P_bbar     = _P_bbar_I_gds.back();
			const auto L_eps = _modelsToEstimate.getBaseLikelihoods(d_ij);
			for (auto a = Base::min; a < Base::max; ++a) {
				const auto g_aa = genotype(a, a);
				const auto P_aa = PmdModels.massFunction(a, d_ij, L_eps);
				P_bbar[g_aa]    = P_aa[d_ij.base];

				_modelsToEstimate.addToRho(d_ij, _P_g_I_ds[i][g_aa], P_aa);
				if (!_genoDist->isInvariant()) {
					for (auto b = coretools::next(a); b < Base::max; ++b) {
						const auto g_ab = genotype(a, b);
						const auto P_ab = PmdModels.massFunction(g_ab, d_ij, L_eps);
						P_bbar[g_ab]    = P_ab[d_ij.base];

						_modelsToEstimate.addToRho(d_ij, _P_g_I_ds[i][g_ab], P_ab);
					}
				}
			}
		}
	}
	_modelsToEstimate.estimateRho();
}

void TRecalibrationEMEstimator::_updateEpsilon(double deltaLL_LL) {
	using coretools::str::toString;
	logfile().list("optimizing Q_beta using a Newton-Raphson algorithm.");
	const auto nTot = _modelsToEstimate.size();

	for (int i = 0; i < _NewtonRaphsonNumIterations; ++i) {
		logfile().startIndent("Running Newton-Raphson iteration ", i + 1, ":");
		_solveDerivative();
		const double curQ = _modelsToEstimate.Q();
		logfile().list("Current Q_beta = ", curQ);

		double lambda   = 1.0;
		size_t nUpdated = 0;
		double deltaQ   = 0;

		while (nUpdated < nTot && lambda > 1.0E-20) {
			logfile().startIndent("Backtracing with lambda = ", lambda, ":");
			_modelsToEstimate.propose(lambda);
			_calculateQ();

			deltaQ   = _modelsToEstimate.Q() - curQ;
			nUpdated = _modelsToEstimate.acceptOrReject();

			logfile().list("Delta Q = ", deltaQ);
			logfile().list(toString(nUpdated), "/", nTot, " models converged.");
			logfile().endIndent();

			// backtrack
			lambda = lambda / 2.0; // backtrack;
		}

		_modelsToEstimate.adjust();

		if (nUpdated < nTot) {
			logfile().conclude("Some models did not improve even with log2(lambda) = ", std::log2(lambda),
							   ", aborting Newton-Raphson.");
			break;
		}

		const double maxF  = _modelsToEstimate.maxF();
		if (maxF < _NewtonRaphsonMaxF) {
			logfile().conclude("max(F) = ", maxF, " < ", _NewtonRaphsonMaxF, ", ending Newton-Raphson.");
			logfile().endIndent();
			break;
		} 
		logfile().conclude("max(F) = ", maxF);

		if (const auto pdQ = std::abs(deltaQ/curQ); pdQ < deltaLL_LL) {
			logfile().conclude("deltaQ/Q = ", pdQ, " < deltaLL/LL = ", deltaLL_LL, ", ending Newton-Raphson.");
			logfile().endIndent();
			break;
		} 
		logfile().endIndent();
	}
	logfile().endIndent();
};

double TRecalibrationEMEstimator::_calculateLL_updatePg(const PMD::TModels &PmdModels) {
	_P_g_I_ds.clear();

	double LL = 0.0;
	for (auto &s_i : _sites) {
		if (s_i.genotype == Genotype::NN) { // unknown genotype
			_P_g_I_ds.emplace_back(1.); // Start at 1,1,1,1,1,1,1,1
			auto &P_g = _P_g_I_ds.back();
			for (auto &d_ij : s_i) {
				const auto P_eps = _modelsToEstimate.getBaseLikelihoods(d_ij);
				const auto P_D   = PmdModels.baseLikelihoods(d_ij, P_eps);
				P_g *= _genoDist->getGenotypeLikelihoods(P_D);
			}
			LL += log(_genoDist->normalize_add(P_g));
		} else { // known genotype.
			_P_g_I_ds.emplace_back(0.); 
			_P_g_I_ds.back()[s_i.genotype] = 1; // Probability of correct genotype is 1
			double P_g = 1.;
			for (auto &d_ij : s_i) {
				const auto L_eps = _modelsToEstimate.getBaseLikelihoods(d_ij);
				const auto L_D   = PmdModels.baseLikelihoods(d_ij, L_eps);
				P_g *= _genoDist->getGenotypeLikelihood(L_D, s_i.genotype);
			}
			LL += log(P_g);
		}
	}
	return LL;
};

void TRecalibrationEMEstimator::_runEM(const std::string &outputName, const PMD::TModels &PmdModels) {
	using coretools::str::toString;
	// run EM
	logfile().startNumbering("Running EM algorithm:");
	logfile().conclude("Initial pi: ", _genoDist->definition());
	logfile().conclude("Initial rho: ", _modelsToEstimate.getRhoDefinition());
	logfile().conclude("Initial model: ", _modelsToEstimate.getModelsDefinition());

	// calculate initial LL
	double oldLL   = _calculateLL_updatePg(PmdModels);
	double deltaLL = abs(oldLL);
	logfile().conclude("Initial log Likelihood = ", oldLL);

	// running iterations
	for (int i = 0; i < _numEMIterations; ++i) {
		logfile().number("EM Iteration:");
		logfile().addIndent();

		logfile().list("Updating pi");
		_genoDist->estimate();

		logfile().list("Updating rho");
		_estimateRho_updatePbbar(PmdModels);

		logfile().startIndent("Updating epsilon");
		_updateEpsilon(std::abs(deltaLL / oldLL));

		logfile().conclude("Current pi: ", _genoDist->definition());
		logfile().conclude("Current rho: ", _modelsToEstimate.getRhoDefinition());
		logfile().conclude("Current epsilon: ", _modelsToEstimate.getModelsDefinition());

		// calculate LL
		const double LL = _calculateLL_updatePg(PmdModels);
		logfile().conclude("Current Log Likelihood = ", LL);

		// check if we break based on LL
		deltaLL = LL - oldLL;
		logfile().conclude("delta LL = ", deltaLL);
		if (i > 0 && deltaLL < _minDeltaLL) {
			logfile().conclude("EM has converged (delta LL < ", _minDeltaLL, ")");
			break;
		}
		oldLL = LL;

		// write current estimates to file
		if (_writeTmpTables) {
			std::string filename = outputName + "_recalibrationEM_Loop" + toString(i) + ".txt";
			logfile().listFlush("Writing current estimates to file '", filename, "' ...");
			_modelsToEstimate.writeRecalFile(*_readGroups, filename);
			logfile().done();
		}

		// end loop
		logfile().endIndent();
		if (i == _numEMIterations - 1) logfile().warning("EM has not converged after maximum number of iterations!");
	}

	// finalize
	logfile().endNumbering();
};

void TRecalibrationEMEstimator::calcLL(TModels &SequencingErrorModels, const PMD::TModels &PmdModels) {
	_initializeModels(SequencingErrorModels);

	logfile().startIndent("Recal Model:");
	logfile().conclude("Rho: ",_modelsToEstimate.getRhoDefinition());
	logfile().conclude("Epsilon: ",_modelsToEstimate.getModelsDefinition());
	logfile().endIndent();

	logfile().startIndent("Calculating log likelihood:");
	const double LL = _calculateLL_updatePg(PmdModels);
	logfile().conclude("Log Likelihood = ", LL);

}

}; // namespace SequencingError

}; // end namespace GenotypeLikelihoods
