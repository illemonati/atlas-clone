/*
 * TRecalibrationEM.cpp
 *
 */
#include "TErrorEstimator.h"

#include <algorithm>
#include <armadillo>
#include <exception>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <stdexcept>

#include "PMD/TModels.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Strings/toString.h"
#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "oldPMD/TModels.h"
#include "TSequencedBase.h"
#include "SequencingError/TModel.h"
#include "SequencingError/TModels.h"
#include "coretools/algorithms.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenotypeLikelihoods {

using coretools::instances::logfile;
using coretools::instances::parameters;
using genometools::Genotype;
using genometools::Base;

namespace impl {
BAM::TReadGroupMap makeRGMap(const BAM::TReadGroups &ReadGroups) {
	if (parameters().parameterExists("poolReadGroups")) {
		std::string poolReadGroupsFile = parameters().getParameter<std::string>("poolReadGroups");
		logfile().startIndent("Will pool read groups (parameter 'poolReadGroups'):");
		return {ReadGroups, poolReadGroupsFile};
		logfile().endIndent();
	} else {
		logfile().list("Will estimate recalibration parameters for each read group. (use 'poolReadGroups' to pool)");
		return {ReadGroups};
	}
}
} // namespace impl


//---------------------------------------------------------------
// TErrorEstimator
//---------------------------------------------------------------
	TErrorEstimator::TErrorEstimator(const BAM::TReadGroups &ReadGroups) : _readGroups(&ReadGroups), _readGroupMap(impl::makeRGMap(ReadGroups)) {

	// genotype distribution: currently only allow for haploid
	const auto dist = parameters().getParameterWithDefault("genoDist", THKY85::name);
	if (dist == THaploidDistribution::name) {
		_genoDist = std::make_unique<THaploidDistribution>();
	} else if (dist == TDiploidDistribution::name) {
		_genoDist = std::make_unique<TDiploidDistribution>();
	} else if (dist == THKY85::name) {
		_genoDist = std::make_unique<THKY85>();
	} else {
		UERROR("Genotype distribution ", dist, " does not exist. Use '", THaploidDistribution::name, "', '", TDiploidDistribution::name, "' or '", THKY85::name, "'!");
	}
	logfile().list("Will use a ", _genoDist->typeString(), " genotype distribution.");

	_recal.initialize(_readGroupMap.size(), parameters().getParameterWithDefault("recalModel", "intercept;quality;position;context;fragmentLength;mappingQuality;"));
	_pmd.initialize(_readGroupMap.size(), parameters().getParameterWithDefault("pmdModel", "CT5;CT3"));

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

	logfile().endIndent();
};

//-----------------------------------------------
// Function to initialize models used for estimation
//-----------------------------------------------
size_t TErrorEstimator::_numSitesDepthTwoOrMore() {
	size_t _numSites = 0;
	for (auto &s : _sites) {
		if (s.depth() > 1) ++_numSites;
	}
	return _numSites;
};

void TErrorEstimator::_initializeModels() {
	using coretools::str::toString;
	using BAM::Mate;
	// count data available for recal
	logfile().listFlush("Counting data available for recal ...");
	// Note: data tables pool read groups!

	RecalEstimatorTools::TRecalDataTables dataTables(*_readGroups, _readGroupMap);
	dataTables.add(_sites);
	logfile().done();

	logfile().conclude("Number of sites with data: ", _sites.size());
	logfile().conclude("Number of sites with depth > 1: ", dataTables.nSites_g1());
	logfile().conclude("Number of bases: ", dataTables.size());
	if (dataTables.nSites_g1() < 100) UERROR("Less than 100 sites with depth >= 2 available - aborting estimation!");

	// identify models with data that can be estimated
	logfile().startIndent("Identifying models to estimate:");
	_recal.pool(_readGroupMap);
	_pmd.pool(_readGroupMap);
	for (auto rg : _readGroupMap.readGroupsInUse()) {
			for (Mate mate = Mate::min; mate < Mate::max; ++mate) {
			const auto& table = dataTables[rg][mate];
			if (table.size() > 0) {
				auto& recal = _recal.RGModel(rg)[mate];
				if (!recal->recalibrates()) UERROR("Cannot estimate readgroup ", rg, ", mate ", mate, "!");
				recal->epsilon()->init(table);
				_epsilons.push_back(recal->epsilon());
				_rhos.push_back(recal->rho());
			} else {
				_recal.reset(rg, mate);
			}
		}
	}
	logfile().endIndent();
};


void TErrorEstimator::_estimateRho_updatePbbar() {
	using genometools::genotype;
	_P_bbarEdij_I_gdijs.clear();
	for (size_t i = 0; i < _sites.size(); ++i) {
		for (const auto &d_ij : _sites[i]) {
			_P_bbarEdij_I_gdijs.emplace_back(0.);
			auto &P_bbarEdij_I_gdij = _P_bbarEdij_I_gdijs.back();
			const auto P_dij_I_bbar = _recal.P_dij(d_ij);
			for (auto a = Base::min; a < Base::max; ++a) {
				const auto aa              = genotype(a, a);
				const auto P_bbar_I_aa_dij = _pmd.P_bbar(a, d_ij, P_dij_I_bbar);
				P_bbarEdij_I_gdij[aa]      = P_bbar_I_aa_dij[d_ij.base];

				_recal.model(d_ij).rho()->add(d_ij.base, _P_g_I_dis[i][aa], P_bbar_I_aa_dij);
				if (!_genoDist->isInvariant()) {
					for (auto b = coretools::next(a); b < Base::max; ++b) {
						const auto ab              = genotype(a, b);
						const auto P_bbar_I_ab_dij = _pmd.P_bbar(ab, d_ij, P_dij_I_bbar);
						P_bbarEdij_I_gdij[ab]      = P_bbar_I_ab_dij[d_ij.base];

						_recal.model(d_ij).rho()->add(d_ij.base, _P_g_I_dis[i][ab], P_bbar_I_ab_dij);
					}
				}
			}
		}
	}
	for (auto& rho: _rhos) rho->estimate();
}

void TErrorEstimator::_updateEpsilon(double deltaLL_LL) {
	using coretools::str::toString;
	logfile().list("optimizing Q_beta using a Newton-Raphson algorithm.");
	const auto nTot = _epsilons.size();

	for (int i = 0; i < _NewtonRaphsonNumIterations; ++i) {
		logfile().startIndent("Running Newton-Raphson iteration ", i + 1, ":");
		_solveDerivative();
		double oldQ = 0.;
		for (const auto& e: _epsilons) oldQ += e->Q();
		logfile().list("Current Q_beta = ", oldQ);


		logfile().startIndent("Setting new epsilon.");
		for (auto& e: _epsilons) e->propose(1.);
		_calculateQ();
		size_t nUpdated = 0;
		for (auto& e: _epsilons) nUpdated += e->acceptOrReject();
		logfile().list(toString(nUpdated), "/", nTot, " models converged.");
		logfile().endIndent();

		double lambda   = 0.5;
		while(nUpdated < nTot) {
			logfile().startIndent("Backtracing with lambda = ", lambda, ":");
			for (auto &e : _epsilons) e->propose(lambda);
			_calculateQ();
			nUpdated = 0;
			for (auto &e : _epsilons) nUpdated += e->acceptOrReject();
			logfile().list(toString(nUpdated), "/", nTot, " models converged.");
			logfile().endIndent();

			lambda /= 2;
			if (lambda < 1e-20) break;
		}
		for (auto &e : _epsilons) e->adjust();

		if (nUpdated < nTot) {
			logfile().conclude(nTot - nUpdated, " models did not improve even with log2(lambda) = ", std::log2(lambda),
							   ", aborting Newton-Raphson.");
			break;
		}

		double maxF  = 0.;
		for (auto &e : _epsilons) maxF = std::max(e->maxF(), maxF);

		if (maxF < _NewtonRaphsonMaxF) {
			logfile().conclude("max(F) = ", maxF, " < ", _NewtonRaphsonMaxF, ", ending Newton-Raphson.");
			logfile().endIndent();
			break;
		} 
		logfile().conclude("max(F) = ", maxF);

		double newQ = 0.;
		for (const auto& e: _epsilons) newQ += e->Q();
		if (const auto pdQ = std::abs(newQ/oldQ); pdQ < deltaLL_LL) {
			logfile().conclude("deltaQ/Q = ", pdQ, " < deltaLL/LL = ", deltaLL_LL, ", ending Newton-Raphson.");
			logfile().endIndent();
			break;
		} 
		logfile().endIndent();
	}
	logfile().endIndent();
};

double TErrorEstimator::_calculateLL_updatePg() {
	_P_g_I_dis.clear();

	double LL = 0.0;
	for (auto &s_i : _sites) {
		if (s_i.genotype == Genotype::NN) { // unknown genotype
			const auto ref = s_i.refBase;
			_P_g_I_dis.emplace_back(1.); // Start at 1,1,1,1,1,1,1,1
			auto &P_g_I_di = _P_g_I_dis.back();
			for (auto &d_ij : s_i) {
				const auto P_dij_I_bbar = _recal.P_dij(d_ij);
				const auto P_dij_I_b    = _pmd.P_dij(d_ij, P_dij_I_bbar);
				const auto P_dij_I_g    = _genoDist->P_dij(P_dij_I_b);
				P_g_I_di *= P_dij_I_g;
			}
			LL += log(_genoDist->normalize_add(P_g_I_di, ref));
		} else { // known genotype.
			_P_g_I_dis.emplace_back(0.); 
			_P_g_I_dis.back()[s_i.genotype] = 1; // Probability of correct genotype is 1
			double P_g = 1.;
			for (auto &d_ij : s_i) {
				const auto L_eps = _recal.P_dij(d_ij);
				const auto L_D   = _pmd.P_dij(d_ij, L_eps);
				P_g *= _genoDist->getGenotypeLikelihood(L_D, s_i.genotype);
			}
			LL += log(P_g);
		}
	}
	return LL;
};

void TErrorEstimator::_runEM() {
	using coretools::str::toString;
	// run EM
	logfile().startNumbering("Running EM algorithm:");
	logfile().conclude("Initial pi: ", _genoDist->definition());
	logfile().startIndent("Initial rho");
	for (const auto& r: _rhos) logfile().list(r->definition());
	logfile().endIndent();
	logfile().startIndent("Initial epsilon");
	for (const auto& e: _epsilons) logfile().list(e->definition());
	logfile().endIndent();

	// calculate initial LL
	double oldLL   = _calculateLL_updatePg();
	double deltaLL = abs(oldLL);
	logfile().conclude("Initial log Likelihood = ", oldLL);

	// running iterations
	for (int i = 0; i < _numEMIterations; ++i) {
		logfile().number("EM Iteration:");
		logfile().addIndent();

		logfile().list("Updating pi");
		_genoDist->estimate();

		logfile().list("Updating rho");
		_estimateRho_updatePbbar();

		logfile().startIndent("Updating epsilon");
		_updateEpsilon(std::abs(deltaLL / oldLL));

		logfile().conclude("Current pi: ", _genoDist->definition());
		logfile().startIndent("Current rho");
		for (const auto &r : _rhos) logfile().list(r->definition());
		logfile().endIndent();
		logfile().startIndent("Current epsilon");
		for (const auto &e : _epsilons) logfile().list(e->definition());
		logfile().endIndent();

		// calculate LL
		const double LL = _calculateLL_updatePg();
		logfile().conclude("Current Log Likelihood = ", LL);

		// check if we break based on LL
		deltaLL = LL - oldLL;
		logfile().conclude("delta LL = ", deltaLL);
		if (i > 0 && deltaLL < _minDeltaLL) {
			logfile().conclude("EM has converged (delta LL < ", _minDeltaLL, ")");
			break;
		}
		oldLL = LL;

		// end loop
		logfile().endIndent();
		if (i == _numEMIterations - 1) logfile().warning("EM has not converged after maximum number of iterations!");
	}

	// finalize
	logfile().endNumbering();
};

void TErrorEstimator::estimate(std::string_view outputName) {
	// initialize models
	_initializeModels();

	// run EM
	_runEM();

	// writing final estimates
	BAM::RGInfo::TReadGroupInfo r(*_readGroups);
	_recal.addToRGInfo(r);
	_pmd.addToRGInfo(r);
	r.write(std::string(outputName) + ".json");
};

void TErrorEstimator::calcLL() {
	_initializeModels();

	logfile().conclude("pi: ", _genoDist->definition());
	logfile().startIndent("rho");
	for (const auto &r : _rhos) logfile().list(r->definition());
	logfile().endIndent();
	logfile().startIndent("epsilon");
	for (const auto &e : _epsilons) logfile().list(e->definition());
	logfile().endIndent();

	logfile().startIndent("Calculating log likelihood:");
	const double LL = _calculateLL_updatePg();
	logfile().conclude("Log Likelihood = ", LL);

}

}; // end namespace GenotypeLikelihoods
