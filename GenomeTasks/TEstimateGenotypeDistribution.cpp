#include "TEstimateGenotypeDistribution.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/devtools.h"
#include <memory>


namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;
using GenotypeLikelihoods::TGenotypeLikelihoods;

double TEstimateGenotypeDistribution::_LL() {
	coretools::TSumLogProbability LL{};
	for (auto &site : _sites) {
		const auto ref = site.refBase;
		TGenotypeLikelihoods P_g_I_di(1.);
		size_t counter = 0;
		for (auto &d_ij : site) {
			const auto P_dij_I_bbar = _genotypeLikelihoodCalculator.sequencingErrorModels().P_dij(d_ij);
			const auto P_dij_I_b    = _genotypeLikelihoodCalculator.postMortemDamageModels().P_dij(d_ij, P_dij_I_bbar);
			const auto P_dij_I_g    = _genoDist->P_dij(P_dij_I_b);
			P_g_I_di *= P_dij_I_g;
			if (++counter > 10) {
				LL.add(coretools::normalize(P_g_I_di));
				counter = 0;
			}
		}
		LL.add(_genoDist->normalize_add(P_g_I_di, ref));
	}
	if (!std::isfinite(LL.getSum())) UERROR("LL = ", LL.getSum(), ", you may need to pool your readgroups!");
	return LL.getSum();
}

void TEstimateGenotypeDistribution::_runEM() {
	using coretools::str::toString;
	// run EM
	logfile().startNumbering("Running EM algorithm:");
	logfile().startIndent("Initial model:");
	_genoDist->log();
	logfile().endIndent();

	// calculate initial LL
	double oldLL   = _LL();
	double deltaLL = abs(oldLL);
	logfile().conclude("Initial log Likelihood = ", oldLL);

	// running iterations
	for (size_t i = 0; i < _numEMIterations; ++i) {
		logfile().number("EM Iteration:");
		logfile().addIndent();

		logfile().list("Updating pi");
		_genoDist->estimate();

		const double LL = _LL();
		deltaLL         = LL - oldLL;
		logfile().startIndent("Current model:");
		_genoDist->log();
		logfile().endIndent();

		logfile().conclude("Current Log Likelihood = ", LL);
		logfile().conclude("delta LL = ", deltaLL);

		// check if we break based on LL
		if (i > 0 && deltaLL < _minDeltaLL) {
			if (deltaLL < 0) logfile().warning("Negative LL!");
			else logfile().conclude("EM has converged (delta LL < ", _minDeltaLL, ")");
			break;
		}
		oldLL = LL;

		// end loop
		logfile().endIndent();
		if (i == _numEMIterations - 1) logfile().warning("EM has not converged after maximum number of iterations!");
	}

	// finalize
	logfile().endNumbering();
}

void TEstimateGenotypeDistribution::_handleWindow() {
	_sites.clear();
	for (const auto &s : _window) {
		if (s.empty() || s.refBase == genometools::Base::N) continue;
		_sites.emplace_back(s);
	}

	_runEM();
	
}

void TEstimateGenotypeDistribution::run() {
	_traverseBAMWindows();
}


TEstimateGenotypeDistribution::TEstimateGenotypeDistribution() {
	_openReference(true);
	_numEMIterations = parameters().get<int>("iterations", 200);
	_minDeltaLL      = parameters().get<double>("minDeltaLL", 1e-6);
	if (parameters().exists("HKY85")) {
		_genoDist = std::make_unique<GenotypeLikelihoods::THKY85>();
		logfile().list("Using HKY85 Genotype Distribution. (parameter 'HKY85')");
	} else {
		_genoDist = std::make_unique<GenotypeLikelihoods::TDiploidDistribution>();
		logfile().list("Estimating 10 Genotypes. (use 'HKY85' for HKY85 model)");
	}
}
}
