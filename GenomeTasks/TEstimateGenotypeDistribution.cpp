#include "TEstimateGenotypeDistribution.h"

#include "TGenotypeDistribution.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSumLog.h"
#include "genometools/GenotypeTypes.h"
#include <iterator>


namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::P;
using GenotypeLikelihoods::TGenotypeLikelihoods;
using genometools::Base;
using genometools::Genotype;

double TEstimateGenotypeDistribution::_LL(const std::vector<GenotypeLikelihoods::TSite> &Sites) {
	coretools::TSumLogProbability LL{};
	for (const auto &site : Sites) {
		const auto ref = site.refBase;
		if (site.empty() || ref == genometools::Base::N) continue;
		TGenotypeLikelihoods P_g_I_di(P(1.));
		for (auto &d_ij : site) {
			const auto P_dij_I_bbar = _genome.errorModels().sequencingErrorModels().P_dij(d_ij);
			const auto P_dij_I_b    = _genome.errorModels().postMortemDamageModels().P_dij(d_ij, P_dij_I_bbar);

			if (_genoDist->isInvariant()) {
				double max = 0.;
				for (auto b = Base::min; b < Base::max; ++b) {
					const auto g = genometools::genotype(b, b);
					P_g_I_di[g] *= P_dij_I_b[b];
					max = std::max<double>(max, P_g_I_di[g]);
				}
				if (max < 1e-2) {
					LL.add(max);
					for (auto b = Base::min; b < Base::max; ++b) {
						const auto g = genometools::genotype(b, b);
						P_g_I_di[g].scale(max);
					}
				}
			} else {
				const auto P_dij_I_g = GenotypeLikelihoods::base2genotype<genometools::Ploidy::diploid>(P_dij_I_b);
				double max           = 0.;
				for (auto g = Genotype::min; g < Genotype::max; ++g) {
					P_g_I_di[g] *= P_dij_I_g[g];
					max = std::max<double>(max, P_g_I_di[g]);
				}
				if (max < 1e-2) {
					LL.add(max);
					for (auto &p : P_g_I_di) p.scale(max);
				}
			}
		}
		LL.add(_genoDist->normalize_add(P_g_I_di, ref));
	}
	if (!std::isfinite(LL.getSum())) UERROR("LL = ", LL.getSum(), "!");
	return LL.getSum();
}

double TEstimateGenotypeDistribution::_runEM(const std::vector<GenotypeLikelihoods::TSite>& Sites) {
	using coretools::str::toString;
	// run EM
	logfile().startIndent("Initial values:");
	_genoDist->log();

	// calculate initial LL
	double oldLL   = _LL(Sites);
	double deltaLL = abs(oldLL);
	logfile().list("log Likelihood = ", oldLL);
	logfile().endIndent();

	// running iterations
	logfile().startNumbering("Running EM algorithm:");
	for (size_t i = 0; i < _numEMIterations; ++i) {
		logfile().addIndent();
		logfile().number("EM Iteration:");
		_genoDist->estimate();

		const double LL = _LL(Sites);
		deltaLL         = LL - oldLL;

		_genoDist->log();
		logfile().list("Log Likelihood = ", LL);
		logfile().list("delta LL = ", deltaLL);

		// check if we break based on LL
		logfile().endIndent();
		oldLL = LL;
		if (i > 0 && deltaLL < _minDeltaLL) {
			if (deltaLL < 0) logfile().warning("Negative LL!");
			else logfile().conclude("EM has converged (delta LL < ", _minDeltaLL, ")");
			break;
		}

		// end loop
		if (i == _numEMIterations - 1) logfile().warning("EM has not converged after maximum number of iterations!");
	}

	// finalize
	logfile().endNumbering();
	return oldLL;
}

void TEstimateGenotypeDistribution::_handleWindow(GenotypeLikelihoods::TWindow& window) {
	if (_genomeWide) {
		_totMaskedSites += window.numMaskedSites();
		_sites.reserve(_sites.size() + window.size());
		std::copy(window.begin(), window.end(), std::back_inserter(_sites));
	} else {
		_out.write(window.chrName(), window.from().position(), window.to().position(), window.depth(),
		           window.numSites(), window.numSitesWithData(), window.fracMissing());

		const auto LL = _runEM(window.sites());
		_genoDist->write(_out);
		_out.writeln(LL);
	}
}

void TEstimateGenotypeDistribution::run() {
	_traverseBAMWindows();
	if (_genomeWide) {
		if (_windows.considerRegions()) {
			_out.write("regions");
		} else {
			_out.write("genome-wide");
		}

		const auto NSites = _sites.size() - _totMaskedSites;
		size_t NData = 0;
		size_t NMissing = 0;
		for (const auto& s: _sites) {
			NData += s.depth();
			if (s.empty()) ++NMissing;
		}

		_out.write("-", "-", NData / NSites, NSites, _sites.size() - NMissing, double(NMissing - _totMaskedSites) / NSites);

		const auto LL = _runEM(_sites);
		_genoDist->write(_out);
		_out.writeln(LL);
	}
}


TEstimateGenotypeDistribution::TEstimateGenotypeDistribution() {
	_windows.requireReference();
	_numEMIterations = parameters().get<int>("iterations", 200);
	_minDeltaLL      = parameters().get<double>("minDeltaLL", 1e-6);
	const size_t ploidy = parameters().get("ploidy", 2);
	if (ploidy == 2) {
		_genoDist = std::make_unique<GenotypeLikelihoods::THKY85>();
	} else if (ploidy == 1){
		_genoDist = std::make_unique<GenotypeLikelihoods::THKY85_mono>();
	} else {
		UERROR("Cannot estimate a HK85 model with a ploidy of ", ploidy, "!");
	}
	logfile().list("Estimating HK85 assuming a ploidy of ", ploidy, ". (parameter 'ploidy')");

	_genomeWide = parameters().exists("genomeWide");
	if (_genomeWide) {
		logfile().list("Will estimating genotype Distribution genome-wide. (parameter 'genomeWide')");
	} else {
		logfile().list("Will estimating genotype Distribution per window. (use 'genomeWide' for genome-wide estimation)");
	}

	std::vector<std::string> header{"Chr", "Start", "End", "depth", "numSites", "numSitesData", "fracMissing"};
	_genoDist->addHeader(header);
	header.push_back("LL");
	_out.open(_genome.outputName() + ".txt.gz");
	_out.writeHeader(header);
}
}
