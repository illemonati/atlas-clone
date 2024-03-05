#include "TEstimateGenotypeDistribution.h"

#include "TGenotypeDistribution.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Strings/toString.h"
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
	const auto isInvariant = _genoDist->isInvariant();
	coretools::TSumLogProbability LL{};
	for (const auto &site : Sites) {
		const auto ref = site.refBase;
		if (site.empty() || ref == genometools::Base::N) continue;
		TGenotypeLikelihoods P_g_I_di(P(1.));
		for (auto &d_ij : site) {
			const auto P_dij_I_bbar = _genome.errorModels().sequencingErrorModels().P_dij(d_ij);
			const auto P_dij_I_b    = _genome.errorModels().postMortemDamageModels().P_dij(d_ij, P_dij_I_bbar);

			if (isInvariant) {
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
	_genoDist->reset();
	logfile().startIndent("Estimating ", _genoDist->typeString(), ":");
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
	_genoDist->write(_out);
	logfile().endIndent();
	return oldLL;
}

void TEstimateGenotypeDistribution::_handleWindow(GenotypeLikelihoods::TWindow& window) {
	if (_genomeWide) {
		_totSites       += window.size();
		_totMaskedSites += window.numMaskedSites();

		// full P
		auto &sites = _sites_P.front();
		auto &stat  = _stats.front();
		for (const auto &site : window) {
			stat.NData += site.depth();
			if (site.empty() || site.refBase == Base::N) {
				++stat.NMissing;
			} else {
				sites.push_back(site);
			}
		}

		// downsample
		for (size_t i = 1; i < _sites_P.size(); ++i) {
			const auto p = _probs[i - 1];
			auto &sites  = _sites_P[i];
			auto &stat   = _stats[i];

			logfile().list("Downsampling reads to probability ", p, ".");
			GenotypeLikelihoods::TWindow downsampled(window, _windows.uptoDepth(), p);
			_windows.filter(downsampled);

			for (const auto &site : downsampled) {
				stat.NData += site.depth();
				if (site.empty() || site.refBase == Base::N) {
					++stat.NMissing;
				} else {
					sites.push_back(site);
				}
			}
		}
	} else {
		// full P
		_out.write(window.chrName(), window.from().position(), window.to().position(), window.depth(),
		           window.numSites(), window.numSitesWithData(), window.fracMissing());

		logfile().list("Using full data.");
		const auto LL = _runEM(window.sites());
		_out.write(LL);

		// downsample
		for (const auto p: _probs) {
			logfile().list("Downsampling reads to probability ", p, ".");
			GenotypeLikelihoods::TWindow downsampled(window, _windows.uptoDepth(), p);
			_windows.filter(downsampled);

			_out.write(downsampled.depth(), downsampled.numSites(), downsampled.numSitesWithData(), downsampled.fracMissing());
			const auto LL = _runEM(downsampled.sites());
			_out.write(LL);
		}

		_out.endln();
	}
}

void TEstimateGenotypeDistribution::run() {
	_traverseBAMWindows();
	if (_genomeWide) {
		if (_windows.considerRegions()) {
			_out.write("regions", "-", "-");
		} else {
			_out.write("genome-wide", "-", "-");
		}
		const auto NSites = _totSites - _totMaskedSites;
		for (size_t i = 0; i < _sites_P.size(); ++i) {
			const auto& sites = _sites_P[i];
			const auto& stat  = _stats[i];
			_out.write(stat.NData / NSites, NSites, _totSites - stat.NMissing, double(stat.NMissing - _totMaskedSites) / NSites);

			const auto LL = _runEM(sites);
			_out.write(LL);
		}
		_out.endln();
	}
}

TEstimateGenotypeDistribution::TEstimateGenotypeDistribution() {
	using coretools::str::toString;
	_windows.requireReference();

	_numEMIterations  = parameters().get<int>("iterations", 200);
	_minDeltaLL       = parameters().get<double>("minDeltaLL", 1e-6);
	const auto ploidy = parameters().get("ploidy", 2);

	if (ploidy == 2) {
		_genoDist = std::make_unique<GenotypeLikelihoods::THKY85>();
	} else if (ploidy == 1) {
		_genoDist = std::make_unique<GenotypeLikelihoods::THKY85_mono>();
	} else {
		UERROR("Cannot estimate a HK85 model with a ploidy of ", ploidy, "!");
	} 
	logfile().list("Estimating HK85 assuming a ploidy of ", ploidy, ". (parameter 'ploidy')");

	// Downsample?
	if (parameters().exists("prob")) {
		parameters().fill("prob", _probs);
		if (_probs.front() == 1.) {
			// we will do the full probability anyway
			_probs.erase(_probs.begin());
		}
		logfile().list("Will do downsampling experiment using all data, and downsampling with probabilities ", _probs, ". (parameter 'prob')");
	} else {
		logfile().list("Will not do downsampling. (use parameter 'prob')");
	}

	// genomewide?
	_genomeWide = parameters().exists("genomeWide");
	if (_genomeWide) {
		logfile().list("Will estimating genotype Distribution genome-wide. (parameter 'genomeWide')");
		_sites_P.resize(_probs.size() + 1); // full P and downsample
		_stats.resize(_probs.size() + 1);
	} else {
		logfile().list("Will estimating genotype Distribution per window. (use 'genomeWide' for genome-wide estimation)");
	}

	// Output file
	std::vector<std::string> header{"Chr", "Start", "End"};
	const auto sp = _probs.empty() ? "": "_p1.0";
	header.insert(header.end(), {toString("depth", sp), toString("numSites", sp), toString("numSitesData", sp), toString("fracMissing", sp)});
	_genoDist->addHeader(header); //, sp);
	header.push_back(toString("LL", sp));

	for (const auto p : _probs) {
		const auto sp = toString("_p", p);
		header.insert(header.end(), {toString("depth", sp), toString("numSites", sp), toString("numSitesData", sp),
									 toString("fracMissing", sp)});
		_genoDist->addHeader(header); //, sp);
		header.push_back(toString("LL", sp));
	}

	_out.open(_genome.outputName() + ".txt.gz");
	_out.writeHeader(header);
}
}
