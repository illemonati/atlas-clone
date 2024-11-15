#include "TEstimateHKY85.h"

#include "PMD/TPsi.h"
#include "TGenotypeDistribution.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include <memory>

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::P;
using genometools::TGenotypeLikelihoods;
using genometools::Base;
using genometools::Genotype;

double TEstimateHKY85::_LL(const std::vector<GenotypeLikelihoods::TSite> &Sites) {
	const auto isInvariant = _genoDist->isInvariant();

	const auto PgI_init = [isInvariant]() {
		TGenotypeLikelihoods Ps(P(1.));
		if (isInvariant) {
			Ps[Genotype::AC] = P(0.);
			Ps[Genotype::AG] = P(0.);
			Ps[Genotype::AT] = P(0.);
			Ps[Genotype::CG] = P(0.);
			Ps[Genotype::CT] = P(0.);
			Ps[Genotype::GT] = P(0.);
		}
		return Ps;
	}();

	coretools::TSumLogProbability LL{};
	for (const auto &site : Sites) {
		const auto ref = site.refBase;
		if (site.empty() || ref == genometools::Base::N) continue;
		TGenotypeLikelihoods P_g_I_di = PgI_init;
		double sum                    = 1.;

		for (auto &d_ij : site) {
			const auto P_dij_I_bbar = _genome.errorModels().sequencingErrorModels().P_dij(d_ij);
			const auto P_dij_I_b    = _pmd ? _pmd->P_dij(d_ij, P_dij_I_bbar) : _genome.errorModels().postMortemDamageModels().P_dij(d_ij, P_dij_I_bbar);

			LL.add(sum);
			const double sum_inv = 1. / sum;
			sum                  = 0.;
			for (auto k = Base::min; k < Base::max; ++k) {
				const auto kk = genometools::genotype(k, k);
				P_g_I_di[kk] *= P(P_dij_I_b[k] * sum_inv);
				sum          += P_g_I_di[kk];
			}
			if (!isInvariant) {
				for (const auto kl :
					 {Genotype::AC, Genotype::AG, Genotype::AT, Genotype::CG, Genotype::CT, Genotype::GT}) {
					const auto k  = genometools::first(kl);
					const auto l  = genometools::second(kl);
					P_g_I_di[kl] *= P(0.5 * (P_dij_I_b[k] + P_dij_I_b[l]) * sum_inv);
					sum          += P_g_I_di[kl];
				}
			}
		}
		LL.add(_genoDist->normalize_add(P_g_I_di, ref));
		if (_pmd) {
			for (auto &d_ij : site) {
				const auto P_dij_I_bbar = _genome.errorModels().sequencingErrorModels().P_dij(d_ij);
				_pmd->psi()->addCT(d_ij, P_g_I_di[Genotype::CC], _pmd->P_b_bbar(Genotype::CC, d_ij, P_dij_I_bbar));
				_pmd->psi()->addGA(d_ij, P_g_I_di[Genotype::GG], _pmd->P_b_bbar(Genotype::GG, d_ij, P_dij_I_bbar));

				if (!isInvariant)
					for (auto g : {Genotype::AC, Genotype::CG, Genotype::CT}) {
						_pmd->psi()->addCT(d_ij, P_g_I_di[g], _pmd->P_b_bbar(g, d_ij, P_dij_I_bbar));
					}
				if (!isInvariant)
					for (auto g : {Genotype::AG, Genotype::CG, Genotype::GT}) {
						_pmd->psi()->addGA(d_ij, P_g_I_di[g], _pmd->P_b_bbar(g, d_ij, P_dij_I_bbar));
					}
			}
		}
	}
	if (!std::isfinite(LL.getSum())) UERROR("LL = ", LL.getSum(), "!");
	return LL.getSum();
}

double TEstimateHKY85::_runEM(const std::vector<GenotypeLikelihoods::TSite>& Sites) {
	using coretools::str::toString;
	// run EM
	_genoDist->reset();
	if (_pmd) _pmd->psi()->reset(_genome.rgInfo()[0][BAM::RGInfo::InfoType::pmd]);

	const auto s = _pmd ? " and PMD:" : ":";
	logfile().startIndent("Estimating ", _genoDist->typeString(), s);
	logfile().startIndent("Initial values:");
	_genoDist->log();
	if (_pmd) _pmd->psi()->log();

	// calculate initial LL
	double oldLL   = _LL(Sites);
	double deltaLL = std::abs(oldLL);
	logfile().list("log Likelihood = ", oldLL);
	logfile().endIndent();

	// running iterations
	logfile().startNumbering("Running EM algorithm:");
	for (size_t i = 0; i < _numEMIterations; ++i) {
		logfile().addIndent();
		logfile().number("EM Iteration:");
		_genoDist->estimate();
		if (_pmd) _pmd->psi()->estimate();

		const double LL = _LL(Sites);
		deltaLL         = LL - oldLL;

		_genoDist->log();
		if (_pmd) _pmd->psi()->log();
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
	}

	if (deltaLL > _minDeltaLL) {
		logfile().warning("EM has not converged after maximum number of iterations!");
		oldLL = -std::numeric_limits<double>::infinity();
	}

	// finalize
	logfile().endNumbering();
	logfile().endIndent();
	return oldLL;
}

void TEstimateHKY85::_handleGenomeWide(GenotypeLikelihoods::TWindow &window) {
	_totSites += window.size();
	_totMaskedSites += window.numMaskedSites();

	// full P
	for (const auto &site : window) {
		_stats_full.NData += site.depth();
		if (site.empty() || site.refBase == Base::N) {
			++_stats_full.NMissing;
		} else {
			_sites_full.push_back(site);
		}
	}

	// downsample
	if (_sample == Sample::reads) {
		for (size_t r = 0; r < _nRounds; ++r) {
			for (size_t i = 0; i < _sites_P[r].size(); ++i) {
				const auto p = P(_depthOrProbs[i]);
				auto &sites  = _sites_P[r][i];
				auto &stat   = _stats_P[r][i];

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
		}
	}
}

void TEstimateHKY85::_handlePerWindow(GenotypeLikelihoods::TWindow &window) {
	using BAM::End;
	// full P

	logfile().list("Using full data.");
	const auto LL = _runEM(window.sites());

	_out.write(window.chrName(), window.from().position(), window.to().position(), window.depth(), window.numSites(),
			   window.numSitesWithData(), window.fracMissing(), _genoDist->pis(), LL);
	if (_pmd) _out.write(_pmd->psi()->vals(End::from5).front(), _pmd->psi()->vals(End::from3).front());

	// downsample
	const auto nIT = _numEMIterations;
	for (const auto dOrP : _depthOrProbs) {

		_numEMIterations = std::min<size_t>(10 * nIT, nIT / dOrP); // may need a bit longer

		if (_sample == Sample::reads) {
			const coretools::Probability p = P(dOrP);
			logfile().list("Downsampling reads to probability ", p, ".");
			GenotypeLikelihoods::TWindow downsampled(window, _windows.uptoDepth(), p);
			_windows.filter(downsampled);

			const auto LL_p  = _runEM(downsampled.sites());
			_out.write(downsampled.depth(), downsampled.numSites(), downsampled.numSitesWithData(),
			           downsampled.fracMissing(), _genoDist->pis(), LL_p);
		} else {
			auto sites      = window.sites();
			size_t depth    = 0;
			size_t withData = 0;
			for (auto &site : sites) {
				if (_sample == Sample::sites) {
					const coretools::Probability p = P(dOrP);
					logfile().list("Downsampling sites to probability ", p, ".");
					site.downsample(p);
				} else /* _sample == Sample::upTo */ {
					const int d = int(dOrP);
					logfile().list("Downsampling sites a maximum depth of", d, ".");
					site.downsample(d);
				}
				depth += site.depth();
				if (!site.empty()) ++withData;
			}
			const auto LL_p  = _runEM(sites);
			_out.write(double(depth)/window.numSites(), window.numSites(), withData, double(window.numSites() - withData) / window.numSites(),
					   _genoDist->pis(), LL_p);
		}
		if (_pmd) _out.write(_pmd->psi()->vals(End::from5).front(), _pmd->psi()->vals(End::from3).front());
	}
	_numEMIterations = nIT;

	_out.endln();
}

void TEstimateHKY85::_handleWindow(GenotypeLikelihoods::TWindow& window) {
	if (_genomeWide) {
		_handleGenomeWide(window);
	} else {
		_handlePerWindow(window);
	}
}

void TEstimateHKY85::run() {
	using BAM::End;
	_traverseBAMWindows();
	if (_genomeWide) {
		_openFile();
		const auto NSites = _totSites - _totMaskedSites;
		//full

		logfile().list("Using full data.");
		const auto LL  = _runEM(_sites_full);
		const auto pis = _genoDist->pis();

		const auto nIT = _numEMIterations;
		for (size_t r = 0; r < _nRounds; ++r) {
			logfile().list("Downsampling round ", r + 1, ".");
			if (_windows.considerRegions()) {
				_out.write("regions", "Round", r+1);
			} else {
				_out.write("genome-wide", "Round", r+1);
			}
			_out.write(_stats_full.NData / NSites, NSites, _totSites - _stats_full.NMissing,
					   double(_stats_full.NMissing - _totMaskedSites) / NSites, pis, LL);
			if (_pmd) _out.write(_pmd->psi()->vals(End::from5).front(), _pmd->psi()->vals(End::from3).front());

			// downsampled
			for (size_t i = 0; i < _sites_P[r].size(); ++i) {
				if (_sample == Sample::upToDepth) {
					_numEMIterations = std::min<size_t>(10*nIT, nIT*_depthOrProbs[0]/_depthOrProbs[i]); // may need a bit longer
				} else {
					_numEMIterations = std::min<size_t>(10*nIT, nIT/_depthOrProbs[i]); // may need a bit longer
				}

				if (_sample == Sample::reads) {
					const auto p = P(_depthOrProbs[i]);
					logfile().list("Downsampling reads to probability ", p, ".");
					const auto &sites = _sites_P[r][i];
					const auto &stat  = _stats_P[r][i];

					const auto LL_i = _runEM(sites);
					_out.write(stat.NData / NSites, NSites, _totSites - stat.NMissing,
					           double(stat.NMissing - _totMaskedSites) / NSites, _genoDist->pis(), LL_i);
				} else {
					auto sites      = _sites_full;
					size_t depth    = 0;
					size_t withData = 0;
					for (auto &site : sites) {
						if (_sample == Sample::sites) {
							const auto p = P(_depthOrProbs[i]);
							logfile().list("Downsampling sites to probability ", p, ".");
							site.downsample(p);
						}else /* _sample == Sample::upTo */ {
							const int d = int(_depthOrProbs[i]);
							logfile().list("Downsampling sites a maximum depth of", d, ".");
							site.downsample(d);
						}
						depth += site.depth();
						if (!site.empty()) ++withData;
					}
					const auto LL_p = _runEM(sites);
					_out.write(double(depth)/NSites, NSites, withData, double(NSites - withData) / NSites, _genoDist->pis(), LL_p);
				}
				if (_pmd) _out.write(_pmd->psi()->vals(End::from5).front(), _pmd->psi()->vals(End::from3).front());
			}
			_out.endln();
		}
	}
}

TEstimateHKY85::TEstimateHKY85() {
	using BAM::RGInfo::InfoType;
	_windows.requireReference();

	if (parameters().exists("doPMD")) {
		logfile().list("Will re-estimate pmd per Window. (parameter 'doPMD')");
		_pmd = std::make_unique<GenotypeLikelihoods::PMD::TWithPMD>();
		_pmd->psi()->log();
	} else {
		logfile().list("Not re-estimating PMD. (use 'doPMD')");
	}

	_numEMIterations  = parameters().get<int>("iterations", 200);
	_minDeltaLL       = parameters().get<double>("minDeltaLL", 1e-6);

	logfile().list("Will perform at max ", _numEMIterations, " EM iterations. (parameter 'iterations')");
	logfile().list("Will stop EM when deltaLL < ", _minDeltaLL, ". (parameter 'minDeltaLL')");

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
	const auto type = parameters().get("sample", "reads");
	if (parameters().exists("prob")) {
		parameters().fill("prob", _depthOrProbs);
		if (_depthOrProbs.front() == 1.) {
			// we will do the full probability anyway
			_depthOrProbs.erase(_depthOrProbs.begin());
		}
		if (type == "reads") {
			_sample = Sample::reads;
		} else if (type == "sites") {
			_sample = Sample::sites;
		} else if (type == "upToDepth") {
			UERROR("Cannot combine ", type, " with prob");
		} else {
			UERROR("Downsampling type ", type, " does not exist");
		}
		logfile().list("Will do downsampling experiment using all data, and downsampling ", type,
					   " with probabilities ", _depthOrProbs, ". (parameter 'prob')");
	} else if (type == "upToDepth") {
		const auto depths = parameters().get<std::vector<int>>("depth");
		for (const auto depth : depths) {
			if (depth < 1) UERROR("Must sample at least depth=1");
			_depthOrProbs.push_back(double(depth));
		}
		_sample = Sample::upToDepth;
	} else {
		logfile().list("Will not do downsampling. (use parameter 'prob')");
	}

	// genomewide?
	_genomeWide = parameters().exists("genomeWide");
	if (_genomeWide) {
		logfile().list("Will estimating genotype Distribution genome-wide. (parameter 'genomeWide')");
		if (parameters().get("genomeWide").empty()) {
			_nRounds = 1;
			if (!_depthOrProbs.empty()) {
				logfile().list("Will do one round of downsampling. (use 'genomeWide' to specify the number)");
			}
		} else {
			if (_depthOrProbs.empty()) {
				UERROR("Cannot du several rounds of downsampling without downsampling probabilities! Use 'prob' to specify.");
			}
			parameters().fill("genomeWide", _nRounds);
			logfile().list("Will do ", _nRounds, " rounds of downsampling. (parameter 'genomeWide')");
		}
		for (size_t r = 0; r < _nRounds; ++r) {
			_sites_P.emplace_back();
			_stats_P.emplace_back();
			_sites_P.back().resize(_depthOrProbs.size());
			_stats_P.back().resize(_depthOrProbs.size());
		}
	} else {
		logfile().list("Will estimating genotype Distribution per window. (use 'genomeWide' for genome-wide estimation)");
		_openFile(); // will be written every Window
	}
}
void TEstimateHKY85::_openFile() {
	using coretools::str::toString;
	std::vector<std::string> header{"Chr", "Start", "End"};
	const auto sp = _depthOrProbs.empty() || _sample == Sample::upToDepth ? "": "p1.0_";
	header.insert(header.end(), {toString(sp, "depth"), toString(sp, "numSites"), toString(sp, "numSitesData"), toString(sp, "fracMissing")});
	_genoDist->addHeader(header, sp); //, sp);
	header.push_back(toString(sp, "LL"));
	if (_pmd) header.insert(header.end(), {"PMD5", "PMD3"});

	for (const auto p : _depthOrProbs) {
		const auto sp = toString("p", p, "_");
		header.insert(header.end(), {toString(sp, "depth"), toString(sp, "numSites"), toString(sp, "numSitesData"),
									 toString(sp, "fracMissing")});
		_genoDist->addHeader(header, sp); //, sp);
		header.push_back(toString(sp, "LL"));
		if (_pmd) header.insert(header.end(), {toString(sp, "PMD5"), toString(sp, "PMD3")});
	}

	_out.open(_genome.outputName() + ".txt.gz");
	_out.writeHeader(header);
}
}
