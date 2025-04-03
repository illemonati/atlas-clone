#include "TEstimateHKY85.h"

#include "PMD/TPsi.h"
#include "TGenotypeDistribution.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Ploidy.h"
#include <iterator>
#include <memory>

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::P;
using genometools::TGenotypeLikelihoods;
using genometools::Base;
using genometools::Genotype;

void TEstimateHKY85::_addSites(const std::vector<GenotypeLikelihoods::TSite>& Sites) {
	const auto isInvariant = _genoDist->ploidy() == genometools::Ploidy::haploid;
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

	_refBases.reserve(_refBases.size() + Sites.size());
	_P_g_I_ds.reserve(_P_g_I_ds.size() + Sites.size());
	for (const auto &site : Sites) {
		if (site.empty() || site.refBase == genometools::Base::N) continue;

		_refBases.push_back(site.refBase);
		TGenotypeLikelihoods P_g_I_di = PgI_init;
		double sum                    = 1.;

		for (auto &d_ij : site) {
			const auto P_dij_I_bbar = _genome.errorModels().sequencingErrorModels().P_dij(d_ij);
			const auto P_dij_I_b    = _genome.errorModels().postMortemDamageModels().P_dij(d_ij, P_dij_I_bbar);

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
		_P_g_I_ds.push_back(P_g_I_di);
	}
}

double TEstimateHKY85::_LL() {
	coretools::TSumLogProbability LL{};
	for (size_t i = 0; i < _refBases.size(); ++i) {
		LL.add(_genoDist->add(_P_g_I_ds[i], _refBases[i]));
	}
	if (!std::isfinite(LL.getSum())) UERROR("LL = ", LL.getSum(), "!");
	return LL.getSum();
}

double TEstimateHKY85::_runEM() {
	using coretools::str::toString;
	// run EM
	_genoDist->reset();

	logfile().startIndent("Estimating ", _genoDist->typeString(), ":");
	logfile().startIndent("Initial values:");
	_genoDist->log();

	// calculate initial LL
	double oldLL   = _LL();
	double deltaLL = std::abs(oldLL);
	logfile().list("log Likelihood = ", oldLL);
	logfile().endIndent();

	// running iterations
	logfile().startNumbering("Running EM algorithm:");
	for (size_t i = 0; i < _numEMIterations; ++i) {
		logfile().addIndent();
		logfile().number("EM Iteration:");
		_genoDist->estimate();

		const double LL = _LL();
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
	_addSites(window.sites());

	// full P
	for (size_t i = 0; i < window.size(); ++i) {
		const auto &site = window[i];
		_stats.NData += site.depth();
		if (site.empty() || site.refBase == Base::N) {
			++_stats.NMissing;
		} else if (_downSample()) {
			_sites.push_back(site);
			if (_downSample()) {
				_readIDs.push_back();
				for (const auto id : window.readIDs()[i]) { _readIDs.push_back(id + _lastReadID); }
			}
		}
	}
	_lastReadID += window.numReadsInWindow();
}

void TEstimateHKY85::_handlePerWindow(GenotypeLikelihoods::TWindow &window) {
	using BAM::End;
	using GenotypeLikelihoods::PMD::Type;
	// full P

	logfile().list("Using full data.");
	_addSites(window.sites());
	const auto LL = _runEM();

	_out.write(window.chrName(), window.from().position(), window.to().position(), window.depth(), window.numSites(),
			   window.numSitesWithData(), window.fracMissing(), _genoDist->pis(), LL);

	// downsample
	const auto nIT = _numEMIterations;
	for (const auto dOrP : _depthOrProbs) {
		constexpr coretools::TStrongArray<std::string_view, Sample> sdepthOrProb{{"probability", "probability", "maximum depth"}};
		logfile().list("Downsampling reads to a ", sdepthOrProb[_sample], " ", dOrP, ".");

		_numEMIterations = std::min<size_t>(10 * nIT, nIT / dOrP); // may need a bit longer

		if (_sample == Sample::reads) {
			const coretools::Probability p = P(dOrP);
			const auto downsampled = window.downsampleReads(p);

			_refBases.clear();
			_P_g_I_ds.clear();
			_addSites(downsampled.sites());
			const auto LL_p  = _runEM();
			_out.write(downsampled.depth(), downsampled.numSites(), downsampled.numSitesWithData(),
			           downsampled.fracMissing(), _genoDist->pis(), LL_p);
		} else {
			auto sites      = window.sites();
			size_t depth    = 0;
			size_t withData = 0;
			for (auto &site : sites) {
				if (_sample == Sample::sites) {
					const coretools::Probability p = P(dOrP);
					site.downsample(p);
				} else /* _sample == Sample::upTo */ {
					const int d = int(dOrP);
					site.limitDepth(d);
				}
				depth += site.depth();
				if (!site.empty()) ++withData;
			}
			_refBases.clear();
			_P_g_I_ds.clear();
			_addSites(sites);
			const auto LL_p  = _runEM();
			_out.write(double(depth)/window.numSites(), window.numSites(), withData, double(window.numSites() - withData) / window.numSites(),
					   _genoDist->pis(), LL_p);
		}
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

void TEstimateHKY85::_handleGenomeWide() {
	using BAM::End;
	using coretools::instances::randomGenerator;
	using GenotypeLikelihoods::PMD::Type;
	using std::back_inserter;

	_openFile();
	const auto NSites = _totSites - _totMaskedSites;

	// full
	logfile().list("Using full data.");
	const auto LL  = _runEM();
	const auto pis = _genoDist->pis();
	const auto nIT = _numEMIterations;

	std::vector<bool> keep;
	std::vector<GenotypeLikelihoods::TSite> sites_p;

	if (_downSample()) {
		sites_p.resize(_sites.size());
		for (size_t i = 0; i < _sites.size(); ++i) { sites_p[i].refBase = _sites[i].refBase; }
	}

	for (size_t r = 0; r < _nRounds; ++r) {
		logfile().list("Downsampling round ", r + 1, ".");
		if (_windows.considerRegions()) {
			_out.write("regions", "Round", r + 1);
		} else {
			_out.write("genome-wide", "Round", r + 1);
		}
		_out.write(_stats.NData / NSites, NSites, _totSites - _stats.NMissing,
		           double(_stats.NMissing - _totMaskedSites) / NSites, pis, LL);


		// downsampled
		for (size_t p = 0; p < _depthOrProbs.size(); ++p) {
			const auto prob = P(_depthOrProbs[p]);

			constexpr coretools::TStrongArray<std::string_view, Sample> sdepthOrProb{
			    {"probability", "probability", "maximum depth"}};
			logfile().list("Downsampling reads to a ", sdepthOrProb[_sample], " ", _depthOrProbs[p], ".");

			if (_sample == Sample::upToDepth) {
				_numEMIterations =
				    std::min<size_t>(10 * nIT, nIT * _depthOrProbs[0] / _depthOrProbs[p]); // may need a bit longer
			} else {
				_numEMIterations = std::min<size_t>(10 * nIT, nIT / _depthOrProbs[p]); // may need a bit longer
			}

			size_t depth    = 0;
			size_t withData = 0;
			if (_sample == Sample::reads) {
				keep.assign(_lastReadID, false);

				for (size_t k = 0; k < keep.size(); ++k) {
					if (randomGenerator().getRand() < prob) { keep[k] = true; }
				}
			}
			for (size_t s = 0; s < _sites.size(); ++s) {
				sites_p[s].data().clear(); // don't remove refBase
				switch (_sample) {
				case Sample::reads: {
					for (size_t j = 0; j < _sites[s].depth(); ++j) {
						if (keep[_readIDs[s][j]]) { sites_p[s].add(_sites[s][j]); }
					}
				} break;
				case Sample::sites: {
					for (size_t j = 0; j < _sites[s].depth(); ++j) {
						if (randomGenerator().getRand() < prob) { sites_p[s].add(_sites[s][j]); }
					}
				} break;
				default: { // Sample::upTodepth
					const auto len = std::min<size_t>(prob, _sites[s].depth());
					sites_p[s].data().reserve(len);
					std::copy(_sites[s].begin(), _sites[s].begin() + len, std::back_inserter(sites_p[s].data()));
				} break;
				}
				depth += sites_p[s].depth();
				if (!sites_p[s].empty()) ++withData;
			}

			_refBases.clear();
			_P_g_I_ds.clear();
			_addSites(sites_p);
			const auto LL_i = _runEM();
			_out.write(double(depth) / NSites, NSites, withData, double(NSites - withData) / NSites, _genoDist->pis(),
			           LL_i);
		}
		_out.endln();
	}
}

void TEstimateHKY85::run() {
	_traverseBAMWindows();
	if (_genomeWide) _handleGenomeWide();
}

TEstimateHKY85::TEstimateHKY85() {
	using BAM::RGInfo::InfoType;
	_windows.requireReference();

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
			_nRounds = parameters().get<size_t>("genomeWide");
			logfile().list("Will do ", _nRounds, " rounds of downsampling. (parameter 'genomeWide')");
		}
	} else {
		logfile().list("Will estimating genotype Distribution per window. (use 'genomeWide' for genome-wide estimation)");
		_openFile(); // will be written every Window
	}
	if (!_depthOrProbs.empty()) {
		_windows.allowDownsampling(true);
	}
}
void TEstimateHKY85::_openFile() {
	using coretools::str::toString;
	std::vector<std::string> header{"Chr", "Start", "End"};
	const auto sp = _depthOrProbs.empty() || _sample == Sample::upToDepth ? "": "p1.0_";
	header.insert(header.end(), {toString(sp, "depth"), toString(sp, "numSites"), toString(sp, "numSitesData"), toString(sp, "fracMissing")});
	_genoDist->addHeader(header, sp); //, sp);
	header.push_back(toString(sp, "LL"));

	for (const auto p : _depthOrProbs) {
		const auto sp = toString("p", p, "_");
		header.insert(header.end(), {toString(sp, "depth"), toString(sp, "numSites"), toString(sp, "numSitesData"),
									 toString(sp, "fracMissing")});
		_genoDist->addHeader(header, sp); //, sp);
		header.push_back(toString(sp, "LL"));
	}

	_out.open(_genome.outputName() + ".txt.gz");
	_out.writeHeader(header);
}
}
