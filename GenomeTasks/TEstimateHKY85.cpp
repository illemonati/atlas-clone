#include "TEstimateHKY85.h"

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
#include <memory>

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::P;
using coretools::user_assert;
using genometools::TGenotypeLikelihoods;
using genometools::Base;

void TEstimateHKY85::_addSites(const GenotypeLikelihoods::TWindow &Window) {
	_refBases.reserve(_refBases.size() + Window.sites().size());
	_P_g_I_ds.reserve(_P_g_I_ds.size() + Window.sites().size());
	for (size_t s = 0; s < Window.sites().size(); ++s) {
		const auto &site = Window.sites()[s];
		if (site.empty() || site.refBase == genometools::Base::N) continue;

		_refBases.push_back(site.refBase);
		_P_d_I_b.push_back();

		TGenotypeLikelihoods P_g_I_di(P(1.));
		double sum = 1.;

		for (auto &d_ij : site) {
			const auto P_dij_I_bbar = _genome.errorModels().sequencingErrorModels().P_dij(d_ij);
			const auto P_dij_I_b    = _genome.errorModels().postMortemDamageModels().P_dij(d_ij, P_dij_I_bbar);
			if (_downSample()) {
				_P_d_I_b.push_back(Standarray::standardize(P_dij_I_b));
			}

			sum = _addToPg(P_g_I_di, P_dij_I_b, sum);
		}
		_P_g_I_ds.push_back(P_g_I_di);

		if (_downSample() && _sample == Sample::reads) {
			_readIDs.push_back();
			for (const auto id : Window.readIDs()[s]) { _readIDs.push_back(id + _lastReadID); }
		}
	}
}

std::pair<size_t, size_t> TEstimateHKY85::_downsampeSites(double ProbOrDepth) {
	using coretools::instances::randomGenerator;

	std::vector<bool> keep;
	if (_sample == Sample::reads) {
		keep.assign(_lastReadID, false);
		for (size_t k = 0; k < keep.size(); ++k) {
			if (randomGenerator().getRand() < ProbOrDepth) { keep[k] = true; }
		}
	}

	_P_g_I_ds.clear();
	size_t depth    = 0;
	size_t withData = 0;
	for (size_t s = 0; s < _P_d_I_b.size(); ++s) {
		TGenotypeLikelihoods P_g_I_di(P(1.));
		double sum = 1.;
		size_t depth_s = 0;
		switch (_sample) {
		case Sample::reads: {
			for (size_t j = 0; j < _P_d_I_b[s].size(); ++j) {
				if (keep[_readIDs[s][j]]) {
					++depth_s;
					sum = _addToPg(P_g_I_di, _P_d_I_b[s][j], sum);
				}
			}
		} break;
		case Sample::sites: {
			for (size_t j = 0; j < _P_d_I_b[s].size(); ++j) {
				if (randomGenerator().getRand() < ProbOrDepth) {
					++depth_s;
					sum = _addToPg(P_g_I_di, _P_d_I_b[s][j], sum);
				}
			}
		} break;
		default: { // Sample::upTodepth

			depth_s = std::min<size_t>(ProbOrDepth, _P_d_I_b[s].size());
			for (size_t j = 0; j < depth_s; ++j) { sum = _addToPg(P_g_I_di, _P_d_I_b[s][j], sum); }
		} break;
		}

		_P_g_I_ds.push_back(P_g_I_di);
		if (depth_s > 0) {
			depth += depth_s;
			++withData;
		}
	}
	return {depth, withData};
}

double TEstimateHKY85::_LL() {
	coretools::TSumLogProbability LL{};
	for (size_t i = 0; i < _refBases.size(); ++i) {
		LL.add(_genoDist->add(_P_g_I_ds[i], _refBases[i]));
	}
	const auto sum = LL.getSum();
	user_assert(std::isfinite(sum), "LL = ", LL.getSum(), "!");
	return sum;
}

double TEstimateHKY85::_runEM() {
	_genoDist->reset();
	logfile().startIndent("Estimating ", _genoDist->typeString(), ":");

	logfile().startIndent("Initial values:");
	_genoDist->log();

	// calculate initial LL
	double oldLL   = _LL();
	double deltaLL = std::abs(oldLL);
	logfile().list("log Likelihood = ", oldLL);
	logfile().endIndent(); // Initial values:

	size_t i = 0;
	for (i = 0; i < _numEMIterations; ++i) {
		_genoDist->estimate();

		const double LL = _LL();
		deltaLL         = LL - oldLL;
		oldLL           = LL;
		if (i > 0 && deltaLL < _minDeltaLL) { break; }
	}

	logfile().startIndent("Final values after ", i + 1, " iterations:");
	_genoDist->log();
	logfile().list("Log Likelihood = ", oldLL);

	if (deltaLL < 0) {
		logfile().warning("Negative LL!");
		oldLL = -std::numeric_limits<double>::infinity();
	} else if (deltaLL > _minDeltaLL) {
		logfile().warning("EM has not converged after maximum number of iterations!");
	}

	logfile().endIndent(); // Final values

	logfile().list("delta LL = ", deltaLL);

	logfile().endIndent(); // Estimating...
	return oldLL;
}

void TEstimateHKY85::_handleGenomeWide(GenotypeLikelihoods::TWindow &Window) {
	_totSites += Window.size();
	_totMaskedSites += Window.numMaskedSites();

	_addSites(Window);
	_lastReadID += Window.numReadsInWindow();

	// full P
	for (size_t i = 0; i < Window.size(); ++i) {
		const auto &site = Window[i];
		_stats.NData += site.depth();
		if (site.empty() || site.refBase == Base::N) {
			++_stats.NMissing;
		}
	}
}

void TEstimateHKY85::_handlePerWindow(GenotypeLikelihoods::TWindow &Window) {
	// full P
	logfile().startIndent("Using full data:");

	_refBases.clear();
	_P_g_I_ds.clear();
	_P_d_I_b.clear();
	_readIDs.clear();
	_lastReadID = 0;

	_addSites(Window);
	_lastReadID = Window.numReadsInWindow();

	const auto LL = _runEM();

	_out.write(Window.chrName(), Window.from().position(), Window.to().position(), Window.depth(), Window.numSites(),
			   Window.numSitesWithData(), Window.fracMissing(), _genoDist->pis(), LL);
	logfile().endIndent(); // Using full data

	// downsample
	const auto nIT    = _numEMIterations;
	const auto NSites = Window.numSites();
	for (const auto dOrP : _depthOrProbs) {
		constexpr coretools::TStrongArray<std::string_view, Sample> sdepthOrProb{
			{"probability", "probability", "maximum depth"}};
		logfile().startIndent("Downsampling reads to a ", sdepthOrProb[_sample], " ", dOrP, ":");

		if (_sample == Sample::upToDepth) {
			_numEMIterations = std::min<size_t>(10 * nIT, nIT * _depthOrProbs[0] / dOrP); // may need a bit longer
		} else {
			_numEMIterations = std::min<size_t>(10 * nIT, nIT / dOrP); // may need a bit longer
		}

		const auto [depth, withData] = _downsampeSites(dOrP);
		const auto LL_i              = _runEM();
		_out.write(double(depth) / NSites, NSites, withData, double(NSites - withData) / NSites, _genoDist->pis(),
				   LL_i);
		logfile().endIndent(); // Downsampling...
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
	_openFile();
	const auto NSites = _totSites - _totMaskedSites;

	// full
	logfile().startIndent("Using full data:");
	const auto LL  = _runEM();
	const auto pis = _genoDist->pis();
	const auto nIT = _numEMIterations;
	logfile().endIndent(); // Using full data

	for (size_t r = 0; r < _nRounds; ++r) {
		logfile().startIndent("Downsampling round ", r + 1, ":");
		if (_windows.considerRegions()) {
			_out.write("regions", "Round", r + 1);
		} else {
			_out.write("genome-wide", "Round", r + 1);
		}
		_out.write(_stats.NData / NSites, NSites, _totSites - _stats.NMissing,
				   double(_stats.NMissing - _totMaskedSites) / NSites, pis, LL);

		// downsampled
		for (size_t p = 0; p < _depthOrProbs.size(); ++p) {
			constexpr coretools::TStrongArray<std::string_view, Sample> sdepthOrProb{
				{"probability", "probability", "maximum depth"}};
			logfile().startIndent("Downsampling reads to a ", sdepthOrProb[_sample], " ", _depthOrProbs[p], ":");

			if (_sample == Sample::upToDepth) {
				_numEMIterations =
					std::min<size_t>(10 * nIT, nIT * _depthOrProbs[0] / _depthOrProbs[p]); // may need a bit longer
			} else {
				_numEMIterations = std::min<size_t>(10 * nIT, nIT / _depthOrProbs[p]); // may need a bit longer
			}

			const auto [depth, withData] = _downsampeSites(_depthOrProbs[p]);
			const auto LL_i              = _runEM();
			_out.write(double(depth) / NSites, NSites, withData, double(NSites - withData) / NSites, _genoDist->pis(),
					   LL_i);
			logfile().endIndent(); // Downsampling reads...
		}
		_out.endln();
		logfile().endIndent(); //Downsampling round...
	}
}

void TEstimateHKY85::run() {
	_traverseBAMWindows();
	if (_genomeWide) _handleGenomeWide();
}

TEstimateHKY85::TEstimateHKY85() {
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
		throw coretools::TUserError("Cannot estimate a HK85 model with a ploidy of ", ploidy, "!");
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
			throw coretools::TUserError("Cannot combine ", type, " with prob");
		} else {
			throw coretools::TUserError("Downsampling type ", type, " does not exist");
		}
		logfile().list("Will do downsampling experiment using all data, and downsampling ", type,
					   " with probabilities ", _depthOrProbs, ". (parameter 'prob')");
	} else if (type == "upToDepth") {
		const auto depths = parameters().get<std::vector<int>>("depth");
		for (const auto depth : depths) {
			user_assert(depth > 0, "Must sample at least depth=1");
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
			user_assert(!_depthOrProbs.empty(), "Cannot du several rounds of downsampling without downsampling probabilities! Use 'prob' to specify.");
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
