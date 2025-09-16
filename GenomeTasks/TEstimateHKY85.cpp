#include "TEstimateHKY85.h"

#include <bitset>
#include <filesystem>
#include <memory>

#include "TBamWindow.h"
#include "TGenotypeDistribution.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/stringProperties.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Genotype.h"

#include "coretools/devtools.h"

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::P;
using coretools::user_assert;
using genometools::TGenotypeLikelihoods;
using genometools::Base;
using genometools::Genotype;

void TEstimateHKY85::_addSite(const GenotypeLikelihoods::TSite &Site) {
	_refBases.push_back(Site.refBase);
	_baseLikelihoods.push_back();

	TGenotypeLikelihoods P_g_I_di(P(1.));
	double sum = 1.;

	for (auto &d_ij : Site) {
		const auto P_dij_I_bbar = _siteTraverser.errorModels().sequencingErrorModels().P_dij(d_ij);
		const auto P_dij_I_b    = _siteTraverser.errorModels().postMortemDamageModels().P_dij(d_ij, P_dij_I_bbar);
		if (_downSample()) { _baseLikelihoods.push_back(Standarray::standardize(P_dij_I_b)); }

		sum = _addToPg(P_g_I_di, P_dij_I_b, sum);
	}
	_siteLikelihoods.push_back(P_g_I_di);
}

std::pair<size_t, size_t> TEstimateHKY85::_downsampeSites(double ProbOrDepth) {
	using coretools::instances::randomGenerator;

	std::bitset<BAM::TBamWindow::maxReadID> keep;
	if (_sample == Sample::reads) {
		for (size_t k = 0; k < keep.size(); ++k) {
			if (randomGenerator().getRand() < ProbOrDepth) { keep.set(k); }
		}
	}

	_siteLikelihoods.clear();
	size_t depth    = 0;
	size_t withData = 0;
	for (size_t s = 0; s < _baseLikelihoods.size(); ++s) {
		TGenotypeLikelihoods P_g_I_di(P(1.));
		double sum = 1.;
		size_t depth_s = 0;
		switch (_sample) {
		case Sample::reads: {
			for (size_t j = 0; j < _baseLikelihoods[s].size(); ++j) {
				if (keep[_readIDs[s][j]]) {
					++depth_s;
					sum = _addToPg(P_g_I_di, _baseLikelihoods[s][j], sum);
				}
			}
		} break;
		case Sample::sites: {
			for (size_t j = 0; j < _baseLikelihoods[s].size(); ++j) {
				if (randomGenerator().getRand() < ProbOrDepth) {
					++depth_s;
					sum = _addToPg(P_g_I_di, _baseLikelihoods[s][j], sum);
				}
			}
		} break;
		default: { // Sample::upTodepth

			depth_s = std::min<size_t>(ProbOrDepth, _baseLikelihoods[s].size());
			for (size_t j = 0; j < depth_s; ++j) { sum = _addToPg(P_g_I_di, _baseLikelihoods[s][j], sum); }
		} break;
		}

		_siteLikelihoods.push_back(P_g_I_di);
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
		LL.add(_genoDist->add(_siteLikelihoods[i], _refBases[i]));
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

void TEstimateHKY85::_handlePosterior() {
	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
			const auto& site = _siteTraverser.site();
			if (site.empty() || site.refBase == Base::N) continue;

			_out.write(_siteTraverser.curChr().name(),_siteTraverser.position().position() + 1, site.refBase);

			std::string sBases, sRecal, sCT, sGA, sPos;
			for (const auto &d : site) {
				sBases += genometools::base2char(d.base);
				sRecal +=
				    coretools::str::toString(P(_siteTraverser.errorModels().sequencingErrorModels().phredInt(d)), ",");
				sPos += coretools::str::toString(d.distFrom5.linear(), ",");
				const auto &model = _siteTraverser.errorModels().postMortemDamageModels().model(d);
				sCT += coretools::str::toString(model.prob(GenotypeLikelihoods::PMD::Type::CT, d), ",");
				sGA += coretools::str::toString(model.prob(GenotypeLikelihoods::PMD::Type::GA, d), ",");
			}
			sRecal.pop_back();
			sPos.pop_back();
			sCT.pop_back();
			sGA.pop_back();

			_out.write(sBases, sRecal, sCT, sGA, sPos);
			auto lik = _siteTraverser.errorModels().calculateGenotypeLikelihoods(site);
			_out.write(lik);
			_genoDist->normalize_add(lik, site.refBase);
			const auto pHet = 1 - lik[Genotype::AA] - lik[Genotype::CC] - lik[Genotype::GG] - lik[Genotype::TT];
			_out.writeln(lik, pHet);
		}
	}
}

void TEstimateHKY85::_handlePerWindow() {
	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextWindow()) {
			_refBases.clear();
			_siteLikelihoods.clear();
			_baseLikelihoods.clear();
			_readIDs.clear();

			auto &window = _siteTraverser.window();
			for (; !window.endOfSites(); window.nextSite()) {
				_addSite(window.site());
				if (_downSample() && _sample == Sample::reads) _readIDs.push_back(window.readIDs());
			}
			_out.write(_siteTraverser.curChr().name(), window.from().position(), window.to().position());

			if (_fullDepth) {
			    logfile().startIndent("Using full data:");
			    const auto LL = _runEM();

				_out.write(window.depth(), window.numSites(), window.numSitesWithData(), window.fracMissing(),
						   _genoDist->pis(), LL);
				logfile().endIndent(); // Using full data
			}

			// downsample
			const auto nIT    = _numEMIterations;
			const auto NSites = window.numSites();
			for (const auto dOrP : _depthOrProbs) {
			    constexpr coretools::TStrongArray<std::string_view, Sample> sdepthOrProb{
			        {"probability", "probability", "maximum depth"}};
			    logfile().startIndent("Downsampling reads to a ", sdepthOrProb[_sample], " ", dOrP, ":");

			    if (_sample == Sample::upToDepth) {
			        _numEMIterations = std::min<size_t>(10 * nIT, nIT * _depthOrProbs[0] / dOrP); // may need a bit longer
				} else { _numEMIterations = std::min<size_t>(10 * nIT, nIT / dOrP); // may need a bit longer
			    }

			    const auto [depth, withData] = _downsampeSites(dOrP);
			    const auto LL_i              = _runEM();
			    _out.write(double(depth) / NSites, NSites, withData, double(NSites - withData) / NSites,
			_genoDist->pis(), LL_i); logfile().endIndent(); // Downsampling...
			}
			_numEMIterations = nIT;
			_out.endln();
		}
	}
}

void TEstimateHKY85::_handleGenomeWide() {
	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
			const auto &site = _siteTraverser.site();
			_addSite(site);
			if (_downSample() && _sample == Sample::reads) {
			    _readIDs.push_back(_siteTraverser.readIDs());
			}
		}
	}

	_openFile();
	const auto NBases = _siteTraverser.numBases();
	const auto NSites = _siteTraverser.numSites();
	const auto NUsed  = _siteTraverser.numSitesWithData();
	const auto nIT    = _numEMIterations;

	// full
	std::vector<double> pis;
	double LL;
	if (_fullDepth) {
		logfile().startIndent("Using full data:");
		LL  = _runEM();
		pis = _genoDist->pis();
		logfile().endIndent(); // Using full data
	}

	for (size_t r = 0; r < _nRounds; ++r) {
		logfile().startIndent("Downsampling round ", r + 1, ":");
		_out.write("genome-wide", "Round", r + 1);
		if (_fullDepth) {
			_out.write((double)NBases / NSites, NSites, NUsed, double(NSites - NUsed) / NSites, pis, LL);
		}

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

			const auto [NBasesDown, NUsedDown] = _downsampeSites(_depthOrProbs[p]);
			const auto LL_i              = _runEM();
			_out.write((double)NBasesDown / NSites, NSites, NUsedDown, double(NSites - NUsedDown) / NSites,
					   _genoDist->pis(), LL_i);
			logfile().endIndent(); // Downsampling reads...
		}
		_out.endln();
		logfile().endIndent(); //Downsampling round...
	}
}

void TEstimateHKY85::run() {
	if (_runType == RunType::genomeWide) {
		_handleGenomeWide();
	} else if (_runType == RunType::windows) {
		_handlePerWindow();
	} else { //posterior
		_handlePosterior();
	}
}

void TEstimateHKY85::_initPosterior() {
	using coretools::str::fromString;
	_runType = RunType::posterior;

	double mu, thetaR, thetaG;
	const auto post = parameters().get("posterior");
	if (std::filesystem::exists(post)) {
		coretools::TInputFile iFile(post, coretools::FileType::Header);
		size_t iMu     = -1;
		size_t iThetaG = -1;
		size_t iThetaR = -1;
		for (size_t i = 0; i < iFile.header().size(); ++i) {
			std::string_view k = iFile.header()[i];

			if (k.size() >= 2 && k.substr(k.size() - 2) == "mu") iMu = std::min(i, iMu);
			if (k.size() >= 7 && k.substr(k.size() - 7) == "theta_r") iThetaR = std::min(i, iThetaR);
			if (k.size() >= 7 && k.substr(k.size() - 7) == "theta_g") iThetaG = std::min(i, iThetaG);
		}

		if (iFile.get("Start") == "Round") {
			mu     = fromString<double, true>(iFile.get(iMu));
			thetaR = fromString<double, true>(iFile.get(iThetaR));
			thetaG = fromString<double, true>(iFile.get(iThetaG));
		} else {
			logfile().list("Calculating posterior probabilities with average HKY85 model over windows.");
			size_t N         = 0;
			double muSum     = 0.;
			double thetaRSum = 0.;
			double thetaGSum = 0.;
			for (; !iFile.empty(); iFile.popFront()) {
				++N;
				muSum += fromString<double, true>(iFile.get(iMu));
				thetaRSum += fromString<double, true>(iFile.get(iThetaR));
				thetaGSum += fromString<double, true>(iFile.get(iThetaG));
			}
			mu     = muSum / N;
			thetaR = thetaRSum / N;
			thetaG = thetaGSum / N;
		}
	} else {
		const auto params = coretools::str::fromString<std::array<double, 3>, true>(post);
		mu                = params[0];
		thetaR            = params[1];
		thetaG            = params[2];
	}
	logfile().list("Calculating posterior probabilities with unique HKY85 model: mu = ", mu, ", theta_r = ", thetaR,
	               ", theta_g = ", thetaG);
	_genoDist = std::make_unique<GenotypeLikelihoods::THKY85>(mu, thetaR, thetaG);

	std::vector<std::string> header{"chr", "pos", "ref", "bases", "recal", "p(CT)", "p(GA)", "dist5'"};
	for (auto g = Genotype::min; g < Genotype::max; ++g) header.push_back("P(D|" + genometools::toString(g) + ")");
	for (auto g = Genotype::min; g < Genotype::max; ++g) header.push_back("P(" + genometools::toString(g) + "|D)");
	header.push_back("P(het)");

	_out.open(_siteTraverser.outputName() + "_post.txt.gz", header);
}

void TEstimateHKY85::_initEstimation() {
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
	if (type == "reads") {
		_sample = Sample::reads;
	} else if (type == "sites") {
		_sample = Sample::sites;
	} else if (type == "upToDepth") {
		_sample = Sample::upToDepth;
	} else {
		throw coretools::TUserError("Cannot understand downsampling strategy '", type, "'!");
	}

	if (parameters().exists("prob")) {
		_fullDepth = false;
		user_assert(type != "upToDepth", "Cannot downsample 'upToDepth with with probability!");
		user_assert(!parameters().exists("depth"), "Cannot use arguments 'prob' and 'depth' at the same time!");

		parameters().fill("prob", _depthOrProbs);
		user_assert(!_depthOrProbs.empty(), "Cannot read probabilities given: '", parameters().get("prob"), "'!");
		std::sort(_depthOrProbs.begin(), _depthOrProbs.end(), std::greater<>());
		user_assert(_depthOrProbs.front() <= 1, "Probability must be <= 1!");
		user_assert(_depthOrProbs.back() > 0, "Probability must be positive!");

		while (!_depthOrProbs.empty() && _depthOrProbs.front() == 1.) {
			_fullDepth = true;
			_depthOrProbs.erase(_depthOrProbs.begin()); 
		}
		if (_depthOrProbs.empty()) {
			_fullDepth = true;
			logfile().warning("You did not specify any downsampling probabilty < 1 , will not do any downsampling!");
		} else {
			logfile().list("Will do downsampling experiment using all data, and downsampling ", type,
						   " with probabilities ", _depthOrProbs, ". (parameter 'prob')");
		}
	} else if ((parameters().exists("depth"))) {
		_fullDepth   = false;
		auto depths = parameters().get<std::vector<double>>("depth");
		user_assert(!depths.empty(), "Cannot read depth given: '", parameters().get("depth"), "'!");
		std::sort(depths.begin(), depths.end(), std::greater<>());
		user_assert(depths.back() > 0, "Depth must be > 0!");

		logfile().list("Will do downsampling experiment using all data, and downsampling ", type, " to depth ", depths,
		               ". (parameter 'depth')");

		if (_sample == Sample::upToDepth) {
			for (const auto depth : depths) {
				user_assert(depth >= 1, "Must sample at least depth=1");
				_depthOrProbs.push_back(depth);
			}
		} else {
			double averageDepth;
			if (parameters().exists("averageDepth")) {
				averageDepth = parameters().get<double>("averageDepth");
			} else {
				logfile().list("No averageDepth given, will calculate it. Use 'averageDepth' to safe time!");
				averageDepth = _siteTraverser.bamFile().averageDepth();

				logfile().list("Average depth estimated to ", averageDepth);
			}

			for (const auto d : depths) {
				if (d >= averageDepth) {
					if (!_fullDepth) {
						logfile().list("Cannot downsample to depth ", d,
						               " as this is >= the average depth. Will use full depth instead.");
						_fullDepth = true;
					} else {
						logfile().list("Will ignore depth ", d,
						               ", this is >= the average depth and full depth will already be done for depth ",
						               depths.front(), ".");
					}
				} else {
					_depthOrProbs.emplace_back(d / averageDepth);
				}
			}
			if (_depthOrProbs.empty()) {
				_fullDepth = true;
				logfile().warning("You did not specify any downsampling depth < ", averageDepth,
				                  ", will not do any downsampling!");
			} else {
				logfile().list("The resulting downsampling probabilities are: ", _depthOrProbs, ".");
			}
		}
	} else {
		_fullDepth = true;
		if (parameters().exists("type")) {
			logfile().warning("Cannot downsample ", parameters().exists("type"),
							  " without argument 'prob' or 'depth'!");
		}
	}
	if (!_downSample()) {
		logfile().list("Will not do downsampling. (use parameter 'prob')");
	}

	// genomewide?
	_runType = parameters().exists("genomeWide") ? RunType::genomeWide : RunType::windows;
	if (_runType == RunType::genomeWide) {
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
	if (!_depthOrProbs.empty() && _sample == Sample::reads) {
		_siteTraverser.requireReadIDs();
	}

	DEV_ASSERT(_fullDepth || _downSample()); // otherwise, there is nothing to do
}

TEstimateHKY85::TEstimateHKY85() {
	_siteTraverser.requireSingleBAM();
	_siteTraverser.requireReference();
	_siteTraverser.filterRefN();
	_siteTraverser.skipEmpty();

	if (parameters().exists("posterior")) {
		_initPosterior();
	} else {
		_initEstimation();
	}
}

void TEstimateHKY85::_openFile() {
	using coretools::str::toString;
	std::vector<std::string> header{"Chr", "Start", "End"};
	if (_fullDepth) {
		const auto sp = !_downSample() || _sample == Sample::upToDepth ? "" : "p1.0_";
		header.insert(header.end(), {toString(sp, "depth"), toString(sp, "numSites"), toString(sp, "numSitesData"),
									 toString(sp, "fracMissing")});
		_genoDist->addHeader(header, sp); //, sp);
		header.push_back(toString(sp, "LL"));
	}

	for (const auto p : _depthOrProbs) {
		const auto sp = toString("p", p, "_");
		header.insert(header.end(), {toString(sp, "depth"), toString(sp, "numSites"), toString(sp, "numSitesData"),
									 toString(sp, "fracMissing")});
		_genoDist->addHeader(header, sp); //, sp);
		header.push_back(toString(sp, "LL"));
	}

	_out.open(_siteTraverser.outputName() + ".txt.gz");
	_out.writeHeader(header);
}
}
