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
#include <iterator>
#include <math.h>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "PMD/TModels.h"
#include "TReadGroupInfo.h"
#include "TReadGroups.h"
#include "TSite.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TSumLog.h"
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

//---------------------------------------------------------------
// TErrorEstimator
//---------------------------------------------------------------
TErrorEstimator::TErrorEstimator()
	: _rgMap(_genome.bamFile().readGroups(), parameters().get<std::string>("pool", "")), _dataTables(_rgMap),
	  _onlyLL(parameters().exists("onlyLL")) {
	_windows.requireReference();
	std::vector<size_t> ploidies;
	parameters().fill("ploidy", ploidies, {2});

	if (parameters().exists("bed")) {
		std::vector<std::string> beds;
		parameters().fill("bed", beds);
		if (ploidies.size() == 1) ploidies.assign(beds.size(), ploidies.front());
		if (ploidies.size() != beds.size()) UERROR("You need to give as many ploidies as chromosomes, or only one ploidy!");

		_regionSites.resize(beds.size());
		for (size_t i = 0; i < beds.size(); ++i) {
			const auto& bedFile   = beds[i];
			_regions.emplace_back(bedFile, _genome.bamFile().chromosomes());

			if (ploidies[i] == 1) {
				_genoDist.push_back(std::make_unique<THKY85_mono>());
			}  else {
				_genoDist.push_back(std::make_unique<THKY85>());
			}
		}
	} else if (parameters().exists("chr")) {
		std::vector<std::string> chrs;
		parameters().fill("chr", chrs);
		if (ploidies.size() == 1) ploidies.assign(chrs.size(), ploidies.front());
		if (ploidies.size() != chrs.size()) UERROR("You need to give as many ploidies as chromosomes, or only one ploidy!");

		_regionSites.resize(chrs.size());
		for (size_t i = 0; i < chrs.size(); ++i) {
			const auto& chr   = chrs[i];
			_refIDs.push_back(_genome.bamFile().chromosomes().refID(chr));

			if (ploidies[i] == 1) {
				_genoDist.push_back(std::make_unique<THKY85_mono>());
			}  else {
				_genoDist.push_back(std::make_unique<THKY85>());
			}
		}
	} else {
		if (ploidies.size() > 1) UERROR("You did not define any regions or chromosomes, only one ploidy can be given!");
		_regionSites.emplace_back();
		if (ploidies.front() == 1) {
			_genoDist.push_back(std::make_unique<THKY85_mono>());
		} else {
			_genoDist.push_back(std::make_unique<THKY85>());
		}
	}

	const auto recalModel = parameters().get("recalModel", "intercept;quality;position;context;fragmentLength;mappingQuality;");
	const auto pmdModel   = parameters().get("pmdModel", "CT5;GA5;CT3;GA3");

	logfile().list("Initial recal model: ", recalModel);
	logfile().list("Initial pmd model: ", pmdModel);
	_recal.initialize(_rgMap.size(), recalModel);
	_recal.pool(_rgMap);
	_pmd.initialize(_rgMap.size(), pmdModel);
	_pmd.pool(_rgMap);


	// estimation parameters
	logfile().startIndent("Settings regarding the EM algorithm:");

	_numEMIterations            = parameters().get<int>("iterations", 200);
	_minDeltaLL                 = parameters().get<double>("minDeltaLL", 1e-6);
	_NewtonRaphsonNumIterations = parameters().get<int>("NRiterations", 20);
	_NewtonRaphsonMaxF          = parameters().get<double>("maxF", 1e-6);
	logfile().list("Will perform at max ", _numEMIterations, " EM iterations. (parameter 'iterations')");
	logfile().list("Will stop EM when deltaLL < ", _minDeltaLL, ". (parameter 'minDeltaLL')");
	logfile().list("Will conduct at max ", _NewtonRaphsonNumIterations, " Newton-Raphson iterations. (parameter NRiterations)");
	logfile().list("Will stop Newton-Raphson when F < ", _NewtonRaphsonMaxF, ". (parameter maxF)");

	// booleans
	_writeRestart = parameters().exists("writeRestart");
	_noPi         = parameters().exists("noPi");
	_noRho        = parameters().exists("noRho");
	_noPsi        = parameters().exists("noPsi");
	_noEpsilon    = parameters().exists("noEpsilon");

	std::string dos;
	std::string nos;
	if (_noPi) nos.append("pi, ");
	else dos.append("pi, ");
	if (_noRho) nos.append("rho, ");
	else dos.append("rho, ");
	if (_noPsi) nos.append("psi, ");
	else dos.append("psi, ");
	if (_noEpsilon) nos.append("epsilon, ");
	else dos.append("epsilon, ");
	if (!dos.empty()) logfile().list("Will estimate the following parameters: ", dos.substr(0, dos.size() - 2), ".");
	if (!nos.empty()) logfile().list("Will not estimate the following parameters: ", nos.substr(0, nos.size() - 2), ".");

	logfile().endIndent();
}

void TErrorEstimator::_initializeModels() {
	_dataTables.write(_genome.outputName());
	using coretools::str::toString;
	using BAM::Mate;
	// count data available for recal

	logfile().listFlush("Counting data available for recal ...");
	logfile().done();

	logfile().startIndent("Number of sites with data:");
	for (size_t i = 0; i < _regionSites.size(); ++i) logfile().list("Region ", i + 1, ": ", _regionSites[i].size());
	logfile().endIndent();

	logfile().conclude("Number of sites with depth > 1: ", _dataTables.nSites_g1());
	logfile().conclude("Number of bases: ", _dataTables.size());
	if (_dataTables.nSites_g1() < 100) UERROR("Less than 100 sites with depth >= 2 available - aborting estimation!");
	_P_g_I_dis.reserve(std::accumulate(_regionSites.begin(), _regionSites.end(), 0, [](auto x1, auto x2){return x1 + x2.size();}));
	_P_bbarEdij_I_gdijs.reserve(_dataTables.size());

	// identify models with data that can be estimated
	logfile().startIndent("Identifying models to estimate:");
	for (auto rg : _rgMap.readGroupsInUse()) {
		if (_dataTables[rg][Mate::first].size() == 0 && _dataTables[rg][Mate::second].size() > 0) UERROR("Second mate data but no first mate data!");

		const auto& pooledWith = _rgMap.readGroupsPooledWith(rg);
		logfile().startIndent("Readgroup ", rg, ":");
		if (pooledWith.size() > 1) logfile().list("Pooled with: ", _rgMap.readGroupsPooledWith(rg), ".");

		for (Mate mate = Mate::min; mate < Mate::max; ++mate) {
			constexpr coretools::TStrongArray<std::string_view, Mate> sMates{{"First", "Second"}};
			const auto &table = _dataTables[rg][mate];
			logfile().list(sMates[mate], " mate: ", table.size(), " bases.");
			if (table.size() > 0) {
				auto &recal = _recal.RGModel(rg)[mate];
				if (!recal->recalibrates()) UERROR("Cannot estimate recal for readgroup ", rg, ", mate ", mate, "!");

				recal->epsilon()->init(table);
				_epsilons.push_back(recal->epsilon());
				_rhos.push_back(recal->rho());
			} else {
				_recal.reset(rg, mate);
			}
		}
		if (_dataTables[rg][Mate::second].size() == 0) logfile().list("Assuming single-ended read.");

		auto &pmd = _pmd.model(rg);
		if (!pmd.hasPMD()) UERROR("Cannot estimate PMD for readgroup ", rg, "!");

		_psis.push_back(pmd.psi());
		_psis.back()->estimateInit();
		logfile().endIndent();
	}
	logfile().endIndent();
}

void TErrorEstimator::_updatePbbar() {
	using genometools::genotype;
	_P_bbarEdij_I_gdijs.clear();
	size_t i     = 0;
	for (size_t r = 0; r < _regionSites.size(); ++r) {
		const auto &sites      = _regionSites[r];
		const auto isInvariant = _genoDist[r]->isInvariant();
		for (const auto& site: sites) {
			const auto &P_g_I_di = _P_g_I_dis[i++];
			for (const auto &d_ij : site) {
				_P_bbarEdij_I_gdijs.emplace_back(0.);
				const auto P_dij_I_bbar = _recal.P_dij(d_ij);

				// PMD
				_pmd.model(d_ij).psi()->addCT(d_ij, P_g_I_di[Genotype::CC], _pmd.P_b_bar(Genotype::CC, d_ij, P_dij_I_bbar));
				_pmd.model(d_ij).psi()->addGA(d_ij, P_g_I_di[Genotype::GG], _pmd.P_b_bar(Genotype::GG, d_ij, P_dij_I_bbar));

				if (!isInvariant) for (auto g : {Genotype::AC, Genotype::CG, Genotype::CT}) {
						_pmd.model(d_ij).psi()->addCT(d_ij, P_g_I_di[g], _pmd.P_b_bar(g, d_ij, P_dij_I_bbar));
				}
				if (!isInvariant) for (auto g : {Genotype::AG, Genotype::CG, Genotype::GT}) {
						_pmd.model(d_ij).psi()->addGA(d_ij, P_g_I_di[g], _pmd.P_b_bar(g, d_ij, P_dij_I_bbar));
				}

				// Rho
				auto &P_bbarEdij_I_gdij = _P_bbarEdij_I_gdijs.back();
				for (auto a = Base::min; a < Base::max; ++a) {
					const auto aa              = genotype(a, a);
					const auto P_bbar_I_aa_dij = _pmd.P_bbar(a, d_ij, P_dij_I_bbar);
					P_bbarEdij_I_gdij[aa]      = P_bbar_I_aa_dij[d_ij.base];

					_recal.model(d_ij).rho()->add(d_ij.base, P_g_I_di[aa], P_bbar_I_aa_dij);

					if (!isInvariant) for (auto b = coretools::next(a); b < Base::max; ++b) {
						const auto ab              = genotype(a, b);
						const auto P_bbar_I_ab_dij = _pmd.P_bbar(ab, d_ij, P_dij_I_bbar);
						P_bbarEdij_I_gdij[ab]      = P_bbar_I_ab_dij[d_ij.base];

						_recal.model(d_ij).rho()->add(d_ij.base, P_g_I_di[ab], P_bbar_I_ab_dij);
					}
				}
			}
		}
	}
}

void TErrorEstimator::_updateEpsilon(double deltaLL) {
	using coretools::str::toString;
	logfile().list("Optimizing Q_beta using a Newton-Raphson algorithm.");
	const auto nTot = _epsilons.size();

	for (size_t i = 0; i < _NewtonRaphsonNumIterations; ++i) {
		logfile().startIndent("Running Newton-Raphson iteration ", i + 1, ":");
		_solveDerivative();
		double oldQ = 0.;
		for (const auto& e: _epsilons) oldQ += e->Q();
		logfile().list("Current Q_beta = ", oldQ);

		logfile().startIndent("Setting new epsilon.");
		for (auto& e: _epsilons) e->propose(1.);
		auto nUpdated = _calculateQ();
		logfile().list(toString(nUpdated), "/", nTot, " models converged.");
		logfile().endIndent();

		int shift = 1;
		while(nUpdated < nTot) {
			logfile().startIndent("Backtracing with lambda = 1/2^", shift, " :");
			for (auto &e : _epsilons) e->propose(1./(1 << shift));
			
			nUpdated = _calculateQ();
			logfile().list(toString(nUpdated), "/", nTot, " models converged.");
			logfile().endIndent();

			++shift;
			if (shift > 30) break;
		}

		double newQ = 0.;
		double maxF  = 0.;
		for (const auto& e: _epsilons) {
			newQ += e->Q();
			maxF = std::max(e->maxF(), maxF);
			e->adjust();
		}
		const auto deltaQ = newQ - oldQ;

		if (nUpdated < nTot) {
			logfile().conclude(nTot - nUpdated, " models did not improve even with log2(lambda) = ", shift - 1,
							   ", aborting Newton-Raphson.");
			logfile().endIndent();
			break;
		}

		if (maxF < _NewtonRaphsonMaxF) {
			logfile().conclude("delta Q = ", deltaQ);
			logfile().conclude("max(F) = ", maxF, " < ", _NewtonRaphsonMaxF, ", ending Newton-Raphson.");
			logfile().endIndent();
			break;
		} else {
			logfile().conclude("max(F) = ", maxF);
		}

		if (deltaQ*100 < deltaLL) {
			logfile().conclude("deltaQ = ", deltaQ, ", < deltaLL/100 = ", 0.01*deltaLL, ", ending Newton-Raphson.");
			logfile().endIndent();
			break;
		} else {
			logfile().conclude("delta Q = ", deltaQ);
		}
		logfile().endIndent();
	}
}

double TErrorEstimator::_calculateLL_updatePg() {
	_P_g_I_dis.clear();

	coretools::TSumLogProbability LL{};
	for (size_t r = 0; r < _regionSites.size(); ++r) {
		const auto &sites = _regionSites[r];
		auto &genoDist    = _genoDist[r];
		for (auto &site : sites) {
			if (site.genotype == Genotype::NN) { // unknown genotype
				const auto ref = site.refBase;
				_P_g_I_dis.emplace_back(1.); // Start at 1,1,1,1,1,1,1,1
				auto &P_g_I_di = _P_g_I_dis.back();
				size_t counter = 0;
				for (auto &d_ij : site) {
					const auto P_dij_I_bbar = _recal.P_dij(d_ij);
					const auto P_dij_I_b    = _pmd.P_dij(d_ij, P_dij_I_bbar);
					const auto P_dij_I_g    = genoDist->P_dij(P_dij_I_b);
					P_g_I_di *= P_dij_I_g;
					if (++counter > 10) {
						LL.add(coretools::normalize(P_g_I_di));
						counter = 0;
					}
				}
				LL.add(genoDist->normalize_add(P_g_I_di, ref));
			} else { // known genotype.
				_P_g_I_dis.emplace_back(0.);
				_P_g_I_dis.back()[site.genotype] = 1; // Probability of correct genotype is 1
				double P_g                       = 1.;
				for (auto &d_ij : site) {
					const auto L_eps = _recal.P_dij(d_ij);
					const auto L_D   = _pmd.P_dij(d_ij, L_eps);
					P_g *= genoDist->getGenotypeLikelihood(L_D, site.genotype);
				}
				LL.add(P_g);
			}
		}
	}
	if (!std::isfinite(LL.getSum())) UERROR("LL = ", LL.getSum(), ", you may need to pool your readgroups!");
	return LL.getSum();
}

void TErrorEstimator::_writeModels(std::string_view Intro) {
	logfile().startIndent(Intro, " pi");
	for (const auto& g: _genoDist) g->log();
	logfile().endIndent();

	logfile().startIndent(Intro, " rho");
	for (const auto& r: _rhos) r->log();
	logfile().endIndent();

	logfile().startIndent(Intro, " psi");
	for (const auto& p: _psis) p->log();
	logfile().endIndent();

	logfile().startIndent(Intro, " epsilon");
	for (const auto& e: _epsilons) e->log();
	logfile().endIndent();
}

void TErrorEstimator::_runEM() {
	using coretools::str::toString;
	// run EM
	logfile().startNumbering("Running EM algorithm:");
	_writeModels("Initial");

	// calculate initial LL
	double oldLL   = _calculateLL_updatePg();
	double deltaLL = abs(oldLL);
	logfile().conclude("Initial log Likelihood = ", oldLL);

	// running iterations
	size_t i = 0;
	for (; i < _numEMIterations; ++i) {
		logfile().number("EM Iteration:");
		logfile().addIndent();

		if (!_noPi) {
			logfile().list("Estimating pi");
			for (auto &g : _genoDist) g->estimate();
		}

		_updatePbbar();

		if (!_noEpsilon) {
			logfile().startIndent("Estimating epsilon:");
			_updateEpsilon(deltaLL);
			logfile().endIndent();
		}

		if (!_noRho) {
			logfile().list("Estimating rho");
			for (auto &rho : _rhos) rho->estimate();
		}

		if (!_noPsi) {
			logfile().list("Estimating psi");
			for (auto &psi : _psis) psi->estimate();
		}

		const double LL = _calculateLL_updatePg();
		deltaLL         = LL - oldLL;
		_writeModels("Current");

		if (_writeRestart) {
			logfile().list("Writing restart file");
			_recal.addToRGInfo(_genome.rgInfo());
			_pmd.addToRGInfo(_genome.rgInfo());
			_genome.rgInfo().write(_genome.outputName() + "_restart.json");
		}

		logfile().conclude("Current Log Likelihood = ", LL);
		logfile().conclude("delta LL = ", deltaLL);
		logfile().endIndent();

		// check if we break based on LL
		if (i > 0 && deltaLL < _minDeltaLL) {
			if (deltaLL < 0) logfile().warning("Negative LL!");
			else logfile().conclude("EM has converged (delta LL < ", _minDeltaLL, ")");
			break;
		}
		oldLL = LL;
	}
	if (i == _numEMIterations - 1) logfile().warning("EM has not converged after maximum number of iterations!");

	// finalize
	logfile().endNumbering();
}

void TErrorEstimator::estimate() {
	// initialize models
	_initializeModels();

	// run EM
	_runEM();

	// writing final estimates
	_recal.addToRGInfo(_genome.rgInfo());
	_pmd.addToRGInfo(_genome.rgInfo());
}

void TErrorEstimator::calcLL() {
	logfile().startIndent("pi");
	_genoDist.front()->log();
	logfile().endIndent();
	logfile().startIndent("rho");
	_recal.RGModel(0)[BAM::Mate::first]->rho()->log();
	logfile().endIndent();
	logfile().startIndent("epsilon");
	_recal.RGModel(0)[BAM::Mate::first]->epsilon()->log();
	logfile().endIndent();
	logfile().startIndent("psi");
	_pmd.model(0).psi()->log();
	logfile().endIndent();

	logfile().startIndent("Calculating log likelihood:");
	const double LL = _calculateLL_updatePg();
	logfile().conclude("Log Likelihood = ", LL);
}
void TErrorEstimator::_handleSite(const TSite &Site, size_t Region) {
	if (Site.empty() || Site.refBase == Base::N) return;

	_regionSites[Region].emplace_back(Site);
	_dataTables.add(Site);
	for (const auto &data : Site) _pmd.model(data).psi()->add(data, Site.refBase);
}

void TErrorEstimator::_handleWindow(GenotypeLikelihoods::TWindow& Window) {
	if (!_regions.empty()) { // Either sites
		for (size_t r = 0; r < _regions.size(); ++r) {
			auto &region = _regions[r];
			for (auto lb = region.begin(Window); lb != region.end() && Window.overlaps(*lb); ++lb) {
				logfile().list("Window overlaps with region ", r + 1, ": [", lb->from().position(), ", ", lb->to().position(), "]");
				const size_t pStart = std::max(lb->from().position(), Window.from().position()) - Window.from().position();
				const size_t pStop  = std::min(lb->to().position(), Window.to().position()) - Window.from().position();
				for (auto p = pStart; p < pStop; ++p) {
					const auto s = Window[p];
					_handleSite(Window[p], r);
				}
			}
		}
	} else { // or chromosomes
		size_t region = 0;
		if (!_refIDs.empty()) {
			const auto rIt = std::find(_refIDs.begin(), _refIDs.end(), Window.refID());
			if (rIt == _refIDs.end()) return;

			region = std::distance(_refIDs.begin(), rIt);
		}
		for (const auto &s : Window) {
			_handleSite(s, region);
		}
	}
}

void TErrorEstimator::run() {
	// read data
	_traverseBAMWindows();

	if (_onlyLL)
		calcLL();
	else
		estimate();
}

} // end namespace GenotypeLikelihoods
