/*
 * TErrorEstimator.h
 *
 */

#ifndef TERRORESTIMATOR_H_
#define TERRORESTIMATOR_H_

#include <memory>
#include <vector>

#include "coretools/Math/TSumLog.h"

#include "genometools/TBed.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Containers.h"

#include "PMD/TModels.h"
#include "SequencingError/TEpsilon.h"
#include "SequencingError/TRecalDataTables.h"
#include "TBamWindowTraverser.h"
#include "TGenotypeDistribution.h"
#include "TSite.h"

namespace GenotypeLikelihoods {

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TErrorEstimator final : public GenomeTasks::TBamWindowTraverser<GenomeTasks::WindowType::SingleBam> {
private:
	// per region
	std::vector<size_t> _refIDs;
	std::vector<genometools::TBed> _regions;
	std::vector<std::vector<TSite>> _regionSites;
	std::vector<std::unique_ptr<TGenotypeDistribution>> _genoDist;

	// per site
	std::vector<genometools::TGenotypeLikelihoods> _P_g_I_dis;

	// per read
	std::vector<genometools::TGenotypeLikelihoods> _P_bbarEdij_I_gdijs;

	BAM::TReadGroupMap _recalMap;
	BAM::TReadGroupMap _pmdMap;
	RecalEstimatorTools::TRecalDataTables _dataTables;

	SequencingError::TModels _recal;
	PMD::TModels _pmd;

	std::vector<SequencingError::TEpsilon*> _epsilons;
	std::vector<SequencingError::TRho*> _rhos;
	std::vector<PMD::TPsi*> _psis;

	// variables for estimation
	size_t _numEMIterations;
	double _minDeltaLL;
	size_t _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;

	bool _writeRestart = false;
	bool _onlyLL       = false;

	size_t _nPi      = -1;
	size_t _nRho     = -1;
	size_t _nPsi     = -1;
	size_t _nEpsilon = -1;

	template<bool UpdateJF> void _calculateQ() {
		size_t i  = 0;
		size_t ij = 0;
		for (size_t r = 0; r < _regionSites.size(); ++r) {
			const auto &sites      = _regionSites[r];
			const auto isInvariant = _genoDist[r]->isInvariant();
			for (const auto &site : sites) {
				const auto &P_g_I_di = _P_g_I_dis[i++];
				for (auto &d_ij : site) {
					const auto &P_bbar_I_gdij = _P_bbarEdij_I_gdijs[ij++];
					if (isInvariant)
						_recal.model(d_ij).epsilon()->add<UpdateJF, true>(d_ij, P_g_I_di, P_bbar_I_gdij);
					else
						_recal.model(d_ij).epsilon()->add<UpdateJF, false>(d_ij, P_g_I_di, P_bbar_I_gdij);
				}
			}
		}
	}

	void _solveDerivative() {
		_calculateQ<true>(); 
		for (auto& e: _epsilons) e->solveJxF();
	}

	size_t _calculateQ() {
		_calculateQ<false>(); 
		size_t nUpdated = 0;
		for (auto& e: _epsilons) nUpdated += e->acceptOrReject();
		return nUpdated;
	}

	void _updateEpsilon(double deltaDeltaLL);
	double _calculateLL_updatePg();

	template<genometools::Ploidy Pl>
	double _calculateLL_updatePg(const std::vector<TSite> &sites, TGenotypeDistribution *genoDist) {
		using genometools::Genotype;
		using genometools::Base;
		using coretools::P;

		constexpr auto PgI_init = []() {
			genometools::TGenotypeLikelihoods Ps(P(1.));
			if constexpr (Pl == genometools::Ploidy::haploid) {
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
		for (const auto &site : sites) {
			if (site.genotype == Genotype::NN) { // unknown genotype
				const auto ref = site.refBase;
				genometools::TGenotypeLikelihoods P_g_I_di = PgI_init;
				double sum                    = 1.;
				for (const auto &d_ij : site) {
					const auto P_dij_I_bbar = _recal.P_dij(d_ij);
					const auto P_dij_I_b    = _pmd.P_dij(d_ij, P_dij_I_bbar);

					LL.add(sum);
					const double sum_inv = 1./sum;
					sum = 0.;
					for (auto k = Base::min; k < Base::max; ++k) {
						const auto kk = genometools::genotype(k, k);
						P_g_I_di[kk] *= P(P_dij_I_b[k]*sum_inv);
						sum      += P_g_I_di[kk];
					}
					if constexpr (Pl == genometools::Ploidy::diploid) {
						for (const auto kl: {Genotype::AC, Genotype::AG, Genotype::AT, Genotype::CG, Genotype::CT, Genotype::GT}) {
							const auto k = genometools::first(kl);
							const auto l = genometools::second(kl);
							P_g_I_di[kl] *= P(0.5*(P_dij_I_b[k] + P_dij_I_b[l])*sum_inv);
							sum     += P_g_I_di[kl];
						}
					}
				}
				LL.add(genoDist->normalize_add(P_g_I_di, ref));
				_P_g_I_dis.push_back(P_g_I_di);
			} else { // known genotype.
				_P_g_I_dis.emplace_back(P(0.));
				_P_g_I_dis.back()[site.genotype] = P(1.); // Probability of correct genotype is 1
				double P_g                       = 1.;
				for (auto &d_ij : site) {
					const auto L_eps = _recal.P_dij(d_ij);
					const auto L_D   = _pmd.P_dij(d_ij, L_eps);
					P_g *= genoDist->getGenotypeLikelihood(L_D, site.genotype);
				}
				LL.add(P_g);
			}
		}
		return LL.getSum();
	}

	void _identifyModels();
	void _runEM();
	void _updatePbbar();

	void _writeModels(std::string_view Intro);
	void _handleSite(const TSite& Site, size_t Region);

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}

public:
	TErrorEstimator();
	void estimate();
	void calcLL();
	void run();
};

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
