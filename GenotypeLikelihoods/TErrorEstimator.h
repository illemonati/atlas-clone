/*
 * TErrorEstimator.h
 *
 */

#ifndef TERRORESTIMATOR_H_
#define TERRORESTIMATOR_H_

#include <memory>
#include <vector>

#include "genometools/Genotypes/Ploidy.h"
#include "genometools/TBed.h"
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
	std::vector<SequencingError::TGenotypeFloats> _P_bbarEdij_I_gdijs;

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
	double _QLL;
	size_t _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	size_t _minData;

	bool _writeIts  = false;
	bool _onlyLL    = false;
	bool _noPi      = false;
	bool _noRho     = false;
	bool _noPsi     = false;
	bool _noEpsilon = false;

	template<bool UpdateJF> void _calculateQ() {
		size_t i  = 0;
		size_t ij = 0;
		for (size_t r = 0; r < _regionSites.size(); ++r) {
			const auto &sites      = _regionSites[r];
			const auto isInvariant = _genoDist[r]->ploidy() == genometools::Ploidy::haploid;
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

	size_t _calculateQ();
	double _updateEpsilon(double deltaDeltaLL);
	double _calculateLL_updatePg();
	double _calculateLL_updatePg(const std::vector<TSite> &sites, TGenotypeDistribution *genoDist, genometools::Ploidy Pl);

	void _identifyModels();
	void _runEM();
	void _updatePbbar(bool DoEps);

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
