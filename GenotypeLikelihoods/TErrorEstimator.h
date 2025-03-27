/*
 * TErrorEstimator.h
 *
 */

#ifndef TERRORESTIMATOR_H_
#define TERRORESTIMATOR_H_

#include <memory>
#include <vector>

#include "TSequencedData.h"
#include "coretools/Containers/TNestedVector.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Ploidy.h"
#include "genometools/TBed.h"

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
	std::vector<std::unique_ptr<TGenotypeDistribution>> _genoDist;

	// per site
	std::vector<std::vector<SequencingError::TGenotypeFloats>> _P_g_I_dis;
	std::vector<std::vector<genometools::Base>> _refBases;

	// per data
	std::vector<coretools::TNestedVector<BAM::TSequencedData, SequencingError::TGenotypeFloats>> _data;

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

	size_t _NRegions() const noexcept {return _data.size();}
	size_t _NSites(size_t Region) const noexcept {return _data[Region].size();}
	size_t _NData(size_t Region, size_t Site) const noexcept {return _data[Region][Site].size();}

	template<bool UpdateJF> void _calculateQ() {
		for (size_t r = 0; r < _NRegions(); ++r) {
			const auto isInvariant = _genoDist[r]->ploidy() == genometools::Ploidy::haploid;
			for (size_t s = 0; s < _NSites(r); ++s) {
				const auto &P_g_I_di = _P_g_I_dis[r][s];
				for (size_t d = 0; d < _NData(r, s); ++d) {
					const auto &d_ij              = _data[r][s].get<0>()[d];
					const auto &P_bbarEdij_I_gdij = _data[r][s].get<1>()[d];
					if (isInvariant)
						_recal.model(d_ij).epsilon()->add<UpdateJF, true>(d_ij, P_g_I_di, P_bbarEdij_I_gdij);
					else
						_recal.model(d_ij).epsilon()->add<UpdateJF, false>(d_ij, P_g_I_di, P_bbarEdij_I_gdij);
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
	double _calculateLL_updatePg(size_t R, const std::vector<TSite> &sites, TGenotypeDistribution *genoDist, genometools::Ploidy Pl);

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
