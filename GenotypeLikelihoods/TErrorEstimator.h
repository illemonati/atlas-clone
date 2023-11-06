/*
 * TErrorEstimator.h
 *
 */

#ifndef TERRORESTIMATOR_H_
#define TERRORESTIMATOR_H_

#include <array>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "PMD/TModels.h"
#include "SequencingError/TEpsilon.h"
#include "TGenome.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TReadGroupInfo.h"
#include "TReadGroups.h"
#include "TSite.h"
#include "coretools/Types/probability.h"
#include "genometools/BED/TBed.h"
#include "genometools/GenotypeTypes.h"

namespace BAM {class TSequencedBase;}
namespace GenotypeLikelihoods::SequencingError {
class TEpsilon;
class TRho;
}

namespace GenotypeLikelihoods {

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TErrorEstimator final : public  GenomeTasks::TGenome_windows {
private:
	// per region
	std::vector<size_t> _refIDs;
	std::vector<genometools::TBed> _regions;
	std::vector<std::vector<TSite>> _regionSites;
	std::vector<std::unique_ptr<TGenotypeDistribution>> _genoDist;
	// per site
	std::vector<TGenotypeLikelihoods> _P_g_I_dis;
	// per read
	std::vector<TGenotypeLikelihoods> _P_bbarEdij_I_gdijs;

	BAM::TReadGroupMap _rgMap;
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
	bool _noPi         = true;
	bool _noRho        = true;
	bool _noPsi        = true;
	bool _noEpsilon    = true;

	void _initializeModels();
	void _runEM();

	// functions to estimate theta_epsilon (sequencing error rates)
	void _updatePbbar();

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
	void _writeModels(std::string_view Intro);

	void _handleWindow() override;
	void _handleAlignment() override {}
	void _handleSite(const TSite& Site, size_t Region);

public:
	TErrorEstimator();
	void estimate();
	void calcLL();
	void run();
};

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
