/*
 * TRecalibrationEM.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMESTIMATOR_H_
#define TRECALIBRATIONEMESTIMATOR_H_

#include <array>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "SequencingError/TEpsilon.h"
#include "genometools/GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TReadGroups.h"
#include "TSite.h"
#include "coretools/Types/probability.h"

namespace BAM {
class TSequencedBase;
}
namespace GenotypeLikelihoods::oldPMD {
class TModels;
}
namespace GenotypeLikelihoods {
namespace SequencingError {
class TWithRecal;
}
} // namespace GenotypeLikelihoods
namespace GenotypeLikelihoods {
namespace SequencingError {
class TModels;
}
} // namespace GenotypeLikelihoods

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator {
private:
	std::vector<TSite> _sites;
	std::vector<TGenotypeLikelihoods> _P_g_I_dis;
	std::vector<TGenotypeLikelihoods> _P_bbarEdij_I_gdijs;

	const BAM::TReadGroupMap *_readGroupMap;
	const BAM::TReadGroups *_readGroups;

	TModels* _recal;
	oldPMD::TModels* _pmd;
	std::unique_ptr<TGenotypeDistribution> _genoDist;

	std::vector<TEpsilon*> _epsilons;
	std::vector<TRho*> _rhos;

	RecalEstimatorTools::TRecalDataTables _dataTables;

	// variables for estimation
	int _numEMIterations;
	double _minDeltaLL;
	int _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	unsigned int _minRequiredObservations;
	std::string _recalFile; // file name in case a file with model is provided

	size_t _numSitesDepthTwoOrMore();
	void _initializeModels();
	void _runEM();

	// functions to estimate theta_epsilon (sequencing error rates)
	void _estimateRho_updatePbbar();

	template<bool updateJF, bool isInvariant> void _calculateQ() {
		size_t ij = 0;
		for (size_t i = 0; i < _sites.size(); ++i) {
			const auto &P_g_I_di = _P_g_I_dis[i];
			for (auto &d_ij : _sites[i]) {
				if (!d_ij) continue;

				const auto &P_bbar_I_gdij = _P_bbarEdij_I_gdijs[ij++];
				_recal->model(d_ij).epsilon()->add<updateJF, isInvariant>(d_ij, P_g_I_di, P_bbar_I_gdij);
			}
		}
	}

	void _solveDerivative() {
		if (_genoDist->isInvariant()) _calculateQ<true, true>();
		else _calculateQ<true, false>(); 
		for (auto& e: _epsilons) e->solveJxF();
	}

	void _calculateQ() {
		if (_genoDist->isInvariant()) _calculateQ<false, true>();
		else _calculateQ<false, false>(); 
	}

	void _updateEpsilon(double deltaDeltaLL);
	double _calculateLL_updatePg();

public:
	TRecalibrationEMEstimator(const BAM::TReadGroups *ReadGroups, const BAM::TReadGroupMap *ReadGroupMap);

	void addSite(const TSite &site) {if (!site.empty()) _sites.emplace_back(site);}

	// function to estimate
	void performEstimation(std::string_view outputName, TModels &SequencingErrorModels, oldPMD::TModels &PmdModels);

	void calcLL(TModels &SequencingErrorModels, oldPMD::TModels &PmdModels);
};

}; // namespace SequencingError

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
