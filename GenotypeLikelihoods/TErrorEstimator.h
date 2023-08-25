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
#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TReadGroups.h"
#include "TSite.h"
#include "coretools/Types/probability.h"

namespace BAM {class TSequencedBase;}
namespace GenotypeLikelihoods::SequencingError {
class TEpsilon;
class TRho;
}

namespace GenotypeLikelihoods {

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TErrorEstimator {
private:
	std::vector<TSite> _sites;
	std::vector<TGenotypeLikelihoods> _P_g_I_dis;
	std::vector<TGenotypeLikelihoods> _P_bbarEdij_I_gdijs;

	const BAM::TReadGroups *_readGroups;
	BAM::TReadGroupMap _readGroupMap;

	SequencingError::TModels _recal;
	PMD::TModels _pmd;
	std::unique_ptr<TGenotypeDistribution> _genoDist;

	std::vector<SequencingError::TEpsilon*> _epsilons;
	std::vector<SequencingError::TRho*> _rhos;
	std::vector<PMD::TPsi*> _psis;


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
	void _estimatePMD_Rho_updatePbbar();

	template<bool updateJF, bool isInvariant> void _calculateQ() {
		size_t ij = 0;
		for (size_t i = 0; i < _sites.size(); ++i) {
			const auto &P_g_I_di = _P_g_I_dis[i];
			for (auto &d_ij : _sites[i]) {
				if (!d_ij) continue;

				const auto &P_bbar_I_gdij = _P_bbarEdij_I_gdijs[ij++];
				_recal.model(d_ij).epsilon()->add<updateJF, isInvariant>(d_ij, P_g_I_di, P_bbar_I_gdij);
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
	TErrorEstimator(const BAM::TReadGroups &ReadGroups);
	void addSite(const TSite &site) {if (!site.empty()) _sites.emplace_back(site);}

	// function to estimate
	void estimate(std::string_view outputName);
	void calcLL();
};

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
