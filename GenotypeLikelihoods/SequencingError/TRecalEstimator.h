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
namespace GenotypeLikelihoods::PMD {
class TModels;
}
namespace GenotypeLikelihoods {
namespace SequencingError {
class TModelRecal;
}
} // namespace GenotypeLikelihoods
namespace GenotypeLikelihoods {
namespace SequencingError {
class TModels;
}
} // namespace GenotypeLikelihoods

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
// store pointers to models so that they can be estimated
//--------------------------------------------------------------------------

class TModelVectorForEstimation {
private:
	// vector of pointers to models that require estimation
	std::vector<TModelRecal *> _models;                    // non-owning
	std::vector<std::array<TModelRecal *, 2>> _modelIndex; // non-owning
public:
	TModelVectorForEstimation() = default;
	void reset(TModels &SequencingErrorModels, const RecalEstimatorTools::TRecalDataTables &DataTables,
			   const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
			   size_t MinRequiredObservations);

	size_t size() const { return _models.size(); };
	TBaseLikelihoods P_dij(const BAM::TSequencedBase &data) const;
	TModelRecal* model(const BAM::TSequencedBase &data) noexcept {return _modelIndex[data.readGroupID][data.isSecondMate()];}
	TModelRecal* model(const BAM::TSequencedBase &data) const noexcept {return _modelIndex[data.readGroupID][data.isSecondMate()];}

	// functions to estimate rho
	void addToRho(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d);
	void estimateRho();

	// functions to estimate beta
	template<bool updateJF>
	void addToEpsilon(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d,
					  coretools::Probability P_bbar_I_gd) {
		model(data)->addToEpsilon<updateJF>(data, P_g_I_d, P_bbar_I_gd);
	}
	void solveJxF();
	void propose(double lambda=1.);
	void scaleParameters();
	size_t acceptOrReject();
	void adjust();

	double Q() const;
	double maxF() const;
	void writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const;

	std::string getModelsDefinition();
	std::string getRhoDefinition();
};

//--------------------------------------------------------------------
// TRecalibrationEMEstimator
//--------------------------------------------------------------------
class TRecalibrationEMEstimator {
private:
	std::vector<TSite> _sites;
	std::unique_ptr<TGenotypeDistribution> _genoDist;
	std::vector<TGenotypeLikelihoods> _P_g_I_dis;
	std::vector<TGenotypeLikelihoods> _P_bbarEdij_I_gdijs;
	const BAM::TReadGroupMap *_readGroupMap;
	const BAM::TReadGroups *_readGroups;

	TModelVectorForEstimation _modelsToEstimate;
	RecalEstimatorTools::TRecalDataTables _dataTables;

	// variables for estimation
	bool equalBaseFrequencies;
	int _numEMIterations;
	double _minDeltaLL;
	int _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	unsigned int _minRequiredObservations;
	bool _writeTmpTables;
	std::string _recalFile; // file name in case a file with model is provided

	size_t _numSitesDepthTwoOrMore();
	void _initializeModels(TModels &SequencingErrorModels);
	void _runEM(const std::string &outputName, const PMD::TModels &PmdModels);

	// functions to estimate theta_epsilon (sequencing error rates)
	void _estimateRho_updatePbbar(const PMD::TModels &PmdModels);
	template<bool updateJF, bool isInvariant>
	void _calculateQ() {
		size_t ij = 0;
		for (size_t i = 0; i < _sites.size(); ++i) {
			const auto &P_g_I_di = _P_g_I_dis[i];
			for (auto &d_ij : _sites[i]) {
				if (!d_ij) continue;
				auto m = _modelsToEstimate.model(d_ij);

				const auto &P_bbar_I_gdij = _P_bbarEdij_I_gdijs[ij++];
				m->addToEpsilon<updateJF, isInvariant>(d_ij, P_g_I_di, P_bbar_I_gdij);
			}
		}
	}
	void _solveDerivative() {
		if (_genoDist->isInvariant()) _calculateQ<true, true>();
		else _calculateQ<true, false>(); 
		_modelsToEstimate.solveJxF();
	}
	void _calculateQ() {
		if (_genoDist->isInvariant()) _calculateQ<false, true>();
		else _calculateQ<false, false>(); 
	}
	void _updateEpsilon(double deltaDeltaLL);
	double _calculateLL_updatePg(const PMD::TModels &PmdModels);

public:
	TRecalibrationEMEstimator(const BAM::TReadGroups *ReadGroups, const BAM::TReadGroupMap *ReadGroupMap);

	// functions to add data
	void addSite(const TSite &site);

	// function to estimate
	void performEstimation(const std::string &outputName, TModels &SequencingErrorModels,
						   const PMD::TModels &PmdModels);

	void calcLL(TModels &SequencingErrorModels, const PMD::TModels &PmdModels);
};

}; // namespace SequencingError

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
