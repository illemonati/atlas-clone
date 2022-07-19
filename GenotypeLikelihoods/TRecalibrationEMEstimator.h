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

#include "GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TGenotypeData.h"
#include "TGenotypeDistribution.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TReadGroups.h"
#include "TSite.h"
#include "probability.h"

namespace BAM {
class TSequencedBase;
}
namespace GenotypeLikelihoods {
class TPostMortemDamage;
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
			   uint32_t MinRequiredObservations);

	size_t size() const { return _models.size(); };
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &data) const;

	// functions to estimate rho
	void addToRho(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d);
	void estimateRho();

	// functions to estimate beta
	void resetQ();
	void addToQFJ(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, coretools::Probability P_bbar_I_gd, bool updateJF = false);
	double curQ();
	void solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();

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
	std::vector<TGenotypeLikelihoods> _P_g_I_ds;
	std::vector<TGenotypeLikelihoods> _P_bbar_I_gds;
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
	void _runEM(const std::string &outputName, const TPostMortemDamage &PmdModels);

	// functions to estimate theta_epsilon (sequencing error rates)
	void _estimateRho_updatePbbar(const TPostMortemDamage &PmdModels);
	void _calculateQ_updateJF(bool updateJF);
	void _solveDerivative() {
		_calculateQ_updateJF(true);
		_modelsToEstimate.solveJxF();
	}
	void _calculateQ() { _calculateQ_updateJF(false); }
	void _updateEpsilon(const TPostMortemDamage &PmdModels, double deltaDeltaLL);
	double _calculateLL_updatePg(const TPostMortemDamage &PmdModels);
	void _writeCurrentEstimates(const std::string &filename);

public:
	TRecalibrationEMEstimator(const BAM::TReadGroups *ReadGroups, const BAM::TReadGroupMap *ReadGroupMap);

	// functions to add data
	void addSite(const TSite &site);

	// function to estimate
	void performEstimation(const std::string &outputName, TModels &SequencingErrorModels,
						   const TPostMortemDamage &PmdModels);

	void calcLL(TModels &SequencingErrorModels, const TPostMortemDamage &PmdModels);
};

}; // namespace SequencingError

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
