/*
 * TErrorEstimator.h
 *
 */

#ifndef TERRORESTIMATOR_H_
#define TERRORESTIMATOR_H_

#include <memory>
#include <vector>

#include "TSequencedData.h"
#include "TSiteTraverser.h"
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
class TErrorEstimator  {
private:
	BAM::TSiteTraverser _siteTraverser;
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

	double _calculateLL_updatePg();
	double _updateEpsilon(double deltaDeltaLL);
	size_t _updateModels();
	void _updatePbbar();
	void _solveDerivative();
	void _calculateQ(bool UpdateJF);

	void _identifyModels();
	void _runEM();

	void _writeModels(std::string_view Intro);
	void _addSite(const TSite& Site, size_t Region);

	void _traverseSites();
public:
	TErrorEstimator();
	void estimate();
	void calcLL();
	void run();
};

}; // end namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMESTIMATOR_H_ */
