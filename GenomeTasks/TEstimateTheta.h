/*
 * TEstimateTheta.h
 *
 *  Created on: Jun 4, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TESTIMATETHETA_H_
#define GENOMETASKS_TESTIMATETHETA_H_

#include <vector>

#include "coretools/Types/probability.h"
#include "genometools/BED/TBed.h"

#include "TBamWindowTraverser.h"
#include "TThetaEstimator.h"
#include "TWindow.h"

namespace GenotypeLikelihoods {
class TThetaEstimatorData;
}

namespace GenomeTasks {

//-----------------------------------
// TEstimateThetaLLSurface
//-----------------------------------
class TEstimateThetaLLSurface final : public TBamWindowTraverser {
private:
	GenotypeLikelihoods::TThetaEstimator _thetaEstimator;
	size_t _steps;

	void _bootstrapThetaEstimation();
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _handleChromosome(const genometools::TChromosome&) override {}

public:
	TEstimateThetaLLSurface();
	void run();
};

//-----------------------------------
// TEstimateThetaDownsamplingQC
//-----------------------------------
class TEstimateTheta final : public TBamWindowTraverser {
private:
	GenotypeLikelihoods::TThetaEstimator _thetaEstimator;
	GenotypeLikelihoods::TThetaOutputFile _thetaOut;
	std::vector<coretools::Probability> downSampleProbVector;
	std::vector<GenotypeLikelihoods::TThetaEstimator> estimators;

	bool _printFullData     = false;
	bool _printAll          = false;
	bool _genomeWide        = false;

	bool _onlyBootstraps    = false;
	size_t _numBootstraps = 0;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _handleChromosome(const genometools::TChromosome&) override {}

	void _addSites(GenotypeLikelihoods::TWindow &window, GenotypeLikelihoods::TThetaEstimator &thetaEstimator);

	void _bootstrapThetaEstimation();
public:
	TEstimateTheta();
	void run();
};

//-----------------------------------
// TEstimateThetaRatio
//-----------------------------------
class TEstimateThetaRatio final : public TBamWindowTraverser {
private:
	GenotypeLikelihoods::TThetaEstimatorRatio _thetaEstimatorRatio;
	genometools::TBed _region1;
	genometools::TBed _region2;

	void _initializeRegion(genometools::TBed &region, int num);
	void _addSites(const GenotypeLikelihoods::TWindow &window, GenotypeLikelihoods::TThetaEstimatorData &data, const genometools::TBed &regions);
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _handleChromosome(const genometools::TChromosome&) override {}
public:
	TEstimateThetaRatio();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATETHETA_H_ */
