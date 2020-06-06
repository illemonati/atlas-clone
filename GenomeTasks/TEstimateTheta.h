/*
 * TEstimateTheta.h
 *
 *  Created on: Jun 4, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TESTIMATETHETA_H_
#define GENOMETASKS_TESTIMATETHETA_H_

#include "TGenome.h"
#include "TThetaEstimator.h"

namespace GenomeTasks{

//-----------------------------------
// TEstimateTheta_base
//-----------------------------------
class TEstimateTheta_base:public TGenome_windows{
protected:
	GenotypeLikelihoods::TThetaEstimator _thetaEstimator;
	GenotypeLikelihoods::TThetaOutputFile _thetaOut;

	void _addSites(TWindow_base & window, GenotypeLikelihoods::TThetaEstimator & thetaEstimator);
	void _addSites();

public:
	TEstimateTheta_base(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//-----------------------------------
// TEstimateTheta
//-----------------------------------
class TEstimateTheta:public TEstimateTheta_base{
private:
	bool _printAll;

	bool _estimateTheta();
	void _handleWindow();

public:
	TEstimateTheta(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateTheta();
};


//-----------------------------------
// TEstimateThetaGenomeWide
//-----------------------------------
class TEstimateThetaGenomeWide:public TEstimateTheta_base{
private:
	uint32_t _numBootstraps;
	bool _onlyBootstraps;

	void _bootstrapThetaEstimation();
	void _handleWindow();
public:
	TEstimateThetaGenomeWide(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateThetaGenomeWide();
};

//-----------------------------------
// TEstimateThetaLLSurface
//-----------------------------------
class TEstimateThetaLLSurface:public TEstimateTheta_base{
private:
	uint32_t _steps;

	void _bootstrapThetaEstimation();
	void _handleWindow();
public:
	TEstimateThetaLLSurface(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateThetaLLSurface();
};

//-----------------------------------
// TEstimateThetaDownsamplingQC
//-----------------------------------
class TEstimateThetaDownsamplingQC:public TEstimateTheta_base{
private:
	std::vector<double> downSampleProbVector;
	std::vector<GenotypeLikelihoods::TThetaEstimator> estimators;
	bool _printFullData;

	//tmp
	TWindow_base destination;

	void _handleWindow();
public:
	TEstimateThetaDownsamplingQC(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void runQC();
};

//-----------------------------------
// TEstimateThetaRatio
//-----------------------------------
class TEstimateThetaRatio:public TGenome_windows{
private:
	GenotypeLikelihoods::TThetaEstimatorRatio _thetaEstimatorRatio;
	BAM::TBedReaderWindows* _region1;
	BAM::TBedReaderWindows* _region2;


	void _initializeRegion(TParameters & Params, BAM::TBedReaderWindows* region, const char num);
	void _addSites(GenotypeLikelihoods::TThetaEstimatorData & data, BAM::TBedReaderWindows & regions);
	void _handleWindow();

public:
	TEstimateThetaRatio(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateThetaRation();
};

}; // end namespace

#endif /* GENOMETASKS_TESTIMATETHETA_H_ */

