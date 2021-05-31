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
#include "TTask.h"

namespace GenomeTasks{

//-----------------------------------
// TEstimateTheta_base
//-----------------------------------
class TEstimateTheta_base:public TGenome_windows{
protected:
	GenotypeLikelihoods::TThetaEstimator _thetaEstimator;
	GenotypeLikelihoods::TThetaOutputFile _thetaOut;

	void _addSites(GenotypeLikelihoods::TWindow_base & window, GenotypeLikelihoods::TThetaEstimator & thetaEstimator);
	void _addSites();

public:
	TEstimateTheta_base(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
	TEstimateTheta(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
	TEstimateThetaGenomeWide(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
	TEstimateThetaLLSurface(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateThetaLLSurface();
};

//-----------------------------------
// TEstimateThetaDownsamplingQC
//-----------------------------------
class TEstimateThetaDownsamplingQC:public TEstimateTheta_base{
private:
	std::vector<Probability> downSampleProbVector;
	std::vector<GenotypeLikelihoods::TThetaEstimator> estimators;
	bool _printFullData;

	//tmp
	GenotypeLikelihoods::TWindow_base destination;

	void _handleWindow();
public:
	TEstimateThetaDownsamplingQC(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void runQC();
};

//-----------------------------------
// TEstimateThetaRatio
//-----------------------------------
class TEstimateThetaRatio:public TGenome_windows{
private:
	GenotypeLikelihoods::TThetaEstimatorRatio _thetaEstimatorRatio;
	BAM::TBed _region1;
	BAM::TBed _region2;

	void _initializeRegion(TParameters & Parameters, BAM::TBed & region, const char num);
	void _addSites(GenotypeLikelihoods::TThetaEstimatorData & data, BAM::TBed & regions);
	void _handleWindow();

public:
	TEstimateThetaRatio(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateThetaRation();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TThetaTask:public TTask{
public:
	TThetaTask(){ _citations.insert("Kousathanas et al. (2017) Genetics"); };
};

class TTask_estimateTheta:public TThetaTask{
public:
	TTask_estimateTheta(){ _explanation = "Estimating heterozygosity (theta) per window"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TEstimateTheta estimator(Parameters, Logfile, _randomGenerator);
		estimator.estimateTheta();
	};
};

class TTask_estimateThetaGenomeWide:public TThetaTask{
public:
	TTask_estimateThetaGenomeWide(){ _explanation = "Estimating heterozygosity (theta) genome-wide"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TEstimateThetaGenomeWide estimator(Parameters, Logfile, _randomGenerator);
		estimator.estimateThetaGenomeWide();
	};
};

class TTask_thetaLLSurface:public TThetaTask{
public:
	TTask_thetaLLSurface(){	_explanation = "Calculating the theta LL surface for each window"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TEstimateThetaLLSurface estimator(Parameters, Logfile, _randomGenerator);
		estimator.estimateThetaLLSurface();
	};
};

class TTask_downsamplingThetaQC:public TThetaTask{
public:
	TTask_downsamplingThetaQC(){ _explanation = "QC recalibration by estimating theta on downsampled data for each window"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TEstimateThetaDownsamplingQC estimator(Parameters, Logfile, _randomGenerator);
		estimator.runQC();
	};
};

class TTask_estimateThetaRatio:public TThetaTask{
public:
	TTask_estimateThetaRatio(){ _explanation = "Estimate the ratio in heterozygosity (theta) between genomic regions"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TEstimateThetaRatio estimator(Parameters, Logfile, _randomGenerator);
		estimator.estimateThetaRation();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TESTIMATETHETA_H_ */

