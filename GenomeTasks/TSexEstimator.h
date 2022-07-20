/*
 * TDepthWriter.h
 *
 *  Created on: Jul 18, 2022
 *      Author: reckel
 */

#ifndef GENOMETASKS_TSEXESTIMATOR_H_
#define GENOMETASKS_TSEXESTIMATOR_H_

#include <string>

#include "TFile.h"
#include "TGenome.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks{

//----------------------------------------
// TSexEstimator
//----------------------------------------
class TSexEstimator:public TGenome_windows{
private:
	coretools::TOutputFile _out;
	coretools::TCountDistribution<> _distPerSite;

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	void writeDepth();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimateSex:public coretools::TTask{
public:
	TTask_estimateSex(){ _explanation = "Estimating the distribution of depth among sites and writing depth per window"; };

	void run(){
		TSexEstimator sexEstimator;
		//sexEstimator.writeDepth();
	};
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
