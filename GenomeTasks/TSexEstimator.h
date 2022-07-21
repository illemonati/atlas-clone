/*
 * TSexEstimator.h
 *
 *  Created on: Jul 20, 2022
 *      Author: raphael
 */

#ifndef GENOMETASKS_TSEXESTIMATOR_H_
#define GENOMETASKS_TSEXESTIMATOR_H_



#include "TGenome.h"
#include "TTask.h"


namespace GenomeTasks{

//----------------------------------------
// TSexEstimator
//----------------------------------------
class TSexEstimator : public TGenome_windows {
private:
	void _handleWindow() override;
	void _handleAlignment() override {};

public:
	TSexEstimator();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimateSex:public coretools::TTask{
public:
	TTask_estimateSex(){ _explanation = "Estimating the distribution of depth among sites and writing depth per window"; };

	void run(){
		TSexEstimator sexEstimator;
	};
};

}; // end namespace


#endif /* GENOMETASKS_TSEXESTIMATOR_H_ */
