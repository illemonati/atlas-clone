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
#include "TBedReaderWindows.h"


namespace GenomeTasks{

//----------------------------------------
// TSexEstimator
//----------------------------------------
class TSexEstimator : public TGenome_windows {
private:
	coretools::TOutputFile _out1;
	coretools::TOutputFile _out2;
	coretools::TCountDistribution<> _distPerSite1;
	coretools::TCountDistribution<> _distPerSite2;
	std::unique_ptr<BAM::TBedReaderWindows> _region1;
	std::unique_ptr<BAM::TBedReaderWindows> _region2;
	void _initializeRegion(std::unique_ptr<BAM::TBedReaderWindows> &region, const int num);
	void _handleWindow() override;
	void _handleAlignment() override {};
	void _considerRegion(std::unique_ptr<BAM::TBedReaderWindows> &region, coretools::TCountDistribution<> &distPerSite, coretools::TOutputFile &out);
	void _writeDepthPerWindow(coretools::TOutputFile &out, const int num);
	void _writeHistogram(coretools::TCountDistribution<> &distPerSite, const int num);


public:
	TSexEstimator();
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
		sexEstimator.writeDepth();
	};
};

}; // end namespace


#endif /* GENOMETASKS_TSEXESTIMATOR_H_ */
