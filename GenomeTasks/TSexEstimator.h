/*
 * TSexEstimator.h
 *
 *  Created on: Jul 20, 2022
 *      Author: raphael
 */

#ifndef GENOMETASKS_TSEXESTIMATOR_H_
#define GENOMETASKS_TSEXESTIMATOR_H_



#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "TBedReaderWindows.h"


namespace GenomeTasks{

//----------------------------------------
// TSexEstimator
//----------------------------------------
class TSexEstimator : public TGenome_windows {
private:
	std::vector<coretools::TCountDistribution<>> _distPerSites;
	std::vector<std::unique_ptr<BAM::TBedReaderWindows>> _regions;
	uint32_t _siteLimit;
	uint16_t _regionNum;
	bool _adaptRegions = false;
	bool _wholeGenome = false;

	void _initializeRegion(std::string regionsFile, const int regionNum);
	void _handleWindow() override;
	void _handleAlignment() override {};
	void _considerRegion(uint16_t regionNum);
	void _writeDepthPerWindow(coretools::TOutputFile &out, const int num);
	void _writeHistogram(uint16_t regionNum);
	void _writeDepthPerChromosome(uint16_t regionNum);


public:
	TSexEstimator();
	void run();
};

} // end namespace


#endif /* GENOMETASKS_TSEXESTIMATOR_H_ */
