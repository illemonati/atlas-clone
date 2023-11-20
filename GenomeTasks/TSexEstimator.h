/*
 * TSexEstimator.h
 *
 *  Created on: Jul 20, 2022
 *      Author: raphael
 */

#ifndef GENOMETASKS_TSEXESTIMATOR_H_
#define GENOMETASKS_TSEXESTIMATOR_H_



#include "TGenome_OLD.h"
#include "coretools/Main/TTask.h"
#include "TBedReaderWindows.h"


namespace GenomeTasks{

//----------------------------------------
// TSexEstimator
//----------------------------------------
class TSexEstimator : public old::TGenome_windows {
private:
	std::vector<coretools::TCountDistribution<>> _distPerSites;
	std::vector<std::unique_ptr<BAM::TBedReaderWindows>> _regions;
	size_t _siteLimit  = 0;
	size_t _regionNum  = 0;
	bool _adaptRegions = false;
	bool _wholeGenome  = false;

	void _initializeRegion(std::string_view regionsFile, int regionNum);
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _handleAlignment(BAM::TAlignment&) override {};
	void _writeDepthPerWindow(coretools::TOutputFile &out, int num);
	void _writeHistogram(size_t regionNum);
	void _writeDepthPerChromosome(size_t regionNum);


public:
	TSexEstimator();
	void run();
};

} // end namespace


#endif /* GENOMETASKS_TSEXESTIMATOR_H_ */
