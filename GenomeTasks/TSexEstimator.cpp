/*
 * TSexEstimator.cpp
 *
 *  Created on: Jul 20, 2022
 *      Author: raphael
 */


#include "TSexEstimator.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::str::toString;

TSexEstimator::TSexEstimator():TGenome_windows() {
	// read the two regions to be used
	logfile().startIndent("Reading regions:");
	_initializeRegion(_region1, 1);
	_initializeRegion(_region2, 2);
}

void TSexEstimator::_initializeRegion(std::unique_ptr<BAM::TBedReaderWindows> &region, const int num) {
	logfile().startIndent((std::string) "Region " + std::to_string(num) + ":");
	std::string regionsFile = parameters().getParameter<std::string>("region" + std::to_string(num));
	logfile().list((std::string) "Reading regions " + std::to_string(num) + " from file '" + regionsFile + " (parameter 'region" + std::to_string(num) +
				   "') ...");

	//regularly exceeds sitelimit as it is only checked at the beginning of every line and not for every site
	uint32_t siteLimit = parameters().getParameterWithDefault<uint32_t>("siteLimit", 0);


	region = std::make_unique<BAM::TBedReaderWindows>(regionsFile, _windowSize, _chromosomes, siteLimit, &logfile());
	logfile().done();
	logfile().conclude("Read " + toString(region->size()) + " sites on " + toString(region->getNumChromosomes()) +
					   " chromosomes.");
	region->print();
}


void TSexEstimator::_handleWindow(){
	_considerRegion(_region1, _distPerSite1, _out1);
	_considerRegion(_region2, _distPerSite2, _out2);
}

void TSexEstimator::_considerRegion(std::unique_ptr<BAM::TBedReaderWindows> &region, coretools::TCountDistribution<> &distPerSite, coretools::TOutputFile &out){
	if(region->containsChromosome(_window.chrName())) {
		//find chromosome of current window in TBedReaderWindows
		auto chromosome = region->findChromosome(_window.chrName());
		uint32_t windowNum = _window.from().position()/_windowSize;
		if (chromosome->windows.count(windowNum)){
			//find Bed window with sites to keep in current window
			auto window = chromosome->windows.find(windowNum)->second;
			genometools::TGenomePosition pos = _window.from();
			auto siteIterator = _window.begin();
			//remove all sites before Bed window begins
			for (; pos.position() < window->start; ++pos){
				siteIterator->clear();
				++siteIterator;
			}
			//it's possible that there are two Bed windows in the current window
			//--> it's possible that sites exist between window->start and window-> end that did not occur in the Bed file
			std::vector<uint32_t>::iterator it = window->positions.begin();
			for (; pos.position() < window->end; ++pos){
				if (pos.position() == *it){
					std::cout << pos.position() << std::endl << *it << std::endl;
					//distPerSite.add(*siteIterator.depth());
					++it;
					++siteIterator;
				} else {
					siteIterator->clear();
					++siteIterator;
				}
			}
			//remove all sites after Bed window has ended
			for (; pos < _window.to(); ++pos){
				siteIterator->clear();
				++siteIterator;
			}
		logfile().listFlush("Writing sequencing depth estimates to file ...");
		out << _window << _window.depth() << std::endl;
		logfile().done();

		logfile().listFlush("Adding per site depth to distribution ...");
		for(auto& s : _window){
			distPerSite.add(s.depth());
		}
		logfile().done();
		}
	}
}

void TSexEstimator::_writeDepthPerWindow(coretools::TOutputFile &out, const int num){
	std::string filename = _outputName + "_depthPerWindow_" + std::to_string(num) + ".txt.gz";
	logfile().list("Writing per window depth estimates to '" + filename + "'.");
	const std::vector<std::string> header = {"window", "depth"};
	out.open(filename, header);
}

void TSexEstimator::_writeHistogram(coretools::TCountDistribution<> &distPerSite, const int num){
	std::string filename = _outputName + "_depthPerSiteHistogram_" + std::to_string(num) + ".txt";
	logfile().list("Writing depth per site distribution to file '" + filename + "' ...");
	distPerSite.write(filename, "depth");
}


void TSexEstimator::writeDepth(){
	_writeDepthPerWindow(_out1, 1);
	_writeDepthPerWindow(_out2, 2);

	_traverseBAMWindows();

	//write distribution
	_writeHistogram(_distPerSite1, 1);
	_writeHistogram(_distPerSite2, 2);
};


} //end namespace GenomeTasks
