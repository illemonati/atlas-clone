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
	//add limit to amount of sites that are processed, no limit if siteLimit==0
	_siteLimit = parameters().getParameterWithDefault<uint32_t>("siteLimit", 0);
	if(_siteLimit == 0){
		logfile().list("Will not use a site limit. (use 'siteLimit' to do so)");
	} else {
		logfile().list("Will only process up to " + std::to_string(_siteLimit) + " sites. (parameter 'siteLimit')");
	}

	// read the two regions to be used
	logfile().startIndent("Reading regions: (parameter 'chromRegions')");
	std::vector<std::string> regions;
	parameters().fillParameterIntoContainer("chromRegions", regions, ',');
	_regionNum = 0;
	for (auto &s: regions){
		_initializeRegion(s, _regionNum);
		_regionNum++;
	}
}

void TSexEstimator::_initializeRegion(std::string regionsFile, int regionNum) {
	logfile().startIndent((std::string) "Region " + std::to_string(regionNum+1) + ":");
	logfile().list((std::string) "Reading region " + std::to_string(regionNum+1) + " from file '" + regionsFile + "'...");

	_regions.push_back(std::make_unique<BAM::TBedReaderWindows>(regionsFile, _windowSize, _chromosomes, _siteLimit, &logfile()));
	coretools::TCountDistribution<> distPerSite;
	_distPerSites.push_back(distPerSite);

	logfile().done();
	if (_regions[regionNum]->getNumChromosomes()==1)
		logfile().conclude("Read " + toString(_regions[regionNum]->size()) + " sites on " + toString(_regions[regionNum]->getNumChromosomes()) + " chromosome.");
	else
		logfile().conclude("Read " + toString(_regions[regionNum]->size()) + " sites on " + toString(_regions[regionNum]->getNumChromosomes()) + " chromosomes.");
	logfile().endIndent();
}

void TSexEstimator::_handleWindow(){
	for (uint16_t i=0; i<_regionNum; i++){
		_considerRegion(i);
	}
}

void TSexEstimator::_considerRegion(uint16_t regionNum){
	if(_regions[regionNum]->containsChromosome(_window.chrName())) {
		//find chromosome of current genome window in TBedReaderWindows
		auto chromosome = _regions[regionNum]->findChromosome(_window.chrName());
		uint32_t windowNum = _window.from().position()/_windowSize;
		if (chromosome->windows.count(windowNum)){
			//find Bed window with sites to keep in current genome window
			auto window = chromosome->windows.find(windowNum)->second;
			//this iterator travels along the genome window, while pointing to the positions of sites
			genometools::TGenomePosition pos = _window.from() + (window->positions.front() - _window.from().position());
			//this iterator also travels along the genome window, but points to TSite objects instead, which can print out the depth at each site
			std::vector<GenotypeLikelihoods::TSite>::iterator siteIterator = _window.begin() + (window->positions.front() - _window.from().position());

			//this iterator that travels along the bed window, which contains the sites that are supposed to be analyzed
			std::vector<uint32_t>::iterator it = window->positions.begin();

			for (; pos.position() <=  window->positions.back(); ++pos){
				//if the position in the genome and bed window are equal, add the depth at this site to the histogram and increment the bed iterator
				if (pos.position() == *it){
					_distPerSites[regionNum].add(siteIterator->depth());
					chromosome->distPerSites.add(siteIterator->depth());
					++it;
					++siteIterator;
				//if the position in the genome and bed window are unequal, don't add the depth
				//and only increase the iterators that travel along the genome window
				} else {
					++siteIterator;
				}
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

void TSexEstimator::_writeDepthPerChromosome(uint16_t regionNum){
	std::string filename = _outputName + "_depthPerChromosome_" + std::to_string(regionNum+1) + ".txt";
	logfile().list("Writing per chromosome depth estimates to '" + filename + "'.");
	const std::vector<std::string> header = {"chromosome", "mean depth"};

	std::vector<std::string> initializedChromosomes;
	_regions[regionNum]->listInitializedChromosomes(initializedChromosomes);
	coretools::TOutputFile out(filename, header);

	for (auto &s: initializedChromosomes){
		out.writeln(s, _regions[regionNum]->findChromosome(s)->distPerSites.mean());
	}
	out.writeln("all", _distPerSites[regionNum].mean());
}

void TSexEstimator::_writeHistogram(uint16_t regionNum){
	std::string filename = _outputName + "_depthPerSiteHistogram_" + std::to_string(regionNum+1) + ".txt";
	logfile().list("Writing depth per site distribution to file '" + filename + "' ...");
	_distPerSites[regionNum].write(filename, "depth");
}


void TSexEstimator::writeDepth(){
	_traverseBAMWindows();

	//write distribution per site and per chromosome
	for (uint16_t i=0; i<_regionNum; i++){
		_writeHistogram(i);
		_writeDepthPerChromosome(i);
	}
	//write ratios
	for (uint16_t i=1; i<_regionNum; i++)
		logfile().list("Ratio of region1_meanDepth/region" + std::to_string(i+1) + "_meanDepth: " + std::to_string(_distPerSites[0].mean()/_distPerSites[i].mean()));
};


} //end namespace GenomeTasks
