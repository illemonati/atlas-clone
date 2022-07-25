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
}


void TSexEstimator::_handleWindow(){
	_considerRegion(_region1, _distPerSite1);
	_considerRegion(_region2, _distPerSite2);
}

void TSexEstimator::_considerRegion(std::unique_ptr<BAM::TBedReaderWindows> &region, coretools::TCountDistribution<> &distPerSite){
	if(region->containsChromosome(_window.chrName())) {
		//find chromosome of current genome window in TBedReaderWindows
		auto chromosome = region->findChromosome(_window.chrName());
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
					std::cout << pos.position() << std::endl << *it << std::endl;
					distPerSite.add(siteIterator->depth());
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

void TSexEstimator::_writeDepthPerChromosome(std::unique_ptr<BAM::TBedReaderWindows> &region, const int num){
	std::string filename = _outputName + "_depthPerChromosome_" + std::to_string(num) + ".txt.gz";
	logfile().list("Writing per chromosome depth estimates to '" + filename + "'.");
	const std::vector<std::string> header = {"chromosome", "mean depth"};

	std::vector<std::string> initializedChromosomes;
	region->listInitializedChromosomes(initializedChromosomes);
	coretools::TOutputFile out(filename, header);

	for (auto &s: initializedChromosomes){
		//out << s << region->findChromosome(s)->distPerSites.mean() << std::endl;
		if(num==2){
			std::cout << s << "\t";
			region->findChromosome(s)->distPerSites.print();
			std::cout << std::endl;}
	}
}

void TSexEstimator::_writeHistogram(coretools::TCountDistribution<> &distPerSite, const int num){
	std::string filename = _outputName + "_depthPerSiteHistogram_" + std::to_string(num) + ".txt";
	logfile().list("Writing depth per site distribution to file '" + filename + "' ...");
	distPerSite.print();
	distPerSite.write(filename, "depth");
}


void TSexEstimator::writeDepth(){
	_traverseBAMWindows();

	//write distribution per site and per chromosome
	_writeHistogram(_distPerSite1, 1);
	_writeHistogram(_distPerSite2, 2);

	_writeDepthPerChromosome(_region1, 1);
	_writeDepthPerChromosome(_region2, 2);
};


} //end namespace GenomeTasks
