/*
 * TThetaEstimatorData.cpp
 *
 *  Created on: Jun 17, 2018
 *      Author: phaentu
 */

#include "TThetaEstimatorData.h"

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <ostream>
#include "TFile.h"
#include "TRandomGenerator.h"
#include "TSite.h"

namespace GenotypeLikelihoods{

//-------------------------------------------------------
//TThetaEstimatorTemporaryFile
//-------------------------------------------------------
TThetaEstimatorTemporaryFile::TThetaEstimatorTemporaryFile(){
	init("");
};

TThetaEstimatorTemporaryFile::TThetaEstimatorTemporaryFile(std::string Filename){
	init(Filename);
};

void TThetaEstimatorTemporaryFile::init(std::string Filename){
	filename = Filename;
	fp = NULL;
	sizeOfData = sizeof(double) * 10;

	isOpen = false;
	isOpenForReading = false;
	isOpenForWriting = false;
	wasWritten = false;
};

void TThetaEstimatorTemporaryFile::openForWriting(){
	if(sizeOfData == 0)
		throw "Can not open temporary data file for theta: file was not initialized!";

	//if file was written, remove it
	clean();

	//now open
	fp = gzopen(filename.c_str(), "wb");
	isOpenForWriting = true;
	isOpenForReading = false;
	wasWritten = true;

	isOpen = true;
};

void TThetaEstimatorTemporaryFile::openForReading(){
	if(!wasWritten)
		throw "Can not parse temporary file: file was never written!";

	//make sure file is closed
	close();

	//now open
	fp = gzopen(filename.c_str(), "rb");
	isOpenForWriting = false;
	isOpenForReading = true;

	isOpen = true;
};

void TThetaEstimatorTemporaryFile::close(){
	if(isOpen){
		gzclose(fp);
		isOpen = false;
		isOpenForWriting = false;
		isOpenForReading = false;
	}
};

void TThetaEstimatorTemporaryFile::clean(){
	close();
	if(wasWritten){
		remove(filename.c_str());
		wasWritten = false;
	}
};

bool TThetaEstimatorTemporaryFile::isEOF(){
	if(!isOpenForReading) return true;
	return gzeof(fp);
}

void TThetaEstimatorTemporaryFile::save(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik){
	if(!isOpenForWriting) throw "Can not add data to '" + filename + "': file is closed!";

	gzwrite(fp, genoLik.pointerToData(), sizeOfData);
};

bool TThetaEstimatorTemporaryFile::read(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik){
	if(!isOpenForReading) throw "Can not read data from '" + filename + "': file is closed!";
	if(gzread(fp, genoLik.pointerToData(), sizeOfData)!=sizeOfData){
		//is end-of-file?
		if(gzeof(fp)) return false;

		//is error
		throw "Failed to read data from temporary file!";
	}
	return true;
};

//-------------------------------------------------------
//TThetaEstimatorData
//-------------------------------------------------------
TThetaEstimatorData::TThetaEstimatorData(){
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded = 0;
	numSitesWithData = 0;
	cumulativeDepth = 0.0;
	sum = 0.0;
	g = 0;

	isBootstrapped = false;
	numBootstrappedSites = false;
	maxKforPoissonPlusOne = 16; //maximum number of bootstrapping copies to consider.
	poissonProb = new double[maxKforPoissonPlusOne];

	readState = false;
	curSite = 0;
	curRep = 0;
	numBootstrapRepsPerEntry = NULL;
	numBootstrapRepsPerEntryInitialized = false;
};

void TThetaEstimatorData::clear(){
	tmpBaseFreq.set(0.0);
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded = 0;
	numSitesWithData = 0;
	cumulativeDepth = 0.0;
	emptyStorage();
	clearBootstrap();
};

void TThetaEstimatorData::add(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genoLik){
	//assumes that emission probabilities were calculated!!
	++totNumSitesAdded;

	//add if site has data
	if(!site.empty()){
		++numSitesWithData;
		cumulativeDepth += site.depth();

		saveSite(genoLik);

		//add site to base frequency estimation
		site.addToBaseFrequencies(tmpBaseFreq);

		//count sites covered >=2
		if(site.depth() > 1)
			++numSitesCoveredTwiceOrMore;
	}
};

TBaseProbabilities TThetaEstimatorData::baseFrequencies(){
	//estimate base frequencies
	tmpBaseFreq.normalize();
	return tmpBaseFreq.asFrequencies();
};

void TThetaEstimatorData::fillP_G(GenotypeLikelihoods::TGenotypeData & P_G, const GenotypeLikelihoods::TGenotypeProbabilities & pGenotype){
	//assumes that pGenotype is set!
	P_G.set(0.0);

	//calculate P_g for each site
	begin();
	do{
		P_g_oneSite.fillBayesian(curGenotypeLikelihoods(), pGenotype);
		P_G += P_g_oneSite;
	} while(next());
};

double TThetaEstimatorData::calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities & pGenotype){
	double LL = 0.0;
	begin();
	do{
		LL += log(curGenotypeLikelihoods().weightedSum(pGenotype));
	} while(next());

	return LL;
};

void TThetaEstimatorData::addToHeader(std::vector<std::string> & header, const std::string prefix){
	header.push_back(prefix + "depth");
	header.push_back(prefix + "fracMissing");
	header.push_back(prefix + "fracTwoOrMore");
};

void TThetaEstimatorData::writeSite(coretools::TOutputFile & out){
	if(isBootstrapped){
		out << "NA";
		out << (double) (totNumSitesAdded - numBootstrappedSites) / (double) totNumSitesAdded;
		out << "NA";
		//out << "NA"; //TODO: check if this NA is needed.
	} else {
		out << cumulativeDepth / (double) totNumSitesAdded;
		out << (double) (totNumSitesAdded - numSitesWithData) / (double) totNumSitesAdded;
		out << (double) numSitesCoveredTwiceOrMore / (double) totNumSitesAdded;
	}
};

void TThetaEstimatorData::fillPoissonForBootstrap(const double lambda){
	double tmp = exp(-lambda);
	poissonProb[0] = tmp;
	for(int k=1; k<maxKforPoissonPlusOne; ++k){
		tmp = tmp * lambda / (double) k;
		poissonProb[k] = poissonProb[k-1] + tmp;
	}

	//make sure last of cumulative is 1.0
	if(poissonProb[maxKforPoissonPlusOne-1] < 0.99999)
		throw "Cumulative Poisson needs more steps: cumulative < 0.99999 at last entry!";
	poissonProb[maxKforPoissonPlusOne-1] = 1.0;
}

void TThetaEstimatorData::bootstrap(coretools::TRandomGenerator & randomGenerator){
	//make sure we start empty
	clearBootstrap();

	if(!numBootstrapRepsPerEntryInitialized){
		numBootstrapRepsPerEntry = new uint8_t[numSitesWithData];
		numBootstrapRepsPerEntryInitialized = true;
	}

	fillPoissonForBootstrap(1.0);

	//now pick among sites with data with replacement and store for each entry how many times it was chosen
	numBootstrappedSites = 0.0;
	for(long l=0; l<numSitesWithData; ++l){
		//do we use this site in the bootstrap?
		numBootstrapRepsPerEntry[l] = randomGenerator.pickOne(maxKforPoissonPlusOne, poissonProb);
		numBootstrappedSites += numBootstrapRepsPerEntry[l];
	}

	//set pointer
	isBootstrapped = true;
};

void TThetaEstimatorData::clearBootstrap(){
	numBootstrappedSites = 0;
	isBootstrapped = false;
};

bool TThetaEstimatorData::begin(){
	curSite = 0; //first site is at index zero
	curRep = 1; //index starts at one

	_begin();

	if(isBootstrapped){
		while(readState && numBootstrapRepsPerEntry[curSite] == 0){
			readNext();
		}
	}

	return(readState);
};

bool TThetaEstimatorData::next(){
	if(!isBootstrapped)
		readNext();
	else {
		if(curRep < numBootstrapRepsPerEntry[curSite]){
			++curRep;
			return true;
		} else {
			curRep = 1;
			readNext();
			while(readState && numBootstrapRepsPerEntry[curSite] == 0){
				readNext();
			}
		}
	}

	return(readState);
};

//-------------------------------------------------------
//TThetaEstimatorDataVector
//-------------------------------------------------------
TThetaEstimatorDataVector::TThetaEstimatorDataVector():TThetaEstimatorData(){};

void TThetaEstimatorDataVector::saveSite(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik){
	//store emission probabilities
	sites.emplace_back(genoLik);
};

void TThetaEstimatorDataVector::emptyStorage(){
	sites.clear();
};

void TThetaEstimatorDataVector::readNext(){
	++siteIt;
	++curSite;

	//readState = (siteIt == sites.end())

	if(siteIt == sites.end())
		readState = false;
}

void TThetaEstimatorDataVector::_begin(){
	siteIt = sites.begin();
	readState = true;
};

bool TThetaEstimatorDataVector::isEnd(){
	return siteIt == sites.end();
};

GenotypeLikelihoods::TGenotypeLikelihoods& TThetaEstimatorDataVector::curGenotypeLikelihoods(){
	return *siteIt;
}

void TThetaEstimatorDataVector::fillP_G(GenotypeLikelihoods::TGenotypeData & P_G, const GenotypeLikelihoods::TGenotypeProbabilities & pGenotype){
	//assumes that pGenotype is set!
	P_G.set(0.0);

	//calculate P_g for each site
	for(auto& it : sites){
		P_g_oneSite.fillBayesian(it, pGenotype);
		P_G += P_g_oneSite;
	}
};

double TThetaEstimatorDataVector::calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities & pGenotype) {
	double LL = 0.0;
	for(auto& it : sites){
		LL += log(it.weightedSum(pGenotype));
	}

	return LL;
};

//-------------------------------------------------------
//TThetaEstimatorDataFile
//-------------------------------------------------------
TThetaEstimatorDataFile::TThetaEstimatorDataFile(std::string TmpFileName):TThetaEstimatorData(){
	dataFileName = TmpFileName;
	sites.init(dataFileName);
	sites.openForWriting();
};

void TThetaEstimatorDataFile::emptyStorage(){
	sites.clean();
};

void TThetaEstimatorDataFile::saveSite(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik){
	sites.save(genoLik);
};

void TThetaEstimatorDataFile::readNext(){
	++curSite;
	if(curSite >= numSitesWithData)
		readState = false;
	else
		readState = sites.read(genotypeLikelihoods);
}

void TThetaEstimatorDataFile::_begin(){
	sites.openForReading();
	--curSite;
	readNext(); //read first! This is required to match begin() of a vector
}

bool TThetaEstimatorDataFile::isEnd(){
	return sites.isEOF();
};

GenotypeLikelihoods::TGenotypeLikelihoods& TThetaEstimatorDataFile::curGenotypeLikelihoods(){
	return genotypeLikelihoods;
};


}; //end namespace
