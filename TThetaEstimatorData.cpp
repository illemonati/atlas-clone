/*
 * TThetaEstimatorData.cpp
 *
 *  Created on: Jun 17, 2018
 *      Author: phaentu
 */

#include "TThetaEstimatorData.h"


//-------------------------------------------------------
//TThetaEstimatorTemporaryFile
//-------------------------------------------------------
TThetaEstimatorTemporaryFile::TThetaEstimatorTemporaryFile(){
	init("", 0);
};

TThetaEstimatorTemporaryFile::TThetaEstimatorTemporaryFile(std::string Filename, int numGenotypes){
	init(Filename, numGenotypes);
};

void TThetaEstimatorTemporaryFile::init(std::string Filename, int numGenotypes){
	filename = Filename;
	fp = NULL;
	sizeOfData = sizeof(double) * numGenotypes;

	isOpen = false;
	isOpenForReading = false;
	isOpenForWriting = false;
	wasWritten = false;
}

void TThetaEstimatorTemporaryFile::openForWriting(){
	if(sizeOfData == 0)
		throw "Can not open temporary data file for theta: file was not initialized!";

	//if file was written, remove it
	clean();

	//now open
	fp = gzopen(filename.c_str(),"wb");
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
	fp = gzopen(filename.c_str(),"rb");
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

void TThetaEstimatorTemporaryFile::save(double* data){
	if(!isOpenForWriting) throw "Can not add data to '" + filename + "': file is closed!";

	gzwrite(fp, data, sizeOfData);
};

bool TThetaEstimatorTemporaryFile::read(double* data){
	if(!isOpenForReading) throw "Can not read data from '" + filename + "': file is closed!";
	if(gzread(fp, data, sizeOfData)!=sizeOfData){
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
TThetaEstimatorData::TThetaEstimatorData(int NumGenotypes){
	numGenotypes = NumGenotypes;
	pointerToData = NULL;
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded = 0;
	numSitesWithData = 0;
	cumulativeDepth = 0.0;
	sum = 0.0;
	g = 0;
	P_g_oneSite = new double[numGenotypes];

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
	initialBaseFreq.clear();
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded = 0;
	numSitesWithData = 0;
	cumulativeDepth = 0.0;
	emptyStorage();
	clearBootstrap();
};

void TThetaEstimatorData::add(TSite & site){
	//assumes that emission probabilities were calculated!!
	++totNumSitesAdded;

	//add if site has data
	if(site.hasData){
		++numSitesWithData;
		cumulativeDepth += site.depth();

		saveSite(site);

		//add site to base frequency estimation
		site.addToBaseFrequencies(initialBaseFreq);

		//count sites covered >=2
		if(site.depth() > 1)
			++numSitesCoveredTwiceOrMore;
	}
};

void TThetaEstimatorData::fillBaseFreq(double* baseFreq){
	initialBaseFreq.normalize();
	for(int i=0; i<4; ++i)
		baseFreq[i] = initialBaseFreq[i];
};

void TThetaEstimatorData::fillP_G(double* & P_G, double* & pGenotype){
	//assumes that pGenotype is set!
	for(g=0; g<numGenotypes; ++g)
		P_G[g] = 0.0;

	//calculate P_g for each site
	double* d;
	begin();
	do{
		sum = 0.0;
		d = curGenotypeLikelihoods();
		for(g=0; g<numGenotypes; ++g){
			P_g_oneSite[g] =  d[g] * pGenotype[g];
			sum += P_g_oneSite[g];
		}
		for(g=0; g<numGenotypes; ++g)
			P_G[g] += P_g_oneSite[g] / sum;

	} while(next());
};

double TThetaEstimatorData::calcLogLikelihood(double* & pGenotype){
	double LL = 0.0;
	double* d;
	begin();
	do{
		sum = 0.0;
		d = curGenotypeLikelihoods();
		for(g=0; g<numGenotypes; ++g)
			sum +=  d[g] * pGenotype[g];
		LL += log(sum);
	} while(next());

	return LL;
};

void TThetaEstimatorData::addToHeader(std::vector<std::string> & header, const std::string prefix){
	header.push_back(prefix + "depth");
	header.push_back(prefix + "fracMissing");
	header.push_back(prefix + "fracTwoOrMore");
};

void TThetaEstimatorData::writeSite(TOutputFile* out){
	if(isBootstrapped){
		*out << "NA";
		*out << (double) (totNumSitesAdded - numBootstrappedSites) / (double) totNumSitesAdded;
		*out << "NA";
		//out << "NA"; //TODO: check if this NA is needed.
	} else {
		*out << cumulativeDepth / (double) totNumSitesAdded;
		*out << (double) (totNumSitesAdded - numSitesWithData) / (double) totNumSitesAdded;
		*out << (double) numSitesCoveredTwiceOrMore / (double) totNumSitesAdded;
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

void TThetaEstimatorData::bootstrap(TRandomGenerator & randomGenerator){
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
TThetaEstimatorDataVector::TThetaEstimatorDataVector(int NumGenotypes):TThetaEstimatorData(NumGenotypes){};

void TThetaEstimatorDataVector::saveSite(TSite & site){
	//store emission probabilities
	sites.push_back(new double[10]);
	pointerToData = *sites.rbegin();

	for(int g=0; g<numGenotypes; ++g){
		pointerToData[g] = site.emissionProbabilities[g];
	}
};


void TThetaEstimatorDataVector::emptyStorage(){
	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt)
		delete[] (*siteIt);
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

double* TThetaEstimatorDataVector::curGenotypeLikelihoods(){
	return *siteIt;
}

void TThetaEstimatorDataVector::fillP_G(double* & P_G, double* & pGenotype){
	//assumes that pGenotype is set!
	for(g=0; g<numGenotypes; ++g)
		P_G[g] = 0.0;

	//calculate P_g for each site
	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt){
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g){
			P_g_oneSite[g] =  (*siteIt)[g] * pGenotype[g];
			sum += P_g_oneSite[g];
		}
		for(g=0; g<numGenotypes; ++g)
			P_G[g] += P_g_oneSite[g] / sum;

	}
};

double TThetaEstimatorDataVector::calcLogLikelihood(double* & pGenotype){
	double LL = 0.0;

	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt){
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g)
			sum +=  (*siteIt)[g] * pGenotype[g];
		LL += log(sum);
	}

	return LL;
};

//-------------------------------------------------------
//TThetaEstimatorDataFile
//-------------------------------------------------------
TThetaEstimatorDataFile::TThetaEstimatorDataFile(int NumGenotypes, std::string TmpFileName):TThetaEstimatorData(NumGenotypes){
	dataFileName = TmpFileName;
	sites.init(dataFileName, numGenotypes);
	sites.openForWriting();

	pointerToData = new double[numGenotypes];
};

void TThetaEstimatorDataFile::emptyStorage(){
	sites.clean();
};

void TThetaEstimatorDataFile::saveSite(TSite & site){
	sites.save(site.emissionProbabilities);
};

void TThetaEstimatorDataFile::readNext(){
	++curSite;
	if(curSite >= numSitesWithData)
		readState = false;
	else
		readState = sites.read(pointerToData);
}

void TThetaEstimatorDataFile::_begin(){
	sites.openForReading();
	--curSite;
	readNext(); //read first! This is required to match begin() of a vector
}

bool TThetaEstimatorDataFile::isEnd(){
	return sites.isEOF();
};

double* TThetaEstimatorDataFile::curGenotypeLikelihoods(){
	return pointerToData;
};



