/*
 * TThetaEstimatorData.h
 *
 *  Created on: Jun 17, 2018
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATORDATA_H_
#define TTHETAESTIMATORDATA_H_

#include "TSite.h"
#include "TFile.h"

using namespace GenotypeLikelihoods;

//---------------------------------------------------------------
//TThetaEstimatorTemporaryFile
//---------------------------------------------------------------
class TThetaEstimatorTemporaryFile{
private:
	std::string filename;
	gzFile fp;
	int sizeOfData;

	bool isOpen;
	bool isOpenForWriting;
	bool isOpenForReading;
	bool wasWritten;

public:
	TThetaEstimatorTemporaryFile();
	TThetaEstimatorTemporaryFile(std::string filename);
	~TThetaEstimatorTemporaryFile(){
		clean();
	};

	void init(std::string filename);

	void openForWriting();
	void openForReading();
	void close();
	void clean();
	bool isEOF();

	void save(TGenotypeLikelihoods & genoLik);
	bool read(TGenotypeLikelihoods & genoLik);
};


//---------------------------------------------------------------
//TThetaEstimatorData
//---------------------------------------------------------------
class TThetaEstimatorData{
protected:
	//counters
	long numSitesCoveredTwiceOrMore;
	long totNumSitesAdded;
	double cumulativeDepth;

	TBaseFrequencies initialBaseFreq;
	bool isBootstrapped;
	long numBootstrappedSites;
	int maxKforPoissonPlusOne;
	double* poissonProb;
	uint8_t* numBootstrapRepsPerEntry;
	bool numBootstrapRepsPerEntryInitialized;

	bool readState;
	long curSite;
	uint8_t curRep;

	//tmp variables
	int g;
	double sum;
	TGenotypePosteriorProbabilities P_g_oneSite;

	virtual void saveSite(TGenotypeLikelihoods & genoLik){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void emptyStorage(){};
	void fillPoissonForBootstrap(const double lambda);
	virtual void _begin(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void readNext(){ throw "Not available in TThetaEstimatorData base class!"; };

public:
	TThetaEstimatorData();
	virtual ~TThetaEstimatorData(){
		clear();
		delete[] poissonProb;
		if(numBootstrapRepsPerEntryInitialized)
			delete[] numBootstrapRepsPerEntry;
	};
	long numSitesWithData;

	void add(const TSite & site, TGenotypeLikelihoods & genoLik);
	void clear();

	void bootstrap(TRandomGenerator & randomGenerator);
	void clearBootstrap();

	virtual bool begin();
	virtual bool next();
	virtual bool isEnd(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual TGenotypeLikelihoods& curGenotypeLikelihoods(){ throw "Not available in TThetaEstimatorData base class!"; };

	long size(){return totNumSitesAdded;};
	long sizeWithData(){
		if(isBootstrapped)
			return numBootstrappedSites;
		return numSitesWithData;
	};

	void addToHeader(std::vector<std::string> & header, const std::string prefix);
	void writeSite(TOutputFile* out);
	void fillBaseFreq(double* baseFreq);
	virtual void fillP_G(TGenotypeData & P_G, const TGenotypeData & pGenotype);
	virtual double calcLogLikelihood(const TGenotypeData & pGenotype);
};

//---------------------------------------------------------------
//TThetaEstimatorDataVector
//---------------------------------------------------------------
class TThetaEstimatorDataVector:public TThetaEstimatorData{
private:
	std::vector<TGenotypeData> sites;
	std::vector<TGenotypeData>::iterator siteIt;

	void saveSite(TGenotypeLikelihoods & genoLik);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataVector();
	virtual ~TThetaEstimatorDataVector(){
		clear();
	};

	void _begin();
	bool isEnd();
	TGenotypeLikelihoods& curGenotypeLikelihoods();
	void fillP_G(TGenotypeData & P_G, const TGenotypeData & pGenotype);
	double calcLogLikelihood(const TGenotypeData & pGenotype);
};

//---------------------------------------------------------------
//TThetaEstimatorDataFile
//---------------------------------------------------------------
class TThetaEstimatorDataFile:public TThetaEstimatorData{
protected:
	std::string dataFileName;

	TThetaEstimatorTemporaryFile sites;
	TGenotypeLikelihoods genotypeLikelihoods;

	void saveSite(TGenotypeLikelihoods & genoLik);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataFile(std::string TmpFileName);
	~TThetaEstimatorDataFile(){
		clear();
	};

	void _begin();
	bool isEnd();
	TGenotypeLikelihoods& curGenotypeLikelihoods();
};





#endif /* TTHETAESTIMATORDATA_H_ */
