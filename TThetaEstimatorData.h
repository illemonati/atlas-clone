/*
 * TThetaEstimatorData.h
 *
 *  Created on: Jun 17, 2018
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATORDATA_H_
#define TTHETAESTIMATORDATA_H_

#include "TSite.h"

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
	TThetaEstimatorTemporaryFile(std::string filename, int numGenotypes);
	~TThetaEstimatorTemporaryFile(){
		clean();
	};

	void init(std::string filename, int numGenotypes);

	void openForWriting();
	void openForReading();
	void close();
	void clean();
	bool isEOF();

	void save(double* data);
	bool read(double* data);
};


//---------------------------------------------------------------
//TThetaEstimatorData
//---------------------------------------------------------------
class TThetaEstimatorData{
protected:
	//counters
	int numGenotypes;
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
	double* pointerToData;

	//tmp variables
	int g;
	double sum;
	double* P_g_oneSite;

	virtual void saveSite(TSite & site){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void emptyStorage(){};
	void fillPoissonForBootstrap(const double lambda);
	virtual void _begin(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void readNext(){ throw "Not available in TThetaEstimatorData base class!"; };

public:
	TThetaEstimatorData(int NumGenotypes);
	virtual ~TThetaEstimatorData(){
		clear();
		delete[] poissonProb;
		delete[] P_g_oneSite;
		if(numBootstrapRepsPerEntryInitialized)
			delete[] numBootstrapRepsPerEntry;
	};
	long numSitesWithData;

	void add(TSite & site);
	void clear();

	void bootstrap(TRandomGenerator & randomGenerator);
	void clearBootstrap();

	virtual bool begin();
	virtual bool next();
	virtual bool isEnd(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual double* curGenotypeLikelihoods(){ throw "Not available in TThetaEstimatorData base class!"; };

	long size(){return totNumSitesAdded;};
	long sizeWithData(){
		if(isBootstrapped)
			return numBootstrappedSites;
		return numSitesWithData;
	};

	void writeHeader(std::ofstream & out);
	void writeSize(std::ofstream & out);
	void fillBaseFreq(double* baseFreq);
	virtual void fillP_G(double* & P_G, double* & pGenotype);
	virtual double calcLogLikelihood(double* & pGenotype);
};

//---------------------------------------------------------------
//TThetaEstimatorDataVector
//---------------------------------------------------------------
class TThetaEstimatorDataVector:public TThetaEstimatorData{
private:
	std::vector<double*> sites;
	std::vector<double*>::iterator siteIt;

	void saveSite(TSite & site);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataVector(int NumGenotypes);
	virtual ~TThetaEstimatorDataVector(){
		clear();
	};

	void _begin();
	bool isEnd();
	double* curGenotypeLikelihoods();
	void fillP_G(double* & P_G, double* & pGenotype);
	double calcLogLikelihood(double* & pGenotype);
};

//---------------------------------------------------------------
//TThetaEstimatorDataFile
//---------------------------------------------------------------
class TThetaEstimatorDataFile:public TThetaEstimatorData{
protected:
	std::string dataFileName;

	TThetaEstimatorTemporaryFile sites;

	void saveSite(TSite & site);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataFile(int NumGenotypes, std::string TmpFileName);
	~TThetaEstimatorDataFile(){
		delete[] pointerToData;
		clear();
	};

	void _begin();
	bool isEnd();
	double* curGenotypeLikelihoods();
};





#endif /* TTHETAESTIMATORDATA_H_ */
