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
#include "TRandomGenerator.h"

namespace GenotypeLikelihoods{

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

	void save(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik);
	bool read(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik);
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

	GenotypeLikelihoods::TBaseData initialBaseFreq;
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
	GenotypeLikelihoods::TGenotypeProbabilities P_g_oneSite;

	virtual void saveSite(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik){ throw "Not available in TThetaEstimatorData base class!"; };
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

	void add(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genoLik);
	void clear();

	void bootstrap(TRandomGenerator & randomGenerator);
	void clearBootstrap();

	virtual bool begin();
	virtual bool next();
	virtual bool isEnd(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual GenotypeLikelihoods::TGenotypeLikelihoods& curGenotypeLikelihoods(){ throw "Not available in TThetaEstimatorData base class!"; };

	long size(){return totNumSitesAdded;};
	long sizeWithData(){
		if(isBootstrapped)
			return numBootstrappedSites;
		return numSitesWithData;
	};

	void addToHeader(std::vector<std::string> & header, const std::string prefix);
	void writeSite(TOutputFile & out);
	void fillBaseFreq(double* baseFreq);
	virtual void fillP_G(GenotypeLikelihoods::TGenotypeData & P_G, const GenotypeLikelihoods::TGenotypeData & pGenotype);
	virtual double calcLogLikelihood(const GenotypeLikelihoods::TGenotypeData & pGenotype);
};

//---------------------------------------------------------------
//TThetaEstimatorDataVector
//---------------------------------------------------------------
class TThetaEstimatorDataVector:public TThetaEstimatorData{
private:
	std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> sites;
	std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator siteIt;

	void saveSite(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataVector();
	virtual ~TThetaEstimatorDataVector(){
		clear();
	};

	void _begin();
	bool isEnd();
	GenotypeLikelihoods::TGenotypeLikelihoods& curGenotypeLikelihoods();
	void fillP_G(GenotypeLikelihoods::TGenotypeData & P_G, const GenotypeLikelihoods::TGenotypeData & pGenotype);
	double calcLogLikelihood(const GenotypeLikelihoods::TGenotypeData & pGenotype);
};

//---------------------------------------------------------------
//TThetaEstimatorDataFile
//---------------------------------------------------------------
class TThetaEstimatorDataFile:public TThetaEstimatorData{
protected:
	std::string dataFileName;

	TThetaEstimatorTemporaryFile sites;
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;

	void saveSite(GenotypeLikelihoods::TGenotypeLikelihoods & genoLik);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataFile(std::string TmpFileName);
	~TThetaEstimatorDataFile(){
		clear();
	};

	void _begin();
	bool isEnd();
	GenotypeLikelihoods::TGenotypeLikelihoods& curGenotypeLikelihoods();
};


}; //end namespace


#endif /* TTHETAESTIMATORDATA_H_ */
