/*
 * TPopulationLikelihoods.h
 *
 *  Created on: Dec 7, 2018
 *      Author: phaentu
 */

#ifndef TPOPULATIONLIKELIHOODS_H_
#define TPOPULATIONLIKELIHOODS_H_

//
// Created by Madleina Caduff on 19.10.18.
//

#ifndef LDLG_TDATA_H
#define LDLG_TDATA_H

#include <iostream>
#include <math.h>
#include "stringFunctions.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TVcfFile.h"
#include "TQualityMap.h"

//------------------------------------------------
//TPopulationSamples
//------------------------------------------------
class TPopulationSamples{
private:
	bool _hasSamples;
	int _numPopulations;
	int _numSamples;
	int* numSamplesPerPop;
	int* startIndexPerPop;
	std::map<std::string, int> populations; // name and pop index
	std::map<std::string, int> samples; //name and pop index
	std::map<std::string, int> sampleOrder; //stores order of samples such that samples of the same population are together
	int* _VCF_order;
	bool _VCF_order_initialized;

	void _init();
	void fillSampleOrder();
public:
	TPopulationSamples();
	TPopulationSamples(std::string filename, TLog* logfile);
	~TPopulationSamples();

	bool hasSamples(){ return _hasSamples; };
	int numSamples(){ return _numSamples; };
	int numPopulations(){ return _numPopulations; };
	std::string getPopulationName(int index);
	int numSamplesInPop(int population){ return numSamplesPerPop[population]; };
	void readSamples(std::string filename, TLog* logfile);
	void readSamplesFromVCFNames(std::vector<std::string> & vcfSampleNames);
	bool sampleIsUsed(const std::string & name);
	int getOrderedSampleIndex(const std::string & name);
	int startIndex(int population){ return startIndexPerPop[population]; };
	std::string getNameFromOrderedIndex(int index);
	void fillVCFOrder(std::vector<std::string> & vcfSampleNames);
	int VCF_order(const int & index){ return _VCF_order[index]; };
	uint8_t* getPointerToDataInPop(uint8_t* data, int population){ return &data[3*startIndexPerPop[population]]; };
	int numSamplesMissingInPop(bool* sampleMissing, int population);
	int numSamplesWithDataInPop(bool* sampleMissing, int population);
};

//------------------------------------------------
// TPopulationLikehoodStorage
//------------------------------------------------
class TPopulationLikehoodSample{
public:
	uint8_t phredLikelihood_0;
	uint8_t phredLikelihood_1;
	uint8_t phredLikelihood_2;

	bool isHaploid;
	bool isMissing;

	uint8_t operator[](int genotype){
		if(genotype == 0) return phredLikelihood_0;
		if(genotype == 1) return phredLikelihood_1;
		if(isHaploid) throw "Genotype has to be 0 or 1 for haploid samples!";
		if(genotype == 2) return phredLikelihood_2;
		throw "Genotype has to be 0, 1 or 2 for diploid samples!";
	};
};

//------------------------------------------------
// TPopulationLikehoodStorage
// class used when reading line by line
//------------------------------------------------
class TPopulationLikehoodStorage{
public:
	int numSamples;
	TPopulationLikehoodSample* samples;

	TPopulationLikehoodStorage(int NumSamples){
		resize(NumSamples);
	};

	~TPopulationLikehoodStorage(){
		clear();
	};

	void clear(){
		if(numSamples > 0)
			delete[] samples;
	};

	void resize(int NumSamples){
		if(NumSamples != numSamples){
			clear();
			numSamples = NumSamples;
			samples = new TPopulationLikehoodSample[numSamples];
		}
	};

	TPopulationLikehoodSample & operator[](int index){
		return samples[index];
	};
};

//------------------------------------------------
// TGenotypeFreqencies
//------------------------------------------------
class TGenotypeFreqencies{
private:

	void setFrequencenciesToZero();
	void normalize();
	void ensureAllFrequenciesAreNonZero();

public:
	double diploidFrequencies[3];
	double haploidFrequencies[2];
	double alleleFrequency;
	double MAF;
	int numDiploidSamples;
	int numHaploidSamples;

	TGenotypeFreqencies();

	void clear();
	bool isMonomorphic();
	void guess(TPopulationLikehoodSample* samples, int numSamples);
	void estimate(TPopulationLikehoodSample* samples, int numSamples, TQualityMap & phredToGTLMap, double epsilonF);
};

//-------------------------------------------------
//TPopulationLikelihoodReader
//-------------------------------------------------
class TPopulationLikelihoodReader{
private:
	TQualityMap phredToGTLMap;
	TVcfFileSingleLine vcfFile;
	bool vcfOpen;
	std::istream* trueFreq;
	bool trueFreqFileOpen = false;

	//settings
	long limitLines;
	int minDepth;
	int minNumSamplesWithData;
	double freqFilter;
	double epsilonF; //F for EM algorithm to estimate allele frequencies
	int minVariantQuality;
	bool estimateGenotypeFrequencies;
	bool storeTrueAlleleFreq;
	long progressFrequency;
	std::string trueAlleleFreqFile;

	//counters
	struct timeval startTime;
	bool vcfParsingStarted;
	long _lineCounter;
	long _notBialleleicCounter;
	long _missingSNPCounter;
	long _lowFreqSNPCounter;
	long _lowVariantQualityCounter;
	long _noPLCounter;
	long _numAcceptedLoci;

	//tmp variables used for reading
	TGenotypeFreqencies genoFrequencies;
	double _trueAlleleFrequency;

	void _init();
	void resetCounters();
	void closeVCF();
	void closeTrueAlleleFreqFile();
    void readDataFromVCF(TParameters & Parameters, TPopulationSamples & samples, TLog* logfile);
    void printProgressFrequencyFiltering(TLog* logfile);

public:
	TPopulationLikelihoodReader();
	TPopulationLikelihoodReader(TParameters & Parameters, TLog* Logfile, bool saveAlleleFreq);
	~TPopulationLikelihoodReader();

	void initialize(TParameters & Parameters, TLog* Logfile, bool saveAlleleFreq);
	void doEstimateGenotypeFrequencies(){ estimateGenotypeFrequencies = true; };
	void doSaveTrueAlleleFrequencies(){ storeTrueAlleleFreq = true; }

	void openVCF(std::string, TLog* logfile);
	void openTrueAlleleFrequenciesFile(std::string filename, bool isZipped);
    bool filterVCF(uint8_t* data, bool* sampleIsMissing, TPopulationSamples & samples, TLog* logfile, std::string & outputName);
    bool readDataFromVCF(TPopulationLikehoodStorage & data, TPopulationSamples & samples, TLog* logfile);
	bool readDataFromVCF(TPopulationLikehoodSample* data, TPopulationSamples & samples, TLog* logfile);
	void concludeFilters(TLog* logfile);

	std::vector<std::string>& getSampleVCFNames(){ return vcfFile.parser.samples; };
	std::string chr(){ return vcfFile.chr(); };
	long position(){ return vcfFile.position(); };
	long numLociParsed(){ return _lineCounter; };
	long numAcceptedLoci(){ return _numAcceptedLoci; };
	TGenotypeFreqencies* genotypeFrequencies(){ return &genoFrequencies; };
	double* diploidGenotypeFrequencies(){ return genoFrequencies.diploidFrequencies; };
	double allelFrequency(){ return genoFrequencies.alleleFrequency; };
	double trueAlleleFrequency(){ return _trueAlleleFrequency; };
	double MAF(){ return genoFrequencies.MAF; };
	int numSamplesWithData();
	int numSamplesWithDataInPopulation(int population);
};

//------------------------------------------------
//TVcfFilter
//------------------------------------------------
class TVcfFilter{
private:
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;
	TLog* logfile;

	// data on individuals
	TPopulationSamples samples;

	//data on loci
	long _numLoci;

//    //looping
//    long curLocusIndex;
//    std::map<int, std::string>::iterator curChrIt;
//    std::map<int, std::string>::iterator nextChrIt;
//    int individualStartIndex;
//    int usedSampleSize;

public:
	TVcfFilter(TParameters & Parameters, TLog* logfile);
    void filterVCF(TParameters & Parameters);


};

//------------------------------------------------
//TAlleleFreqEstimator
//------------------------------------------------
class TAlleleFreqEstimator{
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;
	TLog* logfile;

	// data on individuals
	TPopulationSamples samples;

	//data on loci
	long _numLoci;

public:
	TAlleleFreqEstimator(TParameters & Parameters, TLog* logfile);
    void estimateAlleleFreq(TParameters & Parameters);
};

//-------------------------------------------------
//TPopulationLikelihoods
//-------------------------------------------------
class TPopulationLikelihoods{
private:
	TQualityMap phredToGTLMap;

	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;

	// data on individuals
	TPopulationSamples samples;

	//data on loci
	long _numLoci;
	std::map<int, std::string> chromosomes; //first SNP index and name
	std::vector<long> position;
    std::vector<TPopulationLikehoodSample*> data;
    std::vector<double> alleleFrequencies;
    std::vector<double> trueAlleleFrequencies;
    bool saveAlleleFrequencies;
    bool saveTrueAlleleFrequencies;

    //looping
    long curLocusIndex;
    std::map<int, std::string>::iterator curChrIt;
    std::map<int, std::string>::iterator nextChrIt;
    int individualStartIndex;
    int usedSampleSize;

    // read data from vcf-file
	void init();
    void readDataFromVCF(TParameters & Parameters, TLog* logfile);

public:
    TPopulationLikelihoods();
    TPopulationLikelihoods(TParameters & Parameters, TLog* Logfile);
    ~TPopulationLikelihoods();

    void clean();
    void doSaveAlleleFrequencies(){ saveAlleleFrequencies = true; };
    void doSaveTrueAlleleFrequencies(){ saveTrueAlleleFrequencies = true; };
    std::vector<double> donateAlleleFrequencies(){
    	std::vector<double> alleleFrequenciesCopy = alleleFrequencies;
    	alleleFrequencies.clear();
    	return alleleFrequenciesCopy;
    };
    std::vector<double> donateTrueAlleleFrequencies(){
    	std::vector<double> trueAlleleFrequenciesCopy = trueAlleleFrequencies;
    	trueAlleleFrequencies.clear();
    	return trueAlleleFrequenciesCopy;
    }
    void readData(TParameters & Parameters, TLog* Logfile);
    std::string getVCFName(){ return vcfFilename; };

    //Looping
    void useAllSamples();
    void limitToSinglePopulation(int population);
    void begin();
    void beginAll();
    void beginOnePop(int population);
    bool end();
    void next();
    TPopulationLikehoodSample* curData();
    std::string curSampleName(int index);
    int curSampleSize();
    std::string curChr();
    long curPosition();

    // get main constants (n, L, D, K) and names of environmental variables
    int getNumIndividuals();
    long getNumLoci();
    TPopulationLikehoodSample* getDataAtLocus(long index);
    void print();
};

#endif //LDLG_TDATA_H




#endif /* TPOPULATIONLIKELIHOODS_H_ */
