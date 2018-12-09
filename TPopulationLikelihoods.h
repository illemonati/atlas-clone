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

	void init();
	void fillSampleOrder();
public:
	TPopulationSamples();
	TPopulationSamples(std::string filename, TLog* logfile);
	~TPopulationSamples();

	bool hasSamples(){ return _hasSamples; };
	int numSamples(){ return _numSamples; };
	int numPopulations(){ return _numPopulations; };
	int numSamplesInPop(int population){ return numSamplesPerPop[population]; };
	void readSamples(std::string filename, TLog* logfile);
	void readSamplesFromVCFNames(std::vector<std::string> & vcfSampleNames);
	bool sampleIsUsed(const std::string & name);
	int getOrderedSampleIndex(const std::string & name);
	int startIndex(int population){ return startIndexPerPop[population]; };
	std::string getNameFromOrderedIndex(int index);
	void fillVCFOrder(std::vector<std::string> & vcfSampleNames);
	int VCF_order(const int & index){ return _VCF_order[index]; };
};

//-------------------------------------------------
//TPopulationLikelihoodReader
//-------------------------------------------------
class TPopulationLikelihoodReader{
private:
	TQualityMap phredToGTLMap;
	TVcfFileSingleLine vcfFile;
	bool vcfOpen;

	//settings
	long limitLines;
	int minDepth;
	int minNumSamplesWithData;
	double freqFilter;
	double epsilonF; //F for EM algorithm to estimate allele frequencies
	int minVariantQuality;
	bool estimateGenotypeFrequencies;
	long progressFrequency;

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
	double _genotypeFrequencies[3];
	double _alleleFrequency;
	double _MAF;

	void resetCounters();
	void closeVCF();
    void readDataFromVCF(TParameters & Parameters, TPopulationSamples & samples, TLog* logfile);
    void printProgressFrequencyFiltering(TLog* logfile);

    // EM algorithm for genotype frequencies
    void guessGenotypeFrequencies(unsigned short* phredScores, const int & numSamples);
    void estimateGenotypeFrequenciesNullModel(unsigned short* phredScores, const int & numSamples, double epsilonF);

public:
	TPopulationLikelihoodReader();
	TPopulationLikelihoodReader(TParameters & Parameters, TLog* Logfile);
	~TPopulationLikelihoodReader();

	void initialize(TParameters & Parameters, TLog* Logfile);
	void doEstimateGenotypeFrequencies(){ estimateGenotypeFrequencies = true; };

	void openVCF(std::string, TLog* logfile);
	bool readDataFromVCF(unsigned short* curLocus, TPopulationSamples & samples, TLog* logfile);
	void concludeFilters(TLog* logfile);

	std::vector<std::string>& getSampleVCFNames(){ return vcfFile.parser.samples; };
	std::string chr(){ return vcfFile.chr(); };
	long position(){ return vcfFile.position(); };
	long numLociParsed(){ return _lineCounter; };
	double* genotypeFrequencies(){ return _genotypeFrequencies; };
	double allelFrequency(){ return _alleleFrequency; };
	double MAF(){ return _MAF; };
};

//-------------------------------------------------
//TPopulationLikelihoods
//-------------------------------------------------
class TPopulationLikelihoods{
private:
	TQualityMap phredToGTLMap;

	// about vcf-file
	bool vcfRead;

	// data on individuals
	TPopulationSamples samples;

	//data on loci
	long _numLoci;
	std::map<int, std::string> chromosomes; //first SNP index and name
	std::vector<long> position;
    std::vector<unsigned short*> genotypePhredScores;
    std::vector<double> alleleFrequencies;
    bool saveAlleleFrequencies;

    //looping
    long curLocusIndex;
    std::map<int, std::string>::iterator curChrIt;
    std::map<int, std::string>::iterator nextChrIt;
    int individualStartIndex;
    int usedSampleSize;

    // read data from vcf-file
	void init();
    void openVCF(std::string, TVcfFileSingleLine & vcfFile, TLog* logfile);
    void readDataFromVCF(TParameters & Parameters, TLog* logfile);

public:
    TPopulationLikelihoods();
    TPopulationLikelihoods(TParameters & Parameters, TLog* Logfile);
    ~TPopulationLikelihoods();

    void clean();
    void doSaveAlleleFrequencies(){ saveAlleleFrequencies = true; };
    void readData(TParameters & Parameters, TLog* Logfile);

    //Looping
    void begin();
    void beginOnePop(int population);
    bool end();
    void next();
    unsigned short* curData();
    std::string curSampleName(int index);
    int curSampleSize();
    std::string curChr();
    long curPosition();

    // get main constants (n, L, D, K) and names of environmental variables
    int getNumIndividuals();
    long getNumLoci();
    unsigned short* getDataAtLocus(long index);
    void print();
};

#endif //LDLG_TDATA_H




#endif /* TPOPULATIONLIKELIHOODS_H_ */
