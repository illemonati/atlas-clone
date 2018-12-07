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
	std::map<std::string, int> populations; // name and pop index
	std::map<std::string, int> samples; //name and pop index
	std::map<std::string, int> sampleOrder; //stores order of samples such that samples of the same population are together

	void init();
public:
	TPopulationSamples();
	TPopulationSamples(std::string filename, TLog* logfile);
	~TPopulationSamples();

	bool hasSamples(){ return _hasSamples; };
	int numSamples(){ return _numSamples; };
	int numPopulations(){ return numPopulations; };
	int numSamplesInPop(int pop){ return numSamplesPerPop[pop]; };
	void readSamples(std::string filename, TLog* logfile);
	bool sampleIsUsed(const std::string & name);
	int getOrderedSampleIndex(const std::string & name);
	std::string getNameFromOrderedIndex(int index);
	void fillVCFOrder(int* sampleVcfIndex, std::vector<std::string> & vcfSampleNames);
};

//-------------------------------------------------
//TPopulationLikelihoods
//-------------------------------------------------
class TPopulationLikelihoods{
private:
	TLog* logfile;
	TQualityMap phredToGTLMap;

	// about vcf-file
	bool vcfRead;

	// about env-file
	std::vector<std::string> individualNamesFromEnv;

	//data on loci
	long numAcceptedLoci;
	std::vector<std::string> names;

	// data on individuals
	int numSamples;


    // read data from vcf-file
    void openVCF(std::string, TVcfFileSingleLine & vcfFile);
    void getParametersForReadingVCF(TParameters & Parameters, int & limitLines, int & minDepth, double & maxMissing, double & epsilonF, double & freqFilter, int & progressFrequency);
    void readDataFromVCF(TParameters & Parameters);
    std::string getChrPosString(std::string, long);
    void printProgressFrequencyFiltering(long, long &, long &, long &, struct timeval &);
    void matchIndividualNames(int, std::vector<int> &);
    int findPosOfVcfIndivInEnv(std::string);
    void findExtraIndivInEnvFile(std::vector<std::string> &);

    // EM algorithm for genotype frequencies
    void fillInitialEstimateOfGenotypeFrequencies(double* , unsigned short* );
    void estimateGenotypeFrequenciesNullModel(double* genotypeFrequencies, unsigned short* phredScores, double epsilonF);

    // clean up
	void clean();

public:
	TPopulationLikelihoods(TParameters* Parameters, TLog* Logfile);
	~TPopulationLikelihoods(){
		clean();
    };

	// read all data (vcf and .env)
	void readDataAndInitializeParams();

    // get main constants (n, L, D, K) and names of environmental variables
    int getNumIndiv();
    long getNumLoci();
    int getNumEnvVar();
    int getNumLatFac();
    std::vector<std::string> getEnvVarNames();

    // store genotype likelihoods
    std::vector<unsigned short *> genotypePhredScores;


};

#endif //LDLG_TDATA_H




#endif /* TPOPULATIONLIKELIHOODS_H_ */
