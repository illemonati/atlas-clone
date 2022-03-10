/*
 * TPopulationLikelihoods.h
 *
 *  Created on: Dec 7, 2018
 *      Author: phaentu
 */

#ifndef TPOPULATIONLIKELIHOODS_H_
#define TPOPULATIONLIKELIHOODS_H_

#include <iostream>
#include <math.h>
#include "TPopulationLikelihoodLocus.h"
#include "GenotypeTypes.h"
#include "stringFunctions.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TGenotypeFrequencies.h"
#include "TBed.h"
#include "TVcfFile.h"
#include "TDistances.h"
#include "TPopulation.h"
#include "TSampleLikelihoods.h"

namespace PopulationTools{

using TSampleLikelihoods = genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability>;

//------------------------------------------------
//TVcfFilter
//------------------------------------------------
class TVcfFilter{
private:
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;
	coretools::TLog* logfile;

	// data on individuals
    genometools::TPopulationSamples samples;

	//data on loci
	long _numLoci;

//    //looping
//    long curLocusIndex;
//    std::map<int, std::string>::iterator curChrIt;
//    std::map<int, std::string>::iterator nextChrIt;
//    int individualStartIndex;
//    int usedSampleSize;

public:
	TVcfFilter(coretools::TParameters & Parameters, coretools::TLog* logfile);
    void filterVCF(coretools::TParameters & Parameters);
};

//-------------------------------------------------
//TPopulationLikelihoods
//-------------------------------------------------
class TPopulationLikelihoods{
private:
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;

	// data on individuals
    genometools::TPopulationSamples samples;

	//data on loci
    std::shared_ptr<genometools::TPositionsRaw> _positions;
    std::vector<TSampleLikelihoods*> data;
    std::vector<double> alleleFrequencies;
    std::vector<double> trueAlleleFrequencies;
    bool saveAlleleFrequencies;
    bool saveTrueAlleleFrequencies;

    //looping
    long curLocusIndex;
    int individualStartIndex;
    int usedSampleSize;

    // read data from vcf-file
	void init();
    void readDataFromVCF(coretools::TParameters & Parameters, coretools::TLog* logfile);

public:
    TPopulationLikelihoods();
    TPopulationLikelihoods(coretools::TParameters & Parameters, coretools::TLog* Logfile);
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
    void readData(coretools::TParameters & Parameters, coretools::TLog* Logfile);
    std::string getVCFName(){ return vcfFilename; };

    //Looping
    void useAllSamples();
    void limitToSinglePopulation(int population);
    void begin();
    void beginAll();
    void beginOnePop(int population);
    bool end();
    void next();
    TSampleLikelihoods* curData();
    std::string curSampleName(int index);
    int curSampleSize();
    std::string curChr();
    long curPosition();
    const std::shared_ptr<genometools::TPositionsRaw> & positions() const;

    // get main constants
    size_t getNumIndividuals() const;
    size_t getNumLoci() const;
    TSampleLikelihoods* getDataAtLocus(long index);
    void print();
};

}; //end namespace

#endif /* TPOPULATIONLIKELIHOODS_H_ */
