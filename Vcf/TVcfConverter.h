//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include "../stringFunctions.h"
#include "../TLog.h"
#include "../TParameters.h"
#include "TVcfFile.h"
#include <TPopulationLikelihoodLocus.h>
#include "TGenotypeFrequencies.h"
#include "TGLF.h"

struct TFilterCounters {
    TFilterCounters();

    unsigned long lineCounter;
    unsigned long notBialleleicCounter;
    unsigned long missingSNPCounter;
    unsigned long lowFreqSNPCounter;
    unsigned long lowVariantQualityCounter;
    unsigned long noPLCounter;
    unsigned long notOnChrCounter;
    unsigned long numAcceptedLoci;
};

struct TFilterParameters {
    TFilterParameters();

    long limitLines;
    int minDepth;
    int minNumSamplesWithData;
    double freqFilter;
    double epsilonF; //F for EM algorithm to estimate allele frequencies
    int minVariantQuality;
    bool estimateGenotypeFrequencies;
    std::vector<std::string> chromosomesToKeep;
};

class TVcfFilters{
private:
    void specifyChromosomesToKeep(TParameters & Parameters);
    void readFilteringOptions(TParameters & Parameters);
    TLog * logfile;

public:
    TFilterParameters filterParameters;
    TFilterCounters filterCounters;

    long progressFrequency;
    void printProgressFrequencyFiltering(const struct timeval & startTime);
    int filterOnDepth(TSampleLikelihoods* data, TVcfFileSingleLine & vcfFile, unsigned int numSamples);
    void concludeFilters(const struct timeval & startTime);

    TVcfFilters();
    void init(TParameters & Parameters, TLog* logfile);
};

class TVcfReader {
private:
    TLog* logfile;
    TGlfConverter glfConverter;
    TGenotypeFrequencies genoFrequencies;

public:
    TVcfFileSingleLine vcfFile;
    // info about vcf file
    std::string filename;
    bool isZipped;
    unsigned int numSamples;

    TVcfReader();

    void init(TParameters & Params, TLog* Logfile);
    void openVCF(std::string name);
    bool readOneLineVcf(TVcfFilters & vcfFilters, TSampleLikelihoods * data, std::string & locusName);
};

class TVcfConverter {
protected:
    std::string outname;
    TLog * logfile;
    TVcfReader vcfReader;
    TVcfFilters vcfFilters;
    TGlfConverter glfConverter;
    static int baseToNumber(char base, const std::string & marker);

public:
    TVcfConverter(TLog * Logfile, TParameters & Params);
    ~TVcfConverter();

    void readOutputName(TParameters & Params);
    void readVcfAndWriteFile();
    virtual void writeData(TSampleLikelihoods * data, const std::string & locusName);

    };

class TVcfToBeagle : protected TVcfConverter {
private:
    TOutputFileZipped * beagleFile;
    // beagle
    void writeBeagleHeader();
    void writeData(TSampleLikelihoods * data, const std::string & locusName) override ;

public:
    TVcfToBeagle(TParameters &Params, TLog *Logfile);
    void vcfToBeagle(TParameters & Params);
};

class TVcfToLFMM : protected TVcfConverter {
private:
    TOutputFilePlain * lfmmFile;
    TOutputFilePlain * lociNamesFile;
    // lfmm
    void writeLFMMHeader();
    void writeData(TSampleLikelihoods * data, const std::string & locusName) override ;
    void writeLFMM();
    void storePosteriorGenotypes(TSampleLikelihoods * data);
    double computePosteriorGenotype(TSampleLikelihoods * data, int i);
    void storeLocusNames(const std::string & locusName);
    void writeLociNames();

    std::vector<double *> post_genotypes;
    std::vector<std::string> loci_names;

public:
    TVcfToLFMM(TParameters &Params, TLog *Logfile);
    ~TVcfToLFMM();
    void vcfToLFMM(TParameters & Params);
};


#endif //ATLAS_TVCFCONVERTER_H
