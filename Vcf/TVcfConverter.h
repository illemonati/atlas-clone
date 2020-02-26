//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include "../stringFunctions.h"
#include "../TLog.h"
#include "../TParameters.h"
#include "TVcfFile.h"

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
    int filterOnDepth(unsigned int numSamples, TVcfFileSingleLine & vcfFile);
    void concludeFilters(const struct timeval & startTime);

    TVcfFilters(TParameters & Parameters, TLog* logfile);
};

class TVcfReader {
private:
    TLog* logfile;
    static inline double phredToError(double phred);

public:
    TVcfFileSingleLine vcfFile;
    // info about vcf file
    std::string filename;
    bool isZipped;
    unsigned int numSamples;

    TVcfReader();

    void init(TParameters & Params, TLog* Logfile);
    void openVCF(std::string name);
    bool readOneLineVcf(TVcfFilters & vcfFilters, double * data);
};

class TVcfConverter {
protected:
    std::string outname;
    TLog * logfile;
    TVcfReader vcfReader;
    void readFilters(TParameters & Params);
    static int baseToNumber(char base, std::string & marker);

public:
    TVcfConverter(TLog * Logfile, TParameters & Params);
    ~TVcfConverter();

    void readOutputName(TParameters & Params);
    bool readVcf(TVcfFilters & vcfFilters);
    virtual void writeData(double * data);

    };

class TVcfToBeagle : protected TVcfConverter {
private:
    TOutputFileZipped * beagleFile;
    // beagle
    void writeBeagleHeader();
    void writeData(double * data) override ;

public:
    TVcfToBeagle(TParameters &Params, TLog *Logfile);
    void vcfToBeagle(TParameters & Params);
};


#endif //ATLAS_TVCFCONVERTER_H
