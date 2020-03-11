//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include "../stringFunctions.h"
#include "../TLog.h"
#include "../TParameters.h"
#include "TVcfFile.h"
#include "TPopulationLikelihoods.h"
#include "TGenotypeFrequencies.h"
#include "TGLF.h"

class TVcfConverter {
protected:
    std::string _outname;
    TLog * logfile;
    TGlfConverter glfConverter;
    TPopulationLikelihoodReader * reader;
    TPopulationSamples samples;

public:
    TVcfConverter(TLog * Logfile, TParameters & Params);
    ~TVcfConverter();

    void readOutputName(TParameters & Params);
    void readVcfAndWriteFile(TParameters & Params);
    virtual void writeHeader();
    virtual void writeData(TPopulationLikehoodLocus & data);
};

class TVcfToBeagle : protected TVcfConverter {
private:
    TOutputFileZipped * beagleFile;
    std::map<char, int> baseToNumber;
    // beagle
    void writeHeader() override;
    void writeRefAndAlt();
    void writeData(TPopulationLikehoodLocus & data) override;
    void writePosition();

public:
    TVcfToBeagle(TParameters &Params, TLog *Logfile);
    ~TVcfToBeagle();
    void vcfToBeagle(TParameters & Params);
};

class TVcfToLFMM : protected TVcfConverter {
protected:
    TOutputFilePlain * lfmmFile;
    TOutputFilePlain * lociNamesFile;

    void writeHeader() override;
    void storeLocusNames();
    void writeLociNames();
    void writeLFMM();

    std::vector<float *> genotypes;
    std::vector<std::string> loci_names;

public:
    TVcfToLFMM(TLog *Logfile, TParameters &Params);
    ~TVcfToLFMM();
    void vcfToLFMM(TParameters & Params);
};

class TVcfToLFMMCalledGeno : public TVcfToLFMM {
private:
    void writeData(TPopulationLikehoodLocus & data) override ;
    void storeCalledGenotypes();
public:
    TVcfToLFMMCalledGeno(TParameters &Params, TLog *Logfile);
};

class TVcfToLFMMPostGeno : public TVcfToLFMM {
private:
    // lfmm
    void writeData(TPopulationLikehoodLocus & data) override ;
    void storePosteriorGenotypes(TPopulationLikehoodLocus & data);
    float computePosteriorGenotype(TPopulationLikehoodLocus & data, int i);

public:
    TVcfToLFMMPostGeno(TParameters &Params, TLog *Logfile);
};


#endif //ATLAS_TVCFCONVERTER_H
