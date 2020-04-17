//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include "stringFunctions.h"
#include "TLog.h"
#include "TParameters.h"
#include "TVcfFile.h"
#include "../PopulationTools/TPopulationLikelihoods.h"
#include "../PopulationTools/TGenotypeFrequencies.h"
#include "../GLF/TGLF.h"

class TVcfConverter {
protected:
    std::string _outname;
    TLog * logfile;
    TGlfConverter glfConverter;
    TPopulationLikelihoodReaderLocus * reader;
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
    std::vector<std::string> loci_names;

    void writeHeader() override;
    void storeLocusNames();
    void writeLociNames();

    template <class T>
    void writeLFMM(T genotypes) {
        int numLoci = genotypes.size();
        for (uint32_t i = 0; i < samples.numSamples(); i++){
            for (int l = 0; l < numLoci; l++){
                *(lfmmFile) << static_cast<float>(genotypes[l][i]);
            }
            lfmmFile->endLine();
        }
    }

    void prepareAndReadVcf(TParameters & Params);

public:
    TVcfToLFMM(TLog *Logfile, TParameters &Params);
    ~TVcfToLFMM();
};

class TVcfToLFMMCalledGeno : public TVcfToLFMM {
private:
    void writeData(TPopulationLikehoodLocus & data) override ;
    void storeCalledGenotypes();
    std::vector<uint8_t *> genotypes;

public:
    TVcfToLFMMCalledGeno(TParameters &Params, TLog *Logfile);
    ~TVcfToLFMMCalledGeno();
    void vcfToLFMM(TParameters & Params);
};

class TVcfToLFMMPostGeno : public TVcfToLFMM {
private:
    void writeData(TPopulationLikehoodLocus & data) override ;
    void storePosteriorGenotypes(TPopulationLikehoodLocus & data);
    float computePosteriorGenotype(TPopulationLikehoodLocus & data, uint32_t i);
    std::vector<float *> genotypes;

public:
    TVcfToLFMMPostGeno(TParameters &Params, TLog *Logfile);
    ~TVcfToLFMMPostGeno();
    void vcfToLFMM(TParameters & Params);

};

class TVcfToPosFile : public TVcfConverter {
private:
    TOutputFilePlain * posFile;
    // beagle
    void writeHeader() override;
    void writeRefAndAlt();
    void writeData(TPopulationLikehoodLocus & data) override;
    void writePosition();

public:
    TVcfToPosFile(TParameters &Params, TLog *Logfile);
    ~TVcfToPosFile();
    void vcfToPosFile(TParameters & Params);
};

class TVcfToBedFile : public TVcfConverter {
private:
    TOutputFilePlain * bedFile;
    void writeHeader() override;
    void writeData(TPopulationLikehoodLocus & data) override;
    void writePosition();

public:
    TVcfToBedFile(TParameters &Params, TLog *Logfile);
    ~TVcfToBedFile();
    void vcfToBedFile(TParameters & Params);
};

#endif //ATLAS_TVCFCONVERTER_H
