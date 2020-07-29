//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include <memory>
#include <vector>
#include "stringFunctions.h"
#include "TLog.h"
#include "TParameters.h"
#include "TVcfFile.h"
#include "../PopulationTools/TPopulationLikelihoods.h"
#include "../PopulationTools/TGenotypeFrequencies.h"
#include "../GLF/TGLF.h"
#include "../TBed.h"

class TVcfConverter {
protected:
    std::string _outname;
    TLog * logfile;
    TGlfConverter glfConverter;
    TPopulationLikelihoodReaderLocus * reader;
    TPopulationSamples samples;

    virtual void initOutputFiles();
    virtual void writeHeader();
    virtual void writeData(TPopulationLikehoodLocus & data);
    void readOutputName(TParameters & Params);
public:
    TVcfConverter(TLog * Logfile, TParameters & Params);
    ~TVcfConverter();
    void readVcfAndWriteFile(TParameters & Params);
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
    void initOutputFiles() override;

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
    void initOutputFiles() override;

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
    float calculateMeanOfMeanPosteriorGenotypes(TPopulationLikehoodLocus & data, const float * meanPostGenoForOneLocus);
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
    void initOutputFiles() override;

public:
    TVcfToPosFile(TParameters &Params, TLog *Logfile);
    ~TVcfToPosFile();
    void vcfToPosFile(TParameters & Params);
};

class TVcfToGenotypeTruthSetFile : public TVcfConverter {
private:
    TBed ** bedFiles;
    TOutputFilePlain * genFile;

    // posfile?
    std::vector<std::string> positionsInPosFile;

    int minDistanceToPreviousLocus;
    long positionPreviousLocus;
    int numSamplesPerLocus;
    std::string curChr;

    void writeHeader() override;
    void writeData(TPopulationLikehoodLocus & data) override;
    void filterIndividuals(TPopulationLikehoodLocus & data);
    void mapIndividualsToDepth(std::vector<uint32_t> & samplesToKeep);
    void filterIndividualsWithHighestDepth(std::vector<uint32_t> & samplesToKeep, const std::map< double, std::vector<uint32_t>, std::greater<double> > & depthVsSampleIndexMap);
    void writeToGenFile(const std::vector<uint32_t> & samplesToKeep);
    void storeInBedFile(const std::vector<uint32_t> & samplesToKeep);
    void initOutputFiles() override;
    void resetDistance();
    void readPosfileIntoMemory(std::string fileNamePosfile);

public:
    TVcfToGenotypeTruthSetFile(TParameters &Params, TLog *Logfile);
    ~TVcfToGenotypeTruthSetFile();
    void vcfToGenotypeTruthSetFile(TParameters & Params);
};

class TVcfToVcf: public TVcfConverter {
private:
	void writeHeader() override;
	void writeData(TPopulationLikehoodLocus & data) override;
    void initOutputFiles() override;

public:
	TVcfToVcf(TParameters &Params, TLog *Logfile);
	~TVcfToVcf()= default;

	void vcfToVcf(TParameters & Params);
};

class TVcfExtract : protected TVcfConverter {
private:
    TOutputFilePlain * file;

    enum Types {
        genotypes,
        depth,
        refAndAlt
    };
    Types whatToExtract;

    // write
    void writeHeader() override;
    void _writeRefAndAlt();
    void _writeDepth(uint32_t s);
    void _writeGenotype(uint32_t s);
    void writeData(TPopulationLikehoodLocus & data) override;
    void writePosition();
    void initOutputFiles() override;
    void _fillHeaderIndividuals(std::vector <std::string> & header);

public:
    TVcfExtract(TParameters &Params, TLog *Logfile);
    ~TVcfExtract();
    void vcfExtract(TParameters & Params);
};

class TStitchVcfReader {
private:
    gz::igzstream * _vcf;
    bool _vcfOpen;
    std::vector<std::string> _oneLine;
    uint32_t _numSamples;

public:
    TStitchVcfReader();
    ~TStitchVcfReader();

    void openVcf(const std::string& vcfFilename);
    bool read();
    void parseVCFHeader(std::vector<std::string>& sampleNames);
    void genotypes(std::vector<std::string>& genotypes);
    void posteriorGenotypes(std::vector<std::string>& posteriorGenotypes);
    void maxPosteriorGenotypes(std::vector<std::string>& maxPosteriorGenotypes);
    void meanPosteriorGenotypes(std::vector<std::string>& meanPosteriorGenotypes);
    std::string chr();
    std::string pos();
    std::string refAllele();
    std::string altAllele();
    double infoScore();
    double HWE_pVal();
    uint32_t numSamples() const{return _numSamples;};
};

class TStitchVcfToBeagle {
private:
    TOutputFileZipped beagleFile;
    TStitchVcfReader reader;
    std::map<std::string, int> baseToNumber;
    void parseVCF();
    void parseVCFHeader();

public:
    TStitchVcfToBeagle(TParameters & Params, TLog * logfile);
    ~TStitchVcfToBeagle(){};
};

class TStitchVcfToPosfile {
private:
    TOutputFilePlain posfile;
    TStitchVcfReader reader;
    void parseVCF();
    void parseVCFHeader();

    // filters
    double minInfoScore;
    double minHWE_pval;

public:
    TStitchVcfToPosfile(TParameters & Params, TLog * logfile);
    ~TStitchVcfToPosfile(){};
};

class TStitchVcfExtract {
    // extracts genotypes / posterior genotypes / maxPosteriorGenotypes from STITCH output vcf
private:
    enum Types {
        genotypes,
        posteriorGenotypes,
        maxPostGenotype
    };
    Types whatToExtract;

    TOutputFilePlain file;
    TStitchVcfReader reader;
    void parseVCF();
    void parseVCFHeader();
    void _writeGenotypes();
    void _writePosteriorGenotypes();
    void _writeMaxPosteriorGenotypes();

public:
    TStitchVcfExtract(TParameters & Params, TLog * logfile);
    ~TStitchVcfExtract(){};
};

class TStitchVcfToLFMMPostGeno {
    // converts a stitch output VCF to a LFMM with the mean posterior genotypes
private:
    std::vector<std::vector<std::string>> _genotypes;
    std::vector<std::string> loci_names;
    TOutputFilePlain file;
    TOutputFilePlain fileLociNames;
    TStitchVcfReader reader;
    void parseVCF();
    void parseVCFHeader();
    void _write();
    void _writeLociNames();
    void _writeMeanPosteriorGenotypes();

public:
    TStitchVcfToLFMMPostGeno(TParameters & Params, TLog * logfile);
    ~TStitchVcfToLFMMPostGeno(){};
};

#endif //ATLAS_TVCFCONVERTER_H
