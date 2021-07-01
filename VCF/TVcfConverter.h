//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include "stringFunctions.h"
#include "TLog.h"
#include "TParameters.h"
#include "TPopulationLikelihoods.h"
#include "TGenotypeFrequencies.h"
#include "TGLF.h"
#include "TBed.h"
#include "TVcfFile.h"

namespace VCF{

using coretools::TParameters;
using coretools::TLog;
using coretools::TOutputFile;

//------------------------------------------
// TVcfConverter
//------------------------------------------
class TVcfConverter {
protected:
    std::string _outname;
    TLog * logfile;
    PopulationTools::TPopulationLikelihoodReaderLocus * reader;
    PopulationTools::TPopulationSamples samples;

    virtual void initOutputFiles();
    virtual void writeHeader();
    virtual void writeData(PopulationTools::TPopulationLikehoodLocus & data);
    void readOutputName(TParameters & Params);
public:
    TVcfConverter(TLog * Logfile, TParameters & Params);
    virtual ~TVcfConverter();
    void readVcfAndWriteFile(TParameters & Params);
};

//------------------------------------------
// TVcfToBeagle
//------------------------------------------
class TVcfToBeagle : protected TVcfConverter {
private:
    TOutputFile * beagleFile;
    std::map<char, int> baseToNumber;
    // beagle
    void writeHeader() override;
    void writeRefAndAlt();
    void writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void writePosition();
    void initOutputFiles() override;

public:
    TVcfToBeagle(TParameters &Params, TLog *Logfile);
    ~TVcfToBeagle();
    void vcfToBeagle(TParameters & Params);
};

//------------------------------------------
// TVcfToLFMM
//------------------------------------------
class TVcfToLFMM : protected TVcfConverter {
protected:
    TOutputFile * lfmmFile;
    TOutputFile * lociNamesFile;
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

//------------------------------------------
// TVcfToLFMMCalledGeno
//------------------------------------------
class TVcfToLFMMCalledGeno : public TVcfToLFMM {
private:
    void writeData(PopulationTools::TPopulationLikehoodLocus & data) override ;
    void storeCalledGenotypes();
    std::vector<uint8_t *> genotypes;

public:
    TVcfToLFMMCalledGeno(TParameters &Params, TLog *Logfile);
    ~TVcfToLFMMCalledGeno();
    void vcfToLFMM(TParameters & Params);
};

//------------------------------------------
// TVcfToLFMMPostGeno
//------------------------------------------
class TVcfToLFMMPostGeno : public TVcfToLFMM {
private:
    void writeData(PopulationTools::TPopulationLikehoodLocus & data) override ;
    void storePosteriorGenotypes(PopulationTools::TPopulationLikehoodLocus & data);
    float computePosteriorGenotype(PopulationTools::TPopulationLikehoodLocus & data, uint32_t i);
    std::vector<float *> genotypes;

public:
    TVcfToLFMMPostGeno(TParameters &Params, TLog *Logfile);
    ~TVcfToLFMMPostGeno();
    void vcfToLFMM(TParameters & Params);

};

//------------------------------------------
// TVcfToPosFile
//------------------------------------------
class TVcfToPosFile : public TVcfConverter {
private:
    TOutputFile * posFile;
    // beagle
    void writeHeader() override;
    void writeRefAndAlt();
    void writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void writePosition();
    void initOutputFiles() override;

public:
    TVcfToPosFile(TParameters &Params, TLog *Logfile);
    ~TVcfToPosFile();
    void vcfToPosFile(TParameters & Params);
};

//------------------------------------------
// TVcfToGenotypeTruthSetFile
//------------------------------------------
class TVcfToGenotypeTruthSetFile : public TVcfConverter {
private:
    BAM::TBed ** bedFiles;
    TOutputFile * genFile;

    int minDistanceToPreviousLocus;
    long positionPreviousLocus;
    int numSamplesPerLocus;
    std::string curChr;

    void writeHeader() override;
    void writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void filterIndividuals(PopulationTools::TPopulationLikehoodLocus & data);
    void mapIndividualsToDepth(std::vector<uint32_t> & samplesToKeep);
    void filterIndividualsWithHighestDepth(std::vector<uint32_t> & samplesToKeep, const std::map< double, std::vector<uint32_t>, std::greater<double> > & depthVsSampleIndexMap);
    void writeToGenFile(const std::vector<uint32_t> & samplesToKeep);
    void storeInBedFile(const std::vector<uint32_t> & samplesToKeep);
    void initOutputFiles() override;
    void resetDistance();

public:
    TVcfToGenotypeTruthSetFile(TParameters &Params, TLog *Logfile);
    ~TVcfToGenotypeTruthSetFile();
    void vcfToGenotypeTruthSetFile(TParameters & Params);
};

//------------------------------------------
// TVcfToVcf
//------------------------------------------
//TODO: finish VCfToVCf converter for filtering
class TVcfToVcf: public TVcfConverter {
private:
	void writeRefAndAlt();
	void writeHeader() override;
	void writeData(PopulationTools::TPopulationLikehoodLocus & data) override;

public:
	TVcfToVcf(TParameters &Params, TLog *Logfile);
	virtual ~TVcfToVcf(){};
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_VcfConverter:public TTask{
public:
	TTask_VcfConverter(){ _explanation = "Converting a VCF file to other formats"; };

	void run(TParameters & Parameters, TLog* Logfile){
		std::string format = Parameters.getParameter<std::string>("format");

		if(format == "beagle"){
			Logfile->startIndent("Converting a VCF to Beagle format (parameter 'format'):");
			TVcfToBeagle VcfToBeagle(Parameters, Logfile);
			VcfToBeagle.vcfToBeagle(Parameters);
		} else if(format == "LFMM"){
			Logfile->startIndent("Converting a VCF to LFMM format (parameter 'format'):");

			//posterior or calls?
			std::string genoType = Parameters.getParameterWithDefault<std::string>("genotypes", "calls");
			if(genoType == "posterior"){
				TVcfToLFMMPostGeno vcfToLFMMPostGeno(Parameters, Logfile);
				vcfToLFMMPostGeno.vcfToLFMM(Parameters);
			} else if(genoType == "calls"){
				TVcfToLFMMCalledGeno vcfToLFMMCalledGeno(Parameters, Logfile);
				vcfToLFMMCalledGeno.vcfToLFMM(Parameters);
			} else {
				throw "Unknown genotype method '" + genoType + "'! Use either 'calls' or 'posterior'";
			}
		} else if(format == "posFile"){
			Logfile->startIndent("Converting a VCF file to posfile format used by STITCH (parameter 'format'):");
			TVcfToPosFile VcfToPosFile(Parameters, Logfile);
			VcfToPosFile.vcfToPosFile(Parameters);
		} else if(format == "truthSet"){
			Logfile->startIndent("Converting a VCF file to genotype truth set format (parameter 'format'):");
			TVcfToGenotypeTruthSetFile VcfToGenotypeTruthSetFile(Parameters, Logfile);
			VcfToGenotypeTruthSetFile.vcfToGenotypeTruthSetFile(Parameters);
		} else {
			throw "Unknown format '" + format + "'! Use either 'beagle', 'LFMM', 'posFile' or 'truthSet'.";
		}
		Logfile->endIndent();
	};

};

}; //end namespace

#endif //ATLAS_TVCFCONVERTER_H
