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
#include "TVCFReader.h"

namespace VCF{

//------------------------------------------
// TVcfConverter
//------------------------------------------
class TVcfConverter {
protected:
    std::string _outname;
	coretools::TLog * _logfile;
    genometools::TPopulationLikelihoodReaderLocus _reader;
    genometools::TPopulationSamples _samples;

    virtual void _initOutputFiles() = 0;
    virtual void _writeHeader();
    virtual void _writeData(genometools::TPopulationLikehoodLocus & data) = 0;
	void _readOutputName(coretools::TParameters & Params, const std::string & VCFFilename);

public:
    TVcfConverter(coretools::TLog * Logfile);
    virtual ~TVcfConverter() = default;

    void readVcfAndWriteFile(coretools::TParameters & Params);
};

//------------------------------------------
// TVcfToBeagle
//------------------------------------------
class TVcfToBeagle : protected TVcfConverter {
private:
	coretools::TOutputFile _beagleFile;

    // beagle
    void _writeHeader() override;
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void _initOutputFiles() override;

    void _writeRefAndAlt();
    void _writePosition();

public:
    TVcfToBeagle(coretools::TLog *Logfile);
    ~TVcfToBeagle() override = default;

    void vcfToBeagle(coretools::TParameters & Params);
};

//------------------------------------------
// TVcfToGeno
//------------------------------------------
class TVcfToGeno : protected TVcfConverter {
private:
	coretools::TOutputFile _genoFile;
    coretools::TOutputFile _lociNamesFile;

    // geno
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void _initOutputFiles() override;
    void _closeOutputFiles();
    void _writePosition();

public:
    TVcfToGeno(coretools::TLog *Logfile);
    ~TVcfToGeno() override = default;

    void vcfToGeno(coretools::TParameters & Params);
};

//------------------------------------------
// TVcfToLFMM
//------------------------------------------
class TVcfToLFMM : protected TVcfConverter {
protected:
    coretools::TOutputFile _lfmmFile;
    coretools::TOutputFile _lociNamesFile;
    std::vector<std::string> _loci_names;

    void _storeLocusNames();
    void _writeLociNames();
    void _initOutputFiles() override;
    void _closeOutputFiles();

    template <class T>
    void _writeLFMM(T Genotypes) {
        size_t numLoci = Genotypes.size();
        for (size_t i = 0; i < _samples.numSamples(); i++){
            for (size_t l = 0; l < numLoci; l++){
                _lfmmFile << (double) Genotypes[l][i];
            }
            _lfmmFile.endLine();
        }
    }

    void _prepareAndReadVcf(coretools::TParameters & Params);

public:
    TVcfToLFMM(coretools::TLog *Logfile);
    ~TVcfToLFMM() override;
};

//------------------------------------------
// TVcfToLFMMCalledGeno
//------------------------------------------
class TVcfToLFMMCalledGeno : public TVcfToLFMM {
private:
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override ;
    void storeCalledGenotypes();
    std::vector<uint8_t *> _genotypes;

public:
    TVcfToLFMMCalledGeno(coretools::TLog *Logfile);
    ~TVcfToLFMMCalledGeno() override;
    void vcfToLFMM(coretools::TParameters & Params);
};

//------------------------------------------
// TVcfToLFMMPostGeno
//------------------------------------------
class TVcfToLFMMPostGeno : public TVcfToLFMM {
private:
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override ;
    void _storePosteriorGenotypes(PopulationTools::TPopulationLikehoodLocus & data);
    double _computePosteriorGenotype(PopulationTools::TPopulationLikehoodLocus & data, uint32_t i);
    std::vector<double *> _genotypes;

public:
    TVcfToLFMMPostGeno(coretools::TLog *Logfile);
    ~TVcfToLFMMPostGeno() override;
    void vcfToLFMM(coretools::TParameters & Params);
};

//------------------------------------------
// TVcfToPosFile
//------------------------------------------
class TVcfToPosFile : public TVcfConverter {
private:
    coretools::TOutputFile _posFile;
    // beagle
    void _writeHeader() override;
    void _writeRefAndAlt();
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void _writePosition();
    void _initOutputFiles() override;

public:
    TVcfToPosFile(coretools::TLog *Logfile);
    ~TVcfToPosFile() override = default;
    void vcfToPosFile(coretools::TParameters & Params);
};

//------------------------------------------
// TVcfToGenotypeTruthSetFile
//------------------------------------------
class TVcfToGenotypeTruthSetFile : public TVcfConverter {
private:
    BAM::TBed ** _bedFiles;
    coretools::TOutputFile _genFile;

    int _minDistanceToPreviousLocus;
    long _positionPreviousLocus;
    int _numSamplesPerLocus;
    std::string _curChr;

    void _writeHeader() override;
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void _filterIndividuals(PopulationTools::TPopulationLikehoodLocus & data);
    void _mapIndividualsToDepth(std::vector<uint32_t> & samplesToKeep);
    void _filterIndividualsWithHighestDepth(std::vector<uint32_t> & samplesToKeep, const std::map< double, std::vector<uint32_t>, std::greater<> > & depthVsSampleIndexMap);
    void _writeToGenFile(const std::vector<uint32_t> & samplesToKeep);
    void _storeInBedFile(const std::vector<uint32_t> & samplesToKeep);
    void _initOutputFiles() override;
    void _resetDistance();

public:
    TVcfToGenotypeTruthSetFile(coretools::TLog *Logfile);
    ~TVcfToGenotypeTruthSetFile() override;
    void vcfToGenotypeTruthSetFile(coretools::TParameters & Params);
};

//------------------------------------------
// TVcfToVcf
//------------------------------------------
//TODO: finish VCfToVCf converter for filtering
class TVcfToVcf: public TVcfConverter {
private:
	void _writeHeader() override;
	void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;

public:
	TVcfToVcf(coretools::TLog *Logfile);
	~TVcfToVcf() override = default;
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_VcfConverter:public coretools::TTask{
public:
	TTask_VcfConverter(){ _explanation = "Converting a VCF file to other formats"; };

	void run(){
		using namespace coretools::instances;
		std::string format = parameters().getParameter<std::string>("format");

		if(format == "beagle"){
			logfile().startIndent("Converting a VCF to Beagle format (parameter 'format'):");
			TVcfToBeagle VcfToBeagle(&logfile());
			VcfToBeagle.vcfToBeagle(parameters());
		} else if (format == "geno"){
			logfile().startIndent("Converting a VCF to geno format (parameter 'format'):");
			TVcfToGeno vcfToGeno(&logfile());
		    vcfToGeno.vcfToGeno(parameters());
		} else if(format == "LFMM"){
			logfile().startIndent("Converting a VCF to LFMM format (parameter 'format'):");

            //posterior or calls?
            std::string genoType = parameters().getParameterWithDefault<std::string>("genotypes", "calls");
            if(genoType == "posterior"){
		    TVcfToLFMMPostGeno vcfToLFMMPostGeno(&logfile());
                vcfToLFMMPostGeno.vcfToLFMM(parameters());
            } else if(genoType == "calls"){
		    TVcfToLFMMCalledGeno vcfToLFMMCalledGeno(&logfile());
                vcfToLFMMCalledGeno.vcfToLFMM(parameters());
            } else {
                throw "Unknown genotype method '" + genoType + "'! Use either 'calls' or 'posterior'";
            }
		} else if(format == "posFile"){
			logfile().startIndent("Converting a VCF file to posfile format used by STITCH (parameter 'format'):");
			TVcfToPosFile VcfToPosFile(&logfile());
			VcfToPosFile.vcfToPosFile(parameters());
		} else if(format == "truthSet"){
			logfile().startIndent("Converting a VCF file to genotype truth set format (parameter 'format'):");
			TVcfToGenotypeTruthSetFile VcfToGenotypeTruthSetFile(&logfile());
			VcfToGenotypeTruthSetFile.vcfToGenotypeTruthSetFile(parameters());
		} else {
			throw "Unknown format '" + format + "'! Use either 'beagle', 'geno', 'LFMM', 'posFile' or 'truthSet'.";
		}
		logfile().endIndent();
	};

};

}; //end namespace

#endif //ATLAS_TVCFCONVERTER_H
