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
    TLog * _logfile;
    PopulationTools::TPopulationLikelihoodReaderLocus _reader;
    PopulationTools::TPopulationSamples _samples;

    virtual void _initOutputFiles() = 0;
    virtual void _writeHeader();
    virtual void _writeData(PopulationTools::TPopulationLikehoodLocus & data) = 0;
    void _readOutputName(TParameters & Params, const std::string & VCFFilename);

public:
    TVcfConverter(TLog * Logfile);
    virtual ~TVcfConverter() = default;

    void readVcfAndWriteFile(TParameters & Params);
};

//------------------------------------------
// TVcfToBeagle
//------------------------------------------
class TVcfToBeagle : protected TVcfConverter {
private:
    TOutputFile _beagleFile;

    // beagle
    void _writeHeader() override;
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void _initOutputFiles() override;

    void _writeRefAndAlt();
    void _writePosition();

public:
    TVcfToBeagle(TLog *Logfile);
    ~TVcfToBeagle() override = default;

    void vcfToBeagle(TParameters & Params);
};

//------------------------------------------
// TVcfToGeno
//------------------------------------------
class TVcfToGeno : protected TVcfConverter {
private:
    TOutputFile _genoFile;
    TOutputFile _lociNamesFile;

    // geno
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void _initOutputFiles() override;
    void _closeOutputFiles();
    void _writePosition();

public:
    TVcfToGeno(TLog *Logfile);
    ~TVcfToGeno() override = default;

    void vcfToGeno(TParameters & Params);
};

//------------------------------------------
// TVcfToLFMM
//------------------------------------------
class TVcfToLFMM : protected TVcfConverter {
protected:
    TOutputFile _lfmmFile;
    TOutputFile _lociNamesFile;
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

    void _prepareAndReadVcf(TParameters & Params);

public:
    TVcfToLFMM(TLog *Logfile);
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
    TVcfToLFMMCalledGeno(TLog *Logfile);
    ~TVcfToLFMMCalledGeno() override;
    void vcfToLFMM(TParameters & Params);
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
    TVcfToLFMMPostGeno(TLog *Logfile);
    ~TVcfToLFMMPostGeno() override;
    void vcfToLFMM(TParameters & Params);
};

//------------------------------------------
// TVcfToPosFile
//------------------------------------------
class TVcfToPosFile : public TVcfConverter {
private:
    TOutputFile _posFile;
    // beagle
    void _writeHeader() override;
    void _writeRefAndAlt();
    void _writeData(PopulationTools::TPopulationLikehoodLocus & data) override;
    void _writePosition();
    void _initOutputFiles() override;

public:
    TVcfToPosFile(TLog *Logfile);
    ~TVcfToPosFile() override = default;
    void vcfToPosFile(TParameters & Params);
};

//------------------------------------------
// TVcfToGenotypeTruthSetFile
//------------------------------------------
class TVcfToGenotypeTruthSetFile : public TVcfConverter {
private:
    BAM::TBed ** _bedFiles;
    TOutputFile _genFile;

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
    TVcfToGenotypeTruthSetFile(TLog *Logfile);
    ~TVcfToGenotypeTruthSetFile() override;
    void vcfToGenotypeTruthSetFile(TParameters & Params);
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
	TVcfToVcf(TLog *Logfile);
	~TVcfToVcf() override = default;
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_VcfConverter:public coretools::TTask{
public:
	TTask_VcfConverter(){ _explanation = "Converting a VCF file to other formats"; };

	void run(TParameters & Parameters, TLog* Logfile){
		std::string format = Parameters.getParameter<std::string>("format");

		if(format == "beagle"){
			Logfile->startIndent("Converting a VCF to Beagle format (parameter 'format'):");
			TVcfToBeagle VcfToBeagle(Logfile);
			VcfToBeagle.vcfToBeagle(Parameters);
		} else if (format == "geno"){
		    Logfile->startIndent("Converting a VCF to geno format (parameter 'format'):");
		    TVcfToGeno vcfToGeno(Logfile);
		    vcfToGeno.vcfToGeno(Parameters);
		} else if(format == "LFMM"){
            Logfile->startIndent("Converting a VCF to LFMM format (parameter 'format'):");

            //posterior or calls?
            std::string genoType = Parameters.getParameterWithDefault<std::string>("_genotypes", "calls");
            if(genoType == "posterior"){
                TVcfToLFMMPostGeno vcfToLFMMPostGeno(Logfile);
                vcfToLFMMPostGeno.vcfToLFMM(Parameters);
            } else if(genoType == "calls"){
                TVcfToLFMMCalledGeno vcfToLFMMCalledGeno(Logfile);
                vcfToLFMMCalledGeno.vcfToLFMM(Parameters);
            } else {
                throw "Unknown genotype method '" + genoType + "'! Use either 'calls' or 'posterior'";
            }
		} else if(format == "_posFile"){
			Logfile->startIndent("Converting a VCF file to posfile format used by STITCH (parameter 'format'):");
			TVcfToPosFile VcfToPosFile(Logfile);
			VcfToPosFile.vcfToPosFile(Parameters);
		} else if(format == "truthSet"){
			Logfile->startIndent("Converting a VCF file to genotype truth set format (parameter 'format'):");
			TVcfToGenotypeTruthSetFile VcfToGenotypeTruthSetFile(Logfile);
			VcfToGenotypeTruthSetFile.vcfToGenotypeTruthSetFile(Parameters);
		} else {
			throw "Unknown format '" + format + "'! Use either 'beagle', 'LFMM', '_posFile' or 'truthSet'.";
		}
		Logfile->endIndent();
	};

};

}; //end namespace

#endif //ATLAS_TVCFCONVERTER_H
