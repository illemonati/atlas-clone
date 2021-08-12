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
#include <TPopulationLikelihoodLocus.h>
#include "stringFunctions.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TGenotypeFrequencies.h"
#include "TBed.h"
#include "TVcfFile.h"

namespace PopulationTools{

//------------------------------------------------
// TPopulation
//------------------------------------------------
class TPopulation{
private:
	std::string _name;
	std::vector<std::string> _samples;
	uint32_t _firstSampleIndex;

public:
	TPopulation(const std::string & Name);

	std::string name() const{
		return _name;
	};

	bool operator==(const std::string & Name) const{
		return Name == _name;
	};

	uint32_t numSamples() const {
		return _samples.size();
	};

	void addSample(const std::string & Sample);
	bool sampleExists(const std::string & Sample) const;

	void setFirstSampleIndex(const uint32_t & Index){
		_firstSampleIndex = Index;
	};

	uint32_t firstSampleIndex() const{
		return _firstSampleIndex;
	};

	bool sampleIndexExists(const uint32_t & Index) const;
	uint32_t sampleIndex(const std::string & Sample) const;
	std::string sampleName(const uint32_t & Index) const;
	void addSampleNamesToVector(std::vector<std::string> & vec) const;
	void report(coretools::TLog* Logfile) const;

	//the following functions accept arrays spanning ALL samples and perform calculations on samples in population
	uint32_t numSamplesMissing(bool* sampleMissing) const;

	uint32_t numSamplesWithData(bool* sampleMissing) const {
		return numSamples() - numSamplesMissing(sampleMissing);
	};
};

//------------------------------------------------
//TPopulationSamples
//------------------------------------------------
class TPopulationSamples{
private:
	//populations
	std::vector<TPopulation> _populations;
	uint32_t _numSamples;

	//VCF index: maps index in VCF to internal index, which is ordered by population
	struct vcfInfo{
		bool used;
		uint32_t index;

		vcfInfo(){
			used = false;
			index = 0;
		}
	};
	std::vector<vcfInfo> _vcfIndex;
	std::vector<uint32_t> _indexToVCFIndex;
	std::vector<uint32_t> _indexToPopulationIndex;

	void _init();
	void _fillIndexToPopulationIndex();

public:
	TPopulationSamples();
	TPopulationSamples(std::string filename, coretools::TLog* logfile);
	~TPopulationSamples() = default;

	bool hasSamples() const { return _numSamples > 0; };

	size_t numPopulations() const { return _populations.size(); };
	bool populationExists(const std::string & name) const;
	std::string getPopulationName(const uint32_t & index) const;
	uint32_t populationIndex(const std::string & name) const;
	uint32_t populationIndex(const uint32_t & SampleIndex) const { return _indexToPopulationIndex[SampleIndex]; };
	uint32_t numSamplesInPop(const uint32_t & population) const { return _populations[population].numSamples(); };

	uint32_t numSamples() const { return _numSamples; };
	bool sampleExists(const std::string & name) const;
	uint32_t sampleIndex(const std::string & name) const;
	std::string sampleName(const uint32_t & SampleIndex) const;
	void addSampleNamesToVector(std::vector<std::string> & vec) const;

	void readSamples(std::string filename, coretools::TLog* logfile);
	void readSamplesFromVCFNames(std::vector<std::string> & vcfSampleNames);
	void report(coretools::TLog* Logfile);
	void fillVCFOrder(std::vector<std::string> & vcfSampleNames);

	//bool sampleIsUsed(const std::string & name);
	//uint32_t getOrderedSampleIndex(const std::string & name);
	uint32_t startIndex(int population){ return _populations[population].firstSampleIndex(); };
	uint32_t sampleIndexInVCF(const uint32_t & index);

	uint8_t* getPointerToDataInPop(uint8_t* data, uint32_t population) const;
	uint32_t numSamplesMissingInPop(bool* sampleMissing, uint32_t population) const;
	uint32_t numSamplesWithDataInPop(bool* sampleMissing, uint32_t population) const;
};

//-------------------------------------------------
//TPopulationLikelihoodReader
//-------------------------------------------------
class TPopulationLikelihoodReader{
protected:
	bool _initialized;
	coretools::TLog* logfile;
	VCF::TVcfFileSingleLine vcfFile;
	bool vcfOpen;

	//BED file for windows
	BAM::TBed bedFile;
	std::set<BAM::TGenomeWindow>::iterator _curBedWindow;
	uint32_t _curRefId, _previousRefId;
	bool limitToSitesInBed;

	//settings
	bool limitLines;
	bool filterOnChr;
	uint64_t maxLinesToRead;
	uint32_t minDepth;
	uint32_t minNumSamplesWithData;
	double freqFilter;
	double epsilonF; //F for EM algorithm to estimate allele frequencies
	uint32_t minVariantQuality;
    std::vector<std::string> chromosomesToKeep;
	bool estimateGenotypeFrequencies;
	uint64_t progressFrequency;

	//counters
	coretools::TTimer timer;
	bool vcfParsingStarted;
	uint64_t _lineCounter; //lines read in VCF
	uint64_t _notInBedFile; //sites considered (smaller than # lines if BED file is used
	uint64_t _notBialleleicCounter;
	uint64_t _missingSNPCounter;
	uint64_t _lowFreqSNPCounter;
	uint64_t _lowVariantQualityCounter;
	uint64_t _noPLCounter;
	uint64_t _numAcceptedLoci;
    uint64_t _notOnChrCounter;

    //tmp variables used for reading
	genometools::TGenotypeFrequencies genoFrequencies;
	std::string curChr;

	virtual void _init();
	void resetCounters();
	void closeVCF();

    //void readDataFromVCF(TParameters & Parameters, TPopulationSamples & samples);
    void specifyChromosomesToKeep(coretools::TParameters & Parameters, coretools::TLog* logfile);
    void printProgressFrequencyFiltering();
    int filterOnDepth(TSampleLikelihoods* data, TPopulationSamples & samples);
    virtual bool _readNextLineFromVCF();
    bool _filterSite(TSampleLikelihoods* data, TPopulationSamples & samples);
    bool _updateChromosomeInfo();
    bool _jumpToNextChromosome();

public:
	TPopulationLikelihoodReader();
	TPopulationLikelihoodReader(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq);
	virtual ~TPopulationLikelihoodReader();

	virtual void initialize(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq);
	void doEstimateGenotypeFrequencies(){ estimateGenotypeFrequencies = true; };


	void openVCF(std::string);
	void concludeFilters();
	std::vector<std::string>& getSampleVCFNames(){ return vcfFile.parser.samples; };
	long numLociParsed(){ return _lineCounter; };
	long numAcceptedLoci(){ return _numAcceptedLoci; };
	uint32_t getMinNumSamplesWithData(){return minNumSamplesWithData; };
	void writeUnknownHeader(){ vcfFile.writeHeaderVCF_4_0(); };
	int variantQuality(){ return vcfFile.variantQuality(); };
	void setOutStream(std::ostream & os){ vcfFile.setOutStream(os); };
	void writeVCFLine(){ vcfFile.writeLine(); };
};

class TPopulationLikelihoodReaderLocus:public TPopulationLikelihoodReader{
private:
	//true allele frequencies
	std::istream* trueFreq;
	bool trueFreqFileOpen = false;
	double _trueAlleleFrequency;

	//settings
	bool storeTrueAlleleFreq;
	std::string trueAlleleFreqFile;

	void _init();
	void _closeTrueAlleleFreqFile();
	bool _readNextLineFromVCF();

public:
	TPopulationLikelihoodReaderLocus();
	TPopulationLikelihoodReaderLocus(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq);
	~TPopulationLikelihoodReaderLocus();

	virtual void initialize(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq);
	void doSaveTrueAlleleFrequencies(){ storeTrueAlleleFreq = true; }

	std::string chr(){ return curChr; };
	long position(){ return vcfFile.position(); };
    long positionZeroBased(){ return vcfFile.positionZeroBased(); };
    char refAllele(){ return vcfFile.getRefAllele()[0]; };
    char altAllele(){ return vcfFile.getFirstAltAllele()[0]; };
    std::vector<genometools::BiallelicGenotype> biallelicGenotypes(TPopulationSamples & samples) const;
    genometools::BiallelicGenotype biallelicGenotype(TPopulationSamples & samples, const uint32_t & s) const;
    genometools::Genotype genotype(TPopulationSamples & samples, const uint32_t & s) const;
    double depth(TPopulationSamples & samples,uint32_t s);

    bool readDataFromVCF(TPopulationLikehoodLocus & data, TPopulationSamples & samples);
	bool readDataFromVCF(TSampleLikelihoods* data, TPopulationSamples & samples);

	void openTrueAlleleFrequenciesFile(const std::string filename);
	genometools::TGenotypeFrequencies* genotypeFrequencies(){ return &genoFrequencies; };
	std::array<Probability, 3> diploidGenotypeFrequencies(){ return genoFrequencies.diploidFrequencies(); };
	double allelFrequency(){ return genoFrequencies.alleleFrequency(); };
	double trueAlleleFrequency(){ return _trueAlleleFrequency; };
	double MAF(){ return genoFrequencies.MAF(); };
	int numSamplesWithData();
	void writePosition(coretools::TOutputFile & out);
};

class TPopulationLikelihoodReaderWindow:public TPopulationLikelihoodReader{
private:

	bool _readNextLocusAndUpdateChromosome();

public:
	TPopulationLikelihoodReaderWindow();
	TPopulationLikelihoodReaderWindow(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq);
	~TPopulationLikelihoodReaderWindow();

	bool readDataFromVCF(TPopulationLikehoodWindow & region, TPopulationSamples & samples);
	void writeWindow(coretools::TOutputFile & out);
};

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
	TPopulationSamples samples;

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
	TPopulationSamples samples;

	//data on loci
	long _numLoci;
	std::map<int, std::string> chromosomes; //first SNP index and name
	std::vector<long> position;
    std::vector<TSampleLikelihoods*> data;
    std::vector<double> alleleFrequencies;
    std::vector<double> trueAlleleFrequencies;
    bool saveAlleleFrequencies;
    bool saveTrueAlleleFrequencies;

    //looping
    long curLocusIndex;
    std::map<int, std::string>::iterator curChrIt;
    std::map<int, std::string>::iterator nextChrIt;
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

    // get main constants (n, L, D, K) and names of environmental variables
    int getNumIndividuals();
    long getNumLoci();
    TSampleLikelihoods* getDataAtLocus(long index);
    void print();
};

}; //end namespace

#endif /* TPOPULATIONLIKELIHOODS_H_ */
