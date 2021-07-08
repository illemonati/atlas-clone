/*
 * TPopulationLikelihoods.cpp
 *
 *  Created on: Dec 7, 2018
 *      Author: phaentu
 */

//
// Created by Madleina Caduff on 19.10.18.
//

#include "TPopulationLikelihoods.h"

namespace PopulationTools{

using coretools::str::toString;


//------------------------------------------------
// TPopulation
//------------------------------------------------
TPopulation::TPopulation(const std::string & Name){
	_name = Name;
	_firstSampleIndex = 0;
};

void TPopulation::addSample(const std::string & Sample){
	//check for duplicates
	if(std::find(_samples.begin(), _samples.end(), Sample) != _samples.end()){
		throw "Duplicate sample '" + Sample + "' in population '" + _name + "'!";
	}
	_samples.push_back(Sample);
};

bool TPopulation::sampleExists(const std::string & Sample) const {
	return std::find(_samples.begin(), _samples.end(), Sample) != _samples.end();
};

bool TPopulation::sampleIndexExists(const uint32_t & Index) const {
	return Index >= _firstSampleIndex && Index < _firstSampleIndex + _samples.size();
};

uint32_t TPopulation::sampleIndex(const std::string & Sample) const {
	for(uint32_t i = 0; i < _samples.size(); ++i){
		if(_samples[i] == Sample){
			return _firstSampleIndex + i;
		}
	}
	throw std::runtime_error("uint32_t TPopulation::sampleIndex(const std::string & Sample): sample '" + Sample + "' does not exist!");
};

std::string TPopulation::sampleName(const uint32_t & Index) const{
	return _samples[Index - _firstSampleIndex];
};

void TPopulation::addSampleNamesToVector(std::vector<std::string> & vec) const {
	vec.insert(vec.end(), _samples.begin(), _samples.end());
};

void TPopulation::report(coretools::TLog* Logfile) const {
	for(auto& s : _samples){
		Logfile->list(s);
	}
};

//the following functions accept arrays spanning ALL samples and perform calculations on samples in population
uint32_t TPopulation::numSamplesMissing(bool* sampleMissing) const {
	int numMissing = 0;
	for(uint32_t s = 0; s < numSamples(); ++s){
		if(sampleMissing[_firstSampleIndex + s])
			numMissing++;
	}
	return numMissing;
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationSamples                                                                         //
////////////////////////////////////////////////////////////////////////////////////////////////
TPopulationSamples::TPopulationSamples(){
	_init();
};

TPopulationSamples::TPopulationSamples(std::string filename, coretools::TLog* logfile){
	_init();
	readSamples(filename, logfile);
};

void TPopulationSamples::_init(){
	_numSamples = 0;
};

bool TPopulationSamples::populationExists(const std::string & name) const {
	return std::find(_populations.begin(), _populations.end(), name) != _populations.end();
};

std::string TPopulationSamples::getPopulationName(const uint32_t & index) const {
	if(index >= _populations.size()){
		throw "No population with index " + toString(index) + "!";
	}
	return _populations[index].name();
};

uint32_t TPopulationSamples::populationIndex(const std::string & name) const {
	for(uint32_t p = 0; p < _populations.size(); ++p){
		if(_populations[p] == name){
			return p;
		}
	}
	throw "No population with name " + name + "!";
};

bool TPopulationSamples::sampleExists(const std::string & name) const {
	for(auto& p : _populations){
		if(p.sampleExists(name)){
			return true;
		}
	}
	return false;
};

uint32_t TPopulationSamples::sampleIndex(const std::string & name) const {
	for(auto& p : _populations){
		if(p.sampleExists(name)){
			return p.sampleIndex(name);
		}
	}
	throw std::runtime_error("uint32_t TPopulationSamples::getSampleIndex(const std::string & name): Sample '" + name + "' does not exist!");
};

std::string TPopulationSamples::sampleName(const uint32_t & index) const {
	if(index >= _numSamples){
		throw std::runtime_error("std::string TPopulationSamples::getSampleName(const uint32_t & index) const: index >= _numSamples!");
	}
	for(auto& p : _populations){
		if(p.sampleIndexExists(index)){
			return p.sampleName(index);
		}
	}
	throw std::runtime_error("std::string TPopulationSamples::getSampleName(const uint32_t & index) const: index not found!");
};

void TPopulationSamples::addSampleNamesToVector(std::vector<std::string> & vec) const{
	for(auto& p : _populations){
		p.addSampleNamesToVector(vec);
	}
};

//----------------------------

void TPopulationSamples::readSamples(std::string filename, coretools::TLog* logfile){
	logfile->listFlush("Reading samples from file '" + filename + "' ...");

	//open file
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + " for reading!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string line;
	bool hasPopColumn = false;
	std::map<std::string, uint32_t>::iterator it;

	//now parse and store
	while(file.good() && !file.eof()){
		++lineNum;
		std::getline(file, line);
		coretools::str::fillContainerFromStringWhiteSpace(line, vec, true);

		//skip empty lines
		if(vec.size() > 0){
			if(lineNum == 1){
				if(vec.size() == 2)
					hasPopColumn = true;
				else {
					_populations.emplace_back("Population");
				}
			}

			if(!hasPopColumn && vec.size() != 1)
				throw "Wrong number of columns in file '" + filename + "' on line " + toString(lineNum) + "! Expected 1, but found " + toString(vec.size()) + ".";

			if(hasPopColumn && vec.size() != 2)
				throw "Wrong number of columns in file '" + filename + "' on line " + toString(lineNum) + "! Expected 2, but found " + toString(vec.size()) + ".";

			//add sample to correct population
			if(hasPopColumn){
				auto it = std::find(_populations.begin(), _populations.end(), vec[1]);
				if(it == _populations.end()){
					it = _populations.emplace(_populations.end(), vec[1]);
				}
				it->addSample(vec[0]);
			} else {
				_populations[0].addSample(vec[0]);
			}
		}
	}

	//close file
	file.close();
	logfile->done();

	//fill sample order by population
	//fill firstIndex in each pop
	_numSamples = 0;
	_indexToPopulationIndex.clear();
	for(auto& p : _populations){
		p.setFirstSampleIndex(_numSamples);
		_numSamples += p.numSamples();
	}
	_fillIndexToPopulationIndex();

	//report assignment
	logfile->startIndent("Will consider the following populations:");
	report(logfile);
	logfile->endIndent();
};

void TPopulationSamples::readSamplesFromVCFNames(std::vector<std::string> & vcfSampleNames){
	//create a single populations
	_populations.emplace_back("Population");

	//add samples to that population
	for(size_t s=0; s<vcfSampleNames.size(); ++s){
		_populations[0].addSample(vcfSampleNames[s]);
	}
	_numSamples = _populations[0].numSamples();
	_fillIndexToPopulationIndex();

	fillVCFOrder(vcfSampleNames);
};

void TPopulationSamples::report(coretools::TLog* Logfile){
	for(auto& p : _populations){
		Logfile->startIndent("Population '" + p.name() + "':");
		p.report(Logfile);
		Logfile->endIndent();
	}
};

//------------------------------------------------------------
void TPopulationSamples::fillVCFOrder(std::vector<std::string> & vcfSampleNames){
	_vcfIndex.resize(vcfSampleNames.size());
	_indexToVCFIndex.resize(vcfSampleNames.size());

	std::vector<bool> sampleInVCF(_numSamples);

	for(uint32_t i = 0; i < vcfSampleNames.size(); ++i){
		//check if sample is used (i.e. exists in one of the populations)
		if(sampleExists(vcfSampleNames[i])){
			if(_vcfIndex[i].used){
				throw "Duplicate sample name '" + vcfSampleNames[i] + "' in VCf header!";
			}

			uint32_t index = sampleIndex(vcfSampleNames[i]);
			_indexToVCFIndex[index] = i;
			sampleInVCF[index] = true;

			_vcfIndex[i].used = true;
			_vcfIndex[i].index = index;
		}
	}

	//check if we lack samples
	for(uint32_t s = 0; s < _numSamples; ++s){
		if(!sampleInVCF[s])
			throw "Sample '" + sampleName(s) + "' missing in VCF!";
	}
};

void TPopulationSamples::_fillIndexToPopulationIndex(){
	_indexToPopulationIndex.resize(_numSamples);
	for(uint32_t p = 0; p < _populations.size(); ++p){
		for(uint32_t y = 0; y < _populations[p].numSamples(); ++y){
			_indexToPopulationIndex[_populations[p].firstSampleIndex() + y] = p;
		}
	}
};


uint32_t TPopulationSamples::sampleIndexInVCF(const uint32_t & index){
	return _indexToVCFIndex[index];
};

uint8_t* TPopulationSamples::getPointerToDataInPop(uint8_t* data, uint32_t population) const {
	return &data[ 3 * _populations[population].firstSampleIndex() ];
};

uint32_t TPopulationSamples::numSamplesMissingInPop(bool* sampleMissing, uint32_t population) const {
	return _populations[population].numSamplesMissing(sampleMissing);
};

uint32_t TPopulationSamples::numSamplesWithDataInPop(bool* sampleMissing, uint32_t population) const {
	return _populations[population].numSamplesWithData(sampleMissing);
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationLikelihoodReader                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////
TPopulationLikelihoodReader::TPopulationLikelihoodReader(){
	_init();
};

TPopulationLikelihoodReader::TPopulationLikelihoodReader(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq){
	_init();
	initialize(Parameters, Logfile, saveAlleleFreq);
};

TPopulationLikelihoodReader::~TPopulationLikelihoodReader(){
	closeVCF();
};

void TPopulationLikelihoodReader::_init(){
	//set all values to defaults
	_initialized = false;
	vcfOpen = false;

	//BED file
	limitToSitesInBed = false;

	//settings
	limitLines = false;
	maxLinesToRead = 0;
	minDepth = 1;
	minNumSamplesWithData = 0;
	freqFilter = 0.0;
	epsilonF = 0.001; //F for EM algorithm to estimate allele frequencies
	minVariantQuality = 0;
	estimateGenotypeFrequencies = false;
    filterOnChr = false;

	//counters
	resetCounters();
};

void TPopulationLikelihoodReader::initialize(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq){
	logfile = Logfile;

	//read parsing parameters
	// do we limit the lines to read?
	maxLinesToRead = Parameters.getParameterWithDefault("limitLines", 0L);
	if(maxLinesToRead > 0){
		limitLines = true;
		logfile->list("Will limit analysis to the first " + toString(maxLinesToRead) + " lines of the VCF file. (parameter 'limitLines')");
	} else {
		limitLines = false;
	}

	//limit to sites in bed file?
	if(Parameters.parameterExists("window")){
		std::string filename = Parameters.getParameter<std::string>("window");
		logfile->list("Will limit analysis to windows listed in BED file '" + filename + "'.");
		bedFile.add(filename);
		logfile->conclude("Will use " + toString(bedFile.size()) + " windows of cumulative length " + toString(bedFile.length()) + " bp on " + toString(bedFile.numChromosomesWithWindows()) + " chromosomes.");
		limitToSitesInBed = true;
	}

	// do we set a depth filter?
	minDepth = Parameters.getParameterWithDefault<int>("minDepth", 1);
	if(minDepth < 1)
		throw "minDepth must be >= 1!";
	if(minDepth > 1)
		logfile->list("Will filter samples to a minimum depth of " + toString(minDepth) + ". (parameter 'minDepth')");

	// do we set a missingness filter?
	minNumSamplesWithData = Parameters.getParameterWithDefault<int>("minSamplesWithData", 1);
	if(minNumSamplesWithData < 0)
		throw "minNumSamplesWithData must be >= 0!";
	if(minNumSamplesWithData > 1)
		logfile->list("Will remove loci where less than " + toString(minNumSamplesWithData) + " samples have data. (parameter 'minSamplesWithData')");

	// parameters to set a filter on the allele frequency?
	freqFilter = Parameters.getParameterWithDefault("minMAF", 0.0); // MAF = minor allele frequency
	if(freqFilter < 0.0 || freqFilter >= 0.5)
		throw "MAF filter must be within (0.0,0.5)!";
	if(freqFilter > 0.0 || saveAlleleFreq){
		estimateGenotypeFrequencies = true;
		epsilonF = Parameters.getParameterWithDefault("epsF", 0.0000001);
		logfile->list("Will filter on an allele frequency of " + toString(freqFilter) + ". (parameter 'epsF')");
	} else {
		estimateGenotypeFrequencies = false;
	}

	//filter on variant quality?
	minVariantQuality = Parameters.getParameterWithDefault<int>("minVariantQuality", 0);
	if(minVariantQuality < 0) throw "minVariantQuality must be >= 0!";
	if(minVariantQuality > 0){
		logfile->list("Will only keep sites with variant quality >= " + toString(minVariantQuality) + ". (parameter 'minVariantQuality')");
	}

    // filter for specific chromosomes?
    if(Parameters.parameterExists("keepChromosomes")) {
        specifyChromosomesToKeep(Parameters, logfile);
    }

	//set progress frequency
	progressFrequency = Parameters.getParameterWithDefault<int>("reportFreq", 10000);

	_initialized = true;
};

void TPopulationLikelihoodReader::specifyChromosomesToKeep(coretools::TParameters & Parameters, coretools::TLog* logfile){
    std::string argument = Parameters.getParameter<std::string>("keepChromosomes");
    if(coretools::str::stringContains(argument, ".txt")){ // specified as a file name
        logfile->startIndent("Reading chromosomes that should be kept from '" + argument + "'");
        std::ifstream keepChromosomesFile(argument.c_str());
        if(!keepChromosomesFile)
            throw std::runtime_error("Failed to open file '" + argument + "'!");
        while(keepChromosomesFile.good() && !keepChromosomesFile.eof()){
            std::string line;
            std::getline(keepChromosomesFile, line);
            std::vector<std::string> vec;
            coretools::str::fillContainerFromStringWhiteSpace(line, vec, true);
            //skip empty lines
            if(!vec.empty())
                chromosomesToKeep.push_back(vec[0]);
        }
        keepChromosomesFile.close();
    }
    else { // specified as a vector on command line
        logfile->startIndent("Reading chromosomes from command line.");
        coretools::str::fillContainerFromString(Parameters.getParameter<std::string>("keepChromosomes"), chromosomesToKeep, ',');
    }

    // write to logfile
    logfile->startIndent("Will keep the following chromosomes in the output file:");
    for (auto it = chromosomesToKeep.begin(); it < chromosomesToKeep.end(); it ++)
        logfile->list(*it);

    filterOnChr = true;

    logfile->endIndent();
    logfile->endIndent();
};

void TPopulationLikelihoodReader::resetCounters(){
	vcfParsingStarted = false;
	_lineCounter = 0;
	_notInBedFile = 0;
	_notBialleleicCounter = 0;
	_missingSNPCounter = 0;
	_lowFreqSNPCounter = 0;
	_lowVariantQualityCounter = 0;
	_noPLCounter = 0;
	_numAcceptedLoci = 0;
	_notOnChrCounter = 0;
	curChr = "";
};

void TPopulationLikelihoodReader::openVCF(std::string vcfFilename){
	if(!_initialized){
		throw "Can not open VCF: TPopulationLikelihoodReader was never initialized!";
	}

	//open input stream
	logfile->startIndent("Reading vcf from file '" + vcfFilename + "'.");
	vcfFile.openStream(vcfFilename);

	//enable parsers
	vcfFile.enablePositionParsing();
	vcfFile.enableVariantParsing();
	vcfFile.enableVariantQualityParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();
	vcfOpen = true;

	//reset counters
	resetCounters();
};

void TPopulationLikelihoodReader::closeVCF(){
	vcfOpen = false;
};

int TPopulationLikelihoodReader::filterOnDepth(TSampleLikelihoods* data, TPopulationSamples & samples){
	int numIndividualsWithData = 0;
	for(uint32_t s = 0; s < samples.numSamples(); ++s){
		int vcfIndex = samples.sampleIndexInVCF(s);

		// depth filter: if a locus has < minDepth reads, flag locus as missing (set all genotype likelihoods = 1)
		if (vcfFile.sampleDepth(vcfIndex) < minDepth){
			vcfFile.setSampleMissing(vcfIndex);
		} else {
			numIndividualsWithData++;
		}
	}

	return numIndividualsWithData;
};

bool TPopulationLikelihoodReader::_readNextLineFromVCF(){
	++_lineCounter;

	//print progress
	if(_lineCounter % progressFrequency == 0)
		printProgressFrequencyFiltering();

	// limit lines
	if(limitLines && _lineCounter > maxLinesToRead){
		logfile->list("Reached limit of " + toString(maxLinesToRead) + " lines.");
		return false;
	}

	//read next line
	if(!vcfFile.next()) return false;

	//all fine!
	return true;
};

bool TPopulationLikelihoodReader::_filterSite(TSampleLikelihoods* data, TPopulationSamples & samples){
	//skip sites with != 2 alleles
	if(vcfFile.getNumAlleles() != 2){
		_notBialleleicCounter++;
		return false;
	}

    // keep chromosomes
    if (filterOnChr && std::find(chromosomesToKeep.begin(), chromosomesToKeep.end(), vcfFile.chr()) == chromosomesToKeep.end()){
        _notOnChrCounter ++;
        return false;
    }

	//skip sites with too low variant quality
	if(minVariantQuality > 0 && (vcfFile.variantQualityIsMissing() || vcfFile.variantQuality() < minVariantQuality)){
		_lowVariantQualityCounter++;
		return false;
	}

	//check if GL or PL is given
	if(!vcfFile.formatColExists("GL") && !vcfFile.formatColExists("PL")){
		++_noPLCounter;
		return false;
	}

	uint32_t numIndividualsWithData = filterOnDepth(data, samples);
	for(uint32_t s = 0; s < samples.numSamples(); ++s){
		unsigned int vcfIndex = samples.sampleIndexInVCF(s);
		vcfFile.fillGenotypeLikelihoods(data[s], vcfIndex);
	}

	// missingness filter: if > percentMissingPerLocus of individuals per locus have are missing, remove locus
	if (numIndividualsWithData < minNumSamplesWithData){
		_missingSNPCounter++;
		return false;
	}

	//filter in MAF
	if(freqFilter > 0.0 || estimateGenotypeFrequencies){
		genoFrequencies.estimate(data, samples.numSamples(), epsilonF);

		if(genoFrequencies.MAF < freqFilter){
			_lowFreqSNPCounter++;
			return false;
		}
	}

	//all fine!
	return true;
};

bool TPopulationLikelihoodReader::_updateChromosomeInfo(){
	//function called when VCF is on new chromosome
	// first remove previous chromosome from vector chromosomesToKeep (not needed anymore; we stop when vector is empty)
	auto it = std::find(chromosomesToKeep.begin(), chromosomesToKeep.end(), curChr);
	if (it != chromosomesToKeep.end()){
		chromosomesToKeep.erase(it);
	}

	//update curChr
	curChr = vcfFile.chr();

	// Is there any chromosome left that should be kept?
	if (filterOnChr && chromosomesToKeep.empty()){
		logfile->list("Parsed all chromosomes that were defined with parameter 'keepChromosomes'.");
		return false;
	}

	return true;
};

bool TPopulationLikelihoodReader::_jumpToNextChromosome(){
	while(!vcfFile.eof && vcfFile.chr() == curChr){
		vcfFile.next();
		++_lineCounter;
		++_notInBedFile;
	}

	if(vcfFile.eof){
		return false;
	} else {
		return _updateChromosomeInfo();
	}
};

void TPopulationLikelihoodReader::printProgressFrequencyFiltering(){
	logfile->list("Parsed " + toString(_lineCounter) + " lines, retained " + toString(_numAcceptedLoci) + " loci in " + timer.formattedTime() + ".");
};

void TPopulationLikelihoodReader::concludeFilters(){
	printProgressFrequencyFiltering();
	if(_notInBedFile > 0)
		logfile->conclude(toString(_notInBedFile) + " loci were not in considered windows (BED file).");
	if(_notBialleleicCounter > 0)
		logfile->conclude(toString(_notBialleleicCounter) + " loci were not bi-allelic.");
	if(_lowVariantQualityCounter > 0)
		logfile->conclude(toString(_lowVariantQualityCounter) + " loci had variant quality < " + toString(minVariantQuality) + ".");
	if(_noPLCounter > 0)
		logfile->conclude(toString(_noPLCounter) + " loci had no PL or GL field.");
	if(_missingSNPCounter > 0)
		logfile->conclude(toString(_missingSNPCounter) + " loci had < " + toString(minNumSamplesWithData) + " samples with data.");
	if(_lowFreqSNPCounter > 0)
		logfile->conclude(toString(_lowFreqSNPCounter) + " loci had MAF < " + toString(freqFilter) + ".");
    if(_notOnChrCounter > 0)
        logfile->conclude(toString(_notOnChrCounter) + " loci were on other chromosomes than specified.");
};


////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationLikelihoodReaderLocus                                                           //
////////////////////////////////////////////////////////////////////////////////////////////////

TPopulationLikelihoodReaderLocus::TPopulationLikelihoodReaderLocus():TPopulationLikelihoodReader(){
	//settings regarding true allele frequencies
	_init();
};

TPopulationLikelihoodReaderLocus::~TPopulationLikelihoodReaderLocus(){
	_closeTrueAlleleFreqFile();
};

TPopulationLikelihoodReaderLocus::TPopulationLikelihoodReaderLocus(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq):TPopulationLikelihoodReader(Parameters, Logfile, saveAlleleFreq){
	_init();
};

void TPopulationLikelihoodReaderLocus::_init(){
	trueFreqFileOpen = false;
	storeTrueAlleleFreq = false;
	_trueAlleleFrequency = -1.0;
};

void TPopulationLikelihoodReaderLocus::initialize(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq){
	TPopulationLikelihoodReader::initialize(Parameters, Logfile, saveAlleleFreq);

	//open true allele freq file
	if(Parameters.parameterExists("trueAlleleFreq")){
		openTrueAlleleFrequenciesFile(Parameters.getParameter<std::string>("trueAlleleFreq"));
	}
};

void TPopulationLikelihoodReaderLocus::openTrueAlleleFrequenciesFile(const std::string trueAlleleFreqFileName){
	if(trueAlleleFreqFileName.find(".gz") == std::string::npos){
		logfile->startIndent("Reading true allele frequencies from gzipped file '" + trueAlleleFreqFileName + "'.");
		trueFreq = new gz::igzstream(trueAlleleFreqFileName.c_str());
	} else {
		logfile->startIndent("Reading true allele frequencies from file '" + trueAlleleFreqFileName + "'.");
		trueFreq = new std::ifstream(trueAlleleFreqFileName.c_str());
	}
	if(!(*trueFreq) || trueFreq->fail() || !trueFreq->good())
		throw "Failed to open file '" + trueAlleleFreqFileName + "'!";
	trueFreqFileOpen = true;
};

void TPopulationLikelihoodReaderLocus::_closeTrueAlleleFreqFile(){
	if(trueFreqFileOpen){
		delete trueFreq;
		trueFreqFileOpen = false;
	}
};

bool TPopulationLikelihoodReaderLocus::_readNextLineFromVCF(){

	if(!TPopulationLikelihoodReader::_readNextLineFromVCF())
		return false;

	//read true allele frequency
	//file is assumed to have exact same dimension: always read next line!
	if(storeTrueAlleleFreq){
		std::string temp;
		getline(*trueFreq, temp);
		std::vector<std::string> tmp;
		coretools::str::fillContainerFromString(temp, tmp, '\t');
		if(tmp.size() != 3)
			throw "wrong number of columns in true allele frequency file!";
		std::string chr = tmp[0];
		uint64_t pos = coretools::str::convertString<uint64_t>(tmp[1]);
		_trueAlleleFrequency = coretools::str::convertString<double>(tmp[2]);
		//check if positions match (allele file is 0-based)
		while(pos < vcfFile.position() - 1){
			getline(*trueFreq, temp);
			coretools::str::fillContainerFromString(temp, tmp, '\t');
			if(tmp.size() != 3)
				throw "wrong number of columns in true allele frequency file!";
			pos = coretools::str::convertString<uint64_t>(tmp[1]);
		}
		if(pos > vcfFile.position() - 1)
			throw "current vcf pos=" + toString(vcfFile.position()) + " is not equal to current trueAlleleFreq position=" + toString(pos);
	}

	//all fine!
	return true;
};

bool TPopulationLikelihoodReaderLocus::readDataFromVCF(TPopulationLikehoodLocus & data, TPopulationSamples & samples){
	data.resize(samples.numSamples());
	return readDataFromVCF(data.samples(), samples);
};

bool TPopulationLikelihoodReaderLocus::readDataFromVCF(TSampleLikelihoods* data, TPopulationSamples & samples){
	//set time at beginning
	if(!vcfParsingStarted){
		vcfParsingStarted = true;
		timer.start();

		//start bed
		if(limitToSitesInBed){
			_curBedWindow = bedFile.begin();
		}
	}

	//read next
	while(_readNextLineFromVCF()){ // new line in vcf-file (= new locus)
		//update chr
		if(curChr != vcfFile.chr()){
			if(!_updateChromosomeInfo()){
				return false;
			}

			//jump bed to this chromosome
			if(limitToSitesInBed){
				//jump to next chromosome until we are on a chromsome with BED windows
				while(!bedFile.hasWindowsOnChr(curChr)){
					if(!_jumpToNextChromosome()){
						return false;
					}
				}

				//get ref id of this chromosome
				_curRefId = bedFile.getRefID(curChr);

				//start windows on this chromosome
				_curBedWindow = bedFile.begin(_curRefId);
			}
		}

		//check if site is used
		if(limitToSitesInBed){
			//stop parsing if we reached end of bed file
			if(_curBedWindow == bedFile.end()){
				logfile->list("Reached end of last window (BED file).");
				return false;
			}

			//check if VCF is ahead or behind BED window
			BAM::TGenomePosition pos(_curRefId, vcfFile.positionZeroBased());

			//jump to next window if position is past current window
			while(*_curBedWindow < pos){
				++_curBedWindow;
				if(_curBedWindow == bedFile.end()){
					logfile->list("Reached end of last window (BED file).");
					return false;
				}
			}

			//continue to next pos in VCF is window is ahead
			if(pos < *_curBedWindow){
				++_notInBedFile;
				continue;
			}

			//else accept current position
			if(!_curBedWindow->within(pos)){
				throw std::runtime_error("TPopulationLikelihoodReaderLocus::readDataFromVCF(TSampleLikelihoods* data, TPopulationSamples & samples, TGlfConverter & glfConverter): something went wrong!");
			}
		}

		//filter
		if(_filterSite(data, samples)){
     		//SNP is accepted!
			++_numAcceptedLoci;
			return true;
		}
    }

	//return false at end of file
	logfile->list("Reached end of VCF file.");
	return false;
};

void TPopulationLikelihoodReaderLocus::writePosition(coretools::TOutputFile & out){
	out << vcfFile.chr() << vcfFile.position() << vcfFile.getRefAllele() << vcfFile.getFirstAltAllele();
};

std::vector<genometools::BiallelicGenotype> TPopulationLikelihoodReaderLocus::biallelicGenotypes(TPopulationSamples & samples) const {
	std::vector<genometools::BiallelicGenotype> vec(samples.numSamples());
    for(uint32_t s = 0; s < samples.numSamples(); ++s) {
        uint32_t vcfIndex = samples.sampleIndexInVCF(s);
        vec[s] = vcfFile.sampleBiallelicGenotype(vcfIndex);
    }
    return vec;
};

genometools::BiallelicGenotype TPopulationLikelihoodReaderLocus::biallelicGenotype(TPopulationSamples & samples, const uint32_t & s) const{
    uint32_t vcfIndex = samples.sampleIndexInVCF(s);
    return vcfFile.sampleBiallelicGenotype(vcfIndex);
};

genometools::Genotype TPopulationLikelihoodReaderLocus::genotype(TPopulationSamples & samples, const uint32_t & s) const{
	uint32_t vcfIndex = samples.sampleIndexInVCF(s);
	return vcfFile.sampleGenotype(vcfIndex);
};

double TPopulationLikelihoodReaderLocus::depth(TPopulationSamples & samples, uint32_t s){
    uint32_t vcfIndex = samples.sampleIndexInVCF(s);
    return vcfFile.sampleDepth(vcfIndex);
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationLikelihoodReaderWindow                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////

TPopulationLikelihoodReaderWindow::TPopulationLikelihoodReaderWindow():TPopulationLikelihoodReader(){

};

TPopulationLikelihoodReaderWindow::TPopulationLikelihoodReaderWindow(coretools::TParameters & Parameters, coretools::TLog* Logfile, bool saveAlleleFreq):TPopulationLikelihoodReader(Parameters, Logfile, saveAlleleFreq){

};

TPopulationLikelihoodReaderWindow::~TPopulationLikelihoodReaderWindow(){

};


bool TPopulationLikelihoodReaderWindow::_readNextLocusAndUpdateChromosome(){
	//read next locus
	if(!_readNextLineFromVCF()) return false; //reached end of VCF

	//update chr
	if(curChr != vcfFile.chr()){
		if(!_updateChromosomeInfo()){
			return false;
		}

		//jump to next chromosome until we are on a chromsome with BED windows
		while(!bedFile.hasWindowsOnChr(curChr)){
			if(!_jumpToNextChromosome()){
				return false;
			}
		}

		//get ref id of this chromosome
		_curRefId = bedFile.getRefID(curChr);
	}

	return true;
};

bool TPopulationLikelihoodReaderWindow::readDataFromVCF(TPopulationLikehoodWindow & window, TPopulationSamples & samples){
	//ASSUMPTION: all chromosomes are present in VCF with at least on position!
	//only works if BED file is open
	if(!limitToSitesInBed){
		throw "Can not read data from VCF into window: no windows defined!";
	}

	//initialize if at befinning of VCF
	if(!vcfParsingStarted){
		vcfParsingStarted = true;
		timer.start();

		//read first locus on first chromosome with windows
		_readNextLocusAndUpdateChromosome();
	} else {
		//advance window
		++_curBedWindow;
		if(_curBedWindow == bedFile.end()) return false;
	}

	//resize window
	window.resize(_curBedWindow->size(), samples.numSamples());

	//1) make sure we are at right window and place in VCF
	//----------------------------------------------------
	//if VCF is ahead of window on same chromosome, advance VCF
	while(!vcfFile.eof && _curRefId == _curBedWindow->refID() && vcfFile.positionZeroBased() < _curBedWindow->from().position()){
		_readNextLocusAndUpdateChromosome();
	}

	//if VCF is at end, is on next chromosome or past window, return empty
	if(vcfFile.eof || _curRefId != _curBedWindow->refID() || vcfFile.positionZeroBased() >= _curBedWindow->to().position()){
		window.fillAsMissingDiploid();
		return true;
	}

	//2) fill window!
	//---------------
	uint32_t index = 0;
	while(vcfFile.positionZeroBased() < _curBedWindow->to().position()){
		uint32_t curIndex = vcfFile.positionZeroBased() - _curBedWindow->from().position();

		//fill empty until cur position
		for(; index<curIndex; ++index){
			window[index].fillAsMissingDiploid();
		}

		//fill cur position
		if(_filterSite(window[index].samples(), samples)){
     		//SNP is accepted!
			++_numAcceptedLoci;
		} else {
			window[index].fillAsMissingDiploid();
		}
		++index;

		//read next
		if(!_readNextLocusAndUpdateChromosome()){
			//end of file: fill empty until end
			for(; index<window.numLoci(); ++index){
				window[index].fillAsMissingDiploid();
			}
			break;
		}

		//if VCF is on next chromosome, fill empty until end
		if(_curRefId != _curBedWindow->refID()){
			for(; index<window.numLoci(); ++index){
				window[index].fillAsMissingDiploid();
			}
			break;
		}
	}

	//all fine!
	return true;
};

void TPopulationLikelihoodReaderWindow::writeWindow(coretools::TOutputFile & out){
	out << bedFile.getChromosomeName(_curBedWindow->refID()) << _curBedWindow->from().position() << _curBedWindow->to().position();
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TVcfFilter                                                                                 //
////////////////////////////////////////////////////////////////////////////////////////////////
/*
TVcfFilter::TVcfFilter(TParameters & Parameters, TLog* Logfile){
	vcfRead = false;
	_numLoci = 0;
	logfile = Logfile;
}

void TVcfFilter::filterVCF(TParameters & Parameters){
	if(vcfRead)
		throw "VCF already read!";

	//create reader
	bool saveAlleleFrequencies = false;
	TPopulationLikelihoodReader reader(Parameters, logfile, saveAlleleFrequencies);

	// open vcf file
	vcfFilename = Parameters.getParameter<std::string>("vcf");
	logfile->startIndent("Filtering VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename, logfile);

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());


	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	uint8_t* curLocus = new uint8_t[samples.numSamples() * 3];
	bool* sampleIsMissing = new bool[samples.numSamples()];

	//output file name
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp) + "_filtered.vcf.gz";

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    _numLoci = 0;
    reader.filterVCF(curLocus, sampleIsMissing, samples, logfile, outputName);

    //clean up
	vcfRead = true;
	delete[] curLocus;

    //report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};
*/

////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationLikelihoods                                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////
TPopulationLikelihoods::TPopulationLikelihoods(){
	init();
};


TPopulationLikelihoods::TPopulationLikelihoods(coretools::TParameters & Parameters, coretools::TLog* Logfile){
	init();
	readData(Parameters, Logfile);
};

TPopulationLikelihoods::~TPopulationLikelihoods(){
	clean();
};

void TPopulationLikelihoods::init(){
	vcfRead = false;
	_numLoci = 0;
	saveAlleleFrequencies = false;
	saveTrueAlleleFrequencies = false;
};

void TPopulationLikelihoods::clean(){
	for(TSampleLikelihoods* it : data){
		delete it;
	}
	data.clear();
	vcfRead = false;
};

void TPopulationLikelihoods::readData(coretools::TParameters & Parameters, coretools::TLog* Logfile){
	//check if we limit samples
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameter<std::string>("samples"), Logfile);

	//read Data
	readDataFromVCF(Parameters, Logfile);

	//make sure loop goes across all samples
	useAllSamples();
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read data from VCF-file                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////
void TPopulationLikelihoods::readDataFromVCF(coretools::TParameters & Parameters, coretools::TLog* logfile){
	if(vcfRead)
		throw "VCF already read!";

	//create reader
	TPopulationLikelihoodReaderLocus reader(Parameters, logfile, saveAlleleFrequencies);
	if(saveAlleleFrequencies){
		reader.doEstimateGenotypeFrequencies();
	}
	if(saveTrueAlleleFrequencies)
		reader.doSaveTrueAlleleFrequencies();

	// open vcf file
	vcfFilename = Parameters.getParameter<std::string>("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename);

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());


	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	data.emplace_back(new TSampleLikelihoods[samples.numSamples()]);

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    std::string curChr = "";
    while(reader.readDataFromVCF(data.back(), samples)){
		//update chromosome name
		if(reader.chr() != curChr){
			chromosomes.emplace(_numLoci, reader.chr());
			curChr = reader.chr();
		}

		//store SNP info
		position.emplace_back(reader.position());
		if(saveAlleleFrequencies)
			alleleFrequencies.emplace_back(reader.allelFrequency());
		if(saveTrueAlleleFrequencies)
			trueAlleleFrequencies.emplace_back(reader.trueAlleleFrequency());

		//update for next
		data.emplace_back(new TSampleLikelihoods[samples.numSamples()]);
    }

    //remove last, empty container
    delete data.back();
    data.pop_back();

    //clean up
	vcfRead = true;
	_numLoci = data.size();

    //report final status
	logfile->endIndent();
	reader.concludeFilters();
	if(_numLoci < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// get-functions                                                                                //
//////////////////////////////////////////////////////////////////////////////////////////////////

int TPopulationLikelihoods::getNumIndividuals(){
	return samples.numSamples();
};

long TPopulationLikelihoods::getNumLoci(){
	return _numLoci;
};

TSampleLikelihoods* TPopulationLikelihoods::getDataAtLocus(long index){
	return data[index];
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// functions to loop over data                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////
void TPopulationLikelihoods::useAllSamples(){
	individualStartIndex = 0;
	usedSampleSize = samples.numSamples();
};

void TPopulationLikelihoods::limitToSinglePopulation(int population){
	individualStartIndex = samples.startIndex(population);
	usedSampleSize = samples.numSamplesInPop(population);
};

void TPopulationLikelihoods::begin(){
	curLocusIndex = 0;
	curChrIt = chromosomes.begin();
	nextChrIt = std::next(chromosomes.begin(),1);
};

void TPopulationLikelihoods::beginAll(){
	useAllSamples();
	begin();
};

void TPopulationLikelihoods::beginOnePop(int population){
	limitToSinglePopulation(population);
	begin();
};

bool TPopulationLikelihoods::end(){
	return curLocusIndex == _numLoci;
};

void TPopulationLikelihoods::next(){
	++curLocusIndex;

	//are we on new chromosome?
	if(nextChrIt != chromosomes.end() && curLocusIndex == nextChrIt->first){
		++curChrIt;
		++nextChrIt;
	}
};

TSampleLikelihoods* TPopulationLikelihoods::curData(){
	return data[curLocusIndex];
};

std::string TPopulationLikelihoods::curSampleName(int index){
	return samples.sampleName(individualStartIndex + index);
};

int TPopulationLikelihoods::curSampleSize(){
	return usedSampleSize;
};

std::string TPopulationLikelihoods::curChr(){
	return curChrIt->second;
};

long TPopulationLikelihoods::curPosition(){
	return position[curLocusIndex];
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// print data                                                                                   //
//////////////////////////////////////////////////////////////////////////////////////////////////

void TPopulationLikelihoods::print(){
	begin();

	//write header
	std::cout << "Chr\tPos";
	for(int s=0; s<curSampleSize(); ++s){
		std::cout << "\t" << curSampleName(s);
	}
	std::cout << std::endl;

	//print data
	for(; !end(); next()){
		//print chromosome and position
		std::cout << curChr() << "\t" << curPosition();

		//print data
		TSampleLikelihoods* data = curData();
		for(int s=0; s<curSampleSize(); ++s){
			std::cout << "\t";
			data[s].print();
		}

		std::cout << std::endl;
	}
};

}; //end namespace
