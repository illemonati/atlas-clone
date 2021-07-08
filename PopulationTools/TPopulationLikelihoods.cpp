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

////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationSamples                                                                         //
////////////////////////////////////////////////////////////////////////////////////////////////
TPopulationSamples::TPopulationSamples(){
	_init();
};

TPopulationSamples::TPopulationSamples(std::string filename, TLog* logfile){
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

uint32_t TPopulationSamples::getPopulationIndex(const std::string & name) const {
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

uint32_t TPopulationSamples::getSampleIndex(const std::string & name) const {
	for(auto& p : _populations){
		if(p.sampleExists(name)){
			return p.sampleIndex(name);
		}
	}
	throw std::runtime_error("uint32_t TPopulationSamples::getSampleIndex(const std::string & name): Sample '" + name + "' does not exist!");
};

std::string TPopulationSamples::getSampleName(const uint32_t & index) const {
	if(index >= _numSamples){
		throw std::runtime_error("std::string TPopulationSamples::getSampleName(const uint32_t & index) const: index >= _numSamples!");
	}
	for(auto& p : _populations){
		if(p.sampleIndexExists(index)){
			return p.sampleName(index);
		}
	}
};

void TPopulationSamples::addSampleNamesToVector(std::vector<std::string> & vec) const{
	for(auto& p : _populations){
		p.addSampleNamesToVector(vec);
	}
};

//----------------------------

void TPopulationSamples::readSamples(std::string filename, TLog* logfile){
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
		fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);

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
	for(auto& p : _populations){
		p.setFirstSampleIndex(_numSamples);
		_numSamples += p.numSamples();
	}

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

	fillVCFOrder(vcfSampleNames);
};

void TPopulationSamples::report(TLog* Logfile){
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

			uint32_t index = getSampleIndex(vcfSampleNames[i]);
			_indexToVCFIndex[index] = i;
			sampleInVCF[index] = true;

			_vcfIndex[i].used = true;
			_vcfIndex[i].index = index;
		}
	}

	//check if we lack samples
	for(uint32_t s = 0; s < _numSamples; ++s){
		if(!sampleInVCF[s])
			throw "Sample '" + getSampleName(s) + "' missing in VCF!";
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

TPopulationLikelihoodReader::TPopulationLikelihoodReader(TParameters & Parameters, TLog* Logfile, bool saveAlleleFreq){
	_init();
	initialize(Parameters, Logfile, saveAlleleFreq);
};

TPopulationLikelihoodReader::~TPopulationLikelihoodReader(){
	closeVCF();
	if(limitToSitesInBed)
		delete bedFile;
};

void TPopulationLikelihoodReader::_init(){
	//set all values to defaults
	_initialized = false;
	vcfOpen = false;

	//BED file
	bedFile = nullptr;
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

void TPopulationLikelihoodReader::initialize(TParameters & Parameters, TLog* Logfile, bool saveAlleleFreq){
	logfile = Logfile;

	//read parsing parameters
	// do we limit the lines to read?
	maxLinesToRead = Parameters.getParameterLongWithDefault("limitLines", 0);
	if(maxLinesToRead > 0){
		limitLines = true;
		logfile->list("Will limit analysis to the first " + toString(maxLinesToRead) + " lines of the VCF file. (parameter 'limitLines')");
	} else {
		limitLines = false;
	}

	//limit to sites in bed file?
	if(Parameters.parameterExists("window")){
		std::string filename = Parameters.getParameterString("window");
		logfile->list("Will limit analysis to windows listed in BED file '" + filename + "'.");
		bedFile = new TBed(filename);
		logfile->conclude("Will use " + toString(bedFile->size()) + " windows of cumulative length " + toString(bedFile->length()) + " bp on " + toString(bedFile->getNumChromosomes()) + " chromosomes.");
		limitToSitesInBed = true;
	}

	// do we set a depth filter?
	minDepth = Parameters.getParameterIntWithDefault("minDepth", 1);
	if(minDepth < 1)
		throw "minDepth must be >= 1!";
	if(minDepth > 1)
		logfile->list("Will filter samples to a minimum depth of " + toString(minDepth) + ". (parameter 'minDepth')");

	// do we set a missingness filter?
	minNumSamplesWithData = Parameters.getParameterIntWithDefault("minSamplesWithData", 1);
	if(minNumSamplesWithData < 0)
		throw "minNumSamplesWithData must be >= 0!";
	if(minNumSamplesWithData > 1)
		logfile->list("Will remove loci where less than " + toString(minNumSamplesWithData) + " samples have data. (parameter 'minSamplesWithData')");

	// parameters to set a filter on the allele frequency?
	freqFilter = Parameters.getParameterDoubleWithDefault("minMAF", 0.0); // MAF = minor allele frequency
	if(freqFilter < 0.0 || freqFilter >= 0.5)
		throw "MAF filter must be within (0.0,0.5)!";
	if(freqFilter > 0.0 || saveAlleleFreq){
		estimateGenotypeFrequencies = true;
		epsilonF = Parameters.getParameterDoubleWithDefault("epsF", 0.0000001);
		logfile->list("Will filter on an allele frequency of " + toString(freqFilter) + ". (parameter 'epsF')");
	} else {
		estimateGenotypeFrequencies = false;
	}

	//filter on variant quality?
	minVariantQuality = Parameters.getParameterIntWithDefault("minVariantQuality", 0);
	if(minVariantQuality < 0) throw "minVariantQuality must be >= 0!";
	if(minVariantQuality > 0){
		logfile->list("Will only keep sites with variant quality >= " + toString(minVariantQuality) + ". (parameter 'minVariantQuality')");
	}

    // filter for specific chromosomes?
    if(Parameters.parameterExists("keepChromosomes")) {
        specifyChromosomesToKeep(Parameters, logfile);
    }

	//set progress frequency
	progressFrequency = Parameters.getParameterIntWithDefault("reportFreq", 10000);

	_initialized = true;
};

void TPopulationLikelihoodReader::specifyChromosomesToKeep(TParameters & Parameters, TLog* logfile){
    std::string argument = Parameters.getParameterString("keepChromosomes");
    if(stringContains(argument, ".txt")){ // specified as a file name
        logfile->startIndent("Reading chromosomes that should be kept from '" + argument + "'");
        std::ifstream keepChromosomesFile(argument.c_str());
        if(!keepChromosomesFile)
            throw std::runtime_error("Failed to open file '" + argument + "'!");
        while(keepChromosomesFile.good() && !keepChromosomesFile.eof()){
            std::string line;
            std::getline(keepChromosomesFile, line);
            std::vector<std::string> vec;
            fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
            //skip empty lines
            if(!vec.empty())
                chromosomesToKeep.push_back(vec[0]);
        }
        keepChromosomesFile.close();
    }
    else { // specified as a vector on command line
        logfile->startIndent("Reading chromosomes from command line.");
        fillVectorFromString(Parameters.getParameterString("keepChromosomes"), chromosomesToKeep, ',');
    }

    // write to logfile
    logfile->startIndent("Will keep the following chromosomes in the output file:");
    for (auto it = chromosomesToKeep.begin(); it < chromosomesToKeep.end(); it ++)
        logfile->list(*it);

    filterOnChr = true;

    logfile->endIndent();
    logfile->endIndent();
}

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
	bool isZipped = false;
	if(vcfFilename.find(".gz") == std::string::npos){
		logfile->startIndent("Reading vcf from file '" + vcfFilename + "'.");
	} else {
		logfile->startIndent("Reading vcf from gzipped file '" + vcfFilename + "'.");
		isZipped = true;
	}

	vcfFile.openStream(vcfFilename, isZipped);

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
		if (vcfFile.sampleDepth(vcfIndex) < minDepth)
			vcfFile.setSampleMissing(vcfIndex);
		else
			numIndividualsWithData++;

		//store if missing and haploid
		data[s].isMissing = vcfFile.sampleIsMissing(vcfIndex);
		data[s].isHaploid = vcfFile.sampleIsHaploid(vcfIndex);
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

bool TPopulationLikelihoodReader::_filterSite(TSampleLikelihoods* data, TPopulationSamples & samples, TGlfConverter & glfConverter){
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
	uint32_t numIndividualsWithData = 0;
	if(vcfFile.formatColExists("GL")){
		numIndividualsWithData = filterOnDepth(data, samples);
		double tmp[3];
		for(uint32_t s = 0; s < samples.numSamples(); ++s){
			int vcfIndex = samples.sampleIndexInVCF(s);
			vcfFile.fillLog10GenotypeLikelihoods(vcfIndex, tmp[0], tmp[1], tmp[2]);
			data[s].glfLikelihood_0 = glfConverter.log10ToGlfFormat(tmp[0]);
			data[s].glfLikelihood_1 = glfConverter.log10ToGlfFormat(tmp[1]);
			data[s].glfLikelihood_2 = glfConverter.log10ToGlfFormat(tmp[2]);
		}
	} else if(vcfFile.formatColExists("PL")){
		numIndividualsWithData = filterOnDepth(data, samples);
		uint8_t tmp[3];
		for(uint32_t s = 0; s < samples.numSamples(); ++s){
			int vcfIndex = samples.sampleIndexInVCF(s);
			vcfFile.fillPhredScore(vcfIndex, tmp[0], tmp[1], tmp[2]);
			data[s].glfLikelihood_0 = glfConverter.phredToGlfFormat(tmp[0]);
			data[s].glfLikelihood_1 = glfConverter.phredToGlfFormat(tmp[1]);
			data[s].glfLikelihood_2 = glfConverter.phredToGlfFormat(tmp[2]);
		}
	} else {
		++_noPLCounter;
		return false;
	}

	// missingness filter: if > percentMissingPerLocus of individuals per locus have are missing, remove locus
	if (numIndividualsWithData < minNumSamplesWithData){
		_missingSNPCounter++;
		return false;
	}

	//filter in MAF
	if(freqFilter > 0.0 || estimateGenotypeFrequencies){
		genoFrequencies.estimate(data, samples.numSamples(), glfConverter, epsilonF);

		if(genoFrequencies.MAF < freqFilter){
			_lowFreqSNPCounter++;
			return false;
		}
	}

	//all fine!
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
		//update chr
		curChr = vcfFile.chr();

		if(limitToSitesInBed){
			bedFile->setChr(vcfFile.chr());
		}
		return true;
	}
};

void TPopulationLikelihoodReader::printProgressFrequencyFiltering(){
	struct timeval end{};
	gettimeofday(&end, nullptr);
	float runtime = static_cast<double>(end.tv_sec  - startTime.tv_sec)/60.0;
	logfile->list("Parsed " + toString(_lineCounter) + " lines, retained " + toString(_numAcceptedLoci) + " loci in " + toString(runtime) + " min");
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

TPopulationLikelihoodReaderLocus::TPopulationLikelihoodReaderLocus(TParameters & Parameters, TLog* Logfile, bool saveAlleleFreq):TPopulationLikelihoodReader(Parameters, Logfile, saveAlleleFreq){
	_init();
};

void TPopulationLikelihoodReaderLocus::_init(){
	trueFreqFileOpen = false;
	storeTrueAlleleFreq = false;
	_trueAlleleFrequency = -1.0;
};

void TPopulationLikelihoodReaderLocus::initialize(TParameters & Parameters, TLog* Logfile, bool saveAlleleFreq){
	TPopulationLikelihoodReader::initialize(Parameters, Logfile, saveAlleleFreq);

	//open true allele freq file
	if(Parameters.parameterExists("trueAlleleFreq")){
		openTrueAlleleFrequenciesFile(Parameters.getParameterString("trueAlleleFreq"));
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
		fillVectorFromString(temp, tmp, '\t');
		if(tmp.size() != 3)
			throw "wrong number of columns in true allele frequency file!";
		std::string chr = tmp[0];
		uint32_t pos = convertString<uint32_t>(tmp[1]);
		_trueAlleleFrequency = convertString<double>(tmp[2]);
		//check if positions match (allele file is 0-based)
		while(pos < vcfFile.position() - 1){
			getline(*trueFreq, temp);
			fillVectorFromString(temp, tmp, '\t');
			if(tmp.size() != 3)
				throw "wrong number of columns in true allele frequency file!";
			pos = convertString<uint32_t>(tmp[1]);
		}
		if(pos > vcfFile.position() - 1)
			throw "current vcf pos=" + toString(vcfFile.position()) + " is not equal to current trueAlleleFreq position=" + toString(pos);
	}

	//all fine!
	return true;
};

bool TPopulationLikelihoodReaderLocus::readDataFromVCF(TPopulationLikehoodLocus & data, TPopulationSamples & samples, TGlfConverter & glfConverter){
	data.resize(samples.numSamples());
	return readDataFromVCF(data.samples(), samples, glfConverter);
};

bool TPopulationLikelihoodReaderLocus::readDataFromVCF(TSampleLikelihoods* data, TPopulationSamples & samples, TGlfConverter & glfConverter){
	//set time at beginning
	if(!vcfParsingStarted){
		vcfParsingStarted = true;
		gettimeofday(&startTime, NULL);
	}

	//read next
	while(_readNextLineFromVCF()){ // new line in vcf-file (= new locus)
		//update chr
		if(curChr != vcfFile.chr()){
            // first remove previous chromosome from vector chromosomesToKeep (not needed anymore; we stop when vector is empty)
            if (filterOnChr && std::find(chromosomesToKeep.begin(), chromosomesToKeep.end(), curChr) != chromosomesToKeep.end()){
                chromosomesToKeep.erase(std::remove(chromosomesToKeep.begin(), chromosomesToKeep.end(), curChr), chromosomesToKeep.end());
            }
			curChr = vcfFile.chr();

			if(limitToSitesInBed){
				bedFile->setChr(vcfFile.chr());
			}
		}

		//check if site is used
		if(limitToSitesInBed){
			//stop parsing if we reached end of bed file
			if(bedFile->reachedEnd()){
				logfile->list("Reached end of last window (BED file).");
				return false;
			}

			//jump to next window if position is past current window
			while(!bedFile->reachedEndOfChr() && vcfFile.positionZeroBased() >= bedFile->curWindowEnd()){
				bedFile->nextWindow();
			}

			if(bedFile->reachedEndOfChr()){
				++_notInBedFile;
				continue;
			}

			//skip if position ahead of current window
			if(vcfFile.positionZeroBased() < bedFile->curWindowStart()){
				++_notInBedFile;
				continue;
			}

			//else accept current position
		}

		// Is there any chromosome left that should be kept?
		if (filterOnChr && chromosomesToKeep.empty()){
		    logfile->list("Parsed all chromosomes that were defined with parameter 'keepChromosomes'.");
		    return false;
		}

		//filter
		if(_filterSite(data, samples, glfConverter)){
     		//SNP is accepted!
			++_numAcceptedLoci;
			return true;
		}
    }

	//return false at end of file
	logfile->list("Reached end of VCF file.");
	return false;
};

void TPopulationLikelihoodReaderLocus::writePosition(TOutputFile & out){
	out << vcfFile.chr() << vcfFile.position() << vcfFile.getRefAllele() << vcfFile.getFirstAltAllele();
};

void TPopulationLikelihoodReaderLocus::fillGenotypes(TPopulationSamples & samples, u_int8_t * genotypes){
    for(uint32_t s = 0; s < samples.numSamples(); ++s) {
        uint32_t vcfIndex = samples.sampleIndexInVCF(s);
        genotypes[s] = vcfFile.sampleGenotype(vcfIndex);
    }
}

uint8_t TPopulationLikelihoodReaderLocus::genotype(TPopulationSamples & samples, uint32_t s){
    uint32_t vcfIndex = samples.sampleIndexInVCF(s);
    return vcfFile.sampleGenotype(vcfIndex);
}

double TPopulationLikelihoodReaderLocus::depth(TPopulationSamples & samples, uint32_t s){
    uint32_t vcfIndex = samples.sampleIndexInVCF(s);
    return vcfFile.sampleDepth(vcfIndex);
}

////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationLikelihoodReaderWindow                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////

TPopulationLikelihoodReaderWindow::TPopulationLikelihoodReaderWindow():TPopulationLikelihoodReader(){

};

TPopulationLikelihoodReaderWindow::TPopulationLikelihoodReaderWindow(TParameters & Parameters, TLog* Logfile, bool saveAlleleFreq):TPopulationLikelihoodReader(Parameters, Logfile, saveAlleleFreq){

};

TPopulationLikelihoodReaderWindow::~TPopulationLikelihoodReaderWindow(){

};

bool TPopulationLikelihoodReaderWindow::readDataFromVCF(TPopulationLikehoodWindow & window, TPopulationSamples & samples, TGlfConverter & glfConverter){
	//set time at beginning
	if(!vcfParsingStarted){
		vcfParsingStarted = true;
		gettimeofday(&startTime, NULL);

		//read first locus
		if(!_readNextLineFromVCF()) return false; //reached end of VCF
	}

	//only works if BED file is open
	if(!limitToSitesInBed){
		throw "Can not read data from VCF into window: no windows defined!";
	}

	//1) make sure we are at right window and place in VCF
	//----------------------------------------------------
	//go to next window on current chromosome
	bedFile->nextWindow();
	if(bedFile->reachedEnd() || vcfFile.eof) return false;

	//if bed is at end of chr, jump to next chr
	while(bedFile->reachedEndOfChr()){
		if(!_jumpToNextChromosome()){
			//return false at end of file
			logfile->list("Reached end of VCF file.");
			return false;
		}
		bedFile->setChr(curChr);

		if(bedFile->reachedEnd()) return false;
	}

	//if VCF is on same chromosome but behind, advance VCF
	while(!vcfFile.eof && vcfFile.chr() == bedFile->curChr() && vcfFile.positionZeroBased() < bedFile->curWindowStart()){
		_readNextLineFromVCF();
	}

	//resize window
	window.resize(bedFile->curWindowSize(), samples.numSamples());

	//if VCF is at end, is on next chromosome or past window, return empty
	if(vcfFile.eof || vcfFile.chr() != bedFile->curChr() || vcfFile.positionZeroBased() > bedFile->curWindowEnd()){
		window.fillAsMissing();
		return true;
	}

	//2) fill window!
	//---------------
	uint32_t index = 0;
	while(vcfFile.positionZeroBased() < bedFile->curWindowEnd()){
		uint32_t curIndex = vcfFile.positionZeroBased() - bedFile->curWindowStart();

		//fill empty until cur position
		for(; index<curIndex; ++index){
			window[index].fillAsMissing();
		}

		//fill cur position
		if(_filterSite(window[index].samples(), samples, glfConverter)){
     		//SNP is accepted!
			++_numAcceptedLoci;
		} else {
			window[index].fillAsMissing();
		}
		++index;

		//read next
		if(!_readNextLineFromVCF()){
			//end of file: fill empty until end
			for(; index<window.numLoci(); ++index){
				window[index].fillAsMissing();
			}
			break;
		}

		//if VCF is on next chromosome, fill empty until end
		if(vcfFile.chr() != bedFile->curChr()){
			for(; index<window.numLoci(); ++index){
				window[index].fillAsMissing();
			}
			break;
		}
	}

	//all fine!
	return true;
};

void TPopulationLikelihoodReaderWindow::writeWindow(TOutputFile & out){
	out << bedFile->curChr() << bedFile->curWindowStart() + 1 << bedFile->curWindowEnd() + 1;
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
	vcfFilename = Parameters.getParameterString("vcf");
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
	std::string outputName = Parameters.getParameterStringWithDefault("out", tmp) + "_filtered.vcf.gz";

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


TPopulationLikelihoods::TPopulationLikelihoods(TParameters & Parameters, TLog* Logfile){
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

void TPopulationLikelihoods::readData(TParameters & Parameters, TLog* Logfile){
	//check if we limit samples
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameterString("samples"), Logfile);

	//read Data
	readDataFromVCF(Parameters, Logfile);

	//make sure loop goes across all samples
	useAllSamples();
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read data from VCF-file                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////
void TPopulationLikelihoods::readDataFromVCF(TParameters & Parameters, TLog* logfile){
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
	vcfFilename = Parameters.getParameterString("vcf");
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
    while(reader.readDataFromVCF(data.back(), samples, glfConverter)){
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
	return samples.getSampleName(individualStartIndex + index);
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
			std::cout << "\t" << toString(data[s].glfLikelihood_0) << "," << toString(data[s].glfLikelihood_1) << "," << toString(data[s].glfLikelihood_2);
		}

		std::cout << std::endl;
	}
};

