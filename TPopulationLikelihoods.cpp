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
	init();
};

TPopulationSamples::TPopulationSamples(std::string filename, TLog* logfile){
	init();
	readSamples(filename, logfile);
};

TPopulationSamples::~TPopulationSamples(){
	if(_hasSamples){
		delete[] numSamplesPerPop;
		delete[] startIndexPerPop;
	}
	if(_VCF_order_initialized)
		delete[] _VCF_order;
};

void TPopulationSamples::init(){
	_numPopulations = 0;
	_numSamples = 0;
	numSamplesPerPop = NULL;
	_hasSamples = false;
	_VCF_order = NULL;
	_VCF_order_initialized = false;
};

void TPopulationSamples::readSamples(std::string filename, TLog* logfile){
	logfile->startIndent("Reading samples from file '" + filename + "':");

	//open file
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + " for reading!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string line;
	bool hasPopColumn = false;
	_numPopulations = 0;

	//now parse and store
	while(file.good() && !file.eof()){
		++lineNum;
		std::getline(file, line);
		fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(lineNum == 1 && vec.size() == 2) hasPopColumn = true;

			if(!hasPopColumn && vec.size() != 1)
				throw "Wrong number of columns in file '" + filename + "' on line " + toString(lineNum) + "! Expected 1, but found " + toString(vec.size()) + ".";

			if(hasPopColumn && vec.size() != 2)
				throw "Wrong number of columns in file '" + filename + "' on line " + toString(lineNum) + "! Expected 2, but found " + toString(vec.size()) + ".";

			//check if population exists
			std::map<std::string, int>::iterator it = populations.find(vec[1]);
			if(it == populations.end()){
				populations.emplace(vec[2], _numPopulations);
				it = populations.find(vec[1]);
				_numPopulations++;
			}

			//store sample
			if(samples.find(vec[0]) != samples.end())
				throw "Duplicate sample name '" + vec[0] + "' in file '" + filename + "' on line " + toString(lineNum) + "!";
			samples.emplace(vec[0], it->second);
		}
	}

	//close file
	file.close();

	//fill sample order by population
	fillSampleOrder();
};

void TPopulationSamples::fillSampleOrder(){
	//fill counters
	_numSamples = samples.size();
	if(_numSamples > 0){
		_hasSamples = true;
		numSamplesPerPop = new int[_numPopulations];
		startIndexPerPop = new int[_numPopulations];
		//put samples in order by populations
		int nextIndex = 0;
		for(int p=0; p<_numPopulations; ++p){
			numSamplesPerPop[p] = 0;
			startIndexPerPop[p] = nextIndex;
			for(std::map<std::string, int>::iterator it = samples.begin(); it != samples.end(); ++it){
				if(it->second == p){
					sampleOrder.emplace(it->first, nextIndex);
					++nextIndex;
					++numSamplesPerPop[p];
				}
			}
		}
	}
};

bool TPopulationSamples::sampleIsUsed(const std::string& name){
	if(samples.find(name) == samples.end()) return false;
	return true;
};

int TPopulationSamples::getOrderedSampleIndex(const std::string & name){
	if(!sampleIsUsed(name))
		throw "Sample '" + name + "' does not exist!";
	return sampleOrder.find(name)->second;
};

std::string TPopulationSamples::getNameFromOrderedIndex(int index){
	if(index >= _numSamples)
		throw "Sample index " + toString(index) + " out of range!";

	std::string name;
	for(std::map<std::string, int>::iterator it = sampleOrder.begin(); it != sampleOrder.end(); ++it){
		if(it->second == index)
			name = it->first;
	}
	return name;
};

void TPopulationSamples::readSamplesFromVCFNames(std::vector<std::string> & vcfSampleNames){
	_numSamples = vcfSampleNames.size();
	_VCF_order = new int[_numSamples];
	_VCF_order_initialized = true;

	//set a single populations
	_numPopulations = 1;
	populations.emplace("Population", 0);

	//save samples
	for(size_t s=0; s<vcfSampleNames.size(); ++s){
		samples.emplace(vcfSampleNames[s], 0);
		_VCF_order[s] = s;
	}

	//fill sample order by population
	fillSampleOrder();
};


void TPopulationSamples::fillVCFOrder(std::vector<std::string> & vcfSampleNames){
	_VCF_order = new int[_numSamples];
	_VCF_order_initialized = true;
	bool sampleUsed[_numSamples];
	for(int s=0; s<_numSamples; ++s)
		sampleUsed[s] = false;

	int vcfIndex = 0;
	for(std::vector<std::string>::iterator it = vcfSampleNames.begin(); it != vcfSampleNames.end(); ++it, ++vcfIndex){
		if(sampleIsUsed(*it)){
			int orderedIndex = getOrderedSampleIndex(*it);
			if(sampleUsed[orderedIndex])
				throw "Duplicate sample name '" + *it + "' in VCf header!";
			_VCF_order[orderedIndex] = vcfIndex;
			sampleUsed[orderedIndex] = true;
		}
	}

	//check if we lack samples
	for(int s=0; s<_numSamples; ++s){
		if(!sampleUsed[s])
			throw "Sample '" + getNameFromOrderedIndex(s) + "' missing in VCF!";
	}
};




////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationLikelihoodReader                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////
TPopulationLikelihoodReader::TPopulationLikelihoodReader(){
	vcfOpen = false;

	//settings
	limitLines = 0;
	minDepth = 1;
	maxMissing = 0;
	freqFilter = 0.0;
	epsilonF = 0.001; //F for EM algorithm to estimate allele frequencies
	minVariantQuality = 0;
	estimateGenotypeFrequencies = false,

	//counters
	resetCounters();

	//allele frequency
	_alleleFrequency = 0.0;
	_MAF = 0.0;
};

TPopulationLikelihoodReader::TPopulationLikelihoodReader(TParameters & Parameters, TLog* Logfile){
	initialize(Parameters, Logfile);
};

TPopulationLikelihoodReader::~TPopulationLikelihoodReader(){
	closeVCF();
};

void TPopulationLikelihoodReader::initialize(TParameters & Parameters, TLog* logfile){
	//read parsing parameters
	// do we limit the lines to read?
	limitLines = Parameters.getParameterLongWithDefault("limitLines", -1);
	if(limitLines > 0)
		logfile->list("Will limit analysis to the first " + toString(limitLines) + " lines of the VCF file.");

	// do we set a depth filter?
	minDepth = Parameters.getParameterIntWithDefault("minDepth", 1);
	if(minDepth > 1)
		logfile->list("Will filter samples to a minimum depth of " + toString(minDepth) + ".");

	// do we set a missingness filter?
	maxMissing = Parameters.getParameterIntWithDefault("maxMissing", 1.0);
	if(maxMissing < 0.0 || maxMissing > 1.0)
		throw "maxMissing must be within (0, 1)!";
	if(maxMissing > 0)
		logfile->list("Will remove loci where more than " + toString(maxMissing) + " of samples are missing.");

	// parameters to set a filter on the allele frequency?
	freqFilter = Parameters.getParameterDoubleWithDefault("minMAF", 0.0); // MAF = minor allele frequency
	if(freqFilter < 0.0 || freqFilter >= 0.5)
		throw "MAF filter must be within (0.0,0.5)!";
	if(freqFilter > 0.0){
		estimateGenotypeFrequencies = true;
		epsilonF = Parameters.getParameterDoubleWithDefault("epsF", 0.0001);
		logfile->list("Will filter on an allele frequency of " + toString(freqFilter) + ".");
	}

	//filter on variant quality?
	minVariantQuality = Parameters.getParameterIntWithDefault("minVariantQuality", 0);
	if(minVariantQuality < 0) throw "minVariantQuality must be >= 0!";
	if(minVariantQuality > 0){
		logfile->list("Will only keep sites with variant quality >= " + toString(minVariantQuality) + ".");
	}

	//set progress frequency
	progressFrequency = Parameters.getParameterIntWithDefault("reportFreq", 10000);
};

void TPopulationLikelihoodReader::resetCounters(){
	vcfParsingStarted = false;
	_lineCounter = 0;
	_notBialleleicCounter = 0;
	_missingSNPCounter = 0;
	_lowFreqSNPCounter = 0;
	_lowVariantQualityCounter = 0;
	_noPLCounter = 0;
	_numAcceptedLoci = 0;
};

void TPopulationLikelihoodReader::openVCF(std::string vcfFilename, TLog* logfile){
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

bool TPopulationLikelihoodReader::readDataFromVCF(unsigned short* curLocus, TPopulationSamples & samples, TLog* logfile){
	//set time at beginning
	if(!vcfParsingStarted){
		vcfParsingStarted = true;
		gettimeofday(&startTime, NULL);
	}

	//read next
	while(vcfFile.next()){ // new line in vcf-file (= new locus)
		++_lineCounter;

		//print progress
		if(_lineCounter % progressFrequency == 0)
			printProgressFrequencyFiltering(logfile);

		// limit lines
		if(limitLines > 0 && _lineCounter > limitLines){
			logfile->list("Reached limit of " + toString(limitLines) + " lines.");
			return false;
		}

        //skip sites with != 2 alleles
        if(vcfFile.getNumAlleles() != 2){
        	_notBialleleicCounter++;
        	continue;
        }

        //skip sites with too low variant quality
        if(minVariantQuality > 0 && (vcfFile.variantQualityIsMissing() || vcfFile.variantQuality() < minVariantQuality)){
        	_lowVariantQualityCounter++;
        	continue;
        }

		//check if PL is given
		if(!vcfFile.formatColExists("PL")){
			++_noPLCounter;
			continue;
		}

		// create an array containing the genotype likelihoods of all considered individuals of current locus
        long numIndividualsWithMissingSNP = 0;
		for(int s = 0; s < samples.numSamples(); ++s){
			int vcfIndex = samples.VCF_order(s);

			// depth filter: if a locus has < minDepth reads, flag locus as missing (set all genotype likelihoods = 1)
			if (vcfFile.sampleDepth(vcfIndex) < minDepth){
				vcfFile.setSampleMissing(vcfIndex);
				numIndividualsWithMissingSNP++;
			}
			vcfFile.fillPhredScore(vcfIndex, &curLocus[3 * s]);
		}

		// missingness filter: if > percentMissingPerLocus of individuals per locus have are missing, remove locus
		if ( (double) numIndividualsWithMissingSNP / (double) samples.numSamples() > maxMissing){
			_missingSNPCounter++;
			continue; // skip rest of loop (don't store)
		}

		//filter in MAF
		if(freqFilter > 0.0){
			// estimate allele frequency (EM algorithm)
			estimateGenotypeFrequenciesNullModel(curLocus, samples.numSamples(), epsilonF);
			double f = _genotypeFrequencies[0] + 0.5 * _genotypeFrequencies[1];
			if(f > 0.5) f = 1.0 - f;

			if(f < freqFilter){
				_lowFreqSNPCounter++;
				continue;
			}
		}

		//SNP is accepted!
		++_numAcceptedLoci;
		return true;
    }

	//return false at end of file
	logfile->list("Reached end of VCF file.");
	return false;
};

void TPopulationLikelihoodReader::printProgressFrequencyFiltering(TLog* logfile){
	struct timeval end;
	gettimeofday(&end, NULL);
	float runtime = (end.tv_sec  - startTime.tv_sec)/60.0;
	logfile->list("Parsed " + toString(_lineCounter) + " lines and retained " + toString(_numAcceptedLoci) + " loci in " + toString(runtime) + " min");
};

void TPopulationLikelihoodReader::concludeFilters(TLog* logfile){
	printProgressFrequencyFiltering(logfile);
	if(_notBialleleicCounter > 0)
		logfile->conclude(toString(_notBialleleicCounter) + " loci were not bi-allelic.");
	if(_lowVariantQualityCounter > 0)
		logfile->conclude(toString(_lowVariantQualityCounter) + " loci had variant quality < " + toString(minVariantQuality) + ".");
	if(_noPLCounter > 0)
		logfile->conclude(toString(_noPLCounter) + " loci had no PL field.");
	if(_missingSNPCounter > 0)
		logfile->conclude(toString(_missingSNPCounter) + " loci had > " + toString(maxMissing) + " samples with missing data.");
	if(_lowFreqSNPCounter > 0)
		logfile->conclude(toString(_lowFreqSNPCounter) + " loci had MAF < " + toString(freqFilter) + ".");
};

void TPopulationLikelihoodReader::guessGenotypeFrequencies(unsigned short* phredScores, const int & numSamples){
	//calculate by using MLE genotype for each individual
	_genotypeFrequencies[0] = 0.0; _genotypeFrequencies[1] = 0.0; _genotypeFrequencies[2] = 0.0;
	for(int i = 0 ; i < 3 * numSamples; i += 3){
		if(phredScores[i + 1] < phredScores[i]){
			if(phredScores[i + 2] < phredScores[i + 1]) _genotypeFrequencies[2] += 1.0;
			else _genotypeFrequencies[1] += 1.0;
		} else {
			if(phredScores[i + 2] < phredScores[i]) _genotypeFrequencies[2] += 1.0;
			else _genotypeFrequencies[0] += 1.0;
		}
	}

	double sum = 0.0;
	for(int g = 0; g < 3; ++g){
		_genotypeFrequencies[g] /= (double) numSamples;
		if(_genotypeFrequencies[g] <= 0.0) _genotypeFrequencies[g] = 0.01;
		if(_genotypeFrequencies[g] >= 1.0) _genotypeFrequencies[g] = 0.99;
		sum += _genotypeFrequencies[g];
	}
	for(int g = 0; g < 3; ++g){
		_genotypeFrequencies[g] /= sum;
	}
}

void TPopulationLikelihoodReader::estimateGenotypeFrequenciesNullModel(unsigned short* phredScores, const int & numSamples, double epsilonF){
	//prepare variables
	double sum;
	double weightsNull[3];
	double genoFreq_old[3];

	//estimate initial frequencies from MLEs
	guessGenotypeFrequencies(phredScores, numSamples);

	//run EM for max 1000 steps
	for (int s = 0; s < 1000; ++s){
		//set genofreq and calc P(g|f)
		for(int g = 0; g < 3; ++g){
			genoFreq_old[g] = _genotypeFrequencies[g];
			_genotypeFrequencies[g] = 0.0;
		}

		//estimate new genotype frequencies
		for(int i = 0; i < 3 * numSamples; i += 3){
			sum = 0.0;
			for(int g = 0; g < 3; ++g){
				weightsNull[g] = phredToGTLMap[phredScores[i + g]] * genoFreq_old[g];
				sum += weightsNull[g];
			}
			_genotypeFrequencies[0] += weightsNull[0] / sum;
			_genotypeFrequencies[2] += weightsNull[2] / sum;
		}

		_genotypeFrequencies[0] /= (double) numSamples;
		_genotypeFrequencies[2] /= (double) numSamples;
		_genotypeFrequencies[1] = 1.0 - _genotypeFrequencies[0] - _genotypeFrequencies[2];

		//check if we stop
		if(fabs(_genotypeFrequencies[0] - genoFreq_old[0]) < epsilonF && fabs(_genotypeFrequencies[2] - genoFreq_old[2]) < epsilonF) break;
	}

	//now set allele frequencies
	_alleleFrequency = _genotypeFrequencies[0] + 0.5 * _genotypeFrequencies[1];
	if(_alleleFrequency > 0.5) _MAF = 1.0 - _alleleFrequency;
	else _MAF = _alleleFrequency;
};

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
};

void TPopulationLikelihoods::clean(){
	for(std::vector<unsigned short * >::iterator it=genotypePhredScores.begin(); it < genotypePhredScores.end(); it++)
		delete[] *it;
	vcfRead = false;
};

void TPopulationLikelihoods::readData(TParameters & Parameters, TLog* Logfile){
	//check if we limit samples
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameterString("samples"), Logfile);

	//read Data
	readDataFromVCF(Parameters, Logfile);
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read data from VCF-file                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////
void TPopulationLikelihoods::readDataFromVCF(TParameters & Parameters, TLog* logfile){
	if(vcfRead)
		throw "VCF already read!";

	//create reader
	TPopulationLikelihoodReader reader(Parameters, logfile);
	if(saveAlleleFrequencies)
		reader.doEstimateGenotypeFrequencies();

	// open vcf file
	std::string vcfFilename = Parameters.getParameterString("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename, logfile);

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());


	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
    unsigned short* curLocus = new unsigned short[samples.numSamples() * 3];

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    _numLoci = 0;
    std::string curChr = "";
    while(reader.readDataFromVCF(curLocus, samples, logfile)){
		//update chromosome name
		if(reader.chr() != curChr){
			chromosomes.emplace(_numLoci, reader.chr());
			curChr = reader.chr();
		}

		//store SNP
		genotypePhredScores.emplace_back(curLocus);
		position.emplace_back(reader.position());
		if(saveAlleleFrequencies)
			alleleFrequencies.emplace_back(reader.allelFrequency());

		//update for next
		curLocus = new unsigned short[samples.numSamples() * 3];
		++_numLoci;
    }

    //clean up
	delete[] curLocus;
	vcfRead = true;

    //report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
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

unsigned short * TPopulationLikelihoods::getDataAtLocus(long index){
	return genotypePhredScores[index];
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// functions to loop over data                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////

void TPopulationLikelihoods::begin(){
	curLocusIndex = 0;
	curChrIt = chromosomes.begin();
	nextChrIt = std::next(chromosomes.begin(),1);
	individualStartIndex = 0;
	usedSampleSize = samples.numSamples();
};

void TPopulationLikelihoods::beginOnePop(int population){
	curLocusIndex = 0;
	curChrIt = chromosomes.begin();
	nextChrIt = std::next(chromosomes.begin(),1);
	individualStartIndex = samples.startIndex(population);
	usedSampleSize = samples.numSamplesInPop(population);
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

unsigned short* TPopulationLikelihoods::curData(){
	return &genotypePhredScores[curLocusIndex][individualStartIndex];
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

	for(std::map<int, std::string>::iterator it = chromosomes.begin(); it != chromosomes.end(); ++it){
		std::cout << "CHR '" << it->second << " starts at " << it->first << std::endl;
	}


	for(begin(); !end(); next()){
		//print chromosome and position
		std::cout << curChr() << "\t" << curPosition();

		//print data
		unsigned short* data = curData();
		for(int s=0; s<curSampleSize(); ++s){
			int index = 3*s;
			std::cout << "\t" << toString(data[index]) << "," << toString(data[index+1]) << "," << toString(data[index+2]);
		}

		std::cout << std::endl;
	}
};

