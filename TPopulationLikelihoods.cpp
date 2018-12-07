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
	delete numSamplesPerPop;
};

void TPopulationSamples::init(){
	_numPopulations = 0;
	_numSamples = 0;
	numSamplesPerPop = NULL;
	_hasSamples = false;
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
			if(lineNum == 1 && vec.size == 2) hasPopColumn = true;

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
			if(samples.find(vec[0]))
				throw "Duplicate sample name '" + vec[0] + "' in file '" + filename + "' on line " + toString(lineNum) + "!";
			samples.emplace(vec[0], it->second);
		}
	}

	//close file
	file.close();

	//fill counters
	_numSamples = samples.size();
	if(_numSamples > 0){
		numSamplesPerPop = new int[_numPopulations];

		//put samples in order by populations
		int nextIndex = 0;
		for(int p=0; p<_numPopulations; ++p){
			numSamplesPerPop[p] = 0;
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

void TPopulationSamples::fillVCFOrder(int* sampleVcfIndex, std::vector<std::string> & vcfSampleNames){
	sampleVcfIndex = new int[_numSamples];
	bool sampleUsed[_numSamples];
	for(int s=0; s<_numSamples; ++s)
		sampleUsed[s] = false;

	int vcfIndex = 0;
	for(std::vector<std::string>::iterator it = vcfSampleNames.begin(); it != vcfSampleNames.end(); ++it, ++vcfIndex){
		if(sampleIsUsed(*it)){
			int orderedIndex = getOrderedSampleIndex(*it);
			if(sampleUsed[orderedIndex])
				throw "Duplicate sample name '" + *it + "' in VCf header!";
			sampleVcfIndex[orderedIndex] = vcfIndex;
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
// TPopulationLikelihoods                                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////

TPopulationLikelihoods::TPopulationLikelihoods(TParameters & Parameters, TLog* Logfile){
	logfile = Logfile;

	vcfRead = false;
	relevantNumIndivs = 0;
	numAcceptedLoci = 0;

	//check if we limit samples
	TPopulationSamples samples;
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameterString("samples"), logfile);

	//read Data
	readDataFromVCF(TParameters & Parameters, samples);

};

void TPopulationLikelihoods::clean(){
	for(std::vector<unsigned short * >::iterator it=genotypePhredScores.begin(); it < genotypePhredScores.end(); it++)
		delete[] *it;
};


//////////////////////////////////////////////////////////////////////////////////////////////////
// Read data from VCF-file                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////

void TPopulationLikelihoods::openVCF(std::string vcfFilename, TVcfFileSingleLine & vcfFile){
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
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();
	vcfFile.enableVariantParsing();
};

void TPopulationLikelihoods::getParametersForReadingVCF(TParameters & Parameters, int & limitLines, int & minDepth, double & maxMissing, double & epsilonF, double & freqFilter, int & progressFrequency){

}

void TPopulationLikelihoods::readDataFromVCF(TParameters & Parameters, TPopulationSamples & samples){
	if(vcfRead)
		throw "VCF already read!";

	// open vcf file
	TVcfFileSingleLine vcfFile;
	std::string vcfFilename = Parameters.getParameterString("vcf");
	openVCF(vcfFilename, vcfFile);

	//Match samples
	int* sampleVcfOrder;
	if(samples.hasSamples()){
		samples.fillVCFOrder(sampleVcfOrder, vcfFile.parser.samples);
		numSamples = samples.numSamples();
	} else {
		numSamples = vcfFile.numSamples();
		sampleVcfOrder = new int[numSamples];
		for(int s=0; s<numSamples; ++s)
			sampleVcfOrder[s] = s;
	}

	//read parsing parameters
	// do we limit the lines to read?
	long limitLines = Parameters.getParameterLongWithDefault("limitLines", -1);
	if(limitLines > 0)
	logfile->list("Will limit analysis to the first " + toString(limitLines) + " lines of the VCF file.");

	// do we set a depth filter?
	int minDepth = Parameters.getParameterIntWithDefault("minDepth", 0);
	logfile->list("Will filter samples to a minimum depth of " + toString(minDepth) + ".");

	// do we set a missingness filter?
	double maxMissing = Parameters.getParameterIntWithDefault("maxMissing", 1.0);
	if(maxMissing < 0.0 || maxMissing > 1.0)
		throw "maxMissing must be within (0, 1)!";
	logfile->list("Will remove loci where more than " + toString(maxMissing) + " of samples are missing.");

	// parameters to set a filter on the allele frequency?
	bool doFreqFilter = false;
	double freqFilter = Parameters.getParameterDoubleWithDefault("minMAF", 0.0); // MAF = minor allele frequency
	if(freqFilter < 0.0 || freqFilter >= 0.5)
		throw "MAF filter must be within (0.0,0.5)!";
	double epsilonF = 0.0;
	if(freqFilter > 0.0){
		doFreqFilter = true;
		epsilonF = Parameters.getParameterDoubleWithDefault("epsF", 0.0001);
		logfile->list("Will filter on an allele frequency of " + toString(freqFilter) + ".");
	}

	//filter on variant quality?
	bool doVariantQualityFilter = false;
	int minVariantQuality = Parameters.getParameterIntWithDefault("minVariantQuality", 0);
	if(minVariantQuality < 0) throw "minVariantQuality must be >= 0!";
	if(minVariantQuality > 0){
		doVariantQualityFilter = true;
		logfile->list("Will only keep sites with variant quality >= " + toString(minVariantQuality) + ".");
		vcfFile.enableInfoParsing();
	}

	// progress frequency
	int progressFrequency = Parameters.getParameterIntWithDefault("reportFreq", 10000);


	// initialize variables for vcf-file
	long lineCounter = 0;
	long validSNPCounter = 0;
	long numIndividualsWithMissingSNP = 0;
	long missingSNPCounter = 0;
	long lowFreqSNPCounter = 0;
	std::string curChrName = "";
	struct timeval start;
    unsigned short* curLocus;

	// initialize variables for EM
	double* genotypeFrequencies = new double[3];
	double f = 0;

	// debug
	// double numberloci = parameters->getParameterInt("L");
	// Fil.set_size(numIndividualsFromVcf, numberloci);

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(vcfFile.next()){ // new line in vcf-file (= new locus)
    	numIndividualsWithMissingSNP = 0;
    	++lineCounter;

        //skip sites with != 2 alleles
        if(vcfFile.getNumAlleles() != 2) continue;
        ++validSNPCounter;

        //skip sites with too low variant quality
        if(vcfFile.parser.

		// update chromosome name
		if(vcfFile.chr() != curChrName)
			curChrName = vcfFile.chr();

		// create an array containing the genotype likelihoods of all considered individuals of current locus
		curLocus = new unsigned short[numSamples * 3];
		for(int s = 0; s < numSamples; ++s){
			int vcfIndex = sampleVcfOrder[s];

			// depth filter: if a locus has < minDepth reads, flag locus as missing (set all genotype likelihoods = 1)
			if (vcfFile.sampleDepth(vcfIndex) < minDepth){
				vcfFile.setSampleMissing(vcfIndex);
				numIndividualsWithMissingSNP++;
			}
			vcfFile.fillPhredScore(vcfIndex, &curLocus[3 * s]);
		}

		// missingness filter: if > percentMissingPerLocus of individuals per locus have are missing, remove locus
		if ( (double) numIndividualsWithMissingSNP / (double) numSamples > maxMissing){
			missingSNPCounter += 1;
			continue; // skip rest of loop (don't store)
		}

		// estimate allele frequency (EM algorithm)

		estimateGenotypeFrequenciesNullModel(genotypeFrequencies, curLocus, epsilonF);
		f = genotypeFrequencies[0] + 0.5 * genotypeFrequencies[1];
		if(f > 0.5) f = 1.0 - f;

		// only store SNPs that pass the frequency filter
		if(f >= freqFilter){

			// store the genotype likelihoods of the current locus
			genotypePhredScores.push_back(curLocus);

			// add name of locus to loci
			loci->initializeNames(getChrPosString(curChrName, vcfFile.position()));

			++numAcceptedLoci;
		}
		else {
			lowFreqSNPCounter++;
		}

		//report progress
		if(lineCounter % progressFrequency == 0)
			printProgressFrequencyFiltering(lineCounter, numAcceptedLoci, missingSNPCounter, lowFreqSNPCounter, start);

        // limit lines
		if(limitLines > 0 && lineCounter > limitLines) break;
    }

	printProgressFrequencyFiltering(lineCounter, numAcceptedLoci, missingSNPCounter, lowFreqSNPCounter, start);
	loci->initializeStorage(numAcceptedLoci, numLatFac, numEnvVar);
	logfile->endIndent();
	logfile->endIndent();

    // clean up
	delete[] genotypeFrequencies;

	//set that vcf was read
	vcfRead = true;

	// std::cout << "Fil:\n"<< Fil << std::endl;

    // for debugging loci:

    // loci.printLoci();

    //for debugging genotype likelihoods:
	/*
	std::cout << "Size of genotype likelihoods vector: " << genotypeLikelihoods.size() << std::endl;
    for (std::vector< float *>::iterator it = genotypeLikelihoods.begin(); it < genotypeLikelihoods.end(); it++){
    	std::cout << "\nnew Locus " << std::endl;
    	for (int j = 0; j < relevantNumIndivs * 3; j++){
    		std::cout << (*it)[j] << "   ";
    	}
    }
    */
}

std::string TPopulationLikelihoods::getChrPosString(std::string chromosomeName, long snpPosition){
	return chromosomeName + "/" + toString(snpPosition);
}

void TPopulationLikelihoods::printProgressFrequencyFiltering(long lines, long & numRetainedLoci, long & numLociMissing, long & lowFreqSNPCounter, struct timeval & start){
	struct timeval end;
	gettimeofday(&end, NULL);
	float runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->startIndent("Parsed " + toString(lines) + " lines and retained " + toString(numRetainedLoci) + " loci in " + toString(runtime) + " min");
	logfile->list(toString(numLociMissing) + " loci were filtered out due to missingness.");
	logfile->list(toString(lowFreqSNPCounter) + " loci were filtered out due to low allele frequencies.");
	logfile->endIndent();
}

void TPopulationLikelihoods::matchIndividualNames(int numIndividualsFromVcf, std::vector<int> & orderOfIndividuals){
	// save names of individuals
	std::string sampleName;
    std::vector<std::string> individualNamesFromVcf;
    for (int individual = 0; individual < numIndividualsFromVcf; ++individual){
    	sampleName = vcfFile.sampleName(individual);
    	orderOfIndividuals.push_back(findPosOfVcfIndivInEnv(sampleName));
    	individualNamesFromVcf.push_back(sampleName); // contains all individual names of vcf, also those that are ignored during analysis
	    }
    // count the number of individuals present in vcf AND env file
    relevantNumIndivs = numIndividualsFromVcf - std::count(orderOfIndividuals.begin(), orderOfIndividuals.end(), -1);
    if (relevantNumIndivs < numIndividualsFromEnv){
    	findExtraIndivInEnvFile(individualNamesFromVcf);
    }

    logfile->list("A total number of " + toString(relevantNumIndivs) + " individuals will be considered for analysis.");
    // initialize K
    numLatFac = parameters->getParameterInt("numLatFac");
    if (numLatFac > relevantNumIndivs)
    	throw "Number of latent factors (" + toString(numLatFac) + ") can not exceed number of individuals (" + toString(relevantNumIndivs) + ")!";
    logfile->list("Will use K = " + toString(numLatFac) + " latent factors for analysis.");
}

int TPopulationLikelihoods::findPosOfVcfIndivInEnv(std::string individualNameVcf){
	// find the position of the .vcf-individual in the vector of .env-individual names
	int position = 0;
	int numberOfMatches = 0;
	int matchPos;
	for (std::vector<std::string>::iterator it = individualNamesFromEnv.begin(); it < individualNamesFromEnv.end(); it++, position++){
		if (individualNameVcf == *it){
			matchPos = position;
			numberOfMatches++;
		}
	}
	if (numberOfMatches == 0){
		logfile->list("Individual '" + individualNameVcf + "' from vcf-file does not match any individual from env-file - will ignore it.");
		matchPos = -1;
	}
	else if (numberOfMatches >= 2){
		throw "Individual '" + individualNameVcf + "' from vcf-file occurs " + toString(numberOfMatches) + "x in the env-file! Please change env-file accordingly.";
	}
	// std::cout << "matchPos: " << matchPos << std::endl;
	return matchPos;
}

void TPopulationLikelihoods::findExtraIndivInEnvFile(std::vector<std::string> & individualNamesFromVcf){
	// find the position of the .env-individual in the vector of .vcf-individual names
	for (std::vector<std::string>::iterator nameEnv = individualNamesFromEnv.begin(); nameEnv < individualNamesFromEnv.end(); nameEnv++){
		int numberOfMatches = 0;
		for (std::vector<std::string>::iterator nameVcf = individualNamesFromVcf.begin(); nameVcf < individualNamesFromVcf.end(); nameVcf++){
			if (*nameEnv == *nameVcf){
				numberOfMatches++;
			}
		}
		if (numberOfMatches == 0){
			throw "Individual '" + *nameEnv + "' from env-file does not match any individual from vcf-file! Please change env-file accordingly.";
		}
		else if (numberOfMatches >= 2){
			throw "Individual '" + *nameEnv + "' from env-file occurs " + toString(numberOfMatches) + "x in the vcf-file! Please change vcf-file accordingly.";
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Run EM algorithm to filter the loci on their frequency                                       //
//////////////////////////////////////////////////////////////////////////////////////////////////

void TPopulationLikelihoods::fillInitialEstimateOfGenotypeFrequencies(double* genoFreq, unsigned short* phredScores){
	//calculate by using MLE genotype for each individual
	genoFreq[0] = 0.0; genoFreq[1] = 0.0; genoFreq[2] = 0.0;
	for(int i = 0 ; i < 3 * relevantNumIndivs; i += 3){
		if(phredScores[i + 1] < phredScores[i]){
			if(phredScores[i + 2] < phredScores[i + 1]) genoFreq[2] += 1.0;
			else genoFreq[1] += 1.0;
		} else {
			if(phredScores[i + 2] < phredScores[i]) genoFreq[2] += 1.0;
			else genoFreq[0] += 1.0;
		}
	}

	double sum = 0.0;
	for(int g = 0; g < 3; ++g){
		genoFreq[g] /= (double) relevantNumIndivs;
		if(genoFreq[g] <= 0.0) genoFreq[g] = 0.01;
		if(genoFreq[g] >= 1.0) genoFreq[g] = 0.99;
		sum += genoFreq[g];
	}
	for(int g = 0; g < 3; ++g){
		genoFreq[g] /= sum;
	}
}

void TPopulationLikelihoods::estimateGenotypeFrequenciesNullModel(double* genotypeFrequencies, unsigned short* phredScores, double epsilonF){
	//prepare variables
	double sum;
	int i, g;
	double weightsNull[3];
	double genoFreq_old[3];

	//estimate initial frequencies from MLEs
	fillInitialEstimateOfGenotypeFrequencies(genotypeFrequencies, phredScores);

	//run EM for max 1000 steps
	for (int s = 0; s < 1000; ++s){
		//set genofreq and calc P(g|f)
		for(g = 0; g < 3; ++g){
			genoFreq_old[g] = genotypeFrequencies[g];
			genotypeFrequencies[g] = 0.0;
		}

		//estimate new genotype frequencies
		for(i = 0; i < 3 * relevantNumIndivs; i += 3){
			sum = 0.0;
			for(g = 0; g < 3; ++g){
				weightsNull[g] = phredToGTLMap[phredScores[i + g]] * genoFreq_old[g];
				sum += weightsNull[g];
			}
			genotypeFrequencies[0] += weightsNull[0] / sum;
			genotypeFrequencies[2] += weightsNull[2] / sum;
		}

		genotypeFrequencies[0] /= (double) relevantNumIndivs;
		genotypeFrequencies[2] /= (double) relevantNumIndivs;
		genotypeFrequencies[1] = 1.0 - genotypeFrequencies[0] - genotypeFrequencies[2];

		//check if we stop
		if(fabs(genotypeFrequencies[0] - genoFreq_old[0]) < epsilonF && fabs(genotypeFrequencies[2] - genoFreq_old[2]) < epsilonF) break;
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// get-functions                                                                                //
//////////////////////////////////////////////////////////////////////////////////////////////////

int TPopulationLikelihoods::getNumIndiv(){
	return relevantNumIndivs;
}

long TPopulationLikelihoods::getNumLoci(){
	return numAcceptedLoci;
}

int TPopulationLikelihoods::getNumEnvVar(){
	return numEnvVar;
}

int TPopulationLikelihoods::getNumLatFac(){
	return numLatFac;
}

TLoci* TPopulationLikelihoods::getLoci(){
	return loci;
}

TIndividuals* TPopulationLikelihoods::getIndividuals(){
	return individuals;
}

std::vector<std::string> TPopulationLikelihoods::getEnvVarNames(){
	return envVarNames;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// run all                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////

void TPopulationLikelihoods::readDataAndInitializeParams(){
	readDataFromEnv();
	readDataFromVCF();
}

// ./ldlg task=analyze vcf=simulations.vcf.gz envFile=simulations.env numLatFac=3 verbose
// cat simulations_simParameters.txt | head -21
// ./ldlg task=simulate numIndividuals=5 numLoci=5 numLatComp=3 withoutBX=true withoutUV=true error=0 verbose
// ./ldlg task=analyze vcf=simulations.vcf.gz envFile=simulations.env numLatFac=3 minMAF=0 verbose | grep y_l: | tr -s ' ' | cut -f2-9999 -d' ' > trueGenotypes.r
// echo "a<-read.table('trueGenotypes.r', header=FALSE); print(cov(a))" | R --slave




