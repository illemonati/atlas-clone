/*
 * runSimulations.cpp
 *
 *  Created on: May 19, 2017
 *      Author: wegmannd
 */

#include "runSimulations.h"

void runSimulations(TParameters & params, TLog* logfile){
	//initialize random generator
	TRandomGenerator* randomGenerator;
	logfile->listFlush("Initializing random generator ...");
	if(params.parameterExists("fixedSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator=new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");

	//initialize simulator
	TSimulator simulator(logfile, randomGenerator);

	//read parameters
	logfile->startIndent("Reading simulation parameters:");
	std::string outname = params.getParameterStringWithDefault("out", "ATLAS_simulations");
	logfile->list("Will write output files with tag '" + outname + "'.");
	int sampleSize = params.getParameterIntWithDefault("sampleSize", 1);
	logfile->list("Will simulate data from " + toString(sampleSize) + " individuals.");
	if(sampleSize < 1)
		throw "Sample size needs to be at least 1!";
	double theta = params.getParameterDoubleWithDefault("theta", 0.001);
	logfile->list("Will simulate data with theta = " + toString(theta) + ".");
	double referenceDivergence = params.getParameterDoubleWithDefault("refDiv", 0.01);
	logfile->list("Will simulate data with theta = " + toString(theta) + ".");

	//read length & depth
	float depth = params.getParameterDoubleWithDefault("depth", 10.0);
	simulator.setDepth(depth);
	int readLength = params.getParameterIntWithDefault("readLength", 100);
	simulator.setReadLength(readLength);
	logfile->list("Will simulate reads of length " + toString(readLength) + " to a depth of " + toString(depth) + ".");

	//chromosomes
	std::vector<std::string> string_vec;
	std::vector<long> chrLength;
	params.fillParameterIntoVectorWithDefault("chrLength", string_vec, ',', "1000000");
	repeatIndexes(string_vec, chrLength);
	std::vector<int> ploidy;
	params.fillParameterIntoVectorWithDefault("ploidy", string_vec, ',', "2");
	repeatIndexes(string_vec, ploidy);
	if(ploidy.size() != chrLength.size())
		throw "List fo chromosome lengths and ploidies differ in length!";
	std::vector<bool> haploid;
	for(std::vector<int>::iterator it=ploidy.begin(); it!=ploidy.end(); ++it){
		if(*it == 1) haploid.push_back(true);
		else if(*it == 2) haploid.push_back(false);
		else throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!";
	}

	if(chrLength.size() < 1)
		throw "Issue understanding length of chromosomes!";
	if(chrLength.size() == 1){
		int numChr = params.getParameterIntWithDefault("numChr", 1);
		std::string text = "Will simulate ";
		if(haploid[0]) text += "haploid ";
		else text += "diploid ";
		text += toString(numChr) + " chromosome(s) of length " + toString(chrLength[0]) + " each.";
		logfile->list(text);
		simulator.initializeChromosomes(numChr, chrLength[0], haploid[0]);
	} else {
		logfile->startIndent("Will simulate " + toString(chrLength.size()) + " chromosome(s) of the following length:");
		std::vector<bool>::iterator hIt=haploid.begin();
		std::string text;
		for(std::vector<long>::iterator it=chrLength.begin(); it!=chrLength.end(); ++it, ++hIt){
			text = toString(*it) + " (";
			if(*hIt) text += "haploid)";
			else text += "diploid)";
			logfile->list(text);
		}
		simulator.initializeChromosomes(chrLength, haploid);
	}

	//quality distribution
	double meanQual = params.getParameterDoubleWithDefault("meanQual", 30);
	double sdQual = params.getParameterDoubleWithDefault("sdQual", 10);
	logfile->list("Will simulate normal distributed quality scores with mean = " + toString(meanQual) + " and sd = " + toString(sdQual));
	simulator.setQualityDistribution(meanQual, sdQual);
	logfile->endIndent();

	//quality transformation
	if(params.parameterExists("qualTransform")){
		params.fillParameterIntoVector("qualTransform", string_vec, ',');
		std::vector<double> beta;
		repeatIndexes(string_vec, beta);
		if(beta.size() != 24)
			throw "Wrong number of beta values for quality transformation (" + toString(beta.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";
		simulator.setQualityTransformation(beta);
	} else if(params.parameterExists("recal")){
		std::string filename = params.getParameterString("recal");
		logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");
		std::ifstream file(filename.c_str());
		if(!file) throw "Failed to open file '" + filename + "' for reading!";

		//tmp variables for reading
		std::string tmp;
		std::vector<std::string> vec;
		std::vector<double> beta;

		//skip header
		std::getline(file, tmp);

		//parse file to read details for each read group
		std::getline(file, tmp);

		fillVectorFromString(tmp, vec, "\t");
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 25) throw "Found " + toString(vec.size()) + " instead of 25 columns in '" + filename + "' on line 2!";
			for(int i=1; i < 25; ++i) beta.push_back(stringToFloat(vec[i]));
			logfile->done();
			if(beta.size() != 24)
				throw "Wrong number of beta values for quality transformation (" + toString(beta.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";
			simulator.setQualityTransformation(beta);
		}
	}

	//initialize PMD
	TPMD pmdObject;
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		if(params.parameterExists("pmd")){
			std::string pmdString = params.getParameterString("pmd");
			logfile->list("Initializing PMD for both C->T and G->A with function '" + pmdString +"'.");
			pmdObject.initializeFunction(pmdString, pmdGA, logfile);
			pmdObject.initializeFunction(pmdString, pmdCT, logfile);
			logfile->conclude(pmdObject.getFunctionString(pmdCT));
			if(params.parameterExists("pmdCT")) logfile->warning("Ignoring argument 'pmdCT'!");
			if(params.parameterExists("pmdGA")) logfile->warning("Ignoring argument 'pmdGA'!");
		} else {
			//first C->T
			if(!params.parameterExists("pmdCT")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdCT' has to be provided!";
			std::string pmdStringCT = params.getParameterString("pmdCT");
			logfile->list("Initializing post mortem C->T damage with function '" + pmdStringCT +"'.");
			pmdObject.initializeFunction(pmdStringCT, pmdCT, logfile);
			logfile->conclude(pmdObject.getFunctionString(pmdCT));

			//second G->A
			if(!params.parameterExists("pmdGA")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdGA' has to be provided!";
			std::string pmdStringGA = params.getParameterString("pmdGA");
			logfile->list("Initializing post mortem G->A damage with function '" + pmdStringGA +"'.");
			pmdObject.initializeFunction(pmdStringGA, pmdGA, logfile);
			logfile->conclude(pmdObject.getFunctionString(pmdGA));
		}
		simulator.setPMD(&pmdObject);
	}

	//simulate differently depending on number of individuals
	if(sampleSize == 1){
		simulator.simulateSingleIndividual(theta, referenceDivergence, outname);
	} else {
		//prepare SFS
		logfile->startIndent("Preparing SFS:");
		float theta = params.getParameterDouble("theta");
		logfile->listFlush("Generating SFS for " + toString(sampleSize) + " chromosomes and with theta = " + toString(theta) + " ...");
		SFS sfs(sampleSize, theta);
		logfile->write(" done!");
		std::string filename = outname + "_trueSFS.txt";
		logfile->listFlush("Writing true SFS to '" + filename + "' ...");
		sfs.writeToFile(filename);
		logfile->write(" done!");
		logfile->endIndent();

		throw "NOT YET IMPLEMENTED!!!";
	}

	//clean up
	delete randomGenerator;
}
