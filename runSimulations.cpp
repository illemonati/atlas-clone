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
	int sampleSize = params.getParameterInt("sampleSize");
	logfile->list("Will simulate data from " + toString(sampleSize) + " individuals.");
	if(sampleSize < 1)
		throw "Sample size needs to be at least 1!";
	float theta = params.getParameterDouble("theta");
	logfile->list("Will simulate data with theta = " + toString(theta) + ".");

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

	//read length & depth
	float depth = params.getParameterDoubleWithDefault("depth", 10.0);
	simulator.setDepth(depth);
	int readLength = params.getParameterIntWithDefault("readLength", 100);
	simulator.setReadLength(readLength);
	logfile->list("Will simulate reads of length " + toString(readLength) + " to a depth of " + toString(depth) + ".");

	//chromosomes
	int numChr = params.getParameterIntWithDefault("numChr", 1);
	long chrLength = params.getParameterLongWithDefault("chrLength", 1000000);
	logfile->list("Will simulate " + toString(numChr) + " chromosome(s) of length " + toString(chrLength) + " each.");
	simulator.initializeChromosomes(numChr, chrLength);

	//quality distribution
	double meanQual = params.getParameterDoubleWithDefault("meanQual", 30);
	double sdQual = params.getParameterDoubleWithDefault("meanQual", 10);
	logfile->list("Will simulate normal distributed quality scores with mean = " + toString(meanQual) + " and sd = " + toString(sdQual));
	simulator.setQualityDistribution(meanQual, sdQual);
	logfile->endIndent();

	//simulate differently depending on number of individuals
	if(sampleSize == 1){
		simulator.simulateSingleIndividual(theta, outname);
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
