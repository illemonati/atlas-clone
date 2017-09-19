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

	//output name
	std::string outname = params.getParameterStringWithDefault("out", "ATLAS_simulations");
	logfile->list("Will write output files with tag '" + outname + "'.");

	//initialize simulator
	TSimulator simulator(logfile, randomGenerator, params);

	//read other parameters
	int sampleSize = params.getParameterIntWithDefault("sampleSize", 1);
	if(sampleSize < 1)
		throw "Sample size needs to be at least 1!";

	//simulate differently depending on number of individuals
	if(sampleSize == 1){
		logfile->startIndent("Simulating a single individual:");

		std::vector<std::string> tmp;
		params.fillParameterIntoVectorWithDefault("theta", tmp, ',', "0.001");
		std::vector<double> thetas;
		repeatIndexes(tmp, thetas);
		if(thetas.size() == 1){
			logfile->list("Will simulate data with theta = " + toString(thetas[0]) + ".");
			simulator.simulateSingleIndividual(thetas[0], outname);
		} else {
			logfile->list("Will simulate data with chromosome specific thetas " + concatenateString(thetas, ", "));
			simulator.simulateSingleIndividual(thetas, outname);
		}
	} else if(sampleSize == 2 && params.parameterExists("phi")){
		logfile->startIndent("Simulating two individuals with predefined genetic distance:");
		//simulate according to genetic distance
		//parse distance
		std::vector<std::string> tmp;
		std::vector<double> phis;
		params.fillParameterIntoVector("phi", tmp, ',');
		repeatIndexes(tmp, phis);

		//now run simulations
		simulator.simulateIndividualPair(phis, outname);
	} else {
		logfile->startIndent("Simulating multiple individuals based on an SFS:");
		//prepare SFS
		std::vector<SFS*> sfs;
		if(params.parameterExists("sfs")){
			logfile->startIndent("Reading SFS from files:");

			std::vector<std::string> tmp;
			std::vector<std::string> sfsFileNames;
			params.fillParameterIntoVector("sfs", tmp, ',');
			repeatIndexes(tmp, sfsFileNames);
			bool folded = params.parameterExists("folded");

			simulator.simulatePopulationFromSFS(sfsFileNames, folded, sampleSize, outname);
		} else if(params.parameterExists("theta")){
			//parse theta from command line
			std::vector<std::string> tmp;
			params.fillParameterIntoVector("theta", tmp, ',');
			std::vector<double> thetas;
			repeatIndexes(tmp, thetas);

			if(thetas.size() == 1){
				logfile->list("Will simulate from SFS with theta = " + toString(thetas[0]) + ".");
				simulator.simulatePopulationFromSFS(thetas[0], sampleSize, outname);
			} else {
				logfile->list("Will simulate data from chromosome specific SFS with thetas " + concatenateString(thetas, ", "));
				simulator.simulatePopulationFromSFS(thetas, sampleSize, outname);
			}
			logfile->endIndent();
		}
	}
	logfile->endIndent();

	//clean up
	delete randomGenerator;
}
