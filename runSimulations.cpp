/*
 * runSimulations.cpp
 *
 *  Created on: May 19, 2017
 *      Author: wegmannd
 */

#include "runSimulations.h"

void runSimulations(TParameters & params, TLog* logfile){
	//initialize simulator
	TSimulator* simulator;
	std::string method = params.getParameterStringWithDefault("type", "one");
	if(method == "one")
		simulator = new TSimulatorOneIndividual(logfile, params);
	else if(method == "pair")
		simulator = new TSimulatorPairOfIndividuals(logfile, params);
	else if(method == "SFS")
		simulator = new TSimulatorSFS(logfile, params);
	else if(method == "HW")
		simulator = new TSimulatorHardyWeinberg(logfile, params);
	else throw "Unknown simulation method '" + method + "'!";

	//now run simulations
	simulator->runSimulations();

	//clean up
	delete simulator;
};
