/*
 * TGenotypeLikelihoods.cpp
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#include "TGenotypeLikelihoods.h"

namespace GenotypeLikelihoods{

TGenotypeLikelihoods::TGenotypeLikelihoods(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile){
	logfile = Logfile;
	readGroups = ReadGroups;
	readGroupMap = new TReadGroupMap(ReadGroups);


	//initialize PMD
	if(params.parameterExists("pmdFile")){
		std::string filename = params.getParameterString("pmdFile");
		logfile->listFlush("Initializing Post Mortem Damage (PMD) from file '" + filename + "' ...");
		pmd.initializeFromFile(*readGroups, filename);
		logfile->done();
	} else {
		logfile->list("Assuming there is no PMD in the data.");
	}

	//initialize recalibration
	errorModels = new TSequencingErrorModels(readGroups, readGroupMap, logfile);
	if(params.parameterExists("recal")){
		errorModels->createModels(params.getParameterString("recal"));
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
		errorModels->createEmptyModels();
	}


};



}; //end namespace

