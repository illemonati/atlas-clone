/*
 * TRecalibration.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#include "TRecalibration.h"

//---------------------------------------------------------------
//TRecalibration
//---------------------------------------------------------------

TRecalibration::TRecalibration(){
	_type = "none";
};

double TRecalibration::getErrorRate(TBase & base){
	return base.errorRate;
};

int TRecalibration::getQuality(TBase & base){
	return _qualityMap.errorToQuality(base.errorRate);
};

//---------------------------------------------------------------
//TRecalibrationEM
//---------------------------------------------------------------
TRecalibrationEM::TRecalibrationEM(std::string string, TReadGroups* ReadGroups, TLog* Logfile):TRecalibration(){
	logfile = Logfile;
	_type = "recal";

	//read groups
	readGroups = ReadGroups;

	//models
	models = new TRecalibrationEMModels(readGroups->size(), logfile);

	initializeRecalibrationParameters(string);
};

void TRecalibrationEM::initializeRecalibrationParameters(std::string string){
	//initialize from string or file
	int pos = string.find_first_of('[');
	if(pos == std::string::npos)
		_initializeRecalibrationParametersFromFile(string);
	else
		_initializeRecalibrationParametersFromString(string);
};

void TRecalibrationEM::_initializeRecalibrationParametersFromString(std::string & string){
	//string has format model[param1, param2, param3, ...]
	logfile->startIndent("Initializing recal with string '" + string + "' for all read groups:");

	//read model tag
	size_t pos = string.find_first_of('[');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing '['!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::string modelTag = string.substr(0,pos);
	string.erase(0, pos+1);

	//read parameters: quality, position and context separted by semicolon (;)
	pos = string.find_first_of(']');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing ']'!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::vector<std::string> tmpVec;
	fillVectorFromString(string.substr(0, pos), tmpVec, ";");
	if(tmpVec.size() != 3)
		throw "Failed to understand recal string: wrong number of parameter sets (" + toString(tmpVec.size()) + " instead of 3)!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";

	//initialize model
	models->addSingleModelForAllReadGroups(modelTag, tmpVec, true);

	logfile->endIndent();
};

void TRecalibrationEM::_initializeRecalibrationParametersFromFile(std::string filename){
	//read parameters from file
	logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "' for reading!";

	//skip header
	std::string line;
	std::getline(file, line);

	//tmp variables for reading
	int lineNum = 0;
	std::vector<std::string> vec;

	//parse file to read details for each read group
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 6)
				throw "Wrong number of entries in file '" + filename + "' on line " + toString(lineNum) + ": expected 6 but found " + toString(vec.size()) + "!";

			//check if read group exists
			if(readGroups->readGroupExists(vec[0])){
				//read read group, mate and model
				int rg = readGroups->find(vec[0]);
				bool isSecondMate;
				if(vec[1] == "second"){
					isSecondMate = true;
				} else if(vec[1] == "first"){
					isSecondMate = false;
				} else {
					throw "unknown second mate information!";
				}

				std::string modelTag = vec[2];

				//clean up vec to only contain parameters (remove read group, mate, model and LL)
				vec.erase(vec.begin(), vec.begin() + 3);

				//create model
				models->addModel(rg, isSecondMate, modelTag, vec, false);
			} else {
				logfile->warning("Read group '" + vec[0] + "' does not exist in the BAM header! Are you using the correct recal file?");
			}
		}
	}
	logfile->done();

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	if(models->hasReadGroupsWithoutModel()){
		logfile->warning("Missing read groups in file '" + filename + "'!");
		logfile->startIndent("Will assume no recalibration for the following read groups:");
		models->reportReadGroupsNotUsed(*readGroups);
		models->addNoRecalModelIfMissing();
		logfile->endIndent();
	}
};

