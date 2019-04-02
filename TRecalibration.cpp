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

	models->createModels(string, *ReadGroups);
};


