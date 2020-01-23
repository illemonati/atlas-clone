/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */


#include "TRecalibrationEMModule.h"

void TRecalibrationEMModule::_freeBetas(){
	if(_initialized){
		delete[] _betas;
		delete[] _oldBetas;
	}
};

void TRecalibrationEMModule::_initializeBetas(){
	_freeBetas();
	_betas = new double[_numParameters];
	_oldBetas = new double[_numParameters];
};

void TRecalibrationEMModule::_initializValues(std::vector<std::string> & values){
	if(!values.empty()){
		if(values.size() != _numParameters){
			throw "Failed to initialize recalibration module: wrong number of values (" + toString(values.size()) + " instead of " + toString(_numParameters) + ")!"
		}
	}
};


