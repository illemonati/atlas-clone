/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */


#include "TRecalibrationEMModule.h"

void TRecalibrationEMModule::_freeBetas(){
	if(initialized()){
		delete[] _betas;
		delete[] _oldBetas;
	}
};

void TRecalibrationEMModule::_initializeBeta(){
	_freeBetas();
	_betas = new double[_numParameters];
	_oldBetas = new double[_numParameters];
};


