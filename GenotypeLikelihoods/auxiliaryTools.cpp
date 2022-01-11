/*
 * TRecalibrationEMAuxiliaryTools.cpp
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */


#include "auxiliaryTools.h"

namespace GenotypeLikelihoods{


//--------------------------------------------------------------------
// TRecalibrationEMFirstDerivative
//--------------------------------------------------------------------
TRecalibrationEMFirstDerivatives::TRecalibrationEMFirstDerivatives(size_t Size){
	resize(Size);
};

void TRecalibrationEMFirstDerivatives::resize(size_t Size){
	_derivatives.resize(Size);
};

size_t TRecalibrationEMFirstDerivatives::size() const{
	return _derivatives.size();
};

void  TRecalibrationEMFirstDerivatives::restart(){
	_cur = _derivatives.begin();
};

void TRecalibrationEMFirstDerivatives::add(uint16_t parameterIndex, double derivative){
	_cur->index = parameterIndex;
	_cur->derivative = derivative;
	++_cur;
};

TRecalibrationEMFirstDerivativesIterator TRecalibrationEMFirstDerivatives::begin(){
	return _derivatives.begin();
};

TRecalibrationEMFirstDerivativesIterator TRecalibrationEMFirstDerivatives::end(){
	return _derivatives.end();
};

//--------------------------------------------------------------------
// TRecalibrationEMSecondDerivatives
//--------------------------------------------------------------------
TRecalibrationEMSecondDerivatives::TRecalibrationEMSecondDerivatives(size_t Size){
	resize(Size);
};

void TRecalibrationEMSecondDerivatives::resize(size_t Size){
	_derivatives.resize(Size);
};

size_t TRecalibrationEMSecondDerivatives::size() const{
	return _derivatives.size();
};

void  TRecalibrationEMSecondDerivatives::restart(){
	_cur = _derivatives.begin();
};

void TRecalibrationEMSecondDerivatives::add(uint16_t parameterIndex1, uint16_t parameterIndex2, double derivative){
	_cur->index1 = parameterIndex1;
	_cur->index2 = parameterIndex2;
	_cur->derivative = derivative;
	++_cur;
};

TRecalibrationEMSecondDerivativesIterator TRecalibrationEMSecondDerivatives::begin(){
	return _derivatives.begin();
};

TRecalibrationEMSecondDerivativesIterator TRecalibrationEMSecondDerivatives::end(){
	return _derivatives.end();
};

}; //end namespace
