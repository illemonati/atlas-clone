/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMModel.h"

//---------------------------------------------------------------
//TRecalibrationEMModel_Base
//---------------------------------------------------------------
TRecalibrationEMModel_Base::TRecalibrationEMModel_Base(){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate
	_initialized = false;
	_numSitesAdded = 0;
	_myShift = 0;
	_betas = NULL;
	_oldBetas = NULL;
	_numParameters = 0;
	_name = "base";
};

TRecalibrationEMModel_Base::TRecalibrationEMModel_Base(int Shift){
	_initialized = false;
	_numSitesAdded = 0;
	_myShift = Shift;
	_betas = NULL;
	_oldBetas = NULL;
	_numParameters = 0;
	_name = "base";
};

void TRecalibrationEMModel_Base::_parseParameterString(std::vector<std::string> & vec, std::vector<double>* values){
	std::vector<std::string> tmpVec;

	//quality
	fillVectorFromString(vec[0], tmpVec, ',');
	repeatIndexes(tmpVec, values[0]);

	//position
	fillVectorFromString(vec[1], tmpVec, ',');
	repeatIndexes(tmpVec, values[1]);

	//context
	fillVectorFromString(vec[2], tmpVec, ',');
	repeatIndexes(tmpVec, values[2]);
};

void TRecalibrationEMModel_Base::_allocateBetaMemory(){
	//set initial parameters: all to 0 except first (beta quality) = 1
	_betas = new double[_numParameters];
	_oldBetas = new double[_numParameters];
	_betas[0] = 1.0;
	for(unsigned int i=1; i<_numParameters; ++i)
		_betas[i] = 0.0;

	_initialized = true;
};

void TRecalibrationEMModel_Base::_freeBetaMemory(){
	if(_initialized){
		delete[] _betas;
		delete[] _oldBetas;
		_initialized = false;
	}
	_numSitesAdded = 0;
};

double TRecalibrationEMModel_Base::_calcEpsilon(double & eta){
	if(eta > 23.03) return 0.9999999999;
	if(eta < -23.03) return 0.0000000001;

	eta = exp(eta);
	return eta / (1.0 + eta);
};

void TRecalibrationEMModel_Base::proposeNewParameters(double & lambda, arma::mat & JxF){
	//save old parameters
	for(unsigned int i=0; i<_numParameters; ++i)
		_oldBetas[i] = _betas[i];

	//update new ones
	for(unsigned int i=0; i<_numParameters; ++i)
		_betas[i] = _oldBetas[i] - lambda * JxF(_myShift + i);
};

void TRecalibrationEMModel_Base::rejectProposedParameters(){
	for(unsigned int i=0; i<_numParameters; ++i)
		_betas[i] = _oldBetas[i];
};

void TRecalibrationEMModel_Base::writeParametersToFile(TOutputFilePlain & out){
	out << _name;
	out << getQualityString();
	out << getPositionString();
	out << getContextString();
};

std::string TRecalibrationEMModel_Base::getModelString(){
	return _name + "[" + getQualityString() + ";" + getPositionString() + ";" + getContextString() + "]";
};

//---------------------------------------------------------------
//TRecalibrationEMModel_noRecal
//---------------------------------------------------------------
TRecalibrationEMModel_noRecal::TRecalibrationEMModel_noRecal(int Shift):TRecalibrationEMModel_Base(Shift){
	_numParameters = 0;
	_name = noRecal_name;
};

double TRecalibrationEMModel_noRecal::getErrorRate(TBase & base){
	return base.errorRate;
};

void TRecalibrationEMModel_noRecal::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne){
	//now fill table
	transformedQuality = new int**[MaxQualPlusOne];
	for(int q=0; q<MaxQualPlusOne; ++q){
		transformedQuality[q] = new int*[MaxPosPlusOne];
		for(int p=0; p<MaxPosPlusOne; ++p){
			transformedQuality[q][p] = new int[20];
			for(int c=0; c<20; ++c){
				//no recal!
				transformedQuality[q][p][c] = q;
			}
		}
	}
};

//---------------------------------------------------------------
// TRecalibrationEMModel_qualFuncPosFunc
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(int Shift):TRecalibrationEMModel_Base(Shift){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 1 intercept for all contexts
	// -> in total, 5 variables to estimate
	_numParameters = 5;
	_name = qualfuncPosFunc_name;
	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(std::vector<std::string> & vec, int Shift):TRecalibrationEMModel_qualFuncPosFunc(Shift){
	std::vector<double> values[3];
	_parseParameterString(vec, values);

	//quality: should be two numbers
	if(values[0].size() != 2)
		throw "Wrong number of quality parameters for model " + _name + ": expected 2 but found " + toString(values[0].size()) + "!";

	_betas[0] = values[0][0];
	_betas[1] = values[0][1] / 100.0;  //scale!

	//position
	if(values[1].size() != 2)
		throw "Wrong number of position parameters for model " + _name + ": expected 2 but found " + toString(values[1].size()) + "!";

	_betas[2] = values[1][0] / 100.0; //scale!
	_betas[3] = values[1][1] / 100.0; //scale!

	//context: should be a single value for all contexts
	if(values[2].size() != 1)
		throw "Wrong number of position parameters for model " + _name + ": expected 1 but found " + toString(values[2].size()) + "!";
	_betas[4] = values[2][0];
};

double TRecalibrationEMModel_qualFuncPosFunc::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality, quality squared, position and position squared
	double eta = _qualPosMap.eta[data.quality] * _betas[0];
	eta += _qualPosMap.etaSquared[data.quality] * _betas[1];
	eta += _qualPosMap.position[data.position] * _betas[2];
	eta += _qualPosMap.positionSquared[data.position] * _betas[3];

	//add intercept
	eta += _betas[4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFunc::addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[4];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];
	q[2] = _qualPosMap.position[data.position];
	q[3] = _qualPosMap.positionSquared[data.position];

	//add to F
	//-------------------------------------
	//quality, quality squared, position, position squared: Derivatives are given by the q's
	F(_myShift    ) += weightF * q[0];
	F(_myShift + 1) += weightF * q[1];
	F(_myShift + 2) += weightF * q[2];
	F(_myShift + 3) += weightF * q[3];

	//now intercept: is at position 4 in F!
	F(_myShift + 4) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//all rows except context
	for(int row=0; row<4; ++row){
		for(int col=row; col<4; ++col){
			Jacobian(_myShift + row, _myShift + col) +=  weightJacobian * q[row] * q[col];
		}
	}

	//intercept column
	int tmpIndex = _myShift + 4;
	for(int p=0; p<4; ++p){
		Jacobian(_myShift + p, tmpIndex) += weightJacobian * q[p];
	}
	//intercept x intercept
	Jacobian(tmpIndex, tmpIndex) += weightJacobian;

	++_numSitesAdded;
};

std::string TRecalibrationEMModel_qualFuncPosFunc::getQualityString(){
	return toString(_betas[0]) + "," + toString(_betas[1] * 100.0);
};

std::string TRecalibrationEMModel_qualFuncPosFunc::getPositionString(){
	return toString(_betas[2] * 100.0) + "," + toString(_betas[3] * 100.0);
};

std::string TRecalibrationEMModel_qualFuncPosFunc::getContextString(){
	return toString(_betas[4]);
};

double TRecalibrationEMModel_qualFuncPosFunc::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = _betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += _betas[1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += _betas[2] * (double) base.distFrom5Prime;

	//q[3] is square of position
	eta += _betas[3] * (double) (base.distFrom5Prime * base.distFrom5Prime);

	//add intercept
	eta += _betas[4];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFunc::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne){
	//quality term
	double* qualTermForTransformation = new double[MaxQualPlusOne];
	double tmp;
	tmp = pow(10.0, -(double) 0.0000000001 / 10.0);
	qualTermForTransformation[0] = log(tmp / (1.0 - tmp));

	for(int i=1; i<MaxQualPlusOne; ++i){
		tmp = pow(10.0, -(double) i / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}

	//position term
	double* posTermForTransformation = new double[MaxPosPlusOne];
	for(int i=0; i<MaxPosPlusOne; ++i){
		posTermForTransformation[i] = _betas[2] * i + _betas[3] * i*i;
	}

	//now fill table
	for(int q=0; q<MaxQualPlusOne; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			//error is independent of context!
			double constant = posTermForTransformation[p] + _betas[4] - qualTermForTransformation[q];
			double transQual;
			if(4.0 * _betas[1] * constant > _betas[0] * _betas[0]){
				throw "beta[0]^2 cannot be smaller than 4*beta[1](position + context constants)";
			}
			if(_betas[1] == 0.0){
				transQual = -constant / _betas[0];
			} else {
				tmp = sqrt(_betas[0] * _betas[0] - 4.0 * _betas[1] * constant);
				transQual = (tmp - _betas[0]) / 2.0 / _betas[1];
			}

			transQual = exp(transQual);
			if(transQual == 0) throw "Choose different quality transformation parameters! transQual == 0";

			int newQual = round(-10.0 * log10(transQual / (1.0 + transQual)));

			//now store for each context
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = newQual;
		}
	}

	//clean up
	delete[] qualTermForTransformation;
	delete[] posTermForTransformation;
};

//---------------------------------------------------------------
//TRecalibrationEMModel
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosFuncContext::TRecalibrationEMModel_qualFuncPosFuncContext(int Shift):TRecalibrationEMModel_Base(Shift){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate
	_numParameters = 24;
	_name = qualfuncPosFuncContext_name;
	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosFuncContext::TRecalibrationEMModel_qualFuncPosFuncContext(std::vector<std::string> & vec, int Shift):TRecalibrationEMModel_qualFuncPosFuncContext(Shift){
	std::vector<double> values[3];
	_parseParameterString(vec, values);

	//quality: should be two numbers
	if(values[0].size() != 2)
		throw "Wrong number of quality parameters for model " + _name + ": expected 2 but found " + toString(values[0].size()) + "!";

	_betas[0] = values[0][0];
	_betas[1] = values[0][1] / 100.0;  //scale!

	//position
	if(values[1].size() != 2)
		throw "Wrong number of position parameters for model " + _name + ": expected 2 but found " + toString(values[1].size()) + "!";

	_betas[2] = values[1][0] / 100.0;  //scale!
	_betas[3] = values[1][1] / 100.0;  //scale!

	//context
	if(values[2].size() != 20)
		throw "Wrong number of context parameters for model " + _name + ": expected 20 but found " + toString(values[2].size()) + "!";

	for(int i=0; i<20; i++)
		_betas[4+i] = values[2][i];
};

double TRecalibrationEMModel_qualFuncPosFuncContext::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality, quality squared, position and position squared
	double eta = _qualPosMap.eta[data.quality] * _betas[0];
	eta += _qualPosMap.etaSquared[data.quality] * _betas[1];
	eta += _qualPosMap.position[data.position] * _betas[2];
	eta += _qualPosMap.positionSquared[data.position] * _betas[3];

	//add context
	eta += _betas[data.context + 4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFuncContext::addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[4];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];
	q[2] = _qualPosMap.position[data.position];
	q[3] = _qualPosMap.positionSquared[data.position];

	//add to F
	//-------------------------------------
	//quality, quality squared, position, position squared: Derivatives are given by the q's
	F(_myShift    ) += weightF * q[0];
	F(_myShift + 1) += weightF * q[1];
	F(_myShift + 2) += weightF * q[2];
	F(_myShift + 3) += weightF * q[3];

	//now context: start at position 4 in F!
	F(data.context + 4 + _myShift) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//all rows except context
	for(int row=0; row<4; ++row){
		for(int col=row; col<4; ++col){
			Jacobian(_myShift + row, _myShift + col) +=  weightJacobian * q[row] * q[col];
		}
	}

	//context column
	int tmpIndex = _myShift + data.context + 4;
	for(int p=0; p<4; ++p){
		Jacobian(_myShift + p, tmpIndex) += weightJacobian * q[p];
	}

	//context x context: only add to diagonal, as all others are 0
	Jacobian(tmpIndex, tmpIndex) += weightJacobian;


	++_numSitesAdded;
};

std::string TRecalibrationEMModel_qualFuncPosFuncContext::getQualityString(){
	return toString(_betas[0]) + "," + toString(_betas[1] * 100.0); //scale quadratic effect!
};

std::string TRecalibrationEMModel_qualFuncPosFuncContext::getPositionString(){
	return toString(_betas[2] * 100.0) + "," + toString(_betas[3] * 100.0); //scale!
};

std::string TRecalibrationEMModel_qualFuncPosFuncContext::getContextString(){
	return concatenateString(&_betas[4], 20, ",");
};

double TRecalibrationEMModel_qualFuncPosFuncContext::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = _betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += _betas[1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += _betas[2] * (double) base.distFrom5Prime;

	//q[3] is square of position
	eta += _betas[3] * (double) (base.distFrom5Prime * base.distFrom5Prime);

	//q[4] until q[23] are indicators for the context. Just pick the matching one!
	eta += _betas[base.context + 4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFuncContext::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne){
	//quality term
	double* qualTermForTransformation = new double[MaxQualPlusOne];
	double tmp;
	tmp = pow(10.0, -(double) 0.0000000001 / 10.0);
	qualTermForTransformation[0] = log(tmp / (1.0 - tmp));

	for(int i=1; i<MaxQualPlusOne; ++i){
		tmp = pow(10.0, -(double) i / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}

	//position term
	double* posTermForTransformation = new double[MaxPosPlusOne];
	for(int i=0; i<MaxPosPlusOne; ++i){
		posTermForTransformation[i] = _betas[2] * i + _betas[3] * i*i;
	}

	//now fill table
	for(int q=0; q<MaxQualPlusOne; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c){
				//quality scores
				//now calc transformed quality
				double constant = posTermForTransformation[p] + _betas[c+4] - qualTermForTransformation[q];
				double transQual;

				if(4.0 * _betas[1] * constant > _betas[0] * _betas[0]){
					throw "beta[0]^2 cannot be smaller than 4*beta[1](position + context constants)";
				}
				if(_betas[1] == 0.0){
					transQual = -constant / _betas[0];
				} else {
					tmp = sqrt(_betas[0] * _betas[0] - 4.0 * _betas[1] * constant);
					transQual = (tmp - _betas[0]) / 2.0 / _betas[1];
				}

				transQual = exp(transQual);
				if(transQual == 0) throw "Choose different quality transformation parameters! transQual == 0";
				transformedQuality[q][p][c] = round(-10.0 * log10(transQual / (1.0 + transQual)));
			}
		}
	}

	//clean up
	delete[] qualTermForTransformation;
	delete[] posTermForTransformation;
};

//---------------------------------------------------------------
// TRecalibrationEMModel_qualFuncPosSpecificContext
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosSpecificContext::TRecalibrationEMModel_qualFuncPosSpecificContext(int Shift, int MaxPos):TRecalibrationEMModel_Base(Shift){
	// - transformed quality
	// - square of transformed quality
	// - one parameter per position from 0 to maxPos
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 22 + maxPos + 1 variables to estimate
	_numParamsWithoutPositions = 22;
	_maxPosPlusOne = MaxPos + 1;
	_numParameters = _numParamsWithoutPositions + _maxPosPlusOne;
	_name = qualfuncPosSpecificContext_name;

	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosSpecificContext::TRecalibrationEMModel_qualFuncPosSpecificContext(std::vector<std::string> & vec, int Shift):TRecalibrationEMModel_Base(Shift){
	_numParamsWithoutPositions = 22;
	std::vector<double> values[3];
	_parseParameterString(vec, values);

	//first position so that the memory can be allocated
	if(values[1].size() < 1)
		throw "Missing position values for model " + _name + "!";

	//allocate memory
	_maxPosPlusOne = values[1].size(); //starts at zero!
	_numParameters = _numParamsWithoutPositions + _maxPosPlusOne;

	_allocateBetaMemory();

	//copy position (starts at 22!)
	for(int i=0; i<_maxPosPlusOne; i++)
		_betas[22 + i] = values[1][i];

	//quality: should be two numbers
	if(values[0].size() != 2)
		throw "Wrong number of quality parameters for model " + _name + ": expected 2 but found " + toString(values[0].size()) + "!";

	_betas[0] = values[0][0];
	_betas[1] = values[0][1] / 100.0; //scale!

	//context (starts at 2!)
	if(values[2].size() != 20)
		throw "Wrong number of context parameters for model " + _name + ": expected 20 but found " + toString(values[2].size()) + "!";

	for(int i=0; i<20; i++)
		_betas[2+i] = values[2][i];
};

double TRecalibrationEMModel_qualFuncPosSpecificContext::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality, quality squared
	double eta = _qualPosMap.eta[data.quality] * _betas[0];
	eta += _qualPosMap.etaSquared[data.quality] * _betas[1];

	//add context
	eta += _betas[data.context + 4];

	//add position
	//Note: no check on maxPos! Assuming it was properly initialized for estimation
	eta += _betas[_numParamsWithoutPositions + data.position]; //Position starts at 0

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosSpecificContext::addToFandJacobian(arma::vec & F, arma::mat & Jacobian, const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[4];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];

	//add to F
	//-------------------------------------
	//quality, quality squared: Derivatives are given by the q's
	F(_myShift    ) += weightF * q[0];
	F(_myShift + 1) += weightF * q[1];

	//now context: start at position 2 in F!
	F(data.context + 2 + _myShift) += weightF;

	//now position: start at 22 in F!
	F(_myShift + 22 + data.position) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//quality and quality squared
	Jacobian(_myShift, _myShift) +=  weightJacobian * q[0] * q[0];
	Jacobian(_myShift, _myShift + 1) +=  weightJacobian * q[0] * q[1];
	Jacobian(_myShift + 1, _myShift) +=  weightJacobian * q[1] * q[0];
	Jacobian(_myShift + 1, _myShift + 1) +=  weightJacobian * q[1] * q[1];

	//context x quality
	int tmpIndexContext = _myShift + 2 + data.context;
	Jacobian(_myShift, tmpIndexContext) += weightJacobian * q[0];
	Jacobian(_myShift + 1, tmpIndexContext) += weightJacobian * q[1];

	//context x context: only add to diagonal, as all others are 0
	Jacobian(tmpIndexContext, tmpIndexContext) += weightJacobian;

	//position x quality
	int tmpIndexPos = _myShift + 22 + data.position;
	Jacobian(_myShift, tmpIndexPos) += weightJacobian * q[0];
	Jacobian(_myShift + 1, tmpIndexPos) += weightJacobian * q[1];

	//position x context
	Jacobian(tmpIndexContext, tmpIndexPos) += weightJacobian;

	//position x position
	Jacobian(tmpIndexPos, tmpIndexPos) += weightJacobian;

	++_numSitesAdded;
};

std::string TRecalibrationEMModel_qualFuncPosSpecificContext::getQualityString(){
	return toString(_betas[0]) + "," + toString(_betas[1] * 100.0); //scale quadratic effect
};

std::string TRecalibrationEMModel_qualFuncPosSpecificContext::getPositionString(){
	return concatenateString(&_betas[22], _maxPosPlusOne, ",");
};

std::string TRecalibrationEMModel_qualFuncPosSpecificContext::getContextString(){
	return concatenateString(&_betas[4], 20, ",");
};

double TRecalibrationEMModel_qualFuncPosSpecificContext::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = _betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += _betas[1] * originalErrorRate * originalErrorRate;

	//q[2] until q[21] are indicators for the context. Just pick the matching one!
	eta += _betas[base.context + 2];

	//As of q[22]: position specific effect
	if(base.distFrom5Prime >= _maxPosPlusOne)
		//TODO: give better error. But need read group info for that!
		throw "Position " + toString(base.distFrom5Prime + 1) + " beyond largest position for which recal parameters are available (" + toString(_maxPosPlusOne) + ")!";
	else
		eta += _betas[_numParamsWithoutPositions + base.distFrom5Prime];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosSpecificContext::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne){
	if(MaxPosPlusOne > _maxPosPlusOne)
		throw "Can not fill transformation table for simulations up to position " + toString(MaxPosPlusOne) + ": position specific effects only available up to position " + toString(_maxPosPlusOne) + "!";

	//quality term
	double* qualTermForTransformation = new double[MaxQualPlusOne];
	double tmp;
	tmp = pow(10.0, -(double) 0.0000000001 / 10.0);
	qualTermForTransformation[0] = log(tmp / (1.0 - tmp));

	for(int i=1; i<MaxQualPlusOne; ++i){
		tmp = pow(10.0, -(double) i / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}

	//now fill table
	for(int q=0; q<MaxQualPlusOne; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c){
				//quality scores
				//now calc transformed quality
				double constant = _betas[22 + p] + _betas[c+4] - qualTermForTransformation[q];
				double transQual;

				if(4.0 * _betas[1] * constant > _betas[0] * _betas[0]){
					throw "beta[0]^2 cannot be smaller than 4*beta[1](position + context constants)";
				}
				if(_betas[1] == 0.0){
					transQual = -constant / _betas[0];
				} else {
					tmp = sqrt(_betas[0] * _betas[0] - 4.0 * _betas[1] * constant);
					transQual = (tmp - _betas[0]) / 2.0 / _betas[1];
				}

				transQual = exp(transQual);
				if(transQual == 0) throw "Choose different quality transformation parameters! transQual == 0";
				transformedQuality[q][p][c] = round(-10.0 * log10(transQual / (1.0 + transQual)));
			}
		}
	}

	//clean up
	delete[] qualTermForTransformation;
};

//--------------------------------------------------------------------
// Global functions to create models
//--------------------------------------------------------------------
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string & modelTag, std::vector<std::string> & values, int shift, bool verbose, TLog* logfile){
	trimString(modelTag);

	if(modelTag == qualfuncPosFuncContext_name){
		if(verbose) logfile->list("Will use full model with quality, quality squared, position, position squared and 20 context specific intercepts.");
		return new TRecalibrationEMModel_qualFuncPosFuncContext(values, shift);
	} else if(modelTag == qualfuncPosFunc_name){
		if(verbose) logfile->list("Will use simplified model with only quality, quality squared, position, position squared and one intercept.");
		return new TRecalibrationEMModel_qualFuncPosFunc(values, shift);
	} else if(modelTag == qualfuncPosSpecificContext_name){
		if(verbose) logfile->list("Will use a model with quality, quality squared, 20 context and each position.");
		return new TRecalibrationEMModel_qualFuncPosSpecificContext(values, shift);
	} else if(modelTag == noRecal_name){
		if(verbose) logfile->list("Will use a model that does not recalibrate.");
		return new TRecalibrationEMModel_noRecal(shift);
	} else throw "Unknown recalibration model '" + modelTag + "'!";
};

TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string & modelTag, int maxPos, int shift, bool verbose, TLog* logfile){
	trimString(modelTag);

	if(modelTag == qualfuncPosFuncContext_name){
		return new TRecalibrationEMModel_qualFuncPosFuncContext(shift);
	} else if(modelTag == qualfuncPosFunc_name){
		return new TRecalibrationEMModel_qualFuncPosFunc(shift);
	} else if(modelTag == qualfuncPosSpecificContext_name){
		return new TRecalibrationEMModel_qualFuncPosSpecificContext(shift, maxPos);
	} else if(modelTag == noRecal_name){
		return new TRecalibrationEMModel_noRecal(shift);
	} else throw "Unknown recalibration model '" + modelTag + "'!";
};

//--------------------------------------------------------------------
// TRecalibrationEMModels
//--------------------------------------------------------------------
TRecalibrationEMModels::TRecalibrationEMModels(int numReadGroups, TLog* Logfile){
	totNumParameters = 0;
	readGroupIndex.initialize(numReadGroups);
	logfile = Logfile;
};

TRecalibrationEMModels::~TRecalibrationEMModels(){
	for(TRecalibrationEMModel_Base* it : models)
		delete it;
};

void TRecalibrationEMModels::_addModel(std::string & modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	models.push_back(createTRecalibrationEMModel(modelTag, values, totNumParameters, verbose, logfile));

	totNumParameters += models.back()->numParameters();
};

void TRecalibrationEMModels::addSingleModelForAllReadGroups(std::string modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	//add to read group index
	readGroupIndex.setAllToSingleIndex();

	//add model
	_addModel(modelTag, values, verbose);
};

void TRecalibrationEMModels::addModel(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	//add to read group index
	readGroupIndex.setAsUsed(readGroupId, isSecondMate);

	//add model according to tag
	_addModel(modelTag, values, verbose);
};

void TRecalibrationEMModels::addModel(int readGroupId, bool isSecondMate, std::string modelTag, int maxPos){
	trimString(modelTag);

	//add to read group index
	readGroupIndex.setAsUsed(readGroupId, isSecondMate);

	//add model according to tag
	if(modelTag == qualfuncPosFuncContext_name){
		models.push_back(new TRecalibrationEMModel_qualFuncPosFuncContext(totNumParameters));
	} else if(modelTag == qualfuncPosFunc_name){
		models.push_back(new TRecalibrationEMModel_qualFuncPosFunc(totNumParameters));
	} else if(modelTag == qualfuncPosSpecificContext_name){
		models.push_back(new TRecalibrationEMModel_qualFuncPosSpecificContext(totNumParameters, maxPos));
	} else if(modelTag == noRecal_name){
		models.push_back(new TRecalibrationEMModel_noRecal(totNumParameters));
	} else throw "Unknown recalibration model '" + modelTag + "'!";

	totNumParameters += models.back()->numParameters();
};

bool TRecalibrationEMModels::hasReadGroupsWithoutModel(){
	return readGroupIndex.hasCasesWithoutIndex();
};

void TRecalibrationEMModels::addNoRecalModelIfMissing(){
	//report read groups for which no recal model was given and initialize them as "no_recal" model
	std::pair<int, bool> missingReadGroupInfo;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo))
		addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, noRecal_name, 0);
};

void TRecalibrationEMModels::reportReadGroupsNotUsed(TReadGroups & readGroups){
	readGroupIndex.reportReadGroupsNotUsed(logfile, readGroups);
};

void TRecalibrationEMModels::setEMParamsToZero(){
	Jacobian.resize(totNumParameters, totNumParameters);
	F.resize(totNumParameters);
	JxF.resize(totNumParameters, 1);

	Jacobian.zeros();
	F.zeros();
};

void TRecalibrationEMModels::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	models[ readGroupIndex.index(data) ]->addToFandJacobian(F, Jacobian, data, weightF, weightJacobian);
};

bool TRecalibrationEMModels::solveJxF(const int numSites){
	//Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
	for(unsigned int i=0; i<(totNumParameters-1); ++i){
		for(unsigned int j=i+1; j<totNumParameters; ++j){
			//copy from upper triangle to lower triangle
			Jacobian(j,i) = Jacobian(i,j);
		}
	}

	//scale F and J by 1/#sites
	Jacobian = Jacobian / (double) numSites;
	F = F / (double) numSites;

	//now solve J^-1 x F
	return solve(JxF, Jacobian, F);
};

void TRecalibrationEMModels::proposeNewParameters(double lambda){
	for(TRecalibrationEMModel_Base* it : models)
		it->proposeNewParameters(lambda, JxF);
};

void TRecalibrationEMModels::rejectProposedParameters(){
	for(TRecalibrationEMModel_Base* it : models)
		it->rejectProposedParameters();
};

double TRecalibrationEMModels::getSteepestGradient(){
	double maxF = 0.0;
	for(unsigned int i=0; i<totNumParameters; ++i){
		if(fabs(F(i)) > maxF) maxF = fabs(F(i));
	}
	return maxF;
};

void TRecalibrationEMModels::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << Jacobian << std::endl << std::endl;
};

void TRecalibrationEMModels::printFToStdOut(){
	std::cout << std::endl << std::endl << "F:" << std::endl << F << std::endl << std::endl;
};

void TRecalibrationEMModels::printJxFToStdOut(){
	std::cout << std::endl << std::endl << "JxF:" << std::endl << JxF << std::endl << std::endl;
};

void TRecalibrationEMModels::writeHeader(TOutputFilePlain & out){
	out.writeHeader({"readGroup", "mate", "model", "quality", "position", "context"});
};

void TRecalibrationEMModels::writeParameters(TOutputFilePlain & out, TReadGroups & readGroups){
	for(int r=0; r<readGroupIndex.numReadGroups(); ++r){
		_writeParameters(out, readGroups.getName(r), r, false);
		_writeParameters(out, readGroups.getName(r), r, true);
	}
};

void TRecalibrationEMModels::_writeParameters(TOutputFilePlain & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate){
	if(readGroupIndex.inUse(readGroup, isSecondMate)){
		out << readGroupName;
		if(isSecondMate) out << "second";
		else out << "first";
		models[ readGroupIndex.index(readGroup, isSecondMate) ]->writeParametersToFile(out);
		out << std::endl;
	}
};

