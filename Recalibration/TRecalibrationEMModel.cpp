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
TRecalibrationEMModel_Base::TRecalibrationEMModel_Base(TLog* Logfile){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate
	logfile = Logfile;
	_initialized = false;
	_numSitesAdded = 0;
	_betas = NULL;
	_oldBetas = NULL;
	_numParameters = 0;
	_name = "base";
	_NRStepAccepted = false;
	_NRconverged = false;
	_Q = 0.0;
	_oldQ = 0.0;
};

void TRecalibrationEMModel_Base::_parseParameterString(std::vector<std::string> & vec, std::vector< std::vector<double> > & values){
	std::vector<std::string> tmpVec;

	//quality
	values.emplace_back();
	if(vec.size() > 0){
		fillVectorFromString(vec[0], tmpVec, ',');
		repeatIndexes(tmpVec, *values.rbegin());
	}

	//position
	values.emplace_back();
	if(vec.size() > 1){
		fillVectorFromString(vec[1], tmpVec, ',');
		repeatIndexes(tmpVec, *values.rbegin());
	}

	//context
	values.emplace_back();
	if(vec.size() > 2){
		fillVectorFromString(vec[2], tmpVec, ',');
		repeatIndexes(tmpVec, *values.rbegin());
	}

	//read length
	values.emplace_back();
	if(vec.size() > 3){
		fillVectorFromString(vec[3], tmpVec, ',');
		repeatIndexes(tmpVec, *values.rbegin());
	}
};

void TRecalibrationEMModel_Base::parseQualityParameters(double* betaPointer, std::vector<double> & values){
	//quality: should be two numbers
	if(values.size() != 2)
		throw "Wrong number of quality parameters for model " + _name + ": expected 2 but found " + toString(values.size()) + "!";

	betaPointer[0] = values[0];
	betaPointer[1] = values[1] / 100.0;  //scale!
};

void TRecalibrationEMModel_Base::parsePositionParameters(double* betaPointer, std::vector<double> & values){
	if(values.size() != 2)
		throw "Wrong number of position parameters for model " + _name + ": expected 2 but found " + toString(values.size()) + "!";

	betaPointer[0] = values[0] / 100.0;  //scale!
	betaPointer[1] = values[1] / 10000.0;  //scale!
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
	if(eta > 23.03){
		return 0.9999999999;
	}
	if(eta < -23.03){
		return 0.0000000001;
	}

	eta = exp(eta);
	return eta / (1.0 + eta);
};

void TRecalibrationEMModel_Base::setEMParamsToZero(){
	Jacobian.resize(_numParameters, _numParameters);
	F.resize(_numParameters);
	JxF.resize(_numParameters, 1);

	Jacobian.zeros();
	F.zeros();

	_numSitesAdded = 0;
	_NRconverged = false;
	_NRStepAccepted = false;
};

void TRecalibrationEMModel_Base::setQToZero(){
	if(!_NRconverged){
		_oldQ = _Q;
		_Q = 0.0;
	}
};

double TRecalibrationEMModel_Base::_calcQ(const int & genotype, TRecalibrationEMReadData & data){
	double eps = calcEpsilon(data);
	double B = 1.33333333333333333333 * data.D[genotype] - 1.0;
	double P_d_given_g_beta = B * eps - data.D[genotype] + 1.0;
	return log(P_d_given_g_beta);
};

void TRecalibrationEMModel_Base::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	if(!_NRconverged){
		_Q += _calcQ(knownGenotype, data);
	}
};

void TRecalibrationEMModel_Base::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
	if(!_NRconverged){
		//add this data for all genotypes
		for(int g=0; g<4; ++g){
			_Q += P_g_given_d_oldBeta[g] * _calcQ(g, data);
		}
	}
};

bool TRecalibrationEMModel_Base::solveJxF(){
	if(_NRconverged){
		return true;
	} else {
		//Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
		for(unsigned int i=0; i<(_numParameters-1); ++i){
			for(unsigned int j=i+1; j<_numParameters; ++j){
				//copy from upper triangle to lower triangle
				Jacobian(j,i) = Jacobian(i,j);
			}
		}

		//scale F and J by 1/#sites
		Jacobian = Jacobian / (double) _numSitesAdded;
		F = F / (double) _numSitesAdded;

		//now solve J^-1 x F
		return solve(JxF, Jacobian, F);
	}
};

void TRecalibrationEMModel_Base::proposeNewParameters(double & lambda){
	if(!_NRStepAccepted){
		//save old parameters
		for(unsigned int i=0; i<_numParameters; ++i)
			_oldBetas[i] = _betas[i];

		//update new ones
		for(unsigned int i=0; i<_numParameters; ++i)
			_betas[i] = _oldBetas[i] - lambda * JxF(i);
	}
};

bool TRecalibrationEMModel_Base::acceptProposedParametersBasedOnQ(){
	if(_NRStepAccepted) return true;
	if(_Q > _oldQ){
		_NRStepAccepted = true;
	} else {
		_NRStepAccepted = false;
		rejectProposedParameters();
		_Q = _oldQ;
	}
	return _NRStepAccepted;
};

void TRecalibrationEMModel_Base::rejectProposedParameters(){
	for(unsigned int i=0; i<_numParameters; ++i)
		_betas[i] = _oldBetas[i];
};

double TRecalibrationEMModel_Base::getSteepestGradient(){
	if(_NRStepAccepted) return 0.0;
	double maxF = 0.0;
	for(unsigned int i=0; i<_numParameters; ++i){
		if(fabs(F(i)) > maxF) maxF = fabs(F(i));
	}
	return maxF;
};

void TRecalibrationEMModel_Base::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << Jacobian << std::endl << std::endl;
};

void TRecalibrationEMModel_Base::printFToStdOut(){
	std::cout << std::endl << std::endl << "F:" << std::endl << F << std::endl << std::endl;
};

void TRecalibrationEMModel_Base::printJxFToStdOut(){
	std::cout << std::endl << std::endl << "JxF:" << std::endl << JxF << std::endl << std::endl;
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

std::string TRecalibrationEMModel_Base::getQualityString(){
	return toString(_betas[0]) + "," + toString(_betas[1] * 100.0);
};

std::string TRecalibrationEMModel_Base::getPositionString(){
	return toString(_betas[2] * 100.0) + "," + toString(_betas[3] * 10000.0);
};

std::string TRecalibrationEMModel_Base::getContextString(){
	return concatenateString(&_betas[4], 20, ",");
};

//---------------------------------------------------------------
//TRecalibrationEMModel_noRecal
//---------------------------------------------------------------
TRecalibrationEMModel_noRecal::TRecalibrationEMModel_noRecal(TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	_numParameters = 0;
	_name = noRecal_name;
};

double TRecalibrationEMModel_noRecal::getErrorRate(TBase & base){
	return base.errorRate;
};

void TRecalibrationEMModel_noRecal::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne, int MinQual){
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
TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 1 intercept for all contexts
	// -> in total, 5 variables to estimate
	_numParameters = 5;
	_name = qualFuncPosFunc_name;
	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(std::vector<std::string> & vec, TLog* Logfile):TRecalibrationEMModel_qualFuncPosFunc(Logfile){
	std::vector< std::vector<double> > values;
	_parseParameterString(vec, values);

	//quality and position
	parseQualityParameters(&_betas[0], values[0]);
	parsePositionParameters(&_betas[2], values[1]);

	//context: should be a single value for all contexts
	if(values[2].size() != 1)
		throw "Wrong number of context parameters for model " + _name + ": expected 1 but found " + toString(values[2].size()) + "!";
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

void TRecalibrationEMModel_qualFuncPosFunc::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[4];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];
	q[2] = _qualPosMap.position[data.position];
	q[3] = _qualPosMap.positionSquared[data.position];

	//add to F
	//-------------------------------------
	//quality, quality squared, position, position squared: Derivatives are given by the q's
	F(0) += weightF * q[0];
	F(1) += weightF * q[1];
	F(2) += weightF * q[2];
	F(3) += weightF * q[3];

	//now intercept: is at position 4 in F!
	F(4) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//all rows except context
	for(int row=0; row<4; ++row){
		for(int col=row; col<4; ++col){
			Jacobian(row, col) +=  weightJacobian * q[row] * q[col];
		}
	}

	//intercept column
	int tmpIndex = 4;
	for(int p=0; p<4; ++p){
		Jacobian(p, tmpIndex) += weightJacobian * q[p];
	}
	//intercept x intercept
	Jacobian(tmpIndex, tmpIndex) += weightJacobian;

	++_numSitesAdded;
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

void TRecalibrationEMModel_qualFuncPosFunc::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne, int MinQual){
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

	//fill q < minQual
	for(int q=0; q<MinQual; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = 0;
		}
	}

	//clean up
	delete[] qualTermForTransformation;
	delete[] posTermForTransformation;
};


//---------------------------------------------------------------
//TRecalibrationEMModel_qualFuncPosFuncContext
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosFuncContext::TRecalibrationEMModel_qualFuncPosFuncContext(TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate
	_numParameters = 24;
	_name = qualFuncPosFuncContext_name;
	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosFuncContext::TRecalibrationEMModel_qualFuncPosFuncContext(std::vector<std::string> & vec, TLog* Logfile):TRecalibrationEMModel_qualFuncPosFuncContext(Logfile){
	std::vector< std::vector<double> > values;
	_parseParameterString(vec, values);

	//quality and position
	parseQualityParameters(&_betas[0], values[0]);
	parsePositionParameters(&_betas[2], values[1]);

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

void TRecalibrationEMModel_qualFuncPosFuncContext::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[4];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];
	q[2] = _qualPosMap.position[data.position];
	q[3] = _qualPosMap.positionSquared[data.position];

	//add to F
	//-------------------------------------
	//quality, quality squared, position, position squared: Derivatives are given by the q's
	F(0) += weightF * q[0];
	F(1) += weightF * q[1];
	F(2) += weightF * q[2];
	F(3) += weightF * q[3];

	//now context: start at position 4 in F!
	F(4 + data.context) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//all rows except context
	for(int row=0; row<4; ++row){
		for(int col=row; col<4; ++col){
			Jacobian(row, col) +=  weightJacobian * q[row] * q[col];
		}
	}

	//context column
	int tmpIndex = 4 + data.context;
	for(int p=0; p<4; ++p){
		Jacobian(p, tmpIndex) += weightJacobian * q[p];
	}

	//context x context: only add to diagonal, as all others are 0
	Jacobian(tmpIndex, tmpIndex) += weightJacobian;


	++_numSitesAdded;
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

//	std::cout << "_calcEpsilon(eta) " << _calcEpsilon(eta) << std::endl;
	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFuncContext::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne, int MinQual){
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
	for(int q=MinQual; q<MaxQualPlusOne; ++q){
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

	//fill q < minQual
	for(int q=0; q<MinQual; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = 0;
		}
	}

	//clean up
	delete[] qualTermForTransformation;
	delete[] posTermForTransformation;
};


//---------------------------------------------------------------
//TRecalibrationEMModel_qualFuncPosFuncContextMapQFunc
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext::TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext(TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - transformed read length
	// - square of transformed read length
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 26 variables to estimate
	_numParameters = 26;
	_name = qualFuncPosFuncLengthFuncContext_name;
	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext::TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext(std::vector<std::string> & vec, TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	std::vector< std::vector<double> > values;
	_parseParameterString(vec, values);

	//quality and position
	parseQualityParameters(&_betas[0], values[0]);
	parsePositionParameters(&_betas[2], values[1]);

	//context
	if(values[2].size() != 20)
		throw "Wrong number of context parameters for model " + _name + ": expected 20 but found " + toString(values[2].size()) + "!";

	for(int i=0; i<20; i++)
		_betas[4+i] = values[2][i];
};

double TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality, quality squared, position and position squared
	double eta = _qualPosMap.eta[data.quality] * _betas[0];
	eta += _qualPosMap.etaSquared[data.quality] * _betas[1];
	eta += _qualPosMap.position[data.position] * _betas[2];
	eta += _qualPosMap.positionSquared[data.position] * _betas[3];

	//fragment length
	eta += data.fragmentLength * _betas[4] + _betas[5]; //TODO: think about proper function

	//add context
	eta += _betas[data.context + 4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[6];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];
	q[2] = _qualPosMap.position[data.position];
	q[3] = _qualPosMap.positionSquared[data.position];
	q[4] = data.fragmentLength; //TODO: think about proper model for fragment length
	q[5] = data.fragmentLength;

	//add to F
	//-------------------------------------
	//quality, position and fragment length: Derivatives are given by the q's
	F(0) += weightF * q[0];
	F(1) += weightF * q[1];
	F(2) += weightF * q[2];
	F(3) += weightF * q[3];
	F(4) += weightF * q[4];
	F(5) += weightF * q[5];

	//now context: start at position 4 in F!
	F(6 + data.context) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//all rows except context
	for(int row=0; row<6; ++row){
		for(int col=row; col<6; ++col){
			Jacobian(row, col) +=  weightJacobian * q[row] * q[col];
		}
	}

	//context column
	int tmpIndex = 4 + data.context;
	for(int p=0; p<4; ++p){
		Jacobian(p, tmpIndex) += weightJacobian * q[p];
	}

	//context x context: only add to diagonal, as all others are 0
	Jacobian(tmpIndex, tmpIndex) += weightJacobian;


	++_numSitesAdded;
};

double TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = _betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += _betas[1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += _betas[2] * (double) base.distFrom5Prime;

	//q[3] is square of position
	eta += _betas[3] * (double) (base.distFrom5Prime * base.distFrom5Prime);

	//q[4]: fragment length. TODO: think about model!
	eta += _betas[4] * base.fragmentLength;

	//q[5]: fragment length. TODO: think about model!
	eta += _betas[5] * base.fragmentLength;

	//q[6] until q[25] are indicators for the context. Just pick the matching one!
	eta += _betas[base.context + 6];

//	std::cout << "_calcEpsilon(eta) " << _calcEpsilon(eta) << std::endl;
	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne, int MinQual){
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
	for(int q=MinQual; q<MaxQualPlusOne; ++q){
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

	//fill q < minQual
	for(int q=0; q<MinQual; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = 0;
		}
	}

	//clean up
	delete[] qualTermForTransformation;
	delete[] posTermForTransformation;
};

//---------------------------------------------------------------
// TRecalibrationEMModel_qualFuncPosSpecific
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosSpecific::TRecalibrationEMModel_qualFuncPosSpecific(int MaxPos, TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	// - transformed quality
	// - square of transformed quality
	// - one parameter per position from 0 to maxPos
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 22 + maxPos + 1 variables to estimate
	_numParamsWithoutPositions = 2;
	_maxPosPlusOne = MaxPos + 1;
	_numParameters = _numParamsWithoutPositions + _maxPosPlusOne;
	_name = qualFuncPosSpecific_name;
	lengthWarningPrinted = false;

	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosSpecific::TRecalibrationEMModel_qualFuncPosSpecific(std::vector<std::string> & vec, TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	_numParamsWithoutPositions = 2;
	_name = qualFuncPosSpecific_name;
	std::vector<double> values[2];

	//parse parameter strings
	std::vector<std::string> tmpVec;

	//quality
	fillVectorFromString(vec[0], tmpVec, ',');
	repeatIndexes(tmpVec, values[0]);

	//position
	fillVectorFromString(vec[1], tmpVec, ',');
	repeatIndexes(tmpVec, values[1]);

	//context
	fillVectorFromString(vec[2], tmpVec, ',');


	//first position so that the memory can be allocated
	if(values[1].size() < 1)
		throw "Missing position values for model " + _name + "!";

	//allocate memory
	_maxPosPlusOne = values[1].size(); //starts at zero!
	_numParameters = _numParamsWithoutPositions + _maxPosPlusOne;

	_allocateBetaMemory();

	//quality
	parseQualityParameters(&_betas[0], values[0]);

	//context: should be a single value for all contexts
	if(vec[2] != "-")
		throw "Wrong context parameter for model " + _name + ": expected \"-\" but found " + vec[2] + "!";

	//copy position (starts at 3!)
	for(int i=0; i<_maxPosPlusOne; i++)
		_betas[_numParamsWithoutPositions + i] = values[1][i];

	//other settings
	lengthWarningPrinted = false;
};

void TRecalibrationEMModel_qualFuncPosSpecific::checkParameterRange(std::vector<int> & Qualities, int maxPos){
	if(_maxPosPlusOne != maxPos + 1)
		throw "Largest position in data (" + toString(maxPos + 1) + ") does not match number of position parameters (" + toString(_maxPosPlusOne) + ")!";
};

double TRecalibrationEMModel_qualFuncPosSpecific::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality, quality squared
	double eta = _qualPosMap.eta[data.quality] * _betas[0];
	eta += _qualPosMap.etaSquared[data.quality] * _betas[1];

	//no intercept to add -> is contained in position betas!

	//add position
	//Note: no check on maxPos! Assuming it was properly initialized for estimation
	eta += _betas[_numParamsWithoutPositions + data.position]; //Position starts at 0

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosSpecific::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[2];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];

	//add to F
	//-------------------------------------
	//quality, quality squared: Derivatives are given by the q's
	F(0) += weightF * q[0];
	F(1) += weightF * q[1];

	//no context intercept!

	//now position: start at 2 in F!
	F(_numParamsWithoutPositions + data.position) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//quality and quality squared
	Jacobian(0, 0) +=  weightJacobian * q[0] * q[0];
	Jacobian(0, 1) +=  weightJacobian * q[0] * q[1];
	Jacobian(1, 1) +=  weightJacobian * q[1] * q[1];

	//no context x quality

	//no context x context: only add to diagonal, as all others are 0

	//position x quality
	int tmpIndexPos = _numParamsWithoutPositions + data.position;
	Jacobian(0, tmpIndexPos) += weightJacobian * q[0];
	Jacobian(1, tmpIndexPos) += weightJacobian * q[1];

	//no position x context

	//position x position
	Jacobian(tmpIndexPos, tmpIndexPos) += weightJacobian;

	++_numSitesAdded;
};

std::string TRecalibrationEMModel_qualFuncPosSpecific::getPositionString(){
	return concatenateString(&_betas[_numParamsWithoutPositions], _maxPosPlusOne, ",");
};

std::string TRecalibrationEMModel_qualFuncPosSpecific::getContextString(){
	return "-";
};

double TRecalibrationEMModel_qualFuncPosSpecific::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = _betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += _betas[1] * originalErrorRate * originalErrorRate;

	//no context intercept

	//As of q[2]: position specific effect
	if(base.distFrom5Prime >= _maxPosPlusOne){
		if(!lengthWarningPrinted){
			//TODO: give better error. But need read group info for that!
			logfile->warning("Position " + toString(base.distFrom5Prime + 1) + " is beyond largest position for which recal parameters are available (" + toString(_maxPosPlusOne) + ")! Will use largest position with recal parameters instead. (future warnings suppressed)");
			lengthWarningPrinted = true;
		}
		eta += _betas[_numParamsWithoutPositions + _maxPosPlusOne - 1];
	} else {
		eta += _betas[_numParamsWithoutPositions + base.distFrom5Prime];
	}

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosSpecific::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne, int MinQual){
	//transform correct quality to transformed quality
	if(MaxPosPlusOne > _maxPosPlusOne + 1) //_maxPlosOne starts at zero but MaxPosPlusOne at 1
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
	for(int q=MinQual; q<MaxQualPlusOne; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			//quality scores
			//now calc transformed quality
			double constant = _betas[_numParamsWithoutPositions + p] - qualTermForTransformation[q];
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
			double tmp = round(-10.0 * log10(transQual / (1.0 + transQual)));
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = tmp;
		}
	}

	//fill q < minQual
	for(int q=0; q<MinQual; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = 0;
		}
	}

	//clean up
	delete[] qualTermForTransformation;
};

//---------------------------------------------------------------
// TRecalibrationEMModel_qualFuncPosSpecificContext
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosSpecificContext::TRecalibrationEMModel_qualFuncPosSpecificContext(int MaxPos, TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	// - transformed quality
	// - square of transformed quality
	// - one parameter per position from 0 to maxPos
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 22 + maxPos + 1 variables to estimate
	_numParamsWithoutPositions = 22;
	_maxPosPlusOne = MaxPos + 1;
	_numParameters = _numParamsWithoutPositions + _maxPosPlusOne;
	_name = qualFuncPosSpecificContext_name;
	lengthWarningPrinted = false;

	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosSpecificContext::TRecalibrationEMModel_qualFuncPosSpecificContext(std::vector<std::string> & vec, TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	_numParamsWithoutPositions = 22;
	_name = qualFuncPosSpecificContext_name;
	std::vector< std::vector<double> > values;
	_parseParameterString(vec, values);

	//first position so that the memory can be allocated
	if(values[1].size() < 1)
		throw "Missing position values for model " + _name + "!";

	//allocate memory
	_maxPosPlusOne = values[1].size(); //starts at zero!
	_numParameters = _numParamsWithoutPositions + _maxPosPlusOne;

	_allocateBetaMemory();
	
	//quality
	parseQualityParameters(&_betas[0], values[0]);

	//context (starts at 2!)
	if(values[2].size() != 20)
		throw "Wrong number of context parameters for model " + _name + ": expected 20 but found " + toString(values[2].size()) + "!";

	for(int i=0; i<20; i++)
		_betas[2+i] = values[2][i];

	//copy position (starts at 22!)
	for(int i=0; i<_maxPosPlusOne; i++)
		_betas[22 + i] = values[1][i];

	//other settings
	lengthWarningPrinted = false;
};

void TRecalibrationEMModel_qualFuncPosSpecificContext::checkParameterRange(std::vector<int> & Qualities, int maxPos){
	if(_maxPosPlusOne < maxPos + 1)
		throw "Largest position in data (" + toString(maxPos + 1) + ") does not match number of position parameters (" + toString(_maxPosPlusOne) + ")!";
};

void TRecalibrationEMModel_qualFuncPosSpecificContext::proposeNewParameters(double & lambda){
	//call default of base
	TRecalibrationEMModel_Base::proposeNewParameters(lambda);

	//change betas such that the context are zero on average
	//first regular context and positions 1 through maxPos
	double mean = 0.0;
	for(int c=2; c<18; c++)
		mean += _betas[c];
	mean /= 16.0;
	for(int c=2; c<18; c++)
		_betas[c] -= mean;
	for(int p = _numParamsWithoutPositions + 1; p < _numParamsWithoutPositions + _maxPosPlusOne; p++)
		_betas[p] += mean;

	//now for first position
	mean = 0.0;
	for(int c=18; c<22; c++)
		mean += _betas[c];
	mean /= 4.0;
	for(int c=18; c<22; c++)
		_betas[c] -= mean;
	_betas[_numParamsWithoutPositions] += mean;
};

double TRecalibrationEMModel_qualFuncPosSpecificContext::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality, quality squared
	double eta = _qualPosMap.eta[data.quality] * _betas[0];
	eta += _qualPosMap.etaSquared[data.quality] * _betas[1];

	//add context
	eta += _betas[2 + data.context];

	//add position
	//Note: no check on maxPos! Assuming it was properly initialized for estimation
	eta += _betas[_numParamsWithoutPositions + data.position]; //Position starts at 0

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosSpecificContext::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill q
	double q[2];
	q[0] = _qualPosMap.eta[data.quality];
	q[1] = _qualPosMap.etaSquared[data.quality];

	//add to F
	//-------------------------------------
	//quality, quality squared: Derivatives are given by the q's
	F(0) += weightF * q[0];
	F(1) += weightF * q[1];

	//now context: start at position 2 in F!
	F(2 + data.context) += weightF;

	//now position
	F(_numParamsWithoutPositions + data.position) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//quality and quality squared
	Jacobian(0, 0) +=  weightJacobian * q[0] * q[0];
	Jacobian(0, 1) +=  weightJacobian * q[0] * q[1];
	Jacobian(1, 1) +=  weightJacobian * q[1] * q[1];

	//intercept x quality
	int tmpIndexContext = 2 + data.context;
	Jacobian(0, tmpIndexContext) += weightJacobian * q[0];
	Jacobian(1, tmpIndexContext) += weightJacobian * q[1];

	//intercept x intercept: only add to diagonal, as all others are 0
	Jacobian(tmpIndexContext, tmpIndexContext) += weightJacobian;

	//position x quality
	int tmpIndexPos = _numParamsWithoutPositions + data.position;
	Jacobian(0, tmpIndexPos) += weightJacobian * q[0];
	Jacobian(1, tmpIndexPos) += weightJacobian * q[1];

	//position x intercept
	Jacobian(tmpIndexContext, tmpIndexPos) += weightJacobian;

	//position x position
	Jacobian(tmpIndexPos, tmpIndexPos) += weightJacobian;

	++_numSitesAdded;
};

std::string TRecalibrationEMModel_qualFuncPosSpecificContext::getPositionString(){
	return concatenateString(&_betas[22], _maxPosPlusOne, ",");
};

std::string TRecalibrationEMModel_qualFuncPosSpecificContext::getContextString(){
	return concatenateString(&_betas[2], 20, ",");
};

double TRecalibrationEMModel_qualFuncPosSpecificContext::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = _betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += _betas[1] * originalErrorRate * originalErrorRate;

	//q[2] until q[21] are indicators for the context. Just pick the matching one!
	eta += _betas[2 + base.context];

	//As of q[22]: position specific effect
	if(base.distFrom5Prime >= _maxPosPlusOne){
		if(!lengthWarningPrinted){
			//TODO: give better error. But need read group info for that!
			logfile->warning("Position " + toString(base.distFrom5Prime + 1) + " is beyond largest position for which recal parameters are available (" + toString(_maxPosPlusOne) + ")! Will use largest position with recal parameters instead. (future warnings suppressed)");
			lengthWarningPrinted = true;
		}
		eta += _betas[_numParamsWithoutPositions + _maxPosPlusOne - 1];
	} else {
		eta += _betas[_numParamsWithoutPositions + base.distFrom5Prime];
	}

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosSpecificContext::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne, int MinQual){
	if(MaxPosPlusOne > _maxPosPlusOne + 1) //_maxPlosOne starts at zero but MaxPosPlusOne at 1
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
	std::cout << "MaxQualPlusOne " << MaxQualPlusOne << std::endl;
	for(int q=MinQual; q<MaxQualPlusOne; ++q){
		std::cout << "q " << q << std::endl;
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c){
				//quality scores
				//now calc transformed quality
				double constant = _betas[22 + p] + _betas[2 + c] - qualTermForTransformation[q];
				double transQual;

				if(4.0 * _betas[1] * constant > _betas[0] * _betas[0]){
					std::cout << 4.0 * _betas[1] * constant << ">" << _betas[0] * _betas[0] << std::endl;
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

	//fill q < minQual
	for(int q=0; q<MinQual; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = 0;
		}
	}

	//clean up
	delete[] qualTermForTransformation;
};

//---------------------------------------------------------------
// TRecalibrationEMModel_qualSpecficPosSpecific
//---------------------------------------------------------------
TRecalibrationEMModel_qualSpecficPosSpecific::TRecalibrationEMModel_qualSpecficPosSpecific(std::vector<int> & Qualities, int MaxQual, int MaxPos, TLog* Logfile):TRecalibrationEMModel_Base(Logfile){
	// - one parameter per qualty
	// - one parameter per position from 0 to maxPos
	// -> in total numQual + maxPos + 1 variables to estimate (maxPos is included)
	_maxPosPlusOne = MaxPos + 1;
	_numQualities = Qualities.size();
	_maxQualPlusOne = MaxQual + 1;
	_numParameters = _numQualities + _maxPosPlusOne;
	_name = qualSpecificPosSpecific_name;
	lengthWarningPrinted = false;
	qualityWarningPrinted = false;

	//set which qualities are used
	_qualityIndex = new int[_maxQualPlusOne];
	int index = 0;
	for(int q=0; q<_maxQualPlusOne; ++q){
		if(index < Qualities.size() && Qualities[index] == q){
			_qualityIndex[q] = index;
			++index;
		} else {
			_qualityIndex[q] = -1;
		}
	}

	_allocateBetaMemory();
};

TRecalibrationEMModel_qualSpecficPosSpecific::TRecalibrationEMModel_qualSpecficPosSpecific(std::vector<std::string> & vec, TLog* Logfile):TRecalibrationEMModel_Base(Logfile){

	//DOES NOT WORK!!!

	_name = qualSpecificPosSpecific_name;

	std::vector< std::vector<double> > values;
	_parseParameterString(vec, values);

	//first position so that the memory can be allocated
	if(values[1].size() < 1)
		throw "Missing position values for model " + _name + "!";

	//allocate memory
	_maxPosPlusOne = values[1].size(); //starts at zero!
	_numParameters = _maxPosPlusOne;

	_allocateBetaMemory();

	//quality
	parseQualityParameters(&_betas[0], values[0]);

	//context (starts at 2!)
	if(values[2].size() != 20)
		throw "Wrong number of context parameters for model " + _name + ": expected 20 but found " + toString(values[2].size()) + "!";

	for(int i=0; i<20; i++)
		_betas[2+i] = values[2][i];

	//copy position (starts at 22!)
	for(int i=0; i<_maxPosPlusOne; i++)
		_betas[22 + i] = values[1][i];

	//other settings
	lengthWarningPrinted = false;
	qualityWarningPrinted = false;
};

void TRecalibrationEMModel_qualSpecficPosSpecific::checkParameterRange(std::vector<int> & Qualities, int maxPos){
	if(_maxPosPlusOne < maxPos + 1)
		throw "Largest position in data (" + toString(maxPos + 1) + ") does not match number of position parameters (" + toString(_maxPosPlusOne) + ")!";
	for(const int& q : Qualities){
		if(q >= _maxQualPlusOne || _qualityIndex[q] < 0)
			throw "Quality " + toString(q) + " not within range of quality parameters!";
	}
};

double TRecalibrationEMModel_qualSpecficPosSpecific::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality
	//Note: no check! Assuming it was properly initialized for estimation
	double eta = _betas[_qualityIndex[data.quality]];

	//add position
	//Note: no check on maxPos! Assuming it was properly initialized for estimation
	eta += _betas[_numQualities + data.position]; //Position starts at 0

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualSpecficPosSpecific::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//add to F
	//-------------------------------------
	//quality
	int qIndex = _qualityIndex[data.quality];
	F(qIndex) += weightF;

	//now position
	int pIndex = _numQualities + data.position;
	F(pIndex) += weightF;

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	//quality
	Jacobian(qIndex, qIndex) +=  weightJacobian;

	//position x position
	Jacobian(pIndex, pIndex) += weightJacobian;

	//position x quality
	Jacobian(qIndex, pIndex) += weightJacobian;

	++_numSitesAdded;
};

std::string TRecalibrationEMModel_qualSpecficPosSpecific::getQualityString(){
	return concatenateString(&_betas[0], _maxQualPlusOne, ",");
};

std::string TRecalibrationEMModel_qualSpecficPosSpecific::getPositionString(){
	return concatenateString(&_betas[_numQualities], _maxPosPlusOne, ",");
};

std::string TRecalibrationEMModel_qualSpecficPosSpecific::getContextString(){
	return "-";
};

double TRecalibrationEMModel_qualSpecficPosSpecific::getErrorRate(TBase & base){
	//quality: calculate from error rate
	//TODO: SHOULD BE PHRED INT????
	int q = qualMap.errorToQuality(base.errorRate);

	if(q >= _maxQualPlusOne || _qualityIndex[q] < 0)
		//TODO: give better error. But need read group info for that!
		throw "No recal parameter for quality " + toString(base.distFrom5Prime + 1) + " !";

	double eta = _betas[q];

	//position specific effect
	if(base.distFrom5Prime >= _maxPosPlusOne){
		if(!lengthWarningPrinted){
			//TODO: give better error. But need read group info for that!
			logfile->warning("Position " + toString(base.distFrom5Prime + 1) + " is beyond largest position for which recal parameters are available (" + toString(_maxPosPlusOne) + ")! Will use largest position with recal parameters instead. (future warnings suppressed)");
			lengthWarningPrinted = true;
		}
		eta += _betas[_numQualities + _maxPosPlusOne - 1];
	} else {
		eta += _betas[_numQualities + base.distFrom5Prime];
	}

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualSpecficPosSpecific::fillTransformationTableForSimulation(int*** transformedQuality, int MaxPosPlusOne, int MaxQualPlusOne, int MinQual){
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
	for(int q=MinQual; q<MaxQualPlusOne; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c){
				//quality scores
				//now calc transformed quality
				double constant = _betas[22 + p] + _betas[2 + c] - qualTermForTransformation[q];
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

	//fill q < minQual
	for(int q=0; q<MinQual; ++q){
		for(int p=0; p<MaxPosPlusOne; ++p){
			for(int c=0; c<20; ++c)
				transformedQuality[q][p][c] = 0;
		}
	}

	//clean up
	delete[] qualTermForTransformation;
};


//--------------------------------------------------------------------
// Global functions to create models
//--------------------------------------------------------------------
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, std::vector<std::string> & values, bool verbose, TLog* logfile){
	trimString(modelTag);

	if(modelTag == noRecal_name){
		if(verbose) logfile->list("Will use a model that does not recalibrate.");
		return new TRecalibrationEMModel_noRecal(logfile);
	} else if(modelTag == qualFuncPosFunc_name){
		if(verbose) logfile->list("Will use a model with a quadratic function for quality and position, and one intercept.");
		return new TRecalibrationEMModel_qualFuncPosFunc(values, logfile);
	} else if(modelTag == qualFuncPosFuncContext_name){
		if(verbose) logfile->list("Will use full model with a quadratic function for quality and position and 20 context specific intercepts.");
		return new TRecalibrationEMModel_qualFuncPosFuncContext(values, logfile);
	} else if(modelTag == qualFuncPosFuncLengthFuncContext_name){
		if(verbose) logfile->list("Will use full model with a quadratic function for quality, position and mapping quality, and 20 context specific intercepts.");
		return new TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext(logfile);
	} else if(modelTag == qualFuncPosSpecific_name){
		if(verbose) logfile->list("Will use a model with a quadratic function for quality, and position specific intercepts.");
		return new TRecalibrationEMModel_qualFuncPosSpecific(values, logfile);
	} else if(modelTag == qualFuncPosSpecificContext_name){
		if(verbose) logfile->list("Will use a model with a quadratic function for quality, and an intercept specific for each position and 20 contexts.");
		return new TRecalibrationEMModel_qualFuncPosSpecificContext(values, logfile);
	} else throw "Unknown recalibration model '" + modelTag + "'!";
};

TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, int maxPos, bool verbose, TLog* logfile){
	trimString(modelTag);

	if(modelTag == noRecal_name){
		return new TRecalibrationEMModel_noRecal(logfile);
	} else if(modelTag == qualFuncPosFunc_name){
		return new TRecalibrationEMModel_qualFuncPosFunc(logfile);
	} else if(modelTag == qualFuncPosFuncContext_name){
		return new TRecalibrationEMModel_qualFuncPosFuncContext(logfile);
	} else if(modelTag == qualFuncPosFuncLengthFuncContext_name){
		return new TRecalibrationEMModel_qualFuncPosFuncLengthFuncContext(logfile);
	} else if(modelTag == qualFuncPosSpecific_name){
		return new TRecalibrationEMModel_qualFuncPosSpecific(maxPos, logfile);
	} else if(modelTag == qualFuncPosSpecificContext_name){
		return new TRecalibrationEMModel_qualFuncPosSpecificContext(maxPos, logfile);
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

	models.push_back(createTRecalibrationEMModel(modelTag, values, verbose, logfile));

	totNumParameters += models.back()->numParameters();
};

void TRecalibrationEMModels::addSameModelForAllReadGroups(std::string modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	std::pair<int, bool> missingReadGroupInfo;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo))
		addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, modelTag, values, false);
};

void TRecalibrationEMModels::addModel(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	//add to read group index
	readGroupIndex.setAsUsed(readGroupId, isSecondMate);

	//add model according to tag
	_addModel(modelTag, values, verbose);
};

void TRecalibrationEMModels::addModel(int readGroupId, bool isSecondMate, std::string modelTag, int maxPos){
	//add to read group index
	readGroupIndex.setAsUsed(readGroupId, isSecondMate);

	//create model
	models.push_back(createTRecalibrationEMModel(modelTag, maxPos, false, logfile));
	totNumParameters += models.back()->numParameters();
};

void TRecalibrationEMModels::addModelIfItDoesNotExist(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<int> & Qualities, int maxPos){
	if(readGroupIndex.inUse(readGroupId,isSecondMate)){
		//check model
		models[readGroupIndex.index(readGroupId, isSecondMate)]->checkParameterRange(Qualities, maxPos);
	} else {
		//add to read group index
		readGroupIndex.setAsUsed(readGroupId, isSecondMate);

		//create model
		models.push_back(createTRecalibrationEMModel(modelTag, maxPos, false, logfile));
		totNumParameters += models.back()->numParameters();
	}
};

void TRecalibrationEMModels::removeModel(int readGroupId, bool isSecondMate){
	//get index of model
	int index = readGroupIndex.index(readGroupId, isSecondMate);

	//adjust total number of parameters
	totNumParameters -= models[index]->numParameters();

	//erase
	models.erase(models.begin() + index);

	//update read group index
	readGroupIndex.setAsNotUsed(readGroupId, isSecondMate);
};

void TRecalibrationEMModels::createModels(std::string s, TReadGroups & readGroups, TReadGroupMap & readGroupMap){
	//initialize from string or file
	size_t pos = s.find_first_of('[');
	if(pos == std::string::npos){
		_createModelsFromFile(s, readGroups, readGroupMap);
	} else {
		_createModelsFromString(s, readGroups);
	}
};

void TRecalibrationEMModels::_createModelsFromString(std::string & s, TReadGroups & readGroups){
	//s has format model[param1, param2, param3, ...]
	logfile->startIndent("Initializing recal with string '" + s + "' for all read groups:");

	//read model tag
	size_t pos = s.find_first_of('[');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing '['!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::string modelTag = s.substr(0,pos);
	s.erase(0, pos+1);

	//read parameters: quality, position and context separted by semicolon (;)
	pos = s.find_first_of(']');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing ']'!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::vector<std::string> tmpVec;
	fillVectorFromString(s.substr(0, pos), tmpVec, ";");

	//if(tmpVec.size() != 3)
	//	throw "Failed to understand recal string: wrong number of parameter sets (" + toString(tmpVec.size()) + " instead of 3)!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";

	//initialize model
	addSameModelForAllReadGroups(modelTag, tmpVec, true);

	logfile->endIndent();
};

void TRecalibrationEMModels::_createModelsFromFile(std::string filename, TReadGroups & readGroups, TReadGroupMap & readGroupMap){
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
			if(readGroups.readGroupExists(vec[0])){
				//read read group, mate and model
				int rg = readGroupMap.getIndex(readGroups.find(vec[0]));

				bool isSecondMate;
				if(vec[1] == "second")
					isSecondMate = true;
				else if(vec[1] == "first")
					isSecondMate = false;
				else
					throw "Unknown mate '" + vec[1] + "' in file '" + filename + "' on line " + toString(lineNum) + "!";

				//if read group is pooled. And if so, only create model using the values of the first read group of the pool
				if(!modelExists(rg, isSecondMate)){
					std::string modelTag = vec[2];

					//clean up vec to only contain parameters (remove read group, mate, model and LL)
					vec.erase(vec.begin(), vec.begin() + 3);

					//create model
					addModel(rg, isSecondMate, modelTag, vec, false);
				}
			} else {
				logfile->warning("Read group '" + vec[0] + "' does not exist in the BAM header! Are you using the correct recal file?");
			}
		}
	}
	logfile->done();

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	if(hasReadGroupsWithoutModel()){
		warningForMissingReadGroups(readGroups, readGroupMap);
		logfile->startIndent("Will assume the following read groups to be single end (no recalibration provided for second mate):");
		reportReadGroupsConsideredSingleEnd(readGroups, readGroupMap);
		addNoRecalModelIfMissing();
		logfile->endIndent();
	}
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

void TRecalibrationEMModels::reportReadGroupsNotUsed(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	readGroupIndex.reportReadGroupsNotUsed(logfile, readGroups, ReadGroupMap);
};

void TRecalibrationEMModels::reportReadGroupsNotUsed(TReadGroups & readGroups){
	readGroupIndex.reportReadGroupsNotUsed(logfile, readGroups);
};

void TRecalibrationEMModels::reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	readGroupIndex.reportReadGroupsConsideredSingleEnd(logfile, readGroups, ReadGroupMap);
};

void TRecalibrationEMModels::reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups){
	readGroupIndex.reportReadGroupsConsideredSingleEnd(logfile, readGroups);
};

void TRecalibrationEMModels::warningForMissingReadGroups(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	readGroupIndex.warningForMissingReadGroups(logfile, readGroups, ReadGroupMap);
};

void TRecalibrationEMModels::warningForMissingReadGroups(TReadGroups & readGroups){
	readGroupIndex.warningForMissingReadGroups(logfile, readGroups);
};

void TRecalibrationEMModels::setEMParamsToZero(){
	for(TRecalibrationEMModel_Base* model : models){
		model->setEMParamsToZero();
	}
};

void TRecalibrationEMModels::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	models[ readGroupIndex.index(data) ]->addToFandJacobian(data, weightF, weightJacobian);
};

void TRecalibrationEMModels::setQToZero(){
	for(TRecalibrationEMModel_Base* model : models){
		model->setQToZero();
	}
};

void TRecalibrationEMModels::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
	models[ readGroupIndex.index(data) ]->addToQ(data, P_g_given_d_oldBeta);
};

void TRecalibrationEMModels::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	models[ readGroupIndex.index(data) ]->addToQ(data, knownGenotype);
};

double TRecalibrationEMModels::curQ(){
	double Q = 0.0;
	for(TRecalibrationEMModel_Base* model : models){
		Q += model->curQ();
	}
	return Q;
};

bool TRecalibrationEMModels::solveJxF(){
	bool couldSolve = true;
	for(TRecalibrationEMModel_Base* model : models){
		if(!model->solveJxF()){
			model->printJacobianToStdOut();
			throw "Issue solving JxF! This may be due to a lack of data. Consider adding more sites.";
			couldSolve = false;
		}
	}
	return couldSolve;
};

void TRecalibrationEMModels::proposeNewParameters(double lambda){
	for(TRecalibrationEMModel_Base* it : models){
		it->proposeNewParameters(lambda);
	}
};

unsigned int TRecalibrationEMModels::acceptProposedParametersBasedOnQ(){
	unsigned int numAccepted = 0;
	for(TRecalibrationEMModel_Base* model : models){
		numAccepted += model->acceptProposedParametersBasedOnQ();
	}
	return numAccepted;
};

void TRecalibrationEMModels::rejectProposedParameters(){
	for(TRecalibrationEMModel_Base* it : models)
		it->rejectProposedParameters();
};

double TRecalibrationEMModels::getSteepestGradient(){
	double maxF = 0.0;
	for(TRecalibrationEMModel_Base* model : models){
		double tmp = model->getSteepestGradient();
		if(fabs(tmp) > maxF) maxF = fabs(tmp);
	}
	return maxF;
};

void TRecalibrationEMModels::writeHeader(TOutputFilePlain & out){
	out.writeHeader({"readGroup", "mate", "model", "quality", "position", "context"});
};

void TRecalibrationEMModels::writeParameters(TOutputFilePlain & out, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	for(size_t r=0; r<readGroups.size(); ++r){
		int index = ReadGroupMap[r];
		if(readGroupIndex.inUse(index, false)){
			_writeParameters(out, readGroups.getName(r), index, false);
		}
		if(readGroupIndex.inUse(index, true)){
			_writeParameters(out, readGroups.getName(r), index, true);
		}
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

