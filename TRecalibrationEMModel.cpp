/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMModel.h"

#define qualfuncPosFuncContext_name "qualFuncPosFuncContext"
#define qualfuncPosFunc_name "qualFuncPosFunc"
#define qualfuncPosSpecificContext_name "qualFuncPosSpecificContext"
#define noRecal_name "none"

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

void TRecalibrationEMModel_Base::_allocateBetaMemory(){
	//set initial parameters: all to 0 except first (beta quality) = 1
	_betas = new double[_numParameters];
	_oldBetas = new double[_numParameters];
	_betas[0] = 1.0;
	for(int i=1; i<_numParameters; ++i)
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
	if(eta > 16.11) return 0.9999999;
	if(eta < -16.11) return 0.0000001;

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

void TRecalibrationEMModel_noRecal::writeParametersToFile(TOutputFilePlain & out){
	out << _name << "\t-\t-\t-";
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
	//quality: should be two numbers
	std::vector<double> values;
	fillVectorFromString(vec[0], values, ',');
	if(values.size() != 2)
		throw "Wrong number of quality parameters for model " + _name + ": expected 2 but found " + toString(vec.size()) + "!";

	_betas[0] = values[0];
	_betas[1] = values[1];

	//position
	fillVectorFromString(vec[1], values, ',');
	if(values.size() != 2)
		throw "Wrong number of position parameters for model " + _name + ": expected 2 but found " + toString(vec.size()) + "!";

	_betas[2] = values[0];
	_betas[3] = values[1];

	//context
	fillVectorFromString(vec[2], values, ',');
	if(values.size() != 20)
		throw "Wrong number of context parameters for model " + _name + ": expected 20 but found " + toString(vec.size()) + "!";

	for(int i=0; i<20; i++)
		_betas[4+i] = values[i];
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

void TRecalibrationEMModel_qualFuncPosFuncContext::writeParametersToFile(TOutputFilePlain & out){
	//name
	out << _name;

	//quality
	out << toString(_betas[0]) + "," + toString(_betas[1]);

	//position
	out << toString(_betas[3]) + "," + toString(_betas[3]);

	//context
	out << concatenateString(&_betas[4], 20, ",");
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

//---------------------------------------------------------------
//TRecalibrationEMModelNoContext
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(int Shift):TRecalibrationEMModel_qualFuncPosFuncContext(Shift){
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
	//quality: should be two numbers
	std::vector<double> values;
	fillVectorFromString(vec[0], values, ',');
	if(values.size() != 2)
		throw "Wrong number of quality parameters for model " + _name + ": expected 2 but found " + toString(vec.size()) + "!";

	_betas[0] = values[0];
	_betas[1] = values[1];

	//position
	fillVectorFromString(vec[1], values, ',');
	if(values.size() != 2)
		throw "Wrong number of position parameters for model " + _name + ": expected 2 but found " + toString(vec.size()) + "!";

	_betas[2] = values[0];
	_betas[3] = values[1];

	//context: should be a dash
	if(vec[2] != "-")
		throw "Wrong number of context parameters for model " + _name + ": expected zero but found " + toString(vec.size()) + "!";
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

	//intercept column
	int tmpIndex = _myShift + 4;
	for(int p=0; p<4; ++p){
		Jacobian(_myShift + p, tmpIndex) += weightJacobian * q[p];
	}
	//intercept x intercept
	Jacobian(tmpIndex, tmpIndex) += weightJacobian;

	++_numSitesAdded;
};

void TRecalibrationEMModel_qualFuncPosFunc::writeParametersToFile(TOutputFilePlain & out){
	//name
	out << _name << "\t";

	//quality
	out << toString(_betas[0]) + "," + toString(_betas[1]);

	//position
	out << toString(_betas[2]) + "," + toString(_betas[3]);

	//context
	out << "-";
};

double TRecalibrationEMModel_qualFuncPosFunc::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = _betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += _betas[1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += _betas[2] * (double) base.posInRead;

	//q[3] is square of position
	eta += _betas[3] * (double) (base.posInRead * base.posInRead);

	//add intercept
	eta += _betas[4];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};

//---------------------------------------------------------------
//TRecalibrationEMModelPositionSpecific
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosSpecificContext::TRecalibrationEMModel_qualFuncPosSpecificContext(int Shift, int MaxPos):TRecalibrationEMModel_Base(Shift){
	// - transformed quality
	// - square of transformed quality
	// - one parameter per position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 22 + maxPos variables to estimate
	numParamsWithoutPositions = 22;
	maxPos = MaxPos;
	_numParameters = numParamsWithoutPositions + maxPos;
	_name = qualfuncPosSpecificContext_name;
};

TRecalibrationEMModel_qualFuncPosSpecificContext::TRecalibrationEMModel_qualFuncPosSpecificContext(std::vector<std::string> & vec, int Shift):TRecalibrationEMModel_Base(Shift){
	numParamsWithoutPositions = 22;

	//first position so that the memory can be allocated
	std::vector<double> values;
	fillVectorFromString(vec[1], values, ',');
	if(values.size() < 1)
		throw "Missing position values for model " + _name + "!";

	//allocate memory
	maxPos = vec.size();
	_numParameters = numParamsWithoutPositions + maxPos;
	_allocateBetaMemory();

	//copy position (starts at 22!)
	for(int i=0; i<maxPos; i++)
		_betas[22 + i] = values[i];

	//quality: should be two numbers
	fillVectorFromString(vec[0], values, ',');
	if(values.size() != 2)
		throw "Wrong number of quality parameters for model " + _name + ": expected 2 but found " + toString(vec.size()) + "!";

	_betas[0] = values[0];
	_betas[1] = values[1];

	//context (starts at 2!)
	fillVectorFromString(vec[2], values, ',');
	if(values.size() != 20)
		throw "Wrong number of context parameters for model " + _name + ": expected 20 but found " + toString(vec.size()) + "!";

	for(int i=0; i<20; i++)
		_betas[2+i] = values[i];
};

double TRecalibrationEMModel_qualFuncPosSpecificContext::calcEpsilon(const TRecalibrationEMReadData & data){
	//quality, quality squared
	double eta = _qualPosMap.eta[data.quality] * _betas[0];
	eta += _qualPosMap.etaSquared[data.quality] * _betas[1];

	//add context
	eta += _betas[data.context + 4];

	//add position
	//Note: no check on maxPos! Assuming it was properly initialized for estimation
	eta += _betas[numParamsWithoutPositions + data.position]; //Position starts at 0

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

void TRecalibrationEMModel_qualFuncPosSpecificContext::writeParametersToFile(TOutputFilePlain & out){
	//name
	out << _name;

	//quality
	out << toString(_betas[0]) + "," + toString(_betas[1]);

	//position
	out << concatenateString(&_betas[22], maxPos, ",");

	//context
	out << concatenateString(&_betas[4], 20, ",");
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
	if(base.posInRead > maxPos)
		//TODO: give better error. But need read group info for that!
		throw "Position " + toString(base.posInRead + 1) + " beyond largest position for which recal parameters are available (" + toString(maxPos + 1) + ")!";
	else
		eta += _betas[numParamsWithoutPositions + base.posInRead];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
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

	if(modelTag == qualfuncPosFuncContext_name){
		if(verbose) logfile->list("Will use full model with quality, quality squared, position, position squared and 20 context specific intercepts.");
		models.push_back(new TRecalibrationEMModel_qualFuncPosFuncContext(values, totNumParameters));
	} else if(modelTag == qualfuncPosFunc_name){
		if(verbose) logfile->list("Will use simplified model with only quality, quality squared, position, position squared and one intercept.");
		models.push_back(new TRecalibrationEMModel_qualFuncPosFunc(values, totNumParameters));
	} else if(modelTag == qualfuncPosSpecificContext_name){
		if(verbose) logfile->list("Will use a model with quality, quality squared, 20 context and each position.");
		models.push_back(new TRecalibrationEMModel_qualFuncPosSpecificContext(values, totNumParameters));
	} else if(modelTag == noRecal_name){
		if(verbose) logfile->list("Will use a model that does not recalibrate.");
		models.push_back(new TRecalibrationEMModel_noRecal(totNumParameters));
	} else throw "Unknown recalibration model '" + modelTag + "'!";

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
	bool foundOne;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo))
		addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, "noRecal", 0);
};

void TRecalibrationEMModels::reportReadGroupsNotUsed(std::string* readGroupNames){
	readGroupIndex.reportReadGroupsNotUsed(logfile, readGroupNames);
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

void TRecalibrationEMModels::writeParameters(TOutputFilePlain & out, std::string* readGroupNames){
	for(int r=0; r<readGroupIndex.numReadGroups(); ++r){
		_writeParameters(out, readGroupNames[r], r, false);
		_writeParameters(out, readGroupNames[r], r, true);
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

