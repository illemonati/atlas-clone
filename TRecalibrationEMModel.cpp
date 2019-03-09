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
	numParamsPerReadGroup = NULL;

	//other parameters
	initialized = false;
	EMParamsInitialized = false;
	numSitesAdded = 0;
	numReadGroups = 0;
	totNumParams = 0;
	readGroupShifts = NULL;
	betas = NULL;
	oldBetas = NULL;
};

void TRecalibrationEMModel_Base::_allocateBetaMemory(){
	//set initial parameters: all to 0 except first (beta quality) = 1
	betas = new double*[numReadGroups];
	oldBetas = new double*[numReadGroups];
	readGroupShifts = new int[numReadGroups];
	totNumParams = 0;
	for(int r=0; r<numReadGroups; ++r){
		betas[r] = new double[numParamsPerReadGroup[r]];
		oldBetas[r] = new double[numParamsPerReadGroup[r]];
		for(int i=0; i<numParamsPerReadGroup[r]; ++i)
			betas[r][i] = 0.0;
		betas[r][0] = 1.0;

		totNumParams += numParamsPerReadGroup[r];
	}

	//read group shifts
	readGroupShifts[0] = 0;
	for(int r=1; r<numReadGroups; ++r)
		readGroupShifts[r] = readGroupShifts[r-1] + numParamsPerReadGroup[r-1];

	initialized = true;
};

void TRecalibrationEMModel_Base::_freeBetaMemory(){
	if(initialized){
		delete[] readGroupShifts;
		for(int r=0; r<numReadGroups; ++r){
			delete[] betas[r];
			delete[] oldBetas[r];
		}
		delete[] betas;
		delete[] oldBetas;
		initialized = false;
	}

	totNumParams = 0;
	EMParamsInitialized = false;
	numSitesAdded = 0;
};

void TRecalibrationEMModel_Base::initializeEMParams(){
	//initialize variables for EM
	Jacobian.resize(totNumParams, totNumParams);
	Jacobian.zeros();
	F.resize(totNumParams);
	F.zeros();
	JxF.resize(totNumParams, 1);
	JxF.zeros();
	EMParamsInitialized = true;
}

bool TRecalibrationEMModel_Base::setParams(std::vector<std::string> & vec, int & rg){
	if(!EMParamsInitialized)
		throw "Can not set recal parameters: model not initialized!";

	if(vec.size() != numParamsPerReadGroup[rg])
		return false;

	for(int i=0; i<numParamsPerReadGroup[rg]; ++i)
		betas[rg][i] = stringToDouble(vec[i]);

	return true;
};

void TRecalibrationEMModel_Base::setEMParamsToZero(){
	if(!EMParamsInitialized)
		throw "In TRecalibrationEMModel::setEMParamsToZero(): EM Parameters have never been initialized!";

	Jacobian.zeros();
	F.zeros();
	numSitesAdded = 0;
};

double TRecalibrationEMModel_Base::_calcEpsilon(double & eta){
	if(eta > 16.11) return 0.9999999;
	if(eta < -16.11) return 0.0000001;

	eta = exp(eta);
	return eta / (1.0 + eta);
};

bool TRecalibrationEMModel_Base::solveJxF(){
	//Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
	for(unsigned int i=0; i<(totNumParams-1); ++i){
		for(unsigned int j=i+1; j<totNumParams; ++j){
			//copy from upper triangle to lower triangle
			Jacobian(j,i) = Jacobian(i,j);
		}
	}

	//scale F and J by 1/#sites
	Jacobian = Jacobian / (double) numSitesAdded;
	F = F / (double) numSitesAdded;

	//now solve J^-1 x F
	return solve(JxF, Jacobian, F);
};

void TRecalibrationEMModel_Base::proposeNewParameters(double & lambda){
	//save old parameters
	for(int r=0; r<numReadGroups; ++r){
		for(unsigned int i=0; i<numParamsPerReadGroup[r]; ++i){
			oldBetas[r][i] = betas[r][i];
		}
	}

	//update new ones
	for(int r=0; r<numReadGroups; ++r){
		int tmpIndex = r*numParamsPerReadGroup[r];
		for(unsigned int i=0; i<numParamsPerReadGroup[r]; ++i){
			betas[r][i] = oldBetas[r][i] - lambda * JxF(tmpIndex + i);
		}
	}
}

void TRecalibrationEMModel_Base::rejectProposedParameters(){
	for(int r=0; r<numReadGroups; ++r){
		for(unsigned int i=0; i<numParamsPerReadGroup[r]; ++i){
			betas[r][i] = oldBetas[r][i];
		}
	}
}

double TRecalibrationEMModel_Base::getSteepestGradient(){
	double maxF = 0.0;
	for(unsigned int i=0; i<totNumParams; ++i){
		if(fabs(F(i)) > maxF) maxF = fabs(F(i));
	}
	return maxF;
}

void TRecalibrationEMModel_Base::writeParametersToFile(std::ofstream & out, const uint8_t & readGroup){
	for(unsigned int i=0; i<totNumParams; ++i){
		out << "\t" << betas[readGroup][i];
	}
}

void TRecalibrationEMModel_Base::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << Jacobian << std::endl << std::endl;
}

void TRecalibrationEMModel_Base::printFToStdOut(){
	std::cout << std::endl << std::endl << "F:" << std::endl << F << std::endl << std::endl;
}

void TRecalibrationEMModel_Base::printJxFToStdOut(){
	std::cout << std::endl << std::endl << "JxF:" << std::endl << JxF << std::endl << std::endl;
}



//---------------------------------------------------------------
//TRecalibrationEMModel
//---------------------------------------------------------------
TRecalibrationEMModel::TRecalibrationEMModel(){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate
	numParams = 24;
};

void TRecalibrationEMModel::_initialize(int NumReadGroups){
	numReadGroups = NumReadGroups;

	numParamsPerReadGroup = new uint8_t[numReadGroups];
	for(int r=0; r<numReadGroups; ++r)
		numParamsPerReadGroup[r] = numParams;

	//initialize beta memory
	_allocateBetaMemory();
};

void TRecalibrationEMModel::initializeWithValues(std::vector<std::string> & vec, int NumReadGroups){
	if(vec.size() != numParams)
		throw "Wrong number of recal parameters: expected " + toString(numParams) + " but found " + toString(vec.size()) + "!";

	//initialize objects
	_initialize(NumReadGroups);

	//set parameters
	for(int rg=0; rg<numReadGroups; ++rg)
		setParams(vec, rg);
};

void TRecalibrationEMModel::initializeWithHeader(std::vector<std::string> & vec, int NumReadGroups){
	if(vec.size() != numParams + 2)
		throw "Wrong number of recal parameters: expected " + toString(numParams) + " but found " + toString(vec.size()) + "!";

	//initialize objects
	_initialize(NumReadGroups);
};

void TRecalibrationEMModel::initializeWithDataTable(TRecalibrationEMDataTable & dataTable){
	_initialize(dataTable.numReadGroups);
};

double TRecalibrationEMModel::calcEpsilon(TRecalibrationEMReadData & data){
	//quality, quality squared, position and position squared
	double eta = qualPosMap.eta[data.quality] * betas[data.readGroup][0];
	eta += qualPosMap.etaSquared[data.quality] * betas[data.readGroup][1];
	eta += qualPosMap.position[data.position] * betas[data.readGroup][2];
	eta += qualPosMap.positionSquared[data.position] * betas[data.readGroup][3];

	//add context
	eta += betas[data.readGroup][data.context + 4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel::addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data){
	//loop across reads
	for(int k=0; k<numReads; ++k){
		//fill q
		double q[4];
		q[0] = qualPosMap.eta[data[k].quality];
		q[1] = qualPosMap.etaSquared[data[k].quality];
		q[2] = qualPosMap.position[data[k].position];
		q[3] = qualPosMap.positionSquared[data[k].position];

		//add to F
		//-------------------------------------
		double tmp = P_g_given_d_oldBeta * weights[k];

		//quality, quality squared, position, position squared: Derivatives are given by the q's
		F(readGroupShifts[data[k].readGroup]    ) += tmp * q[0];
		F(readGroupShifts[data[k].readGroup] + 1) += tmp * q[1];
		F(readGroupShifts[data[k].readGroup] + 2) += tmp * q[2];
		F(readGroupShifts[data[k].readGroup] + 3) += tmp * q[3];

		//now context: start at position 4 in F!
		F(data[k].context + 4 + readGroupShifts[data[k].readGroup]) += tmp;

		//add to Jacobian (only upper triangle)
		//-------------------------------------
		tmp = weightsJacobian[k];

		//all rows except context
		for(int row=0; row<4; ++row){
			for(int col=row; col<4; ++col){
				Jacobian(readGroupShifts[data[k].readGroup] + row, readGroupShifts[data[k].readGroup] + col) +=  tmp * q[row] * q[col];
			}
		}

		//context column
		int tmpIndex = readGroupShifts[data[k].readGroup] + data[k].context + 4;
		for(int p=0; p<4; ++p){
			Jacobian(readGroupShifts[data[k].readGroup] + p, tmpIndex) += tmp * q[p];
		}

		//context x context: only add to diagonal, as all others are 0
		Jacobian(tmpIndex, tmpIndex) += tmp;
	}

	++numSitesAdded;
};

void TRecalibrationEMModel::writeHeader(std::ofstream & out){
	out << "quality\tquality^2\tposition\tposition^2";
	TGenotypeMap genoMap;
	for(int i=0; i<genoMap.getNumContext(); ++i)
		out << "\t" << genoMap.getContextString(i);
}

double TRecalibrationEMModel::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = betas[base.readGroup][0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += betas[base.readGroup][1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += betas[base.readGroup][2] * (double) base.distFrom5Prime;

	//q[3] is square of position
	eta += betas[base.readGroup][3] * (double) (base.distFrom5Prime * base.distFrom5Prime);

	//q[4] until q[23] are indicators for the context. Just pick the matching one!
	eta += betas[base.readGroup][base.context + 4];

	return _calcEpsilon(eta);
};

//---------------------------------------------------------------
//TRecalibrationEMModelNoContext
//---------------------------------------------------------------
TRecalibrationEMModelNoContext::TRecalibrationEMModelNoContext(){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 1 intercept for all contexts
	// -> in total, 5 variables to estimate
	numParams = 5;
};

double TRecalibrationEMModelNoContext::calcEpsilon(TRecalibrationEMReadData & data){
	//quality, quality squared, position and position squared
	double eta = qualPosMap.eta[data.quality] * betas[data.readGroup][0];
	eta += qualPosMap.etaSquared[data.quality] * betas[data.readGroup][1];
	eta += qualPosMap.position[data.position] * betas[data.readGroup][2];
	eta += qualPosMap.positionSquared[data.position] * betas[data.readGroup][3];

	//add intercept
	eta += betas[data.readGroup][4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModelNoContext::addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data){
	//loop across reads
	for(int k=0; k<numReads; ++k){
		//fill q
		double q[4];
		q[0] = qualPosMap.eta[data[k].quality];
		q[1] = qualPosMap.etaSquared[data[k].quality];
		q[2] = qualPosMap.position[data[k].position];
		q[3] = qualPosMap.positionSquared[data[k].position];

		//add to F
		//-------------------------------------
		double tmp = P_g_given_d_oldBeta * weights[k];

		//quality, quality squared, position, position squared: Derivatives are given by the q's
		F(readGroupShifts[data[k].readGroup]    ) += tmp * q[0];
		F(readGroupShifts[data[k].readGroup] + 1) += tmp * q[1];
		F(readGroupShifts[data[k].readGroup] + 2) += tmp * q[2];
		F(readGroupShifts[data[k].readGroup] + 3) += tmp * q[3];

		//now context: start at position 4 in F!
		F(data[k].context + 4 + readGroupShifts[data[k].readGroup]) += tmp;

		//add to Jacobian (only upper triangle)
		//-------------------------------------
		tmp = weightsJacobian[k];

		//all rows except context
		for(int row=0; row<4; ++row){
			for(int col=row; col<4; ++col){
				Jacobian(readGroupShifts[data[k].readGroup] + row, readGroupShifts[data[k].readGroup] + col) +=  tmp * q[row] * q[col];
			}
		}

		//intercept column
		int tmpIndex = readGroupShifts[data[k].readGroup] + 4;
		for(int p=0; p<4; ++p){
			Jacobian(readGroupShifts[data[k].readGroup] + p, tmpIndex) += tmp * q[p];
		}
		//intercept x intercept
		Jacobian(tmpIndex, tmpIndex) += tmp;
	}

	++numSitesAdded;
};


void TRecalibrationEMModelNoContext::writeHeader(std::ofstream & out){
	out << "quality\tquality^2\tposition\tposition^2\tintercept";
}

double TRecalibrationEMModelNoContext::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = betas[base.readGroup][0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += betas[base.readGroup][1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += betas[base.readGroup][2] * (double) base.distFrom5Prime;

	//q[3] is square of position
	eta += betas[base.readGroup][3] * (double) (base.distFrom5Prime * base.distFrom5Prime);

	//add intercept
	eta += betas[base.readGroup][4];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
}

//---------------------------------------------------------------
//TRecalibrationEMModelPositionSpecific
//---------------------------------------------------------------
TRecalibrationEMModelPositionSpecific::TRecalibrationEMModelPositionSpecific(){
	// - transformed quality
	// - square of transformed quality
	// - one parameter per position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 22 + maxPos variables to estimate
	numParamsWithoutPositions = 22;
	maxPos = 0;
	maxPosPerReadGroup = NULL;
};

void TRecalibrationEMModelPositionSpecific::_initializeFromVectorSize(std::vector<std::string> & vec, int NumReadGroups){
	//must be at least 23 parameters (1 position)
	if(vec.size() < numParamsWithoutPositions + 1)
		throw "Wrong number of recal parameters: expected at least 23 but found " + toString(vec.size()) + "!";

	numReadGroups = NumReadGroups;

	//figure out many positions are considered and initialize
	int numPos = vec.size() - numParamsWithoutPositions;
	maxPosPerReadGroup = new unsigned int[numReadGroups];
	numParamsPerReadGroup = new uint8_t[numReadGroups];
	for(int r=0; r<numReadGroups; ++r){
		maxPosPerReadGroup[r] = numPos;
		numParamsPerReadGroup[r] = numParamsWithoutPositions + maxPosPerReadGroup[r];
	}

	//initialize beta memory
	_allocateBetaMemory();
};

void TRecalibrationEMModelPositionSpecific::initializeWithValues(std::vector<std::string> & vec, int NumReadGroups){
	_initializeFromVectorSize(vec, NumReadGroups);

	//set parameters
	for(int rg=0; rg<numReadGroups; ++rg)
		setParams(vec, rg);
};

void TRecalibrationEMModelPositionSpecific::initializeWithDataTable(TRecalibrationEMDataTable & dataTable){
	numReadGroups = dataTable.numReadGroups;

	//find max pos for each read group
	maxPosPerReadGroup = new unsigned int[numReadGroups];
	numParamsPerReadGroup = new uint8_t[numReadGroups];
	for(int r=0; r<numReadGroups; ++r){
		maxPosPerReadGroup[r] = dataTable.maxPos[r][0];
		numParamsPerReadGroup[r] = numParamsWithoutPositions + maxPosPerReadGroup[r];
	}

	//initialize beta memory
	_allocateBetaMemory();
};

void TRecalibrationEMModelPositionSpecific::initializeWithHeader(std::vector<std::string> & vec, int NumReadGroups){
	_initializeFromVectorSize(vec, NumReadGroups);
};

double TRecalibrationEMModelPositionSpecific::calcEpsilon(TRecalibrationEMReadData & data){
	//quality, quality squared
	double eta = qualPosMap.eta[data.quality] * betas[data.readGroup][0];
	eta += qualPosMap.etaSquared[data.quality] * betas[data.readGroup][1];

	//add context
	eta += betas[data.readGroup][data.context + 4];

	//add position
	if(data.position < maxPosPerReadGroup[data.readGroup])
		eta += betas[data.readGroup][numParamsWithoutPositions + data.position]; //Position starts at 0
	else
		throw "No recalibration info for position " + toString(data.position + 1) + "!"; //TODO: give better error. But need read group info for that!

	return _calcEpsilon(eta);
};

void TRecalibrationEMModelPositionSpecific::addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data){
	//loop across reads
	for(int k=0; k<numReads; ++k){
		//fill q
		double q[4];
		q[0] = qualPosMap.eta[data[k].quality];
		q[1] = qualPosMap.etaSquared[data[k].quality];

		//add to F
		//-------------------------------------
		double tmp = P_g_given_d_oldBeta * weights[k];

		//quality, quality squared: Derivatives are given by the q's
		F(readGroupShifts[data[k].readGroup]    ) += tmp * q[0];
		F(readGroupShifts[data[k].readGroup] + 1) += tmp * q[1];

		//now context: start at position 2 in F!
		F(data[k].context + 2 + readGroupShifts[data[k].readGroup]) += tmp;

		//now position: start at 22 in F!
		F(readGroupShifts[data[k].readGroup] + 22 + data[k].position) += tmp;

		//add to Jacobian (only upper triangle)
		//-------------------------------------
		tmp = weightsJacobian[k];

		//quality and quality squared
		Jacobian(readGroupShifts[data[k].readGroup], readGroupShifts[data[k].readGroup]) +=  tmp * q[0] * q[0];
		Jacobian(readGroupShifts[data[k].readGroup], readGroupShifts[data[k].readGroup] + 1) +=  tmp * q[0] * q[1];
		Jacobian(readGroupShifts[data[k].readGroup] + 1, readGroupShifts[data[k].readGroup]) +=  tmp * q[1] * q[0];
		Jacobian(readGroupShifts[data[k].readGroup] + 1, readGroupShifts[data[k].readGroup] + 1) +=  tmp * q[1] * q[1];

		//context x quality
		int tmpIndexContext = readGroupShifts[data[k].readGroup] + 2 + data[k].context;
		Jacobian(readGroupShifts[data[k].readGroup], tmpIndexContext) += tmp * q[0];
		Jacobian(readGroupShifts[data[k].readGroup] + 1, tmpIndexContext) += tmp * q[1];

		//context x context: only add to diagonal, as all others are 0
		Jacobian(tmpIndexContext, tmpIndexContext) += tmp;

		//position x quality
		int tmpIndexPos = readGroupShifts[data[k].readGroup] + 22 + data[k].position;
		Jacobian(readGroupShifts[data[k].readGroup], tmpIndexPos) += tmp * q[0];
		Jacobian(readGroupShifts[data[k].readGroup] + 1, tmpIndexPos) += tmp * q[1];

		//position x context
		Jacobian(tmpIndexContext, tmpIndexPos) += tmp;

		//position x position
		Jacobian(tmpIndexPos, tmpIndexPos) += tmp;
	}

	++numSitesAdded;
};

void TRecalibrationEMModelPositionSpecific::writeHeader(std::ofstream & out){
	//q and q^2
	out << "quality\tquality^2";

	//context
	TGenotypeMap genoMap;
		for(int i=0; i<genoMap.getNumContext(); ++i)
			out << "\t" << genoMap.getContextString(i);

		//position
	for(int i=0; i<maxPos; ++i)
		out << "\tposition_" << toString(i+1);
};

double TRecalibrationEMModelPositionSpecific::getErrorRate(int rg, double originalErrorRate, const int & posInRead, const uint8_t & context){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	originalErrorRate = log(originalErrorRate / (1.0 - originalErrorRate));
	double eta = betas[rg][0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += betas[rg][1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += betas[rg][2] * (double) posInRead;

	//q[3] is square of position
	eta += betas[rg][3] * (double) (posInRead * posInRead);

	//add intercept
	eta += betas[rg][4];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};



