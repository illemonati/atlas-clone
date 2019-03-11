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
TRecalibrationEMModel_Base::TRecalibrationEMModel_Base(int Shift){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate

	//other parameters
	initialized = false;
	numSitesAdded = 0;
	myShift = Shift;
	betas = NULL;
	oldBetas = NULL;
	numParameters = 0;
};

void TRecalibrationEMModel_Base::_allocateBetaMemory(){
	//set initial parameters: all to 0 except first (beta quality) = 1
	betas = new double[numParameters];
	oldBetas = new double[numParameters];
	betas[0] = 1.0;
	for(int i=1; i<numParameters; ++i)
		betas[i] = 0.0;

	initialized = true;
};

void TRecalibrationEMModel_Base::_freeBetaMemory(){
	if(initialized){
		delete[] betas;
		delete[] oldBetas;
		initialized = false;
	}
	numSitesAdded = 0;
};

bool TRecalibrationEMModel_Base::setParams(std::vector<std::string> & vec){
	if(vec.size() != numParameters)
		return false;

	for(int i=0; i<numParameters; ++i)
		betas[i] = stringToDouble(vec[i]);

	return true;
};

double TRecalibrationEMModel_Base::_calcEpsilon(double & eta){
	if(eta > 16.11) return 0.9999999;
	if(eta < -16.11) return 0.0000001;

	eta = exp(eta);
	return eta / (1.0 + eta);
};

void TRecalibrationEMModel_Base::proposeNewParameters(double & lambda, arma::mat & JxF){
	//save old parameters
	for(int i=0; i<numParameters; ++i)
		oldBetas[i] = betas[i];

	//update new ones
	for(int i=0; i<numParameters; ++i)
		betas[i] = oldBetas[i] - lambda * JxF(myShift + i);
};

void TRecalibrationEMModel_Base::rejectProposedParameters(){
	for(int i=0; i<numParameters; ++i)
		betas[i] = oldBetas[i];
};

void TRecalibrationEMModel_Base::writeParametersToFile(std::ofstream & out){
	for(unsigned int i=0; i<numParameters; ++i)
		out << "\t" << betas[i];
};

void

//---------------------------------------------------------------
//TRecalibrationEMModel_noRecal
//---------------------------------------------------------------
TRecalibrationEMModel_noRecal::TRecalibrationEMModel_noRecal(int Shift):TRecalibrationEMModel_Base(Shift){
	numParameters = 0;
};

double TRecalibrationEMModel_noRecal::getErrorRate(TBase & base){
	return base.errorRate;
};

//---------------------------------------------------------------
//TRecalibrationEMModel
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(int Shift):TRecalibrationEMModel_Base(Shift){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate
	numParameters = 24;
	_allocateBetaMemory();
};

TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(std::vector<std::string> & vec, int Shift):TRecalibrationEMModel_qualFuncPosFunc(Shift){
	if(vec.size() != numParameters)
		throw "Wrong number of recal parameters: expected " + toString(numParameters) + " but found " + toString(vec.size()) + "!";
	setParams(vec);
};

double TRecalibrationEMModel_qualFuncPosFunc::calcEpsilon(TRecalibrationEMReadData & data){
	//quality, quality squared, position and position squared
	double eta = qualPosMap.eta[data.quality] * betas[0];
	eta += qualPosMap.etaSquared[data.quality] * betas[1];
	eta += qualPosMap.position[data.position] * betas[2];
	eta += qualPosMap.positionSquared[data.position] * betas[3];

	//add context
	eta += betas[data.context + 4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFunc::addToFandJacobian(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMReadData* & data, const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta){
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
		F(myShift    ) += tmp * q[0];
		F(myShift + 1) += tmp * q[1];
		F(myShift + 2) += tmp * q[2];
		F(myShift + 3) += tmp * q[3];

		//now context: start at position 4 in F!
		F(data[k].context + 4 + myShift) += tmp;

		//add to Jacobian (only upper triangle)
		//-------------------------------------
		tmp = weightsJacobian[k];

		//all rows except context
		for(int row=0; row<4; ++row){
			for(int col=row; col<4; ++col){
				Jacobian(myShift + row, myShift + col) +=  tmp * q[row] * q[col];
			}
		}

		//context column
		int tmpIndex = myShift + data[k].context + 4;
		for(int p=0; p<4; ++p){
			Jacobian(myShift + p, tmpIndex) += tmp * q[p];
		}

		//context x context: only add to diagonal, as all others are 0
		Jacobian(tmpIndex, tmpIndex) += tmp;
	}

	++numSitesAdded;
};

void TRecalibrationEMModel_qualFuncPosFunc::writeHeader(std::ofstream & out){
	out << "quality\tquality^2\tposition\tposition^2";
	TGenotypeMap genoMap;
	for(int i=0; i<genoMap.getNumContext(); ++i)
		out << "\t" << genoMap.getContextString(i);
}

double TRecalibrationEMModel_qualFuncPosFunc::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += betas[1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += betas[2] * (double) base.distFrom5Prime;

	//q[3] is square of position
	eta += betas[3] * (double) (base.distFrom5Prime * base.distFrom5Prime);

	//q[4] until q[23] are indicators for the context. Just pick the matching one!
	eta += betas[base.context + 4];

	return _calcEpsilon(eta);
};

//---------------------------------------------------------------
//TRecalibrationEMModelNoContext
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosFuncNoContext::TRecalibrationEMModel_qualFuncPosFuncNoContext(int Shift):TRecalibrationEMModel_qualFuncPosFunc(Shift){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 1 intercept for all contexts
	// -> in total, 5 variables to estimate
	numParameters = 5;
	_allocateBetaMemory();
};


TRecalibrationEMModel_qualFuncPosFuncNoContext::TRecalibrationEMModel_qualFuncPosFuncNoContext(std::vector<std::string> & vec, int Shift):TRecalibrationEMModel_qualFuncPosFuncNoContext(Shift){
	if(vec.size() != numParameters)
		throw "Wrong number of recal parameters: expected " + toString(numParameters) + " but found " + toString(vec.size()) + "!";
	setParams(vec);
};

double TRecalibrationEMModel_qualFuncPosFuncNoContext::calcEpsilon(TRecalibrationEMReadData & data){
	//quality, quality squared, position and position squared
	double eta = qualPosMap.eta[data.quality] * betas[0];
	eta += qualPosMap.etaSquared[data.quality] * betas[1];
	eta += qualPosMap.position[data.position] * betas[2];
	eta += qualPosMap.positionSquared[data.position] * betas[3];

	//add intercept
	eta += betas[4];

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosFuncNoContext::addToFandJacobian(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMReadData* & data, const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta){
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
		F(myShift    ) += tmp * q[0];
		F(myShift + 1) += tmp * q[1];
		F(myShift + 2) += tmp * q[2];
		F(myShift + 3) += tmp * q[3];

		//now context: start at position 4 in F!
		F(data[k].context + 4 + myShift) += tmp;

		//add to Jacobian (only upper triangle)
		//-------------------------------------
		tmp = weightsJacobian[k];

		//all rows except context
		for(int row=0; row<4; ++row){
			for(int col=row; col<4; ++col){
				Jacobian(myShift + row, myShift + col) +=  tmp * q[row] * q[col];
			}
		}

		//intercept column
		int tmpIndex = myShift + 4;
		for(int p=0; p<4; ++p){
			Jacobian(myShift + p, tmpIndex) += tmp * q[p];
		}
		//intercept x intercept
		Jacobian(tmpIndex, tmpIndex) += tmp;
	}

	++numSitesAdded;
};


void TRecalibrationEMModel_qualFuncPosFuncNoContext::writeHeader(std::ofstream & out){
	out << "quality\tquality^2\tposition\tposition^2\tintercept";
};

double TRecalibrationEMModel_qualFuncPosFuncNoContext::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += betas[1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += betas[2] * (double) base.posInRead;

	//q[3] is square of position
	eta += betas[3] * (double) (base.posInRead * base.posInRead);

	//add intercept
	eta += betas[4];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};

//---------------------------------------------------------------
//TRecalibrationEMModelPositionSpecific
//---------------------------------------------------------------
TRecalibrationEMModel_qualFuncPosSpecific::TRecalibrationEMModel_qualFuncPosSpecific(int Shift, int MaxPos):TRecalibrationEMModel_Base(Shift){
	// - transformed quality
	// - square of transformed quality
	// - one parameter per position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 22 + maxPos variables to estimate
	numParamsWithoutPositions = 22;
	maxPos = MaxPos;
	numParameters = numParamsWithoutPositions + maxPos;
};

TRecalibrationEMModel_qualFuncPosSpecific::TRecalibrationEMModel_qualFuncPosSpecific(std::vector<std::string> & vec, int Shift):TRecalibrationEMModel_Base(Shift){
	numParamsWithoutPositions = 22;
	if(vec.size() < numParamsWithoutPositions + 1)
		throw "Wrong number of recal parameters: expected at least " + toString(numParamsWithoutPositions + 1) + " but found " + toString(vec.size()) + "!";

	maxPos = vec.size() - numParamsWithoutPositions;
	numParameters = numParamsWithoutPositions + maxPos;

	setParams(vec);
};

double TRecalibrationEMModel_qualFuncPosSpecific::calcEpsilon(TRecalibrationEMReadData & data){
	//quality, quality squared
	double eta = qualPosMap.eta[data.quality] * betas[0];
	eta += qualPosMap.etaSquared[data.quality] * betas[1];

	//add context
	eta += betas[data.context + 4];

	//add position
	//Note: no check on maxPos! Assuming it was properly initialized for estimation
	eta += betas[numParamsWithoutPositions + data.position]; //Position starts at 0

	return _calcEpsilon(eta);
};

void TRecalibrationEMModel_qualFuncPosSpecific::addToFandJacobian(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMReadData* & data, const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta){
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
		F(myShift    ) += tmp * q[0];
		F(myShift + 1) += tmp * q[1];

		//now context: start at position 2 in F!
		F(data[k].context + 2 + myShift) += tmp;

		//now position: start at 22 in F!
		F(myShift + 22 + data[k].position) += tmp;

		//add to Jacobian (only upper triangle)
		//-------------------------------------
		tmp = weightsJacobian[k];

		//quality and quality squared
		Jacobian(myShift, myShift) +=  tmp * q[0] * q[0];
		Jacobian(myShift, myShift + 1) +=  tmp * q[0] * q[1];
		Jacobian(myShift + 1, myShift) +=  tmp * q[1] * q[0];
		Jacobian(myShift + 1, myShift + 1) +=  tmp * q[1] * q[1];

		//context x quality
		int tmpIndexContext = myShift + 2 + data[k].context;
		Jacobian(myShift, tmpIndexContext) += tmp * q[0];
		Jacobian(myShift + 1, tmpIndexContext) += tmp * q[1];

		//context x context: only add to diagonal, as all others are 0
		Jacobian(tmpIndexContext, tmpIndexContext) += tmp;

		//position x quality
		int tmpIndexPos = myShift + 22 + data[k].position;
		Jacobian(myShift, tmpIndexPos) += tmp * q[0];
		Jacobian(myShift + 1, tmpIndexPos) += tmp * q[1];

		//position x context
		Jacobian(tmpIndexContext, tmpIndexPos) += tmp;

		//position x position
		Jacobian(tmpIndexPos, tmpIndexPos) += tmp;
	}

	++numSitesAdded;
};

void TRecalibrationEMModel_qualFuncPosSpecific::writeHeader(std::ofstream & out){
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

double TRecalibrationEMModel_qualFuncPosSpecific::getErrorRate(TBase & base){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	double originalErrorRate = log(base.errorRate / (1.0 - base.errorRate));
	double eta = betas[0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += betas[1] * originalErrorRate * originalErrorRate;

	//q[2] until q[21] are indicators for the context. Just pick the matching one!
	eta += betas[base.context + 2];

	//As of q[22]: position specific effect
	if(base.posInRead > maxPos)
		//TODO: give better error. But need read group info for that!
		throw "Position " + toString(base.posInRead + 1) + " beyond largest position for which recal parameters are available (" + toString(maxPos + 1) + ")!";
	else
		eta += betas[numParamsWithoutPositions + base.posInRead];

	//now calculate epsilon from eta
	return _calcEpsilon(eta);
};



