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
	mergedInd = false;
	readGroupMap = NULL;
	origNumReadGroups = -1;
	numReadGroups = -1;
	readGroupMapInitialized = false;
}

int TRecalibration::findReadGroupIndex(std::string & name, BamTools::SamReadGroupDictionary & readGroups){
	int i = 0;
	for(std::vector<BamTools::SamReadGroup>::iterator it = readGroups.Begin(); it !=  readGroups.End(); ++it, ++i){
		if(name == it->ID) return i;
	}
	return -1;
}

void TRecalibration::initializeReadGroupMap(BamTools::SamHeader* bamHeader, TParameters & params, TLog* logfile){
	origNumReadGroups = bamHeader->ReadGroups.Size();
	readGroupMap = new int[origNumReadGroups];
	readGroupMapInitialized = true;
	if(params.parameterExists("poolReadGroups")) mergedInd = true;
	else mergedInd = false;


	//construct array from vectors and report
	for(int i=0; i<origNumReadGroups; ++i)	readGroupMap[i] = -1; //map initialized

	if(mergedInd){
		//read read groups and their expected lengths
		std::string filename = params.getParameterString("poolReadGroups");
		if(filename=="") throw "No file specifying read groups to merge provided!";
		logfile->listFlush("Reading read groups to be merged from file '" + filename + "' ...");
		std::vector< std::vector<std::string> > readGroupsToMerge;
		std::vector< std::vector<std::string> >::reverse_iterator rIt;
		std::ifstream file(filename.c_str());
		if(!file) throw "Failed to open file '" + filename + "!";

		//parse file and fill vectors
		int lineNum = 0;
		std::vector<std::string> vec;
		std::string readGroup;
		while(file.good() && !file.eof()){
			++lineNum;
			fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
			if(!vec.empty()){
				if(vec.size() < 2) throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + filename + "'! Read groups cannot be merged with themselves!";
				//add to new header
				//others are those to be merged: find read group in header and store int
				readGroupsToMerge.push_back(std::vector<std::string>());
				rIt = readGroupsToMerge.rbegin();
				for(unsigned int i=0; i<vec.size(); ++i){
					rIt->push_back(vec[i]);
				}
			}
		}
		TReadGroups ReadGroupObject;
		ReadGroupObject.fill(*bamHeader);
		logfile->write(" done!");

		std::vector< std::vector<std::string> >::iterator mergeIt = readGroupsToMerge.begin();
		int oldId;

		for(unsigned int rg = 0; rg < readGroupsToMerge.size(); ++rg, ++mergeIt){
			logfile->startIndent("The following read groups will be combined into one group for recalibration:");
			for(std::vector<std::string>::iterator it = mergeIt->begin(); it != mergeIt->end(); ++it){
				logfile->list(*it);
				oldId = ReadGroupObject.find(*it);
				if(readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
				readGroupMap[oldId] = rg;
			}
			logfile->endIndent();
		}

		numReadGroups = readGroupsToMerge.size();

		//now add read groups that will not be merged
		bool printed = false;
		std::string name;
		for(int i = 0; i < ReadGroupObject.size(); ++i){
			//check if it is mapped, otherwise add
			if(readGroupMap[i] < 0){
				if(!printed){
					logfile->startIndent("The following read groups will be kept as is:");
					printed = true;
				}
				name = ReadGroupObject.getName(i);
				logfile->list(name);
				readGroupMap[i] = numReadGroups;
				++numReadGroups;
			}
		}

		if(printed) logfile->endIndent();
		else logfile->list("All existing read groups will be merged into a new read group.");
	}else{
		numReadGroups = origNumReadGroups;
		for(int i = 0; i < numReadGroups; ++i){
			readGroupMap[i] = i;
		}
	}
}

//---------------------------------------------------------------
//TRecalibrationEMModel
//---------------------------------------------------------------
TRecalibrationEMModel::TRecalibrationEMModel(){
	numParams = 24;
	initialized = false;
	EMParamsInitialized = false;
	numSitesAdded = 0;
	numReadGroups = 0;
	totNumParams = 0;
	readGroupShifts = NULL;
	tmpIndex = 0;
	tmp = 0;
	betas = NULL;
	oldBetas = NULL;
}

TRecalibrationEMModel::TRecalibrationEMModel(int NumReadGroups){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 24 variables to estimate
	numParams = 24;
	initialized = false;
	initialize(NumReadGroups);
}

void TRecalibrationEMModel::initialize(int NumReadGroups){
	EMParamsInitialized = false;
	numSitesAdded = 0;

	numReadGroups = NumReadGroups;
	totNumParams = numParams * numReadGroups;
	readGroupShifts = new int[numReadGroups];
	tmpIndex = 0;
	tmp = 0;

	for(int k=0; k<numReadGroups; ++k)
		readGroupShifts[k] = k * numParams;

	//initialize beta memory
	//set initial parameters: all to 0 except beta_quality = 1
	betas = new double*[numReadGroups];
	oldBetas = new double*[numReadGroups];
	for(int r=0; r<numReadGroups; ++r){
		betas[r] = new double[numParams];
		oldBetas[r] = new double[numParams];
		for(int i=0; i<numParams; ++i)
			betas[r][i] = 0.0;
		betas[r][0] = 1.0;
	}
	initialized = true;
}

bool TRecalibrationEMModel::setParams(std::vector<std::string> & vec, int & rg){
	if(vec.size() < (numParams + 1)) return false;
	std::vector<std::string>::iterator it = vec.begin(); ++it; //skype name of read group in first column
	for(int i=0; i<numParams; ++i, ++it){
		betas[rg][i] = stringToDouble(*it);
	}
	return true;
}

void TRecalibrationEMModel::initializeEMParams(){
	//initialize variables for EM
	Jacobian.resize(totNumParams, totNumParams);
	Jacobian.zeros();
	F.resize(totNumParams);
	F.zeros();
	JxF.resize(totNumParams, 1);
	JxF.zeros();
	EMParamsInitialized = true;
}

void TRecalibrationEMModel::setEMParamsToZero(){
	if(!EMParamsInitialized)
		throw "In TRecalibrationEMModel::setEMParamsToZero(): EM Parameters have never been initialized!";

	Jacobian.zeros();
	F.zeros();
	numSitesAdded = 0;
}

double TRecalibrationEMModel::calcEpsilon(const short & readGroup, float* & q, const short & context){
	//eta = params[readGroup[k]][0];
	tmp = 0.0;
	for(int p=0; p<4; ++p){ //loop over all parameters except context
		tmp += betas[readGroup][p] * q[p];
	}
	//add context
	tmp += betas[readGroup][context + 4];

	if(tmp > 16.11) return 0.9999999;
	if(tmp < -16.11) return 0.0000001;

	tmp = exp(tmp);
	return tmp / (1.0 + tmp);
};

void TRecalibrationEMModel::addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, float** & q, short* & readGroup, short* & context){
	//add to F
	//--------
	int m;
	for(int k=0; k<numReads; ++k){
		tmp = P_g_given_d_oldBeta * weights[k];
		//all 4 covariates except context. Derivatives are given by the q's
		for(m=0; m<4; ++m){ //loop over all parameters except context
			F(m + readGroupShifts[readGroup[k]]) += tmp * q[k][m];
		}
		//now context: start at position 4 in F!
		F(context[k] + 4 + readGroupShifts[readGroup[k]]) += tmp;
	}

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	for(int k=0; k<numReads; ++k){
		tmp = weightsJacobian[k];

		//all rows except context
		for(int row=0; row<4; ++row){
			for(int col=row; col<4; ++col){
				Jacobian(readGroupShifts[readGroup[k]] + row, readGroupShifts[readGroup[k]] + col) +=  tmp * q[k][row] * q[k][col];
			}
		}

		//context column
		tmpIndex = readGroupShifts[readGroup[k]] + context[k] + 4;
		for(int p=0; p<4; ++p){
			Jacobian(readGroupShifts[readGroup[k]] + p, tmpIndex) += tmp * q[k][p];
		}
		//context x context: only add to diagonal, as all others are 0
		Jacobian(tmpIndex, tmpIndex) += tmp;
	}

	++numSitesAdded;
}

bool TRecalibrationEMModel::solveJxF(){
	//Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
	for(int i=0; i<(totNumParams-1); ++i){
		for(int j=i+1; j<totNumParams; ++j){
			//copy from upper triangle to lower triangle
			Jacobian(j,i) = Jacobian(i,j);
		}
	}

	//scale F and J by 1/#sites
	Jacobian = Jacobian / (double) numSitesAdded;
	F = F / (double) numSitesAdded;

	//now solve J^-1 x F
	return solve(JxF, Jacobian, F);
}

void TRecalibrationEMModel::proposeNewParameters(double & lambda){
	//save old parameters
	for(int r=0; r<numReadGroups; ++r){
		for(int i=0; i<numParams; ++i){
			oldBetas[r][i] = betas[r][i];
		}
	}

	//update new ones
	for(int r=0; r<numReadGroups; ++r){
		tmpIndex = r*numParams;
		for(int i=0; i<numParams; ++i){
			betas[r][i] = oldBetas[r][i] - lambda * JxF(tmpIndex + i);
		}
	}
}

void TRecalibrationEMModel::rejectProposedParameters(){
	for(int r=0; r<numReadGroups; ++r){
		for(int i=0; i<numParams; ++i){
			betas[r][i] = oldBetas[r][i];
		}
	}
}

double TRecalibrationEMModel::getSteepestGradient(){
	double maxF = 0.0;
	for(int i=0; i<numParams; ++i){
		if(fabs(F(i)) > maxF) maxF = fabs(F(i));
	}
	return maxF;
}

void TRecalibrationEMModel::writeParametersToFile(std::ofstream & out, const int & readGroup){
	for(int i=0; i<numParams; ++i){
		out << "\t" << betas[readGroup][i];
	}
}

void TRecalibrationEMModel::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << Jacobian.diag() << std::endl << std::endl;
}

double TRecalibrationEMModel::getErrorRate(int rg, double originalErrorRate, const int & posInRead, const int & context){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	originalErrorRate = log(originalErrorRate / (1.0 + originalErrorRate));
	double eta = betas[rg][0] * originalErrorRate;

	//q[1] is square of transformed quality
	eta += betas[rg][1] * originalErrorRate * originalErrorRate;

	//q[2] is position
	eta += betas[rg][2] * (double) posInRead;

	//q[3] is square of position
	eta += betas[rg][3] * (double) (posInRead * posInRead);

	//q[4] until q[23] are indicators for the context. Just pick the matching one!
	eta += betas[rg][context + 4];

	//now calculate epsilon from eta
	if(eta > 22.2) return 0.9999999999;
	if(eta < -23.02685) return 0.0000000001;

	eta = exp(eta);
	return eta / (1.0 + eta);
}

//---------------------------------------------------------------
//TRecalibrationEMModelNoContext
//---------------------------------------------------------------
TRecalibrationEMModelNoContext::TRecalibrationEMModelNoContext(int NumReadGroups){
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 1 intercept for all contexts
	// -> in total, 5 variables to estimate
	numParams = 5;
	initialize(NumReadGroups);
}

 double TRecalibrationEMModelNoContext::calcEpsilon(const short & readGroup, float* & q, const short & context){
	//eta = params[readGroup[k]][0];
	tmp = 0.0;
	for(int p=0; p<4; ++p){ //loop over all parameters except context
		tmp += betas[readGroup][p] * q[p];
	}
	//add intercept
	tmp += betas[readGroup][4];

	if(tmp > 16.11) return 0.9999999;
	if(tmp < -16.11) return 0.0000001;

	tmp = exp(tmp);
	return tmp / (1.0 + tmp);
};

void TRecalibrationEMModelNoContext::addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, float** & q, short* & readGroup, short* & context){
	//add to F
	//--------
	int m;
	for(int k=0; k<numReads; ++k){
		tmp = P_g_given_d_oldBeta * weights[k];
		//all 4 covariates except context. Derivatives are given by the q's
		for(m=0; m<4; ++m){ //loop over all parameters except context
			F(readGroupShifts[readGroup[k]] + m) += tmp * q[k][m];
		}
		//now intercept at position 4 in F!
		F(readGroupShifts[readGroup[k]] + 4) += tmp;
	}

	//add to Jacobian (only upper triangle)
	//-------------------------------------
	for(int k=0; k<numReads; ++k){
		tmp = weightsJacobian[k];

		//all rows except context
		for(int row=0; row<4; ++row){
			for(int col=row; col<4; ++col){
				Jacobian(readGroupShifts[readGroup[k]] + row, readGroupShifts[readGroup[k]] + col) +=  tmp * q[k][row] * q[k][col];
			}
		}

		//intercept column
		tmpIndex = readGroupShifts[readGroup[k]] + 4;
		for(int p=0; p<4; ++p){
			Jacobian(readGroupShifts[readGroup[k]] + p, tmpIndex) += tmp * q[k][p];
		}
		//intercept x intercept
		Jacobian(tmpIndex, tmpIndex) += tmp;
	}

	++numSitesAdded;
}

void TRecalibrationEMModelNoContext::writeParametersToFile(std::ofstream & out, const int & readGroup){
	//write q, q2, p and p2
	for(int i=0; i<4; ++i)
		out << "\t" << betas[readGroup][i];
	//write the same intercept for all contcext
	for(int i=0; i<20; ++i)
		out << "\t" << betas[readGroup][4];
}

double TRecalibrationEMModelNoContext::getErrorRate(int rg, double originalErrorRate, const int & posInRead, const int & context){
	//eta = SUM_i beta[i] * q[i] + beta_c of right context c
	// q[0] is transformed quality
	originalErrorRate = log(originalErrorRate / (1.0 + originalErrorRate));
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
	if(eta > 22.2) return 0.9999999999;
	if(eta < -23.02685) return 0.0000000001;

	eta = exp(eta);
	return eta / (1.0 + eta);
}

//---------------------------------------------------------------
//RecalibrationEMSite
//---------------------------------------------------------------
TRecalibrationEMSite::TRecalibrationEMSite(){
	initialized = false;
	q = NULL;
	context = NULL;
	readGroup = NULL;
	readGroupShifts = NULL;
	D = NULL;
	B = NULL;
	epsilon = NULL;
	P_g_given_d_oldBeta = NULL;
	numReads = 0;
};

TRecalibrationEMSite::TRecalibrationEMSite(TSite & site, int* readGroupMap){
	numReads = site.bases.size();
	q = new float*[numReads];
	D = new float*[4];
	B = new float*[4];
	for(int g=0; g<4; ++g){
		D[g] = new float[numReads];
		B[g] = new float[numReads];
	}

	context = new short[numReads];
	epsilon = new float[numReads];
	readGroup = new short[numReads];
	readGroupShifts = new short[numReads];
	P_g_given_d_oldBeta = new float[4];
	initialized = true;
	int k=0;
	double eps;
	for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it, ++k){
		readGroup[k] = readGroupMap[(*it)->readGroup];
		readGroupShifts[k] = readGroup[k] * 24; //shift by num params per read group!
		q[k] = new float[4];

		//we will work with the following q_ikl:
		// - transformed quality
		// - square of transformed quality
		eps = dePhred((*it)->quality);
		if(eps < 0.0000000001) eps = 0.0000000001;
		else if(eps > 0.9999999999) eps = 0.9999999999;

		q[k][0] = log(eps / (1.0 - eps));
		q[k][1] = q[k][0] * q[k][0];

		// - position
		// - square of position
		q[k][2] = (*it)->posInRead;
		q[k][3] = q[k][2] * q[k][2];

		// - 20 context indicators (either 0.0 or 1.0)
		//only store which is one!
		context[k] = (*it)->context;

		//now also store D: D[ref][read]
		switch((*it)->getBaseAsEnum()){
			case A: D[0][k] = 0.0; //geno = AA
					D[1][k] = 1.0; //geno = CC
					D[2][k] = 1.0 - (*it)->PMD_GA; //geno = GG
					D[3][k] = 1.0; //geno = TT
					break;
			case C: D[0][k] = 1.0; //geno = AA
					D[1][k] = (*it)->PMD_CT; //geno = CC
					D[2][k] = 1.0; //geno = GG
					D[3][k] = 1.0; //geno = TT
					break;
			case G: D[0][k] = 1.0; //geno = AA
					D[1][k] = 1.0; //geno = CC
					D[2][k] = (*it)->PMD_GA; //geno = GG
					D[3][k] = 1.0; //geno = TT
					break;
			case T: D[0][k] = 1.0; //geno = AA
					D[1][k] = 1.0 - (*it)->PMD_CT; //geno = CC
					D[2][k] = 1.0; //geno = GG
					D[3][k] = 0.0; //geno = TT
					break;
			case N:
					D[0][k] = 0.0;
					D[1][k] = 0.0;
					D[2][k] = 0.0;
					D[3][k] = 0.0;
					break;
		}

		//now store B
		for(int g=0; g<4; ++g){
			B[g][k] = 4.0 / 3.0 * D[g][k] - 1.0;
		}
	}
};

TRecalibrationEMSite::~TRecalibrationEMSite(){
	if(initialized){
		for(int i=0; i<4; ++i){
			delete[] D[i];
			delete[] B[i];
		}
		for(int i=0; i<numReads; ++i){
			delete[] q[i];
		}
		delete[] q;
		delete[] D;
		delete[] B;
		delete[] context;
		delete[] readGroup;
		delete[] readGroupShifts;
		delete[] epsilon;
		delete[] P_g_given_d_oldBeta;
	}
}

void TRecalibrationEMSite::calcEpsilon(TRecalibrationEMModel* & model){
	//calc epsilon using parameter estimates provided
	for(int k=0; k<numReads; ++k)
		epsilon[k] = model->calcEpsilon(readGroup[k], q[k], context[k]);
}

double TRecalibrationEMSite::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel* & model, double* freqs){
	calcEpsilon(model);

	//over all genotypes
	double P_g_given_d_theta_denominator = 0.0;
	double tmp;
	for(int g=0; g<4; ++g){
		tmp = 1.0;
		//loop over all reads
		for(int k=0; k<numReads; ++k){
			tmp *= B[g][k] * epsilon[k] - D[g][k] + 1.0;
		}
		P_g_given_d_oldBeta[g] = tmp * freqs[g];
		P_g_given_d_theta_denominator += P_g_given_d_oldBeta[g];
	}

	if(P_g_given_d_theta_denominator < 1.0E-25){
		//do again but in log
		double max = 0.0;
		for(int g=0; g<4; ++g){
			tmp = 0.0;
			//loop over all reads
			for(int k=0; k<numReads; ++k){
				tmp += log(B[g][k] * epsilon[k] - D[g][k] + 1.0);
			}
			P_g_given_d_oldBeta[g] = tmp + log(freqs[g]);
			if(g==0) max = P_g_given_d_oldBeta[g];
			else if(P_g_given_d_oldBeta[g] > max) max = P_g_given_d_oldBeta[g];
		}

		//rescale and delog
		P_g_given_d_theta_denominator = 0.0;
		for(int g=0; g<4; ++g){
			P_g_given_d_oldBeta[g] = exp(P_g_given_d_oldBeta[g] - max);
			P_g_given_d_theta_denominator += P_g_given_d_oldBeta[g];
		}
	}

	//calculate P(g|d, theta)
	for(int g=0; g<4; ++g){
		P_g_given_d_oldBeta[g] = P_g_given_d_oldBeta[g] / P_g_given_d_theta_denominator;
	}

	//return LL = P_g_given_d_theta_denominator
	return log(P_g_given_d_theta_denominator);
}

double TRecalibrationEMSite::calcLL(TRecalibrationEMModel* & model, double* freqs){
	calcEpsilon(model);

	//over all genotypes
	double LL = 0.0;
	double tmp;
	for(int g=0; g<4; ++g){
		tmp = 1.0;
		//loop over all reads
		for(int k=0; k<numReads; ++k){
			tmp *= B[g][k] * epsilon[k] - D[g][k] + 1.0;
		}
		LL += tmp * freqs[g];
	}

	//return LL = P_g_given_d_theta_denominator
	return log(LL);
}

double TRecalibrationEMSite::calcQ(TRecalibrationEMModel* & model){
	calcEpsilon(model);

	//now calculate P(d, g, new params)
	double P_d_given_g_beta;
	double Q = 0.0;

	for(int g=0; g<4; ++g){
		P_d_given_g_beta = 1.0;
		//loop over all reads
		for(int k=0; k<numReads; ++k){
			P_d_given_g_beta *= B[g][k] * epsilon[k] - D[g][k] + 1;
		}

		if(P_d_given_g_beta < 1.0E-50) P_d_given_g_beta = 1.0E-50;
		Q += P_g_given_d_oldBeta[g] * log(P_d_given_g_beta);
	}

	return Q;
}

void TRecalibrationEMSite::addToJacobianAndF(TRecalibrationEMModel* & model){
	//calculate epsilon with current parameters
	calcEpsilon(model);

	//tmp variables
	double* weights = new double[numReads];
	double* eps1MinusEps = new double[numReads];
	double* oneMinus2Eps = new double[numReads];
	double* weightJacobian = new double[numReads];
	for(int k=0; k<numReads; ++k){
		eps1MinusEps[k] = epsilon[k] * (1.0 - epsilon[k]);
		oneMinus2Eps[k] = 1.0 - 2.0 * epsilon[k];
	}

	//fill F and Jacobian
	for(int g=0; g<4; ++g){
		//calc weights
		//------------
		for(int k=0; k<numReads; ++k){
			weights[k] = B[g][k] / (1.0 - D[g][k] + B[g][k] * epsilon[k]) * eps1MinusEps[k];
			weightJacobian[k] = P_g_given_d_oldBeta[g] * weights[k] * (oneMinus2Eps[k] - weights[k]);
		}

		model->addToFandJacobian(numReads, weights, weightJacobian, P_g_given_d_oldBeta[g], q, readGroup, context);
	}

	//delete tmp variables
	delete[] weights;
	delete[] eps1MinusEps;
	delete[] oneMinus2Eps;
	delete[] weightJacobian;
}

//---------------------------------------------------------------
//TRecalibrationEMWindow
//---------------------------------------------------------------
TRecalibrationEMWindow::TRecalibrationEMWindow(TBaseFrequencies* baseFreqs, int* ReadGroupMap){
	freqs = new double[4];
	for(int i=0; i<4; ++i) freqs[i] = (*baseFreqs)[i];
	readGroupMap = ReadGroupMap;
}

void TRecalibrationEMWindow::addSite(TSite & site){
	sites.push_back(new TRecalibrationEMSite(site, readGroupMap));
}

double TRecalibrationEMWindow::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel* & model){
	double LL = 0.0;
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		LL += (*site)->fill_P_g_given_d_beta_AND_calcLL(model, freqs);
	}
	return LL;
}

double TRecalibrationEMWindow::calcLL(TRecalibrationEMModel* & model){
	double LL = 0.0;
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		LL += (*site)->calcLL(model, freqs);
	}
	return LL;
}

double TRecalibrationEMWindow::calcQ(TRecalibrationEMModel* & model){
	double Q = 0.0;
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		Q += (*site)->calcQ(model);
	}
	return Q;
}

void TRecalibrationEMWindow::addToJacobianAndF(TRecalibrationEMModel* & model){
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		(*site)->addToJacobianAndF(model);
	}
}

void TRecalibrationEMWindow::setEuqalBaseFrequencies(){
	for(int i=0; i<4; ++i) freqs[i] = 0.25;
}


//---------------------------------------------------------------
//TRecalibrationEM
//---------------------------------------------------------------
TRecalibrationEM::TRecalibrationEM(BamTools::SamHeader* BamHeader, std::string &name, TParameters & args, TLog* Logfile){
	//read groups and log file
	numSitesAdded = 0;
	equalBaseFrequencies = false;
	bamHeader = BamHeader;
	initializeReadGroupMap(BamHeader, args, Logfile);
	readGroupNames = new std::string[origNumReadGroups];
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<origNumReadGroups; ++r, ++it){
		readGroupNames[r] = it->ID;
	}
	logfile = Logfile;

	//initialize model
	//which model to run?
	std::string modelTag = args.getParameterStringWithDefault("model", "full");
	if(modelTag == "full"){
		logfile->list("Will use full model with quality, quality squared, position, position squared and 20 context specific intercepts.");
		model = new TRecalibrationEMModel(numReadGroups);
	} else if(modelTag == "noContext"){
		logfile->list("Will use full model with quality, quality squared, position, position squared and one intercept.");
		model = new TRecalibrationEMModelNoContext(numReadGroups);
	} else throw "Unknown recalibration model '" + modelTag + "'!";

	//Are the values provided?
	estimatetionRequired = false;
	if(args.parameterExists("recal")){ //ToDo: Super ugly hack.... find better solution.
		//read parameters from file
		std::string filename = name;
		logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");
		std::ifstream file(filename.c_str());
		if(!file) throw "Failed to open file '" + filename + "' for reading!";

		//tmp variables for reading
		std::string tmp;
		int lineNum = 0;
		std::vector<std::string> vec;
		std::vector<std::string>::iterator it;
		int rg;
		bool* rgFound = new bool[numReadGroups];
		for(int r=0; r<numReadGroups; ++r) rgFound[r] = false;

		//skip header
		std::getline(file, tmp);

		//parse file to read details for each read group
		while(file.good() && !file.eof()){
			++lineNum;
			fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
			//skip empty lines
			if(vec.size() > 0){
				//find read group
				it = vec.begin();
				if(!bamHeader->ReadGroups.Contains(*it)) throw "Read group '" + *it + "' does not exist in the BAM header!";
				rg = readGroupMap[findReadGroupIndex(*it, bamHeader->ReadGroups)];
				rgFound[rg] = true;

				//add to model
				if(!model->setParams(vec, rg))
					throw "Issues reading reclibration for readGroup '" + *it + "' on line " + toString(lineNum) + "! Do you use the right model?";
			}
		}

		//check if we miss some read groups
		for(int r=0; r<numReadGroups; ++r){
			if(!rgFound[r]) throw "Read group '" + readGroupNames[r] + "' is missing in file '" + filename + "'!";
		}
		delete[] rgFound;
		logfile->write(" done!");

		//check if we anyway estimate things
		if(args.parameterExists("estimateRecal")) estimatetionRequired = true;
	} else estimatetionRequired = true;

	//read estimation parameters, if required
	if(estimatetionRequired){
		logfile->startIndent("Will run EM to estimate recalibration parameters:");
		numEMIterations = args.getParameterIntWithDefault("iterations", 100);
		logfile->list("Will perform at max " + toString(numEMIterations) + " EM iterations.");
		maxEpsilon = args.getParameterDoubleWithDefault("maxEps", 0.000001);
		logfile->list("Will stop EM when deltaLL < " + toString(maxEpsilon));
		NewtonRaphsonNumIterations = args.getParameterIntWithDefault("NRiterations", 10);
		logfile->list("Will conduct at max " + toString(NewtonRaphsonNumIterations) + " Newton-Raphson iterations");
		NewtonRaphsonMaxF = args.getParameterDoubleWithDefault("maxF", 0.0001);
		logfile->list("Will stop Newton-Raphson when F < " + toString(NewtonRaphsonMaxF));
		equalBaseFrequencies = args.parameterExists("equalBaseFreq");
		if(equalBaseFrequencies) logfile->list("Will assume equal base frequencies {0.25, 0.25, 0.25, 0.25}");
		logfile->endIndent();

		//initialize variables for EM
		model->initializeEMParams();
	} else {
		numEMIterations = -1;
		maxEpsilon = 0.0;
		NewtonRaphsonNumIterations = -1;
		NewtonRaphsonMaxF = 0.0;
		maxCoverage = -1;
	}
}

void TRecalibrationEM::addNewWindow(TBaseFrequencies* freqs){
	windows.push_back(new TRecalibrationEMWindow(freqs, readGroupMap));
	//set iterator
	curWindow = windows.end(); --curWindow;
	if(equalBaseFrequencies) (*curWindow)->setEuqalBaseFrequencies();
}

void TRecalibrationEM::addSite(TSite & site){
	(*curWindow)->addSite(site);
	++numSitesAdded;
}

double TRecalibrationEM::getErrorRate(TBase* base){
	return model->getErrorRate(readGroupMap[base->readGroup], dePhred(base->quality), base->posInRead, base->context);
}

void TRecalibrationEM::runNewtonRaphson(int & maxNewtonRaphsonIteratios, double & maxFThreshold, TLog* logfile, bool & writeTmpTables, std::string debugFilename){
	//variables
	double maxF;
	double lambda; //used in backtracking
	bool acceptMove;
	bool NRconverged = false;

	//calculate Q at current location
	double Q;
	double curQ = 0.0;
	for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
		curQ += (*curWindow)->calcQ(model);
	}

	std::ofstream* myStream = NULL;

	//open debug file
	if(writeTmpTables){
		myStream = new std::ofstream(debugFilename.c_str());
		if(!myStream) throw "Failed to open output file '" + debugFilename + "'!";
		//add header
		*myStream << "iteration";
		for(int i=0; i<model->numParams; ++i) *myStream << "\tbeta'" << i;
		for(int i=0; i<model->numParams; ++i) *myStream << "\tF" << i;
		for(int i=0; i<model->numParams; ++i) *myStream << "\tbeta" << i;
		*myStream << std::endl;
	}

	//run up to maxNewtonRaphsonIteratios iterations, but stop if max(F) < maxFThreshold
	logfile->startIndent("Running Newton-Raphson optimization:");
	for(int i=0; i<maxNewtonRaphsonIteratios; ++i){
		logfile->startIndent("Running iteration " + toString(i+1) + ":");
		logfile->listFlush("Calculating Jacobian and gradient ...");
		if(writeTmpTables) *myStream << i;

		//set to zero
		model->setEMParamsToZero();

		//fill Jacobin and F: loop over all sites
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			(*curWindow)->addToJacobianAndF(model);
		}

		//now solve J^-1 x F
		if(model->solveJxF()){
			logfile->write(" done!");

/*
			std::cout << "----------------------------------------------" << std::endl;
			std::cout << "JxF " << JxF << std::endl;
			std::cout << "----------------------------------------------" << std::endl;
*/

			//update params for each read group using backtracking
			lambda = 1.0;
			acceptMove = false;
			while(!acceptMove){
				logfile->listFlush("Proposing move with lambda = " + toString(lambda) + " ...");
				model->proposeNewParameters(lambda);

				//calculate Q at new location
				Q = 0.0;
				for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
					Q += (*curWindow)->calcQ(model);
				}

				//check if we accept or backtrack
				if(Q > curQ){
					acceptMove = true; //accept
					logfile->write(" accepting move!");
					logfile->conclude("Q was reduced from " + toString(curQ) + " to " + toString(Q));
					curQ = Q;
				} else {
					lambda = lambda / 2.0; //backtrack;
					logfile->write(" rejecting move!");
					model->rejectProposedParameters();
					if(lambda < 0.000000001){
						acceptMove = true; //accept
						NRconverged = true;
						logfile->conclude("No improvement even with lambda = " + toString(lambda) + ", aborting Newton-Raphson.");
					}
				}
			}
		} else {
			model->printJacobianToStdOut();
			throw "Issue solving JxF in TRecalibrationEM::runNewtonRalphson()! This may be due to a lack of data. Consider adding more sites.";
		}

		//get largest gradient (F) to check if we break
		maxF = model->getSteepestGradient();
		logfile->conclude("max(F) = " + toString(maxF));
		logfile->endIndent();
		if(maxF < maxFThreshold || NRconverged) break;
	}
	logfile->endIndent();
}

void TRecalibrationEM::runEM(std::string outputName, bool & writeTmpTables){
	logfile->startNumbering("Running EM algorithm to find MLE recalibration parameters:");
	if(numSitesAdded < 100) throw "Less than 100 sites available for recalibration - aborting estimation!";

	double LL, deltaLL, oldLL = 0.0;
	std::ofstream out;
	std::string filename;

	//running iterations
	for(int iter = 0; iter < numEMIterations; ++iter){
		logfile->number("EM Iteration:"); logfile->addIndent();

		//calculate P(g|d, oldbeta) for all sites and calculate LL
		LL = 0.0;
		logfile->listFlush("Calculating P(g|d, beta') ...");
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			LL += (*curWindow)->fill_P_g_given_d_beta_AND_calcLL(model);
		}
		logfile->write(" done!");
		logfile->conclude("Current Log Likelihood = " + toString(LL));

		//DEBUG--------------------------------------------------------
		//calc Q surface for current old params
		//calcQSurface(outputName + "_Qsurface_EMiteration_" + toString(iter) + ".txt", 21);
		//DEBUG--------------------------------------------------------

		//check if we break based on LL
		if(iter > 0){
			deltaLL = LL - oldLL;
			logfile->conclude("Epsilon = " + toString(deltaLL));
			if(iter > 0 && deltaLL < maxEpsilon){
				logfile->conclude("EM has converged (epsilon < " + toString(maxEpsilon) + ")");
				break;
			} else oldLL = LL;
		} else oldLL = LL;


		//run NewtonRaphson until convergence
		runNewtonRaphson(NewtonRaphsonNumIterations, NewtonRaphsonMaxF, logfile, writeTmpTables, outputName + "_NewtonRaphson_" + toString(iter) + ".txt");

		//write current estimates to file
		if(writeTmpTables){
			filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
			logfile->listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename, LL);
			logfile->write(" done!");
		}

		//end loop
		logfile->endIndent();
	}

	//finalize
	logfile->endNumbering();

	//writing final estimates
	filename = outputName + "_recalibrationEM.txt";
	logfile->listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename, LL);
	logfile->write(" done!");

	//calc LL surface
	//calcLikelihoodSurface(outputName + "_LLsurface.txt", 21);
}

void TRecalibrationEM::writeCurrentEstimates(std::string filename, double & LL){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	writeHeader(out);
	writeParams(out, LL);
	out.close();
}

void TRecalibrationEM::writeHeader(std::ofstream & out){
	out << "readGroup\tquality\tquality^2\tposition\tposition^2";
	TGenotypeMap genoMap;
	for(int i=0; i<genoMap.getNumContext(); ++i)
		out << "\t" << genoMap.getContextString(i);
	out << "\tLL" << std::endl;
}

void TRecalibrationEM::writeParams(std::ofstream & out, double & LL){
	for(int r=0; r<origNumReadGroups; ++r){
		out << readGroupNames[r];
		model->writeParametersToFile(out, readGroupMap[r]);
		out << "\t" << LL;
		out << std::endl;
	}
}

double TRecalibrationEM::calcLL(){
	double LL = 0.0;
	for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
		LL += (*curWindow)->calcLL(model);
	}
	return LL;
}

/*
void TRecalibrationEM::calcLikelihoodSurface(std::string filename, int numMarginalGridPoints){
	double LL;

	//open outputfile
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	out << "beta0\tbeta1\tLL" << std::endl;

	//set min, max and step for each parameter
	double min[5];
	min[0] = -5.0;
	min[1] = -5.0;
	min[2] = -1.0;
	min[3] = -1.0;
	min[4] = -1.0;


	double max[5];
	max[0] = 10.0;
	max[1] = 10.0;
	max[2] = 1.0;
	max[3] = 1.0;
	max[4] = 1.0;

	double step[5];
	for(int i=0; i<5; ++i){
		step[i] = (max[i] - min[i]) / (numMarginalGridPoints - 1.0);
	}

	//without last two
	for(int r=0; r<numReadGroups; ++r){
		params[r][3] = 0.0;
		params[r][4] = 0.0;
	}

	//Loop over parameters
	for(int p1=0; p1<numMarginalGridPoints; ++p1){
		//for(int r=0; r<numReadGroups; ++r) params[r][0] = min[0] + p1 * step[0];
		params[0][0] = min[0] + p1 * step[0];
		for(int p2=0; p2<numMarginalGridPoints; ++p2){
			//for(int r=0; r<numReadGroups; ++r) params[r][1] = min[1] + p2 * step[1];
			params[0][1] = min[1] + p2 * step[1];
			for(int p3=0; p3<numMarginalGridPoints; ++p3){
				//for(int r=0; r<numReadGroups; ++r) params[r][2] = min[2] + p3 * step[2];
				params[0][2] = min[2] + p3 * step[2];
				//for(int p4=0; p4<numMarginalGridPoints; ++p4){
					//for(int r=0; r<numReadGroups; ++r) params[r][3] = min[3] + p4 * step[3];
					//for(int p5=0; p5<numMarginalGridPoints; ++p5){
						//for(int r=0; r<numReadGroups; ++r) params[r][4] = min[4] + p5 * step[4];


						//calculate LL
						LL = 0.0;
						for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
							LL += (*curWindow)->calcLL(params);
						}

						//write to file
						for(int i=0; i<5; ++i) out << params[0][i] << "\t";
						out << LL << std::endl;
					//}
				//}
			}
		}
	}

	//close file
	out.close();
}


void TRecalibrationEM::calcQSurface(std::string filename, int numMarginalGridPoints){
	double Q;

	//open outputfile
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	out << "beta0\tbeta1\tQ" << std::endl;

	//set min, max and step for each parameter
	double min[2];
	min[0] = -5.0;
	min[1] = -5.0;

	double max[2];
	max[0] = 10.0;
	max[1] = 10.0;

	double step[2];
	for(int i=0; i<2; ++i){
		step[i] = (max[i] - min[i]) / (numMarginalGridPoints - 1.0);
	}

	//print old params

	//Loop over parameters
	for(int p1=0; p1<numMarginalGridPoints; ++p1){
		for(int r=0; r<numReadGroups; ++r) newParams[r][0] = min[0] + p1 * step[0];
		for(int p2=0; p2<numMarginalGridPoints; ++p2){
			for(int r=0; r<numReadGroups; ++r) newParams[r][1] = min[1] + p2 * step[1];
			//for(int p3=0; p3<numMarginalGridPoints; ++p3){
				//for(int r=0; r<numReadGroups; ++r) params[r][2] = min[2] + p3 * step[2];
				//for(int p4=0; p4<numMarginalGridPoints; ++p4){
					//for(int r=0; r<numReadGroups; ++r) params[r][3] = min[3] + p4 * step[3];
					//for(int p5=0; p5<numMarginalGridPoints; ++p5){
						//for(int r=0; r<numReadGroups; ++r) params[r][4] = min[4] + p5 * step[4];


						//calculate Q
						Q = 0.0;
						//logfile->listFlush("Calculating Q at {" + toString(params[0][0]) + ", " + toString(params[0][1]) + "} ...");
						for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
							Q += (*curWindow)->calcQ(newParams);
						}
						//logfile->write(" done!");
						//logfile->conclude("Current Q = " + toString(Q));

						//write to file
						for(int i=0; i<2; ++i) out << newParams[0][i] << "\t";
						out << Q << std::endl;
					//}
				//}
			//}
		}
	}

	//close file
	out.close();
}

*/

int TRecalibrationEM::getQuality(TBase* base){
	double q = getErrorRate(base);
	//transform to quality
	return makePhredInt(q);
}
//---------------------------------------------------------------
//TBQSR_cell_base BQSR
//---------------------------------------------------------------
TBQSR_cell_base::TBQSR_cell_base(){
	curEstimate = 0.01;
	estimationConverged = false;
	firstDerivative = 0.0;
	firstDerivativeSave = 0.0;
	secondDerivative = 0.0;
	secondDerivativeSave = 0.0;
	numObservations = 0.0;
	numObservationsTmp = 0.0;
	F = 0.0;
	LL = 0.0;
	myReadGroup = -1;
	store = false;
	batchSize = 100000;
	next = 0;
}

void TBQSR_cell_base::empty(){
	if(!estimationConverged){
		numObservationsTmp = 0;
		firstDerivativeSave = firstDerivative;
		secondDerivativeSave = secondDerivative;
		firstDerivative = 0.0;
		secondDerivative = 0.0;
		LL = 0.0;
	}
}

void TBQSR_cell_base::init(bool Store, int ReadGroup){
	store = Store;
	myReadGroup = ReadGroup;
}

void TBQSR_cell_base::reopenEstimation(){
	estimationConverged = false;
	empty();
}

void TBQSR_cell_base::set(float error, std::string & NumObservations){
	curEstimate = error;
	if(curEstimate <= 0.0) curEstimate = 0.000000001;
	if(curEstimate >= 1.0) curEstimate = 0.9;
	if(NumObservations == "-") numObservations = 0;
	else numObservations = pow(10.0, stringToDouble(NumObservations));
};

float TBQSR_cell_base::getD(TBase* base, Base & RefBase){
	float D = 0.0;
	switch(base->getBaseAsEnum()){
		case A: if(RefBase == A){
					D = 1.0;
					break;
				}
				if(RefBase == G) D = base->PMD_GA;
				break;
		case C: if(RefBase == C) D = 1.0 - base->PMD_CT;
				break;
		case G: if(RefBase == G) D = 1.0 - base->PMD_GA;
				break;
		case T: if(RefBase == C) D = base->PMD_CT;
		        else if(RefBase == T) D = 1.0;
				break;
		case N: throw "Can not add base with unknown reference to BQSR cell!";
	}
	return D;
}

void TBQSR_cell_base::runNewtonRaphson(float & convergenceThreshold){
	if(F != F) throw "F is not a number!";

	curEstimate = curEstimate - firstDerivative / secondDerivative;
	//decide on convergence
	F = fabs(firstDerivative / numObservations);
	if(F < convergenceThreshold) estimationConverged = true;
}


std::string TBQSR_cell_base::getNumObsForPrinting(){
	if(numObservations == 0) return "-";
	else return toString(log10(numObservations));
}

void TBQSR_cell_base::calcLikelihoodSurfaceAt(int numPositions, double* positions, std::string & tag, std::ofstream & out){
	bool estimationConvergedTmp = estimationConverged;
	estimationConverged = false;
	float curEstimateTmp = curEstimate;

	for(int i=0; i<numPositions; ++i){
		curEstimate = positions[i];
		recalculateDerivativesFromDataInMemory();
		recalculateLLFromDataInMemory();
		out << tag << "\t" << positions[i] << "\t" << LL << " \t" << firstDerivative << "\t" << secondDerivative << std::endl;
	}

	curEstimate = curEstimateTmp;
	estimationConverged = estimationConvergedTmp;
}

//---------------------------------------------------------------
//TBQSR_cell BQSR
//---------------------------------------------------------------
TBQSR_cell::TBQSR_cell():TBQSR_cell_base(){
	numMatches = 0;
	pointerToBatch = NULL;
}

void TBQSR_cell::init(float initialError, bool Store, int ReadGroup){
	TBQSR_cell_base::init(Store, ReadGroup);

	curEstimate = initialError;
	if(curEstimate <= 0.0) curEstimate = 0.000000001;
	if(curEstimate >= 1.0) curEstimate = 0.9;

	//storage
	if(store){
		D_storage.push_back(new float[batchSize]);
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
	}
}

void TBQSR_cell::empty(){
	if(!estimationConverged){
		TBQSR_cell_base::empty();
		if(!store) numMatches = 0;
	} else {
		clearStorage();
	}
}

void TBQSR_cell::clearStorage(){
	if(store){
		for(batchIt = D_storage.rbegin(); batchIt != D_storage.rend(); ++batchIt)
			delete[] *batchIt;
		D_storage.clear();
		next = 0;
	}
}

void TBQSR_cell::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		if(store){
			if(next == batchSize){
				//add new batch
				D_storage.push_back(new float[batchSize]);
				batchIt = D_storage.rbegin();
				pointerToBatch = *batchIt;
				next = 0;
			}

			//add D to batch
			pointerToBatch[next] = getD(base, RefBase);
			addToDerivatives(pointerToBatch[next]);
			++next;
		} else {
			float D = getD(base, RefBase);
			addToDerivatives(D);
		}
		++numObservationsTmp;
		if(base->getBaseAsEnum() == RefBase) ++numMatches;
	}
}

void TBQSR_cell::addToDerivatives(float & D){
	float oneMinus4D = 1.0 - 4.0 * D;
	if(curEstimate >= 1.0) curEstimate = 0.999999;
	firstDerivative += oneMinus4D / (-4.0*D*curEstimate + 3.0*D + curEstimate);
	float tmpF = oneMinus4D / (D*(3.0-4.0*curEstimate) + curEstimate);
	secondDerivative -= tmpF * tmpF;
}

void TBQSR_cell::addToLL(float & D){
	LL += log((1.0-D)*curEstimate/3.0 + D*(1.0-curEstimate));
}

void TBQSR_cell::recalculateDerivativesFromDataInMemory(){
	if(!estimationConverged){
		//set to zero
		empty();

		//first the last batch, which is not filled to the end
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
		for(int i=0; i<next; ++i){ //next is set when adding sites
			addToDerivatives(pointerToBatch[i]);
		}

		//and now the other batches
		++batchIt;
		for(;batchIt != D_storage.rend(); ++batchIt){
			pointerToBatch = *batchIt;
			for(int i=0; i<batchSize; ++i){
				addToDerivatives(pointerToBatch[i]);
			}
		}
	}
}

void TBQSR_cell::recalculateLLFromDataInMemory(){
	LL = 0.0;

	//first the last batch, which is not filled to the end
	batchIt = D_storage.rbegin();
	pointerToBatch = *batchIt;
	for(int i=0; i<next; ++i){ //next is set when adding sites
		addToLL(pointerToBatch[i]);
	}

	//and now the other batches
	++batchIt;
	for(;batchIt != D_storage.rend(); ++batchIt){
		pointerToBatch = *batchIt;
		for(int i=0; i<batchSize; ++i){
			addToLL(pointerToBatch[i]);
		}
	}
}

void TBQSR_cell::runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon){
	//need Newton-Raphson to estimate epsilon
	float oldEstimate = curEstimate;

	runNewtonRaphson(convergenceThreshold);

	//check boundaries
	if(curEstimate <= 0.0){
		curEstimate = 0.000000001;
		if(oldEstimate == 0.00000001)
			estimationConverged = true; //if estimate is repeatedly below, accept
	} else if(curEstimate >= 1.0){
		curEstimate = 0.999999999;
		if(oldEstimate == 0.999999999)
			estimationConverged = true; //if estimate is repeatedly above, accept
	}

	//do not allow big jump in quality -> max +/- 10!
	if(curEstimate / oldEstimate > 10.0){
		curEstimate = oldEstimate * 10.0;
	} else if(oldEstimate / curEstimate > 10.0){
		curEstimate = oldEstimate / 10.0;
	}

	//check if quality did not change
	if(abs(makePhred(oldEstimate) - makePhred(curEstimate)) < minEpsilon) estimationConverged = true;
}

bool TBQSR_cell::estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		if(store){
			numObservations = (D_storage.size() - 1) * batchSize + next;
		} else numObservations = numObservationsTmp;

		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
		} else if(numMatches >= numObservations){ //epsilon = 0
			curEstimate = 0.0;
			estimationConverged = true;
		} else if(numMatches < 1.0){ // epsilon = 1.0
			std::cout << "numMatches < 1" << std::endl;
			curEstimate = 1.0;
			estimationConverged = true;
		} else {
			//need Newton-Raphson to estimate epsilon
			runNewtonRaphsonAndCheck(convergenceThreshold, minEpsilon);
		}
	}

	return estimationConverged;
}

void TBQSR_cell::calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out){
	double* positions = new double[numPositions];
	//now fill between 0.000000001 and 0.999999999
	double delta = (0.999999999 - 0.000000001) / (numPositions - 1.0);
	for(int i=0; i<numPositions; ++i){
		positions[i] = 0.000000001 + delta * (double) i;
	}

	//calc surface!
	calcLikelihoodSurfaceAt(numPositions, positions, tag, out);
}

//---------------------------------------------------------------
//TBQSR_cellPosition BQSR
//---------------------------------------------------------------
TBQSR_cellPosition::TBQSR_cellPosition():TBQSR_cell_base(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
	pointerToBatch = NULL;
}

void TBQSR_cellPosition::init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cell_base::init(Store, ReadGroup);
	BQSR_cells_readGroup_quality = gotBQSR_cells_quality_readGroup;
	qualityIndex = QualityIndex;

	//storage
	if(store){
		D_storage.push_back(new BQSRFactorStorage[batchSize]);
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
	}
}

void TBQSR_cellPosition::clearStorage(){
	if(store){
		for(batchIt = D_storage.rbegin(); batchIt != D_storage.rend(); ++batchIt)
			delete[] *batchIt;
		D_storage.clear();
		next = 0;
	}
}

void TBQSR_cellPosition::addToDerivatives(float & D, float & epsilon){
	double epsMinus4Deps = epsilon - 4.0 * D * epsilon;
	firstDerivative += epsMinus4Deps / (-4.0*D*epsilon*curEstimate + 3.0*D + epsilon*curEstimate);
	double tmpF = epsMinus4Deps / (D*(3.0-4.0*epsilon*curEstimate) + epsilon*curEstimate);
	secondDerivative -= tmpF * tmpF;
}

void TBQSR_cellPosition::addToLL(float & D, float & epsilon){
	LL += log((1.0-D)*epsilon*curEstimate/3.0 + D*(1.0-epsilon*curEstimate));
}

float TBQSR_cellPosition::getEpsilon(TBase* base){
	return  BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(base->quality)].curEstimate;
}

void TBQSR_cellPosition::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		if(store){
			if(next == batchSize){
				//add new batch
				D_storage.push_back(new BQSRFactorStorage[batchSize]);
				batchIt = D_storage.rbegin();
				pointerToBatch = *batchIt;
				next = 0;
			}

			//add D to batch
			pointerToBatch[next].D = getD(base, RefBase);
			pointerToBatch[next].epsilon = getEpsilon(base);
			addToDerivatives(pointerToBatch[next].D, pointerToBatch[next].epsilon);
			++next;
		} else {
			float D = getD(base, RefBase);
			float eps = getEpsilon(base);
			addToDerivatives(D, eps);
		}
		++numObservationsTmp;
	}
}

void TBQSR_cellPosition::recalculateDerivativesFromDataInMemory(){
	if(!estimationConverged){
		//set to zero
		empty();

		//first the last batch, which is not filled to the end
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
		for(int i=0; i<next; ++i){
			addToDerivatives(pointerToBatch[i].D, pointerToBatch[i].epsilon);
		}

		//and now the other batches
		++batchIt;
		for(;batchIt != D_storage.rend(); ++batchIt){
			pointerToBatch = *batchIt;
			for(int i=0; i<batchSize; ++i){
				addToDerivatives(pointerToBatch[i].D, pointerToBatch[i].epsilon);
			}
		}
	}
}

void TBQSR_cellPosition::recalculateLLFromDataInMemory(){
	LL = 0.0;

	//first the last batch, which is not filled to the end
	batchIt = D_storage.rbegin();
	pointerToBatch = *batchIt;
	for(int i=0; i<next; ++i){
		addToLL(pointerToBatch[i].D, pointerToBatch[i].epsilon);
	}

	//and now the other batches
	++batchIt;
	for(;batchIt != D_storage.rend(); ++batchIt){
		pointerToBatch = *batchIt;
		for(int i=0; i<batchSize; ++i){
			addToLL(pointerToBatch[i].D, pointerToBatch[i].epsilon);
		}
	}
}

void TBQSR_cellPosition::runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon){
	//need Newton-Raphson to estimate epsilon
	float oldEstimate = curEstimate;
	runNewtonRaphson(convergenceThreshold);

	//check boundaries
	if(curEstimate < 0.0){
		curEstimate = 0.01;
		if(oldEstimate == 0.01)
			estimationConverged = true; //if estimate is repeatedly below, accept
	} else if(curEstimate > 10000.0){
		curEstimate = 100.0;
		if(oldEstimate == 100.0)
			estimationConverged = true; //if estimate is repeatedly above, accept
	}

	//check if quality did not change
	if(fabs(oldEstimate - curEstimate) < minEpsilon) estimationConverged = true;
}


bool TBQSR_cellPosition::estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		if(store){
			numObservations = (D_storage.size() - 1) * batchSize + next;
		} else numObservations = numObservationsTmp;

		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
		} else {
			//need Newton-Raphson to estimate epsilon
			runNewtonRaphsonAndCheck(convergenceThreshold, minEpsilon);
		}
	}

	return estimationConverged;
}

void TBQSR_cellPosition::calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out){
	double* positions = new double[numPositions];
	//now fill between 0.01 and 100 -> use log10
	double delta = 4.0 / (numPositions - 1.0);
	for(int i=0; i<numPositions; ++i){
		positions[i] = pow(10.0, -2.0 + delta * (double) i);
	}

	//calc surface!
	calcLikelihoodSurfaceAt(numPositions, positions, tag, out);
}

//---------------------------------------------------------------
TBQSR_cellPositionRev::TBQSR_cellPositionRev():TBQSR_cellPosition(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
	BQSR_cells_readGroup_position = NULL;
	considerPosition = false;
}

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position = gotBQSR_quality_position;
	considerPosition = true;
}

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position = NULL;
}

float TBQSR_cellPositionRev::getEpsilon(TBase* base){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(base->quality)].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->posInRead].curEstimate;
	return  epsilonAlpha;
}

//---------------------------------------------------------------
TBQSR_cellContext::TBQSR_cellContext():TBQSR_cellPositionRev(){
	BQSR_cells_readGroup_position_rev = NULL;
	considerPositionRev = false;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

float TBQSR_cellContext::getEpsilon(TBase* base){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(base->quality)].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->posInRead].curEstimate;
	if(considerPositionRev) epsilonAlpha *= BQSR_cells_readGroup_position_rev[myReadGroup][base->posInReadRev].curEstimate;
	return  epsilonAlpha;
}

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TLog* Logfile){
	logfile = Logfile;
	bamHeader = BamHeader;
	estimatetionRequired = false;
	estimationConverged = false;
	numContexts = 20;
	qualityIndex = NULL;
	maxPos = 0;

	storeDataInMemory = params.parameterExists("storeInMemory");
	if(storeDataInMemory) logfile->list("Will store D in memory to iterate Newton-Raphson faster");
	//if(mergeReadGroupsRecalibration) logfile->list("Pooling all read groups for BQSR recalibration");
	dataStored = false;

	initializeReadGroupMap(bamHeader, params, logfile);

	//check if BQSR table readGroup x Quality is given, or has to be estimated
	initializeBQSRReadGroupQualityTable(params);

	//Do we also consider the effect of the position in read (cycle)?
	initializeBQSRReadGroupPositionTable(params);
	initializeBQSRReadGroupPositionReverseTable(params);

	//Do we also consider the context (dinucleotide)?
	initializeBQSRReadGroupContextTable(params);

	//read Newton-Raphson arguments from user
	convergenceThreshold_F = params.getParameterDoubleWithDefault("maxF", 0.0000001);
	minEpsilonQuality = params.getParameterDoubleWithDefault("minEpsQuality", 0.000001);
	minEpsilonFactors = params.getParameterDoubleWithDefault("minEpsFactors", 0.0001);
	if(estimatetionRequired){
		logfile->startIndent("Conditions to stop Newton-Raphson algorithm:");
		logfile->list("Stopping Newton-Raphson if F < " + toString(convergenceThreshold_F));
		logfile->list("Stopping Newton-Raphson if the change in quality is < " + toString(minEpsilonQuality));
		logfile->list("Stopping Newton-Raphson if the change in a factor (e.g. position) is < " + toString(minEpsilonFactors));
		logfile->endIndent();
	}

	//get minimal number of observations to conduct estimation
	minObservations = params.getParameterLongWithDefault("minObservations", 32000);

	//do we print LL surfaces?
	numPosLLsurface = params.getParameterIntWithDefault("LLSurface", 0);
	if(numPosLLsurface > 0) printLLSurface = true;
	else printLLSurface = false;
	LLSurfacePrinted = false;
}

void TRecalibrationBQSR::initializeBQSRReadGroupQualityTable(TParameters & params){
	if(params.parameterExists("BQSRQuality")) initializeBQSRReadGroupQualityTableFromFile(params);
	else {
		qualityConverged = false;
		estimateQuality = true;
		estimatetionRequired = true;
		int minQ = params.getParameterIntWithDefault("minQ", 0);
		int maxQ = params.getParameterIntWithDefault("maxQ", 100);
		logfile->list("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
		qualityIndex = new TQualityIndex(minQ, maxQ);

		//initialize BQSR table
		BQSR_cells_readGroup_quality = new TBQSR_cell*[numReadGroups];
		for(int i=0; i<numReadGroups; ++i){
			BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
			for(int q=0; q<qualityIndex->numQ; ++q){
				BQSR_cells_readGroup_quality[i][q].init(dePhred(qualityIndex->getQuality(q)), storeDataInMemory, i);
			}
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupQualityTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRQuality");
	logfile->listFlush("Constructing BQSR readGroup x quality table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_quality = new TBQSR_cell*[origNumReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	int minQ = 100;
	int maxQ = 0;
	int q;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get min and max quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			q = stringToInt(vec[1]);
			if(q > maxQ) maxQ = q;
			if(q < minQ) minQ = q;
		}
	}

	//initialize quality index
	qualityIndex = new TQualityIndex(minQ, maxQ);

	//create corresponding objects
	for(int i=0; i<origNumReadGroups; ++i){
		BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q) BQSR_cells_readGroup_quality[i][q].init(dePhred(qualityIndex->getQuality(q)), storeDataInMemory, i);
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double quality;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				q = stringToInt(vec[1]);
				quality = stringToDouble(vec[3]);
				BQSR_cells_readGroup_quality[readGroup][qualityIndex->getIndex(q)].set(dePhred(quality), vec[4]);
			} else throw "readGroup " + vec[0] + " does not exist in BAM file header!";
		}
	}

	//set that no estimation is not required, unless asked for
	if(params.parameterExists("estimateBQSRQuality")){
		qualityConverged = false;
		estimateQuality = true;
	} else {
		qualityConverged = true;
		estimateQuality = false;
	}

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
}


void TRecalibrationBQSR::initializeBQSRReadGroupPositionTable(TParameters & params){
	if(params.parameterExists("BQSRPosition")) initializeBQSRReadGroupPositionTableFromFile(params);
	else {
		positionConverged = false;
		considerPosition = params.parameterExists("estimateBQSRPosition");
		if(considerPosition){
			estimatetionRequired = true;
			estimatePosition = true;
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions up to " + toString(maxPos));
			BQSR_cells_readGroup_position = new TBQSR_cellPosition*[numReadGroups];
			for(int r=0; r<numReadGroups; ++r){
				BQSR_cells_readGroup_position[r] = new TBQSR_cellPosition[maxPos];
				for(int p=0; p<maxPos; ++p) BQSR_cells_readGroup_position[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
			}
		} else {
			BQSR_cells_readGroup_position = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupPositionTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPosition");
	logfile->listFlush("Constructing BQSR readGroup x position table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position = new TBQSR_cellPosition*[origNumReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get max position
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			p = stringToInt(vec[1]);
			if(p > maxPos) maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[origNumReadGroups];
	for(int r=0; r<origNumReadGroups; ++r){
		BQSR_cells_readGroup_position[r] = new TBQSR_cellPosition[maxPos];
		isListed[r] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				p = stringToInt(vec[1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_position[readGroup][p-1].set(alpha, vec[4]);
				isListed[readGroup][p-1] = true;
			}
		}
	}

	//check if we miss positions
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<origNumReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPosition")){
		positionConverged = false;
		estimatePosition = true;
	} else {
		positionConverged = true;
		estimatePosition = false;
	}
	considerPosition = true;

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering positions up to " + toString(maxPos));
}


//the functions are almost identical to the other position -> put in class!
void TRecalibrationBQSR::initializeBQSRReadGroupPositionReverseTable(TParameters & params){
	if(params.parameterExists("BQSRPositionReverse")) initializeBQSRReadGroupPositionReverseTableFromFile(params);
	else {
		positionReverseConverged = false;
		considerPositionReverse = params.parameterExists("estimateBQSRPositionReverse");
		if(considerPositionReverse){
			estimatePositionReverse = true;
			estimatetionRequired = true;
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions reverse up to " + toString(maxPos));
			BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[numReadGroups];
			for(int r=0; r<numReadGroups; ++r){
				BQSR_cells_readGroup_position_reverse[r] = new TBQSR_cellPositionRev[maxPos];
				for(int p=0; p<maxPos; ++p){
					if(considerPosition) BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
					else BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
				}
			}
		} else {
			BQSR_cells_readGroup_position_reverse = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupPositionReverseTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPositionReverse");
	logfile->listFlush("Constructing BQSR readGroup x position reverse table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position reverse table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[origNumReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get max position
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			p = stringToInt(vec[1]);
			if(p > maxPos) maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[origNumReadGroups];
	for(int r=0; r<origNumReadGroups; ++r){
		BQSR_cells_readGroup_position_reverse[r] = new TBQSR_cellPositionRev[maxPos];
		isListed[r] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			if(considerPosition) BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
			else BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				p = stringToInt(vec[1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_position_reverse[readGroup][p-1].set(alpha, vec[4]);
				isListed[readGroup][p-1] = true;
			}
		}
	}

	//check if we miss positions
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<origNumReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPositionReverse")){
		positionReverseConverged = false;
		estimatePositionReverse = true;
	} else {
		positionReverseConverged = true;
		estimatePositionReverse = false;
	}
	considerPositionReverse = true;

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering positions reverse up to " + toString(maxPos));
}

void TRecalibrationBQSR::initializeBQSRReadGroupContextTable(TParameters & params){
	if(params.parameterExists("BQSRContext")) initializeBQSRReadGroupContextTableFromFile(params);
	else {
		contextConverged = false;
		considerContext = params.parameterExists("estimateBQSRContext");
		if(considerContext){
			estimateContext = true;
			estimatetionRequired = true;
			logfile->list("Considering context");
			BQSR_cells_readGroup_context = new TBQSR_cellContext*[numReadGroups];
			for(int r=0; r<numReadGroups; ++r){
				BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
				for(int c=0; c<numContexts; ++c){
					if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
					else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
					else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
					else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
				}
			}
		} else {
			BQSR_cells_readGroup_context = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupContextTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRContext");
	logfile->listFlush("Constructing BQSR readGroup x context table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x context table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_context = new TBQSR_cellContext*[origNumReadGroups];
	for(int r=0; r<origNumReadGroups; ++r){
		BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
		for(int c=0; c<numContexts; ++c){
			if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
			else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
			else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
			else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
		}
	}

	//create object to check of all contexts have been initialized!
	bool** isListed = new bool*[origNumReadGroups];
	for(int i=0; i<origNumReadGroups; ++i){
		isListed[i] = new bool[numContexts];
		for(int c=0; c<numContexts; ++c){
			isListed[i][c] = false;
		}
	}

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header
	int context;
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				context = genoMap.getContext(vec[1][0], vec[1][1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_context[readGroup][context].set(alpha, vec[4]);
				isListed[readGroup][context] = true;
			}
		}
	}

	//check if we miss contexts
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<origNumReadGroups; ++i, ++it){
		for(int c=0; c<numContexts; ++c){
			if(!isListed[i][c]) throw "Context " + genoMap.getContextString(c) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRContext")){
		contextConverged = false;
		estimateContext = true;
	} else {
		contextConverged = true;
		estimateContext = false;
	}
	considerContext = true;

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering context");
}

void TRecalibrationBQSR::addSite(TSite & site){
	if(site.referenceBase != 'N'){
		Base refBase = site.getRefBaseAsEnum();
		if(!qualityConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
					BQSR_cells_readGroup_quality[readGroupMap[(*it)->readGroup]][qualityIndex->getIndex((*it)->quality)].addBase(*it, refBase);
			}
		}
		else if(considerPosition && !positionConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				if((*it)->posInRead >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position[readGroupMap[(*it)->readGroup]][(*it)->posInRead].addBase(*it, refBase);
			}
		}
		else if(considerPositionReverse && !positionReverseConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				if((*it)->posInReadRev >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position_reverse[readGroupMap[(*it)->readGroup]][(*it)->posInReadRev].addBase(*it, refBase);
			}
		} else if(considerContext && !contextConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				BQSR_cells_readGroup_context[readGroupMap[(*it)->readGroup]][(*it)->context].addBase(*it, refBase);
			}
		}
	}
}

void TRecalibrationBQSR::recalculateDerivativesFromDataInMemory(){
	if(!qualityConverged){
		for(int r=0; r<numReadGroups; ++r){
			for(int j=0; j<qualityIndex->numQ; ++j){
				BQSR_cells_readGroup_quality[r][j].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(considerPosition && !positionConverged){
		for(int r=0; r<numReadGroups; ++r){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(considerPositionReverse && !positionReverseConverged){
		for(int r=0; r<numReadGroups; ++r){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position_reverse[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	} else if(considerContext && !contextConverged){
		for(int r=0; r<numReadGroups; ++r){
			for(int c=0; c<numContexts; ++c){
				BQSR_cells_readGroup_context[r][c].recalculateDerivativesFromDataInMemory();
			}
		}
	}
}

bool TRecalibrationBQSR::estimateEpsilon(std::string filenameTag){
	//recalc derivatives if data is in memory. Otherwise, derivatives were calculated when data was added.
	if(dataStored) recalculateDerivativesFromDataInMemory();

	//estimate epsilon, if not yet done
	estimationConverged = true;
	int numCellsNotConverged = 0;
	double maxF = 0.0;

	//readGroup x quality
	//-------------------------------------------------------
	if(!qualityConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfaceQuality(filenameTag);
			LLSurfacePrinted = true;
		}
		//now do estimation
		logfile->listFlush("Estimating epsilon for readGroup x quality table ...");
		for(int i=0; i<numReadGroups; ++i){
			for(int j=0; j<qualityIndex->numQ; ++j){
				if(!BQSR_cells_readGroup_quality[i][j].estimate(convergenceThreshold_F, minEpsilonQuality, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_quality[i][j].F > maxF) maxF = BQSR_cells_readGroup_quality[i][j].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			qualityConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			qualityConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * qualityIndex->numQ));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!qualityConverged){
			//empty all cells
			for(int i=0; i<numReadGroups; ++i){
				for(int j=0; j<qualityIndex->numQ; ++j){
					BQSR_cells_readGroup_quality[i][j].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writeQualityToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<numReadGroups; ++i){
					for(int j=0; j<qualityIndex->numQ; ++j){
						BQSR_cells_readGroup_quality[i][j].clearStorage();
					}
				}
			}
			dataStored = false;

			//what's next?
			if(!considerPosition && !considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
			LLSurfacePrinted = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for position, if not yet done
	//-------------------------------------------------------
	if(considerPosition && !positionConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfacePosition(filenameTag);
			LLSurfacePrinted = true;
		}

		//now do estimation
		logfile->listFlush("Estimating epsilon for readGroup x position table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position[i][p].estimate(convergenceThreshold_F, minEpsilonFactors, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_position[i][p].F > maxF) maxF = BQSR_cells_readGroup_position[i][p].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			positionConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			positionConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * maxPos));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!positionConverged){
			//empty all cells
			for(int i=0; i<numReadGroups; ++i){
				for(int p=0; p<maxPos; ++p){
					BQSR_cells_readGroup_position[i][p].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writePositionToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<numReadGroups; ++i){
					for(int p=0; p<maxPos; ++p){
						BQSR_cells_readGroup_position[i][p].clearStorage();
					}
				}
			}
			dataStored = false;

			//what's next?
			if(!considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
			LLSurfacePrinted = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for position reverse, if not yet done
	//-------------------------------------------------------
	if(considerPositionReverse && !positionReverseConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfaceReversePosition(filenameTag);
			LLSurfacePrinted = true;
		}

		//now do estimation
		logfile->listFlush("Estimating epsilon for readGroup x position reverse table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position_reverse[i][p].estimate(convergenceThreshold_F, minEpsilonFactors, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_position_reverse[i][p].F > maxF) maxF = BQSR_cells_readGroup_position_reverse[i][p].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			positionReverseConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			positionReverseConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * maxPos));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!positionReverseConverged){
			//empty all cells
			for(int i=0; i<numReadGroups; ++i){
				for(int p=0; p<maxPos; ++p){
					BQSR_cells_readGroup_position_reverse[i][p].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writePositionReverseToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<numReadGroups; ++i){
					for(int p=0; p<maxPos; ++p){
						BQSR_cells_readGroup_position_reverse[i][p].clearStorage();
					}
				}
			}
			dataStored = false;

			//what's next?
			if(!considerContext) estimationConverged = true;
			else estimationConverged = false;
			LLSurfacePrinted = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for context
	//-------------------------------------------------------
	if(considerContext && !contextConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfaceContext(filenameTag);
			LLSurfacePrinted = true;
		}

		//now do estimation
		logfile->listFlush("Estimating epsilon for quality x context table ...");
		for(int r=0; r<numReadGroups; ++r){
			for(int c=0; c<numContexts; ++c){
				if(!BQSR_cells_readGroup_context[r][c].estimate(convergenceThreshold_F, minEpsilonFactors, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_context[r][c].F > maxF) maxF = BQSR_cells_readGroup_context[r][c].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			contextConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			contextConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * numContexts));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!contextConverged){
			//empty all cells
			for(int r=0; r<numReadGroups; ++r){
				for(int c=0; c<numContexts; ++c){
					BQSR_cells_readGroup_context[r][c].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writeContextToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<numReadGroups; ++i){
					for(int c=0; c<numContexts; ++c){
						BQSR_cells_readGroup_context[i][c].clearStorage();
					}
				}
			}
			dataStored = false;
			estimationConverged = true;
		}
		return estimationConverged;
	}

	//return true on final convergence
	return estimationConverged;
}

void TRecalibrationBQSR::writeAllToFile(std::string filenameTag){
	//write readGroup x Quality table
	writeQualityToFile(filenameTag);
	//write readGroup x position table
	if(considerPosition){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	if(considerPositionReverse){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	if(considerContext){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSR::writeCurrentTmpTable(std::string filenameTag){
	//write readGroup x Quality table
	if(!qualityConverged) writeQualityToFile(filenameTag);

	//write readGroup x position table
	else if(considerPosition && !positionConverged){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	else if(considerPositionReverse && !positionReverseConverged){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	else if(considerContext && !contextConverged){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSR::writeQualityToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x quality table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<origNumReadGroups; ++i, ++it){
		for(int q=0; q<qualityIndex->numQ; ++q){
			out << it->ID << "\t" << qualityIndex->getQuality(q) << "\tM\t" << makePhred(BQSR_cells_readGroup_quality[readGroupMap[i]][q].curEstimate) << "\t" << BQSR_cells_readGroup_quality[readGroupMap[i]][q].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_quality[readGroupMap[i]][q].firstDerivativeSave << "\t" << BQSR_cells_readGroup_quality[readGroupMap[i]][q].secondDerivativeSave << "\t" << BQSR_cells_readGroup_quality[readGroupMap[i]][q].F << "\t" << BQSR_cells_readGroup_quality[readGroupMap[i]][q].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->write(" done!");
}

void TRecalibrationBQSR::writePositionToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x position table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tObservations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<origNumReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position[readGroupMap[i]][p].curEstimate << "\t" << BQSR_cells_readGroup_position[readGroupMap[i]][p].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_position[readGroupMap[i]][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position[readGroupMap[i]][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position[readGroupMap[i]][p].F << "\t" << BQSR_cells_readGroup_position[readGroupMap[i]][p].estimationConverged;
			out << "\n";

		}
	}
	out.close();
	logfile->write(" done!");
}

void TRecalibrationBQSR::writePositionReverseToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Reverse_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x position reverse table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tObservations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<origNumReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position_reverse[readGroupMap[i]][p].curEstimate << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMap[i]][p].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMap[i]][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMap[i]][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMap[i]][p].F << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMap[i]][p].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->write(" done!");
}

void TRecalibrationBQSR::writeContextToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x context table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tContext\tEventType\tScaling\tObservations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<origNumReadGroups; ++r, ++it){
		for(int c=0; c<numContexts; ++c){
			out << it->ID << "\t" << genoMap.getContextString(c) << "\tM\t" << BQSR_cells_readGroup_context[readGroupMap[r]][c].curEstimate << "\t" << BQSR_cells_readGroup_context[readGroupMap[r]][c].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_context[readGroupMap[r]][c].firstDerivativeSave << "\t" << BQSR_cells_readGroup_context[readGroupMap[r]][c].secondDerivativeSave << "\t" << BQSR_cells_readGroup_context[readGroupMap[r]][c].F << "\t" << BQSR_cells_readGroup_context[readGroupMap[r]][c].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->write(" done!");
}


void TRecalibrationBQSR::calculateAndPrintLLSurfaceQuality(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x quality and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tQuality\terrorRate\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<numReadGroups; ++r, ++it){
		for(int q=0; q<qualityIndex->numQ; ++q){
			BQSR_cells_readGroup_quality[r][q].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + toString(qualityIndex->getQuality(q)), out);
		}
	}
	out.close();
		logfile->write(" done!");
}

void TRecalibrationBQSR::calculateAndPrintLLSurfacePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tPosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<numReadGroups; ++r, ++it){
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position[r][p].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + toString(p+1), out);
		}
	}
	out.close();
	logfile->write(" done!");
}

void TRecalibrationBQSR::calculateAndPrintLLSurfaceReversePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_ReversePosition_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x reverse position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tReversePosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<numReadGroups; ++r, ++it){
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position_reverse[r][p].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + toString(p+1), out);
		}
	}
	out.close();
	logfile->write(" done!");
}

void TRecalibrationBQSR::calculateAndPrintLLSurfaceContext(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x context and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tContext\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<numReadGroups; ++r, ++it){
		for(int c=0; c<numContexts; ++c){
			BQSR_cells_readGroup_context[r][c].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + genoMap.getContextString(c), out);
		}
	}
	out.close();
	logfile->write(" done!");
}

bool TRecalibrationBQSR::allConverged(){
	if(!qualityConverged) return false;
	if(considerPosition && !positionConverged) return false;
	if(considerPositionReverse && !positionReverseConverged) return false;
	if(considerContext && !contextConverged) return false;
	return true;
}

void TRecalibrationBQSR::reopenEstimation(){
	//resets all cells not to have converged
	if(estimateQuality){
		for(int i=0; i<numReadGroups; ++i){
			for(int q=0; q<qualityIndex->numQ; ++q){
				BQSR_cells_readGroup_quality[i][q].reopenEstimation();
			}
		}
		qualityConverged = false;
	}

	//also for position
	if(considerPosition && estimatePosition){
		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position[i][p].reopenEstimation();
			}
		}
		positionConverged = false;
	}

	//reverse position
	if(considerPositionReverse && estimatePositionReverse){
		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position_reverse[i][p].reopenEstimation();
			}
		}
		positionReverseConverged = false;
	}

	//and context
	if(considerContext && estimateContext){
		for(int r=0; r<numReadGroups; ++r){
			for(int c=0; c<numContexts; ++c){
				BQSR_cells_readGroup_context[r][c].reopenEstimation();
			}
		}
		contextConverged = false;
	}
}

double TRecalibrationBQSR::getErrorRate(TBase* base){
	double q = BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate;
	if(considerPosition) q *= BQSR_cells_readGroup_position[base->readGroup][base->posInRead].curEstimate;
	if(considerPositionReverse) q *= BQSR_cells_readGroup_position_reverse[base->readGroup][base->posInReadRev].curEstimate;
	if(considerContext) q *= BQSR_cells_readGroup_context[base->readGroup][base->context].curEstimate;
	if(q > 1.0) q = 1.0; //make sure the scaling does not lead to errors > 1.0!
	return q;
}

int TRecalibrationBQSR::getQuality(TBase* base){
	double q = getErrorRate(base);
	//transform to quality
	return makePhredInt(q);
}

