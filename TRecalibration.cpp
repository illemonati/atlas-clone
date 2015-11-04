/*
 * TRecalibration.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#include "TRecalibration.h"

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
TRecalibrationEM::TRecalibrationEM(TParameters* arguments, TLog* logfile){
	numParams = 3; //3 beta and 4 gamma
	params = new double[numParams];
	newParams = new double[numParams];
	Jacobian.resize(numParams,numParams);
	Jacobian.zeros();
	F.resize(numParams);
	F.zeros();
	JxF.resize(numParams, numParams);
	JxF.zeros();
	maxF = 0.0;
	numSitesAdded = 0;
	logLikelihood = 0.0;

	//set initial values
	//betas
	logfile->startIndent("Will start recalibration EM with:");
	params[0] = arguments->getParameterDoubleWithDefault("initBeta0", 0);
	logfile->list("beta0 = " + toString(params[0]));
	params[1] = arguments->getParameterDoubleWithDefault("initBeta1", 1);
	logfile->list("beta1 = " + toString(params[1]));
	params[2] = arguments->getParameterDoubleWithDefault("initBeta2", 0);
	logfile->list("beta2 = " + toString(params[2]));

	//gammas -> start all at 0 by default
	params[3] = arguments->getParameterDoubleWithDefault("initGammaA", 0);
	logfile->list("gammaA = " + toString(params[3]));
	params[4] = arguments->getParameterDoubleWithDefault("initGammaC", 0);
	logfile->list("gammaC = " + toString(params[4]));
	params[5] = arguments->getParameterDoubleWithDefault("initGammaG", 0);
	logfile->list("gammaG = " + toString(params[5]));
	params[6] = arguments->getParameterDoubleWithDefault("initGammaT", 0);
	logfile->list("gammaT = " + toString(params[6]));
	logfile->endIndent();
};

TRecalibrationEM::TRecalibrationEM(TLog* logfile){
	numParams = 3; //3 beta and 4 gamma
	params = new double[numParams];
	newParams = new double[numParams];
	Jacobian.resize(numParams,numParams);
	Jacobian.zeros();
	F.resize(numParams);
	F.zeros();
	JxF.resize(numParams, numParams);
	JxF.zeros();
	maxF = 0.0;
	numSitesAdded = 0;
	logLikelihood = 0.0;

	//init params
	for(int i=0; i<numParams; ++i){
		params[i] = 0.0;
		newParams[i] = 0.0;
	}
}

void TRecalibrationEM::setParams(double* Params){
	for(int i=0; i<numParams; ++i){
		params[i] = Params[i];
	}
}

double TRecalibrationEM::calcEta(TBase* base){
	return calcEta(base, params);
}

double TRecalibrationEM::calcEta(TBase* base, double* theseParams){
	//function transform log error

	throw "ERROR in TRecalibrationEM::calcEta!";
	//TODO: find other way to calculate transformedLogError and add back what is commented out below

	double eta = theseParams[0]; // + theseParams[1] * base->transformedLogError + theseParams[2] * (double) base->posInRead;
	//eta += theseParams[base->getBaseAsEnum() + 3];
	return eta;
}

double TRecalibrationEM::calcEpsilon(const double & eta){
	double tmp = exp(eta);
	return tmp / (1.0 + tmp);
}

double TRecalibrationEM::calcEpsilon(TBase* base){
	return calcEpsilon(base, params);
}

double TRecalibrationEM::calcEpsilon(TBase* base, double* theseParams){
	double tmp = calcEta(base, theseParams);
	return calcEpsilon(tmp);
}

void TRecalibrationEM::initEMStep(){
	//fill vector of new params by copying current values
	for(int i=0; i<numParams; ++i){
		newParams[i] = params[i];
	}
}

void TRecalibrationEM::initNetwonRalphsonStep(){
	Jacobian.zeros();
	F.zeros();
	numSitesAdded = 0;
}

void TRecalibrationEM::saveParams(){
	for(int i=0; i<numParams; ++i){
		params[i] = newParams[i];
	}
}

void TRecalibrationEM::addSiteToJacobianAndF(std::vector<TBase*> & bases, TBaseFrequencies* freqs){
	//adds terms to Jacobian for one site (hence a vector of bases that were read)
	//assumes bases to be haploid! -> program does not throw an error if they are not!

	//TODO: find different way to get a log transformed error.

	/*
	double epsilon;
	double epsilonThird;

	//initialize
	double P_d_given_g_theta[4];
	for(int g=0; g<4; ++g){ //over all genotypes A, C, G or T
		P_d_given_g_theta[g] = 1.0;
	}

	//calc P(d|g, theta) for all g
	for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
		//get epsilon
		epsilon  = calcEpsilon(*it); //use current (= old) params for this step
		epsilonThird = epsilon / 3.0;

		for(int g=0; g<4; ++g){ //over all genotypes
			//add to P(d|g, theta)
			if((*it)->getBaseAsEnum() == g)
				P_d_given_g_theta[g] *= (1.0 - epsilon);
			else
				P_d_given_g_theta[g] *= epsilonThird;
		}
	}

	//calculate P(g|d, theta)
	double P_g_given_d_theta_denominator = 0.0;
	for(int g=0; g<4; ++g){
		P_d_given_g_theta[g] *= (*freqs)[g];
		P_g_given_d_theta_denominator += P_d_given_g_theta[g];
	}
	double P_g_given_d_theta[4];
	for(int g=0; g<4; ++g){
		P_g_given_d_theta[g] = P_d_given_g_theta[g] / P_g_given_d_theta_denominator;
	}

	//now add site to Jacobian and F
	//Note: Jacobian is symmetric! But we only add to top triangle -> will copy final numbers down later
	int g, nucIndex;
	double eta, expEta, tmp, onePlusExpEta, weightFcorrect, weightFincorrect;
	for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
		eta = calcEta(*it, newParams); //use new params for this step!
		expEta = exp(eta);
		onePlusExpEta = 1.0 + expEta;
		tmp = expEta / (onePlusExpEta * onePlusExpEta);

		Jacobian(0,0) -= tmp;
		Jacobian(0,1) -= tmp * (*it)->transformedLogError;
		Jacobian(0,2) -= tmp * (double) (*it)->posInRead;
		Jacobian(1,1) -= tmp * (*it)->transformedLogError * (*it)->transformedLogError;
		Jacobian(1,2) -= tmp * (*it)->transformedLogError * (double) (*it)->posInRead;
		Jacobian(2,2) -= tmp * (double) (*it)->posInRead * (double) (*it)->posInRead;

		//now add to those with indicators on the base
		//Note that the derivative is zero for all other bases
		g = (*it)->getBaseAsEnum();
		nucIndex = g + 3;
		Jacobian(nucIndex, nucIndex) -= tmp;
		Jacobian(0, nucIndex) -= tmp;
		Jacobian(1, nucIndex) -= tmp * (*it)->transformedLogError;
		Jacobian(2, nucIndex) -= tmp * (*it)->posInRead;

		//now add to F -> need to add for all genotypes, but with different term
		weightFcorrect = -1.0 / (1.0 + exp(-eta));
		weightFincorrect = 1.0 / onePlusExpEta;
		for(int geno = 0; geno < 4; ++geno){
			//prepare weight
			if(geno == g) tmp = weightFcorrect;
			else tmp = weightFincorrect;
			tmp = P_g_given_d_theta[geno] * tmp;

			//fill for beta0, beta1 and beta2
			F(0) += tmp; //beta0
			F(1) += tmp * (*it)->transformedLogError;
			F(2) += tmp * (double) (*it)->posInRead;
		}
		//derivation is only non-zero for beta3 to beta6 if g corresponds to index
		//F(nucIndex) += P_g_given_d_theta[g] * weightFcorrect;
	}

	++numSitesAdded;
	*/
}

void TRecalibrationEM::runNewtonRalphson(){
	//Need to copy numbers to other triangle in Jacobian, as only upper triangled is filled when parsing sites
	for(int i=0; i<(numParams-1); ++i){
		for(int j=i+1; j<numParams; ++j){
			//copy from upper triangle to lower triangle
			Jacobian(j,i) = Jacobian(i,j);
		}
	}

	//scale F and J by 1/#sites
	Jacobian = Jacobian / (double) numSitesAdded;
	F = F / (double) numSitesAdded;

	std::cout << "JACOBIAN = " << std::endl;
	std::cout << Jacobian << std::endl;

	std::cout << std::endl << "DET = " << det(Jacobian) << std::endl;

	std::cout << "F = " << std::endl;
		std::cout << F << std::endl;

	//get largest gradient (F) to check if we break
	maxF = 0.0;
	for(int i=0; i<(numParams-1); ++i){
		if(F(i) > maxF) maxF = F(i);
	}
	std::cout << "MAX F = " << maxF << std::endl;

	//now calculate J^-1 x F
	std::cout << "New Params:";

	if(solve(JxF, Jacobian, F)){

		std::cout << "JxF = " << std::endl;
		std::cout << JxF << std::endl;

		for(int i=0; i<numParams; ++i){
			newParams[i] -= JxF(i);
			std::cout << "\t" << newParams[i];
		}
	} else throw "Issue solving JxF in TRecalibrationEM::runNewtonRalphson()!";
	std::cout << std::endl;
}

void TRecalibrationEM::writeHeader(std::ofstream & out){
	out << "beta0\tbeta1\tbeta2";
}

void TRecalibrationEM::writeParams(std::ofstream & out){
	out << params[0];
	for(int i=1; i<numParams; ++i){
		out << "\t" << params[i];
	}
}

void TRecalibrationEM::resetLikelihood(){
	logLikelihood = 0.0;
}

void TRecalibrationEM::addSiteToLikelihood(std::vector<TBase*> & bases, TBaseFrequencies* freqs){
	//adds this site (hence a vector of bases that were read) to the log likelihood
	//assumes bases to be haploid! -> program does not throw an error if they are not!
	double epsilon;
	double epsilonThird;
	double* likelihood = new double[4]; //for each genotype
	for(int i=0; i<4; ++i) likelihood[i] = 1.0;

	//calc likelihood term
	for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
		//get epsilon
		epsilon  = calcEpsilon(*it); //use current params for this step
		epsilonThird = epsilon / 3.0;

		for(int g=0; g<4; ++g){ //over all genotypes
			if((*it)->getBaseAsEnum() == g){
				likelihood[g] *= (1.0 - epsilon);
			} else {
				likelihood[g] *= epsilonThird;
			}
		}
	}

	//build weighted sum across genotypes
	double ll = 0.0;
	for(int g=0; g<4; ++g){
		ll += (*freqs)[g] * likelihood[g];
	}
	logLikelihood += log(ll);
}

//---------------------------------------------------------------
//TRecalibrationBQSR_cell BQSR
//---------------------------------------------------------------
TBQSR_cell::TBQSR_cell(){
	curEstimate = 0.0;
	estimationConverged = false;
	firstDerivative = 0.0;
	firstDerivativeSave = 0.0;
	secondDerivative = 0.0;
	secondDerivativeSave = 0.0;
	numObservations = 0.0;
	numObservationsTmp = 0.0;
	numMatches = 0.0;
	F = 0.0;
	LL = 0.0;
}

void TBQSR_cell::init(double initialError){
	curEstimate = initialError;
}

void TBQSR_cell::set(double error, std::string & NumObservations){
	curEstimate = error;
	if(NumObservations == "-") numObservations = 0;
	else numObservations = pow(10.0, stringToDouble(NumObservations));
};

void TBQSR_cell::empty(){
	if(!estimationConverged){
		numObservationsTmp = 0;
		numMatches = 0;
		firstDerivativeSave = firstDerivative;
		secondDerivativeSave = secondDerivative;
		firstDerivative = 0.0;
		secondDerivative = 0.0;
		LL = 0.0;
	}
}

void TBQSR_cell::reopenEstimation(){
	estimationConverged = false;
	empty();
}

double TBQSR_cell::getD(TBase* base, Base & RefBase){
	double D = 0.0;
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

void TBQSR_cell::addToDerivatives(TBase* base, Base & RefBase){
	double D = getD(base, RefBase);

	double oneMinus4D = 1.0 - 4.0 * D;
	firstDerivative += oneMinus4D / (-4.0*D*curEstimate + 3.0*D + curEstimate);
	double tmp = oneMinus4D / (D*(3.0-4.0*curEstimate) + curEstimate);
	secondDerivative -= tmp * tmp;

	++numObservationsTmp;
}


void TBQSR_cell::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		addToDerivatives(base, RefBase);
		if(base->getBaseAsEnum() == RefBase) ++numMatches;
	}
}

void TBQSR_cell::runNewtonRalphson(double & convergenceThreshold){
	curEstimate = curEstimate - firstDerivative / secondDerivative;
	//decide on convergence
	F = fabs(firstDerivative);
	if(F < convergenceThreshold) estimationConverged = true;
}

bool TBQSR_cell::estimate(double & convergenceThreshold, long & minObservations){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		numObservations = numObservationsTmp;

		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
		} else if(numMatches >= numObservations){ //epsilon = 0
			curEstimate = 0.0;
			estimationConverged = true;
		} else if(numMatches < 1.0){ // epsilon = 1.0
			curEstimate = 1.0;
			estimationConverged = true;
		} else {
			//need Newton-Ralphson to estimate epsilon
			double oldEstimate = curEstimate;
			runNewtonRalphson(convergenceThreshold);

			//check boundaries
			if(curEstimate < 0.0){
				curEstimate = 0.000000001;
				if(oldEstimate == 0.00000001)
					estimationConverged = true; //if estimate is repeatedly below, accept
			} else if(curEstimate > 1.0){
				curEstimate = 0.999999999;
				if(oldEstimate == 0.999999999)
					estimationConverged = true; //if estimate is repeatedly above, accept
			}
		}
	}
	return estimationConverged;
}

std::string TBQSR_cell::getNumObsForPrinting(){
	if(numObservations == 0) return "-";
	else return toString(log10(numObservations));
}
//---------------------------------------------------------------
TBQSR_cellPosition::TBQSR_cellPosition():TBQSR_cell(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
}

void TBQSR_cellPosition::init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex){
	BQSR_cells_readGroup_quality = gotBQSR_cells_quality_readGroup;
	qualityIndex = QualityIndex;
}

void TBQSR_cellPosition::addToDerivatives(TBase* base, Base & RefBase, double & epsilon){
	double D = getD(base, RefBase);

	double epsMinus4Deps = epsilon - 4.0 * D * epsilon;
	firstDerivative += epsMinus4Deps / (-4.0*D*epsilon*curEstimate + 3.0*D + epsilon*curEstimate);
	double tmp = epsMinus4Deps / (D*(3.0-4.0*epsilon*curEstimate) + epsilon*curEstimate);
	secondDerivative -= tmp * tmp;

	++numObservationsTmp;

	//LL += log((1.0-D) * epsilon*curEstimate / 3.0 + D * (1.0 - epsilon*curEstimate));
}

void TBQSR_cellPosition::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		addToDerivatives(base, RefBase, BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate);
	}
}

bool TBQSR_cellPosition::estimate(double & convergenceThreshold, long & minObservations){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		numObservations = numObservationsTmp;

		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
			return estimationConverged;
		} else {
			//need Newton-Ralphson to estimate epsilon
			double oldEstimate = curEstimate;
			runNewtonRalphson(convergenceThreshold);

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
		}
	}
	return estimationConverged;
}

//---------------------------------------------------------------
TBQSR_cellPositionRev::TBQSR_cellPositionRev():TBQSR_cellPosition(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
	BQSR_cells_readGroup_position = NULL;
	considerPosition = false;
}

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position = gotBQSR_quality_position;
	considerPosition = true;
}

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position = NULL;
}

void TBQSR_cellPositionRev::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		double epsilonAlpha = BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate;
		if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[base->readGroup][base->posInRead].curEstimate;
		addToDerivatives(base, RefBase, epsilonAlpha);
	}
}

//---------------------------------------------------------------
TBQSR_cellContext::TBQSR_cellContext():TBQSR_cellPositionRev(){
	BQSR_cells_readGroup_position_rev = NULL;
	considerPositionRev = false;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}


void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		double epsilonAlpha = BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate;
		if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[base->readGroup][base->posInRead].curEstimate;
		if(considerPositionRev) epsilonAlpha *= BQSR_cells_readGroup_position_rev[base->readGroup][base->posInReadRev].curEstimate;
		addToDerivatives(base, RefBase, epsilonAlpha);
	}
}

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TLog* Logfile){
	logfile = Logfile;
	bamHeader = BamHeader;
	numReadGroups = bamHeader->ReadGroups.Size();
	estimatetionRequired = false;
	estimationConverged = false;
	numContexts = 20;
	qualityIndex = NULL;
	maxPos = 0;

	//check if BQSR table readGroup x Quality is given, or has to be estimated
	initializeBQSRReadGroupQualityTable(params);

	//Do we also consider the effect of the position in read (cycle)?
	initializeBQSRReadGroupPositionTable(params);
	initializeBQSRReadGroupPositionReverseTable(params);

	//Do we also consider the context (dinucleotide)?
	initializeBQSRReadGroupContextTable(params);

	//read Newton-Ralphson arguments from user
	convergenceThreshold = params.getParameterDoubleWithDefault("maxF", 0.0001);
	if(estimatetionRequired) logfile->list("Stopping Newton-Ralphson if F < " + toString(convergenceThreshold));

	//get minimal number of observations to conduct estimation
	minObservations = params.getParameterLongWithDefault("minObservations", 1000);

	//---------------------------
	/*
	surfaceCalculated = false;
	numLLSurfacePoints = 500;
	LLSurface = new TBQSR_cellPosition[numLLSurfacePoints];
	double delta = 5.0 / (numLLSurfacePoints - 1.0);
	for(int i=0; i<numLLSurfacePoints; ++i){
		LLSurface[i].init(BQSR_cells_readGroup_quality, qualityIndex, pmdObject);
		LLSurface[i].set(0.1 + i*delta);
	}
	*/
	//---------------------------
}

int TRecalibrationBQSR::findReadGroupIndex(std::string & name){
	int i = 0;
	for(std::vector<BamTools::SamReadGroup>::iterator it = bamHeader->ReadGroups.Begin(); it !=  bamHeader->ReadGroups.End(); ++it, ++i){
		if(name == it->ID) return i;
	}
	return -1;
}

void TRecalibrationBQSR::initializeBQSRReadGroupQualityTable(TParameters & params){
	if(params.parameterExists("BQSRQuality")) initializeBQSRReadGroupQualityTableFromFile(params);
	else {
		qualityConverged = false;
		estimatetionRequired = true;
		int minQ = params.getParameterIntWithDefault("minQ", 1);
		int maxQ = params.getParameterIntWithDefault("maxQ", 100);
		logfile->list("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
		qualityIndex = new TQualityIndex(minQ, maxQ);

		//initialize BQSR table
		BQSR_cells_readGroup_quality = new TBQSR_cell*[numReadGroups];
		for(int i=0; i<numReadGroups; ++i){
			BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
			for(int q=0; q<qualityIndex->numQ; ++q) BQSR_cells_readGroup_quality[i][q].init(dePhred(qualityIndex->getQuality(q)));
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupQualityTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRQuality");
	logfile->listFlush("Constructing BQSR readGroup x quality table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_quality = new TBQSR_cell*[numReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	int minQ = 100;
	int maxQ = 0;
	int q;
	std::string tmp;
	std::getline(file, tmp); //skip header

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
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q) BQSR_cells_readGroup_quality[i][q].init(0.01);
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmp); //skip header
	double quality;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0]);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				q = stringToInt(vec[1]);
				quality = stringToDouble(vec[3]);
				BQSR_cells_readGroup_quality[readGroup][qualityIndex->getIndex(q)].set(dePhred(quality), vec[4]);
			}
		}
	}

	//set that no estimation is not required, unless asked for
	if(params.parameterExists("estimateBQSRQuality")) qualityConverged = false;
	else qualityConverged = true;

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
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions up to " + toString(maxPos));
			BQSR_cells_readGroup_position = new TBQSR_cellPosition*[numReadGroups];
			for(int i=0; i<numReadGroups; ++i){
				BQSR_cells_readGroup_position[i] = new TBQSR_cellPosition[maxPos];
				for(int p=0; p<maxPos; ++p) BQSR_cells_readGroup_position[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
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
	BQSR_cells_readGroup_position = new TBQSR_cellPosition*[numReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmp;
	std::getline(file, tmp); //skip header

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
	bool** isListed = new bool*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_readGroup_position[i] = new TBQSR_cellPosition[maxPos];
		isListed[i] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
			isListed[i][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmp); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0]);
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
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPosition")) positionConverged = false;
	else positionConverged = true;
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
			estimatetionRequired = true;
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions reverse up to " + toString(maxPos));
			BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[numReadGroups];
			for(int i=0; i<numReadGroups; ++i){
				BQSR_cells_readGroup_position_reverse[i] = new TBQSR_cellPositionRev[maxPos];
				for(int p=0; p<maxPos; ++p){
					if(considerPosition) BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
					else BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
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
	BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[numReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmp;
	std::getline(file, tmp); //skip header

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
	bool** isListed = new bool*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_readGroup_position_reverse[i] = new TBQSR_cellPositionRev[maxPos];
		isListed[i] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			if(considerPosition) BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
			else BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
			isListed[i][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmp); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0]);
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
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPositionReverse")) positionReverseConverged = false;
	else positionReverseConverged = true;
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
			estimatetionRequired = true;
			logfile->list("Considering context");
			BQSR_cells_readGroup_context = new TBQSR_cellContext*[numReadGroups];
			for(int r=0; r<numReadGroups; ++r){
				BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
				for(int c=0; c<numContexts; ++c){
					if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex);
					else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
					else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex);
					else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex);
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
	BQSR_cells_readGroup_context = new TBQSR_cellContext*[numReadGroups];
	for(int r=0; r<numReadGroups; ++r){
		BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
		for(int c=0; c<numContexts; ++c){
			if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex);
			else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
			else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex);
			else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex);
		}
	}

	//create object to check of all contexts have been initialized!
	bool** isListed = new bool*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		isListed[i] = new bool[numContexts];
		for(int c=0; c<numContexts; ++c){
			isListed[i][c] = false;
		}
	}

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmp;
	std::getline(file, tmp); //skip header
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
			readGroup = findReadGroupIndex(vec[0]);
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
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int c=0; c<numContexts; ++c){
			if(!isListed[i][c]) throw "Context " + genoMap.getContextString(c) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRContext")) contextConverged = false;
	else contextConverged = true;
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
				BQSR_cells_readGroup_quality[(*it)->readGroup][qualityIndex->getIndex((*it)->quality)].addBase(*it, refBase);
			}
		}
		else if(considerPosition && !positionConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				if((*it)->posInRead >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position[(*it)->readGroup][(*it)->posInRead].addBase(*it, refBase);
			}
		}
		else if(considerPositionReverse && !positionReverseConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				if((*it)->posInReadRev >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position_reverse[(*it)->readGroup][(*it)->posInReadRev].addBase(*it, refBase);
			}
		} else if(considerContext && !contextConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				BQSR_cells_readGroup_context[(*it)->readGroup][(*it)->context].addBase(*it, refBase);
			}
		}
	}
}

bool TRecalibrationBQSR::estimateEpsilon(){
	//estimate epsilon for quality and readgroup, if not yet done
	estimationConverged = true;
	int numCellsNotConverged = 0;
	double maxF = 0.0;

	//readGroup x quality
	//-------------------------------------------------------
	if(!qualityConverged){
		logfile->listFlush("Estimating epsilon for readGroup x quality table ...");
		for(int i=0; i<numReadGroups; ++i){
			for(int j=0; j<qualityIndex->numQ; ++j){
				if(!BQSR_cells_readGroup_quality[i][j].estimate(convergenceThreshold, minObservations)){
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
		} else {
			if(!considerPosition && !considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for position, if not yet done
	//-------------------------------------------------------
	if(considerPosition && !positionConverged){
		logfile->listFlush("Estimating epsilon for readGroup x position table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position[i][p].estimate(convergenceThreshold, minObservations)){
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
		} else {
			if(!considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for position reverse, if not yet done
	//-------------------------------------------------------
	if(considerPositionReverse && !positionReverseConverged){
		logfile->listFlush("Estimating epsilon for readGroup x position reverse table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position_reverse[i][p].estimate(convergenceThreshold, minObservations)){
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
		} else {
			if(!considerContext) estimationConverged = true;
			else estimationConverged = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for context
	//-------------------------------------------------------
	if(considerContext && !contextConverged){
		logfile->listFlush("Estimating epsilon for quality x context table ...");
		for(int r=0; r<numReadGroups; ++r){
			for(int c=0; c<numContexts; ++c){
				if(!BQSR_cells_readGroup_context[r][c].estimate(convergenceThreshold, minObservations)){
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
		} else {
			estimationConverged = true;
		}
		return estimationConverged;
	}

	//return true on final convergence
	return estimationConverged;
}

void TRecalibrationBQSR::writeToFile(std::string filenameTag){
	//write readGroup x Quality table
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x quality table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations";
	out << "\tfirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int q=0; q<qualityIndex->numQ; ++q){
			out << it->ID << "\t" << qualityIndex->getQuality(q) << "\tM\t" << makePhred(BQSR_cells_readGroup_quality[i][q].curEstimate) << "\t" << BQSR_cells_readGroup_quality[i][q].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_quality[i][q].firstDerivativeSave << "\t" << BQSR_cells_readGroup_quality[i][q].secondDerivativeSave << "\t" << BQSR_cells_readGroup_quality[i][q].F << "\t" << BQSR_cells_readGroup_quality[i][q].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->write(" done!");

	//write readGroup x position table
	if(considerPosition){
		filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
		logfile->listFlush("Writing BQSR readGroup x position table to '" + filename + "' ...");
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open file '" + filename + "' for writing!";
		out << "ReadGroup\tPosition\tEventType\tScaling\tObservations\n";
		BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
		for(int i=0; i<numReadGroups; ++i, ++it){
			for(int p=0; p<maxPos; ++p){
				out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position[i][p].curEstimate << "\t" << BQSR_cells_readGroup_position[i][p].getNumObsForPrinting();
				//for debugging: also print derivatives, F and whether is has converged
				out << "\t" << BQSR_cells_readGroup_position[i][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position[i][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position[i][p].F << "\t" << BQSR_cells_readGroup_position[i][p].estimationConverged;
				out << "\n";
			}
		}
		out.close();
		logfile->write(" done!");
	}

	//write readGroup x position table
	if(considerPositionReverse){
		filename = filenameTag + "_BQSR_ReadGroup_Position_Reverse_Table.txt";
		logfile->listFlush("Writing BQSR readGroup x position reverse table to '" + filename + "' ...");
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open file '" + filename + "' for writing!";
		out << "ReadGroup\tPosition\tEventType\tScaling\tObservations\n";
		BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
		for(int i=0; i<numReadGroups; ++i, ++it){
			for(int p=0; p<maxPos; ++p){
				out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position_reverse[i][p].curEstimate << "\t" << BQSR_cells_readGroup_position_reverse[i][p].getNumObsForPrinting();
				//for debugging: also print derivatives, F and whether is has converged
				out << "\t" << BQSR_cells_readGroup_position_reverse[i][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[i][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[i][p].F << "\t" << BQSR_cells_readGroup_position_reverse[i][p].estimationConverged;
				out << "\n";
			}
		}
		out.close();
		logfile->write(" done!");
	}


	//write readGroup x context table
	if(considerContext){
		filename = filenameTag + "_BQSR_ReadGroup_Context_Table.txt";
		logfile->listFlush("Writing BQSR readGroup x text table to '" + filename + "' ...");
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open file '" + filename + "' for writing!";
		out << "ReadGroup\tContext\tEventType\tScaling\tObservations\n";
		it = bamHeader->ReadGroups.Begin();
		for(int r=0; r<numReadGroups; ++r, ++it){
			for(int c=0; c<numContexts; ++c){
				out << it->ID << "\t" << genoMap.getContextString(c) << "\tM\t" << BQSR_cells_readGroup_context[r][c].curEstimate << "\t" << BQSR_cells_readGroup_context[r][c].getNumObsForPrinting();
				//for debugging: also print derivatives, F and whether is has converged
				out << "\t" << BQSR_cells_readGroup_context[r][c].firstDerivativeSave << "\t" << BQSR_cells_readGroup_context[r][c].secondDerivativeSave << "\t" << BQSR_cells_readGroup_context[r][c].F << "\t" << BQSR_cells_readGroup_context[r][c].estimationConverged;
				out << "\n";
			}
		}
		out.close();
		logfile->write(" done!");
	}
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
	for(int i=0; i<numReadGroups; ++i){
		for(int q=0; q<qualityIndex->numQ; ++q){
			BQSR_cells_readGroup_quality[i][q].reopenEstimation();
		}
	}

	//also for position
	if(considerPosition){
		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position[i][p].reopenEstimation();
			}
		}
	}

	//and context
	if(considerContext){
		for(int r=0; r<numReadGroups; ++r){
			for(int c=0; c<numContexts; ++c){
				BQSR_cells_readGroup_context[r][c].reopenEstimation();
			}
		}
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

