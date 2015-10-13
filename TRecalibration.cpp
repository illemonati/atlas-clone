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
	/*
	params[3] = arguments->getParameterDoubleWithDefault("initGammaA", 0);
	logfile->list("gammaA = " + toString(params[3]));
	params[4] = arguments->getParameterDoubleWithDefault("initGammaC", 0);
	logfile->list("gammaC = " + toString(params[4]));
	params[5] = arguments->getParameterDoubleWithDefault("initGammaG", 0);
	logfile->list("gammaG = " + toString(params[5]));
	params[6] = arguments->getParameterDoubleWithDefault("initGammaT", 0);
	logfile->list("gammaT = " + toString(params[6]));
	*/
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
	//function transfrom log error
	double eta = theseParams[0] + theseParams[1] * base->transformedLogError + theseParams[2] * (double) base->posInRead;
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
		/*
		Jacobian(nucIndex, nucIndex) -= tmp;
		Jacobian(0, nucIndex) -= tmp;
		Jacobian(1, nucIndex) -= tmp * (*it)->transformedLogError;
		Jacobian(2, nucIndex) -= tmp * (*it)->posInRead;
		 */

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
TRecalibrationBQSR_cell::TRecalibrationBQSR_cell(){
	curEstimate = 0.0;
	estimationConverged = false;
	firstDerivative = 0;
	secondDerivative = 0;
	pmdObject = NULL;
	numObservations = 0;
}

void TRecalibrationBQSR_cell::init(double initialError, TPMD* PmdObject){
	curEstimate = initialError;
	pmdObject = PmdObject;
}


void TRecalibrationBQSR_cell::fillD(){

}

void TRecalibrationBQSR_cell::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		double D = 0.0;

		switch(base->getBaseAsEnum()){
			case A: if(RefBase == A) D = 1.0; break;
					if(RefBase == G) D = pmdObject->getProbGA(base->pos3); break;
			case C: if(RefBase == C) D = pmdObject->getProbGA(base->pos5); break;
			case G: if(RefBase == G) D = 1.0 - pmdObject->getProbGA(base->pos3); break;
			case T: if(RefBase == C) D = pmdObject->getProbGA(base->pos5); break;
					if(RefBase == T) D = 1.0; break;
			case N: throw "Can not add base with unknown reference to BQSR cell!";
		}

		double oneMinus4D = 1.0 - 4.0 * D;
		firstDerivative += oneMinus4D / (-4.0*D*curEstimate + 3.0*D + curEstimate);
		double tmp = oneMinus4D / (D*(3.0-4.0*curEstimate) + curEstimate);
		secondDerivative += tmp * tmp;

		++numObservations;
	}
}

bool TRecalibrationBQSR_cell::estimate(double & convergenceThreshold){
	if(estimationConverged) return estimationConverged;

	//Need Newton-Ralphson to estimate epsilon
	//Make new estimate
	curEstimate = curEstimate + firstDerivative / secondDerivative;

	//decide on convergence
	if(abs(firstDerivative) < convergenceThreshold) estimationConverged = true;
	return estimationConverged;
}

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TPMD* PmdObject, TLog* Logfile){
	logfile = Logfile;
	bamHeader = BamHeader;
	pmdObject = PmdObject;
	numReadGroups = bamHeader->ReadGroups.Size();
	estimationConverged = false;
	qualityConverged = false;
	positionConverged = false;

	//read arguments from user
	minQ = params.getParameterIntWithDefault("minQ", 1);
	maxQ = params.getParameterIntWithDefault("maxQ", 100);
	numQ = maxQ - minQ + 1;
	initialError = params.getParameterDoubleWithDefault("initialError", 0.01);
	convergenceThreshold = params.getParameterDoubleWithDefault("maxF", 0.0001);

	//read which covariates to consider
	considerPosition = params.parameterExists("considerPosition");

	//initialize BQSR table
	BQSR_cells_quality = new TRecalibrationBQSR_cell*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_quality[i] = new TRecalibrationBQSR_cell[numQ];
		for(int q=0; q<numQ; ++q) BQSR_cells_quality[i][q].init(initialError, pmdObject);
	}

	//Do we also consider the effect of the positon in read (cycle)?
	considerPosition = params.parameterExists("considerPosition");
	if(considerPosition){
		maxPos = params.getParameterInt("maxPos");
		if(maxPos < 0) throw "Max position has to be larger than zero!";
		BQSR_cells_position = new TRecalibrationBQSR_cellCovariate*[numQ];
		for(int q=0; q<numQ; ++q){
			BQSR_cells_position[q] = new TRecalibrationBQSR_cellCovariate[maxPos];
			for(int p=0; p<maxPos; ++p) BQSR_cells_position[q][p].init(pmdObject);
		}
	}
}

void TRecalibrationBQSR::addSite(TSite & site){
	if(site.referenceBase != 'N'){
		Base refBase = site.getRefBaseAsEnum();
		for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
			BQSR_cells_quality[(*it)->readGroup][(*it)->quality - minQ].addBase(*it, refBase);
		}
	}
}

bool TRecalibrationBQSR::estimateEpsilon(){
	estimationConverged = true;
	for(int i=0; i<numReadGroups; ++i){
		for(int j=0; j<numQ; ++j){
			if(!BQSR_cells_quality[i][j].estimate(convergenceThreshold))
				estimationConverged = false;
		}
	}
	return estimationConverged;
}

void TRecalibrationBQSR::writeToFile(std::string filename){
	std::ofstream out(filename.c_str());

	//write header
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations\n";

	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int j=0; j<numQ; ++j){
			out << it->ID << "\t" << j + minQ << "\tM\t" << BQSR_cells_quality[i][j].curEstimate << "\t" << BQSR_cells_quality[i][j].numObservations << "\n";
		}
	}
	out.close();
}


/*
TRecalibrationBQSR_cellC::TRecalibrationBQSR_cellC(double Error, int pos, TPMD* PmdObject):TRecalibrationBQSR_cell(Error){
	N_2 = 0;
	N_1 = 0;
	D = PmdObject->getProbCT(pos);
}

TRecalibrationBQSR_cellT::TRecalibrationBQSR_cellT(double Error, int pos, TPMD* PmdObject):TRecalibrationBQSR_cellC(Error, pos, PmdObject){
	N_3 = 0;
}

TRecalibrationBQSR_cellA::TRecalibrationBQSR_cellA(double Error, TPMD* PmdObject, double ConvergenceThreshold):TRecalibrationBQSR_cell(Error){
	firstDerivative = 0;
	secondDerivative = 0;
	pmdObject = PmdObject;
	convergenceThreshold = ConvergenceThreshold;
}
TRecalibrationBQSR_cellG::TRecalibrationBQSR_cellG(double Error, TPMD* PmdObject, double ConvergenceThreshold):TRecalibrationBQSR_cellA(Error, PmdObject, ConvergenceThreshold){}


//----------------------------------------

void TRecalibrationBQSR_cellC::addBase(TBase* base, char & RefBase){
	//assumes base is either A, C, G or T
	if(RefBase == 'C') ++N_2;
	else ++N_1;
}

void TRecalibrationBQSR_cellT::addBase(TBase* base, char & RefBase){
	//assumes base is either A, C, G or T
	if(RefBase == 'C') ++N_2;
	else if(RefBase == 'T') ++N_3;
	else ++N_1;
}


double TRecalibrationBQSR_cellA::getD(TBase* base, char & RefBase){
	if(RefBase == A) return 1.0;
	if(RefBase == G) return pmdObject->getProbGA(base->pos3);
	return 0.0;
}
double TRecalibrationBQSR_cellG::getD(TBase* base, char & RefBase){
	if(RefBase == G) return 1.0 - pmdObject->getProbGA(base->pos3);
	return 0.0;
}
void TRecalibrationBQSR_cellA::addBase(TBase* base, char & RefBase){
	double D = getD(base, RefBase);

	double oneMinus4D = 1.0 - 4.0 * D;
	firstDerivative += oneMinus4D / (-4.0*D*curEpsilon + 3.0*D + curEpsilon);
	double tmp = oneMinus4D / (D*(3.0-4.0*curEpsilon) + curEpsilon);
	secondDerivative += tmp * tmp;
}

//----------------------------------------
bool TRecalibrationBQSR_cellC::estimateEpsilon(){
	double B = 4.0*D/3.0 - 1.0;
	curEpsilon = N_1 * (1.0 - D) / (B*(N_1 + N_2));
	estimationConverged = true;
	return estimationConverged;
}

bool TRecalibrationBQSR_cellT::estimateEpsilon(){
	double A = (1.0 - 4*D)/3.0 - 1.0;
	double a = -A*(N_1 + N_2 + N_3);
	double b = A*(N_1 + N_2) - D*(N_1 + N_3);
	double c = D * N_1;

	//there are two solutions to the quadratic function. Calc second if first does not fit in range.
	curEpsilon = -b + sqrt(b*b - 4.0*a*c)/(2*a);
	if(curEpsilon < 0.0 || curEpsilon > 1.0){
		curEpsilon = -b - sqrt(b*b - 4.0*a*c)/(2*a);
	}

	//return that it converged
	estimationConverged = true;
	return estimationConverged;
}

bool TRecalibrationBQSR_cellA::estimateEpsilon(){
	//Need Newton-Ralphson to estimate epsilon
	//Make new estimate
	curEpsilon = curEpsilon + firstDerivative / secondDerivative;

	//decide on convergence
	if(abs(firstDerivative) < convergenceThreshold) estimationConverged = true;
	return estimationConverged;
}

*/
