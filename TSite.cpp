/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"

//---------------------------------------------------------------
//GenotypeMap
//---------------------------------------------------------------
GenotypeMap::GenotypeMap(){
//set up genotype map
	genotypeMap = new Genotype*[4];
	for(int i=0; i<4; ++i)
		genotypeMap[i] = new Genotype[4];

	int geno = 0;
	for(int i=0; i<4; ++i){
		for(int j=i; j<4; ++j){
			genotypeMap[i][j] = static_cast<Genotype>(geno);
			genotypeMap[j][i] = genotypeMap[i][j];
			++geno;
		}
	}
}

std::string GenotypeMap::getGenotypeString(int num){
	if(num==0) return "AA";
	if(num==1) return "AC";
	if(num==2) return "AG";
	if(num==3) return "AT";
	if(num==4) return "CC";
	if(num==5) return "CG";
	if(num==6) return "CT";
	if(num==7) return "GG";
	if(num==8) return "GT";
	if(num==9) return "TT";
	return "ERROR!";
};

//-------------------------------------------------------
//TBase
//-------------------------------------------------------
void TBase::fillEmissionProbabilities(TPMD & pmdObject){
	fillEmissionProbabilitiesCore(pmdObject, errorRate);
}

void TBase::fillEmissionProbabilitiesRecalibratedError(TPMD & pmdObject, TRecalibration & recal){
	fillEmissionProbabilitiesCore(pmdObject, recal.recalibrate(errorRate));
}

void TBaseDiploidA::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(AA, oneMinusError);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pmd) * errorOneThird + pmd * oneMinusError);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, 0.5 - errorOneThird);
	emissionProbabilities.set(AG, ((1.0 + pmd) * oneMinusError + (1.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, (pmd * oneMinusError + (2.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(CG));
};

void TBaseHaploidA::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(A, oneMinusError);
	emissionProbabilities.set(C, errorOneThird);
	emissionProbabilities.set(G, pmd * oneMinusError + (1.0 - pmd) * errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidC::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(AA , errorOneThird);
	emissionProbabilities.set(CC , (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(GG , errorOneThird);
	emissionProbabilities.set(TT , errorOneThird);
	emissionProbabilities.set(AC , ((1.0 - pmd) * oneMinusError + (1.0 + pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG , errorOneThird);
	emissionProbabilities.set(AT , errorOneThird);
	emissionProbabilities.set(CG , emissionProbabilities.get(AC));
	emissionProbabilities.set(CT , emissionProbabilities.get(AC));
	emissionProbabilities.set(GT , errorOneThird);
};

void TBaseHaploidC::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(G, errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidG::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, errorOneThird);
	emissionProbabilities.set(AG, ((1.0 - pmd) * oneMinusError + (1.0 + pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AG));
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(AG));

};

void TBaseHaploidG::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, errorOneThird);
	emissionProbabilities.set(G, (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidT::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, (1.0 - pmd) * errorOneThird + pmd * oneMinusError);
	emissionProbabilities.set(GG, errorOneThird);
	emissionProbabilities.set(TT, oneMinusError);
	emissionProbabilities.set(AC, (pmd * oneMinusError + (2.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG, errorOneThird);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AC));
	emissionProbabilities.set(CT, ((1.0 + pmd) * oneMinusError + (1.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(GT, 0.5 - errorOneThird);
};

void TBaseHaploidT::fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, pmd * oneMinusError + (1.0 - pmd) * errorOneThird);
	emissionProbabilities.set(G, errorOneThird);
	emissionProbabilities.set(T, oneMinusError);
}
//-------------------------------------------------------
//TPMD
//-------------------------------------------------------
TPMD::TPMD(){
	myFunctions[pmdCT] = NULL;
	myFunctions[pmdGA] = NULL;
	functionsInitialized[pmdCT] = false;
	functionsInitialized[pmdGA] = false;
};

void TPMD::initializeFunction(std::string & pmdString, PMDType type){
	//parse string to get model.  options are
	// none
	// Skoglund[lambda,c]
	// Veeramah[a,b,c]
	std::string example = "Use either Skoglund[p,c] or Veeramah[a,b,c]";

	//check if it is none
	if(pmdString == "none"){
		myFunctions[type] = new TPMDFunction();
	} else {
		std::string::size_type pos = pmdString.find_first_of('[');
		if(pos == std::string::npos) throw "Can not initialize post mortem damage function: wrong format!\n" + example;
		std::string name = pmdString.substr(0,pos);

		//prepare first value
		std::string tmp = pmdString.substr(pos+1, pmdString.length() - pos - 1);
		pos = tmp.find_first_of(',');
		if(pos == std::string::npos) throw "Can not initialize post mortem damage function: wrong format!\n" + example;
		double first = atof(tmp.substr(0, pos).c_str());

		//switch between functions
		if(name == "Skoglund"){
			//get lambda and
			if(first < 0.0) throw "Can not initialize Skoglund function with lambda < 0!";
			if(first > 1.0) throw "Can not initialize Skoglund function with lambda > 1!";
			double c = atof(tmp.substr(pos+1).c_str());
			if(c < 0.0) throw "Can not initialize Skoglund function with c < 0!";
			myFunctions[type] = new TPMDSkoglund(first, c);
		} else if(name == "Veeramah"){
			//get a, b and c
			if(first < 0.0) throw "Can not initialize Veeramah function with a < 0!";

			//get b
			tmp = tmp.substr(pos+1);
			pos = tmp.find_first_of(',');
			if(pos == std::string::npos) throw "Can not initialize post mortem damage function: wrong format!\n" + example;
			double b = atof(tmp.substr(0, pos).c_str());
			if(b < 0.0) throw "Can not initialize Veeramah function with b < 0!";

			//get c
			double c = atof(tmp.substr(pos+1).c_str());
			if(c < 0.0) throw "Can not initialize Veeramah function with c < 0!";

			//test if it can be > 1
			if(first + c > 1) throw "Can not initialize Veeramah function with a + c > 1!";

			//initialze
			myFunctions[type] = new TPMDVeeramah(first, b, c);
		} else throw "Can not initialize post mortem damage function: wrong name!\n" + example;
	}
	functionsInitialized[type] = true;
}

//---------------------------------------------------------------
//Recalibration
//---------------------------------------------------------------
TRecalibration::TRecalibration(std::string recalString){
	if(recalString==""){
		doRecalibration = false;
		a = 0.0;
		b = 0.0;
		c = 1.0;
	} else {
		doRecalibration= true;
		std::string example = "Use '[a,b,c]'";
		std::string::size_type pos = recalString.find_first_of('[');
		if(pos == std::string::npos) throw "Can not initialize recalibration: wrong format! " + example;
		recalString = recalString.substr(pos+1, recalString.length() - pos - 2);
		pos = recalString.find_first_of(',');
		if(pos == std::string::npos) throw "Can not initialize recalibration: wrong format!\n" + example;
		a = stringToDoubleCheck(recalString.substr(0, pos));
		recalString = recalString.substr(pos+1, recalString.length() - pos - 1);
		pos = recalString.find_first_of(',');
		if(pos == std::string::npos) throw "Can not initialize recalibration: wrong format!\n" + example;
		b = stringToDoubleCheck(recalString.substr(0, pos));
		c = stringToDoubleCheck(recalString.substr(pos+1));
	}
}

double TRecalibration::recalibrate(double & error){
	if(!doRecalibration) return error;
	double tmp = log10(error);
	tmp = a * (1.0 - exp(-b * tmp)) + c * tmp;
	if(tmp > -0.1) return 0.8;
	else return pow10(tmp);
}

std::string TRecalibration::getFunctionString(){
	return "log10(error recalibrated) = " + toString(b) + " * (1 - exp(-" + toString(a) + " * log10(error))) + " + toString(c) + " * log10(error)";
}

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
	for(int i=0; i<numParams; ++i){
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
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(){

}



//-------------------------------------------------------
//TSite
//-------------------------------------------------------
void TSite::clear(){
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it)
		delete *it;
	bases.clear();
	hasData = false;
};

double TSite::qualityToLogError(char & quality){
	return -((double)quality - 33.0) / 10.0;
}

void TSiteDiploid::add(char & base, char & quality, int PosInRead, int pos5, int pos3, int & ReadGroup){
	if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
		double logError = qualityToLogError(quality);
		if(logError < 0.0){
			if(base == 'A') bases.push_back(new TBaseDiploidA(logError, PosInRead, pos5, pos3, ReadGroup));
			else if(base == 'C') bases.push_back(new TBaseDiploidC(logError, PosInRead, pos5, pos3, ReadGroup));
			else if(base == 'G') bases.push_back(new TBaseDiploidG(logError, PosInRead, pos5, pos3, ReadGroup));
			else bases.push_back(new TBaseDiploidT(logError, PosInRead, pos5, pos3, ReadGroup));
		}
		hasData = true;
	}
};
void TSiteHaploid::add(char & base, char & quality, int PosInRead, int pos5, int pos3, int & ReadGroup){
	if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
		double logError = qualityToLogError(quality);
		if(logError < 0.0){
			if(base == 'A') bases.push_back(new TBaseHaploidA(logError, PosInRead, pos5, pos3, ReadGroup));
			else if(base == 'C') bases.push_back(new TBaseHaploidC(logError, PosInRead, pos5, pos3, ReadGroup));
			else if(base == 'G') bases.push_back(new TBaseHaploidG(logError, PosInRead, pos5, pos3, ReadGroup));
			else bases.push_back(new TBaseHaploidT(logError, PosInRead, pos5, pos3, ReadGroup));
		}
		hasData = true;
	}
};

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	double weight = 1.0 / bases.size();
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		(*it)->addToBaseFrequencies(frequencies, weight);
	}
}

void TSite::calcEmissionProbabilities(TPMD & pmdObject){
	for(int i=0; i<numGenotypes; ++i){
		emissionProbabilities[i] = 1.0;
	}
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		(*it)->fillEmissionProbabilities(pmdObject);
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] *= (*it)->getEmissionProbability(i);
		}
	}
}

void TSite::calcEmissionProbabilitiesScaledError(TPMD & pmdObject, TRecalibration & recal){
	for(int i=0; i<numGenotypes; ++i){
		emissionProbabilities[i] = 1.0;
	}
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		(*it)->fillEmissionProbabilitiesRecalibratedError(pmdObject, recal);
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] *= (*it)->getEmissionProbability(i);
		}
	}
}

std::string TSite::getBases(){
	if(bases.size()==0) return "-";
	std::string b = "";
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		b += (*it)->getBase();
	}
	return b;
}

std::string TSite::getEmissionProbs(){
	std::string b = toString(emissionProbabilities[0]);
	for(int i=1; i<numGenotypes; ++i){
		b += "\t" + toString(emissionProbabilities[i]);
	}
	return b;
}

void TSite::calculateP_g(double* genotypeProbabilities){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<numGenotypes; ++i){
		P_g[i] =  emissionProbabilities[i] * genotypeProbabilities[i];
		sum += P_g[i];
	}
	for(int i=0; i<10; ++i){
		P_g[i] /= sum;
	}
}

double TSite::calculateWeightedSumOfEmissionProbs(double* weights){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<numGenotypes; ++i){
		sum += emissionProbabilities[i] * weights[i];
	}
	return sum;
}

double TSite::calculateLogLikelihood(double* genotypeProbabilities){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<numGenotypes; ++i){
		sum +=  emissionProbabilities[i] * genotypeProbabilities[i];
	}
	return log(sum);
}

