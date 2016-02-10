/*
 * TPostMortemDamage.cpp
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */


/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#include "TPostMortemDamage.h"

//---------------------------------------------------------------
//TPMDTable
//---------------------------------------------------------------
void TPMDTable::calculateSums(){
	if(!sumsCalculated){
		sums = new long*[maxLength];
		for(int p=0; p<maxLength; ++p){
			sums[p] = new long[4];
			for(int i=0; i<4; ++i){
				sums[p][i] = counts[p][i][0];
				for(int j=1; j<4; ++j){
					sums[p][i] += counts[p][i][j];
				}
			}
		}
		sumsCalculated = true;
	}
};

void TPMDTable::deleteSums(){
	for(int p=0; p<maxLength; ++p){
		delete[] sums[p];
	}
	delete[] sums;
	sumsCalculated = false;
}

TPMDTable::TPMDTable(int MaxLength){
	maxLength = MaxLength;
	sumsCalculated = false;
	sums = NULL;

	//init counts
	counts = new long**[maxLength];
	for(int p=0; p<maxLength; ++p){
		counts[p] = new long*[5];
		for(int i=0; i<5; ++i){
			counts[p][i] = new long[5];
			for(int j=0; j<5; ++j){
				counts[p][i][j] = 0;
			}
		}
	}
};

TPMDTable::~TPMDTable(){
	for(int p=0; p<maxLength; ++p){
		for(int i=0; i<5; ++i){
			delete[] counts[p][i];
		}
		delete[] counts[p];
	}
	delete counts;
	if(sumsCalculated) deleteSums();
};

void TPMDTable::empty(){
	for(int p=0; p<maxLength; ++p){
		for(int i=0; i<5; ++i){
			for(int j=0; j<5; ++j){
				counts[p][i][j] = 0;
			}
		}
	}
};

void TPMDTable::add(int & pos, Base & ref, Base & read){
	if(pos < maxLength)
		++counts[pos][ref][read];
};

void TPMDTable::writeTable(std::ofstream & out, std::string prefix){
	calculateSums();

	//now print out
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			out <<	prefix << genoMap.getBaseAsChar(i) << "->" << genoMap.getBaseAsChar(j);
			for(int p=0; p<maxLength; ++p){
				if(sums[p][i] < 1) out << "\t0.0";
				else out << "\t" << (double) counts[p][i][j] / (double) sums[p][i];
			}
			out << "\n";
		}
	}
};

void TPMDTable::writeTableWithCounts(std::ofstream & out, std::string prefix){
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			out <<	prefix << genoMap.getBaseAsChar(i) << "->" << genoMap.getBaseAsChar(j);
			for(int p=0; p<maxLength; ++p){
				out << "\t" << (double) counts[p][i][j];
			}
			out << "\n";
		}
	}
};

std::string TPMDTable::getPMDStringCT(){
	calculateSums();
	std::string s = "Empiric[";
	double tmpCT, tmpTC;  //tmpRefRead
	for(int p=0; p<maxLength; ++p){
		if(p>0) s += ",";
		if(sums[p][T] < 1 || sums[p][C] < 1) s += "0.0";
		else {
			tmpCT = (double) counts[p][C][T] / (double) sums[p][C];
			tmpTC = (double) counts[p][T][C] / (double) sums[p][T];
			s += toString(std::max(0.0, (tmpCT - tmpTC)/(1.0 - tmpTC)));
		}
	}
	return s + "]";
};

std::string TPMDTable::getPMDStringGA(){
	calculateSums();
	std::string s = "Empiric[";
	double tmpGA, tmpAG;  //tmpRefRead
	for(int p=0; p<maxLength; ++p){
		if(p>0) s += ",";
		if(sums[p][A] < 1 || sums[p][G] < 1) s += "0.0";
		else {
			tmpGA = (double) counts[p][G][A] / (double) sums[p][G];
			tmpAG = (double) counts[p][A][G] / (double) sums[p][A];
			s += toString(std::max(0.0, (tmpGA - tmpAG)/(1.0 - tmpAG)));
		}
	}
	return s + "]";
};

//function to fit an exponential model
void TPMDTable::fillFAndJacobian(arma::vec & F, arma::mat & J, Base & from, Base & to, double* oldParams){
	F.zeros();
	J.zeros();

	double weight, weightJ, tmp;
	double expMinusAlphaP;
	double dExpMinusAlphaP;
	for(int p=0; p<maxLength; ++p){
		//exp
		expMinusAlphaP = exp(-oldParams[2] * p);
		dExpMinusAlphaP = oldParams[1] * expMinusAlphaP;

		//first term
		//----------
		tmp = oldParams[0] + dExpMinusAlphaP;
		weight = counts[p][from][to] / tmp;
		weightJ = weight / tmp;

		//add to F
		F(0) += weight;
		F(1) += weight * expMinusAlphaP;
		F(2) -= weight * p * dExpMinusAlphaP;

		//add to J -> only upper triangle, as it is symmetric
		J(0,0) -= weightJ;
		J(0,1) -= weightJ * expMinusAlphaP;
		J(0,2) += weightJ * p * dExpMinusAlphaP;
		J(1,1) -= weightJ * expMinusAlphaP * expMinusAlphaP;
		J(1,2) -= weightJ * p * oldParams[0] * expMinusAlphaP;
		J(2,2) += weightJ * p * p * oldParams[0] * dExpMinusAlphaP;

		//second term
		//-----------
		tmp = (1.0 - oldParams[0] - dExpMinusAlphaP);
		weight = (sums[p][from] - counts[p][from][to]) / tmp;
		weightJ = weight / tmp;

		//add to F
		F(0) -= weight;
		F(1) -= weight * expMinusAlphaP;
		F(2) += weight * p * dExpMinusAlphaP;

		//add to J -> only upper triangle, as it is symmetric
		J(0,0) -= weightJ;
		J(0,1) -= weightJ * expMinusAlphaP;
		J(0,2) += weightJ * p * dExpMinusAlphaP;
		J(1,1) -= weightJ * expMinusAlphaP * expMinusAlphaP;
		J(1,2) += weightJ * p * (1.0 - oldParams[0]) * expMinusAlphaP;
		J(2,2) -= weightJ * p * p * (1.0 - oldParams[0]) * dExpMinusAlphaP;
	}

	//now fill in lower triangle of J
	J(1,0) = J(0,1);
	J(2,0) = J(0,2);
	J(2,1) = J(1,2);
}

void TPMDTable::fillF(arma::vec & F, Base & from, Base & to, double* oldParams){
	F.zeros();

	double weight;
	double expMinusAlphaP;
	double dExpMinusAlphaP;
	for(int p=0; p<maxLength; ++p){
		//exp
		expMinusAlphaP = exp(-oldParams[2] * (double) p);
		dExpMinusAlphaP = oldParams[1] * expMinusAlphaP;

		//first term
		//----------
		weight = counts[p][from][to] / (oldParams[0] + dExpMinusAlphaP);
		weight -= (sums[p][from] - counts[p][from][to]) / (1.0 - oldParams[0] - dExpMinusAlphaP);

		//add to F
		F(0) += weight;
		F(1) += weight * expMinusAlphaP;
		F(2) -= weight * (double) p * dExpMinusAlphaP;
	}
}

double TPMDTable::calcLL(Base & from, Base & to, double* oldParams){
	double LL = 0.0;
	double dExpMinusAlphaP;
	for(int p=0; p<maxLength; ++p){
		dExpMinusAlphaP = oldParams[1] * exp(-oldParams[2] * p);
		LL += counts[p][from][to] * log(oldParams[0] + dExpMinusAlphaP)
			+ (sums[p][from] - counts[p][from][to]) * log(1.0 - oldParams[0] - dExpMinusAlphaP);
	}
	return LL;
}

std::string TPMDTable::fitExponentialModel(Base from, Base to, int & numNRIterations, double & eps){
	//variables
	arma::mat J(3,3);
	arma::vec F(3);
	arma::vec oldF(3);
	arma::mat JxF;

	//set starting values
	double oldParams[3]; //exponential parameters alpha + beta exp(-gamma * position);
	double newParams[3];
	double LL, oldLL;

	//find last entry with counts
	int lastPositionToConsiderPlusOne = -1;
	for(int p=maxLength-1; p >= 0; ++p){
		if(sums[p][from] > 0){
			lastPositionToConsiderPlusOne = p;
			break;
		}
	}

	//use OLS to find starting point
	//------------------------------
	//some variables
	double gammaStep = 0.01;
	double gammaTmp = -gammaStep + 0.00000001;
	double SSRold = lastPositionToConsiderPlusOne;
	double SSRnew;
	double SSRdiff = -1.0;

	arma::mat X(lastPositionToConsiderPlusOne, 2);
	X.ones();
	arma::mat XtX;
	arma::mat betaHat;

	//fill vector y to fit using OLS
	arma::vec y(lastPositionToConsiderPlusOne);
	double sumYSquared = 0.0;
	for(int p=0; p<lastPositionToConsiderPlusOne; ++p){
		y(p) = (double) counts[p][from][to] / (double) sums[p][from];
		sumYSquared += y(p) * y(p);
	}

	//do until we get a small alpha
	while(fabs(gammaStep) > 0.00000001){
		while(SSRdiff < 0.0){
			//update gamma
			gammaTmp += gammaStep;

			//fill x
			for(int p=0; p<lastPositionToConsiderPlusOne; ++p){
				X(p,1) = exp(-gammaTmp * (double) p);
			}

			betaHat = inv(X.t()*X) * X.t() * y;

			//calc sum of squares
			arma::mat tmp = (betaHat.t() * X.t() * y);
			SSRnew = sumYSquared - tmp(0,0);
			SSRdiff = SSRnew - SSRold;
			SSRold = SSRnew;
		}
		//update alpha step
		gammaStep = -0.1 * gammaStep;
		SSRdiff = -1.0;
	}

	//store params
	oldParams[0] = betaHat(0);
	oldParams[1] = betaHat(1);
	oldParams[2] = gammaTmp;

	//Conduct Newton-Raphson to refine
	//----------------------------------
	oldLL = calcLL(from, to, oldParams);

	for(int i = 0; i<numNRIterations; ++i){
		fillFAndJacobian(F, J, from, to, oldParams);
		if(solve(JxF, J, F)){
			//estimate new params
			for(int x=0; x<3; ++x)
				newParams[x] = oldParams[x] - JxF(x);

				//calculate LL at new location
				LL = calcLL(from, to, newParams);

				//check if we accept or backtrack
				if(LL > oldLL){
					//store new params
					for(int x=0; x<3; ++x)
						oldParams[x] = newParams[x];

					//check if we stop NR
					if(LL - oldLL < eps){
						oldLL = LL;
						break;
					}
					oldLL = LL;
				}
		} else {
			std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << J << std::endl << std::endl;
			throw "Issue solving JxF in TPMDTable::fitExponentialModel!";
		}
	}

	//writing to file -> first need to "decode" the parameters
	//the exponential PMD model is f(C->T) = mu + (1-mu) *[ a*exp(-b * position) + c ]
	//but we fitted f(C->T) = alpha + beta * exp(-gamma * position).
	//Hence we have:
	// mu is estimated from T->C transitions
	// a =  beta / (1 - mu)
	// b = gamma
	// c = (alpha - mu) / (1 - mu)
	double mu = 0.0; double sum = 0.0;
	for(int p=0; p<lastPositionToConsiderPlusOne; ++p){
		mu += counts[p][to][from];
		sum += sums[p][to];
	}
	mu = mu / sum;

	double a = oldParams[1] / (1.0 - mu);
	double b = oldParams[2];
	double c = (oldParams[0] - mu) / (1.0 - mu);

	//calculate variance of estimator
	arma::mat Fisher = inv(J);
	for(int x=0; x<3; ++x)
		std::cout << "Variance " << x << " = " << sqrt(-Fisher.diag()) << std::endl;



	//no return string
	return "Exponential[" + toString(a) + "," + toString(b) + "," + toString(c) + "]";
}

//---------------------------------------------------------------
//TPMDTables
//---------------------------------------------------------------
TPMDTables::TPMDTables(TReadGroups* ReadGroups, int maxLength){
	readGroups = ReadGroups;
	forward = new TPMDTable*[readGroups->numGroups];
	reverse = new TPMDTable*[readGroups->numGroups];
	for(int i=0; i<readGroups->numGroups; ++i){
		forward[i] = new TPMDTable(maxLength);
		reverse[i] = new TPMDTable(maxLength);
	}
};

TPMDTables::~TPMDTables(){
	for(int i=0; i<readGroups->numGroups; ++i){
		delete forward[i];
		delete reverse[i];
	}
	delete[] forward;
	delete[] reverse;
};

void TPMDTables::addForward(int readGroup, int pos, Base & ref, Base & read){
	forward[readGroup]->add(pos, ref, read);
};

void TPMDTables::addReverse(int readGroup, int pos, Base & ref, Base & read){
	reverse[readGroup]->add(pos, ref, read);
};

void TPMDTables::writePMDFile(std::string filename){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups
	for(int i=0; i<readGroups->numGroups; ++i){
		out << readGroups->getName(i) << "\t" << forward[i]->getPMDStringCT() << "\t" << reverse[i]->getPMDStringGA() << "\n";
	}
	out.close();
}

void TPMDTables::writeTable(std::string filename){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups
	for(int i=0; i<readGroups->numGroups; ++i){
		forward[i]->writeTable(out, readGroups->getName(i) + "\tforward\t");
		reverse[i]->writeTable(out, readGroups->getName(i) + "\treverse\t");
	}
	out.close();
};

void TPMDTables::writeTableWithCounts(std::string filename){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups
	for(int i=0; i<readGroups->numGroups; ++i){
		forward[i]->writeTableWithCounts(out, readGroups->getName(i) + "\tforward\t");
		reverse[i]->writeTableWithCounts(out, readGroups->getName(i) + "\treverse\t");
	}
	out.close();
};

void TPMDTables::fitExponentialModel(int numNRIterations, double eps, std::string & filename){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups, fit and write exponential model
	for(int i=0; i<readGroups->numGroups; ++i){
		out << readGroups->getName(i) << "\t" << forward[i]->fitExponentialModel(C, T, numNRIterations, eps) << "\t" << reverse[i]->fitExponentialModel(G, A, numNRIterations, eps) << "\n";
	}
	out.close();
}
//---------------------------------------------------------------
//TPMDFunction
//---------------------------------------------------------------
TPMDSkoglund::TPMDSkoglund(double & Lambda, double & C){
	lambda = Lambda; c = C;
};
double TPMDSkoglund::getProb(int & pos){
	//Note: distance is zero based!
	return c + pow(1.0 - lambda, (double) pos) * lambda;
};
std::string TPMDSkoglund::getString(){
	return "P(pmd|pos) = p * (1 - p)^pos + c = " + toString(lambda) + " * (1 - " + toString(lambda) + ")^pos + " + toString(c);
};

//---------------------------------------------------------------
TPMDExponential::TPMDExponential(double & A, double & B, double & C){
	a = A; b = B; c = C;
}
double TPMDExponential::getProb(int & pos){
	//Note: distance is zero based!
	return a * exp(- (double) pos * b) + c;
};
std::string TPMDExponential::getString(){
	return "P(pmd|pos) = a * exp(- pos * b) + c = " + toString(a) + "* exp(- pos * " + toString(b) + ") + " + toString(c);
};

//---------------------------------------------------------------
TPMDEmpiric::TPMDEmpiric(std::string & values, std::string & example){
	std::vector<double> vec;
	fillVectorFromString(values, vec, ',');
	length = vec.size();
	if(length < 1) throw "Can not initialize post mortem damage function 'empiric[" + values + "]': wrong format!\n" + example;
	probs = new double[length];
	for(int i=0; i<length; ++i){
		probs[i] = vec[i];
	}
	last = probs[length-1];
};

double TPMDEmpiric::getProb(int & pos){
	if(pos < length) return probs[pos];
	else return last;
};

std::string TPMDEmpiric::getString(){
	std::string s = "P(pmd|pos) = {" + toString(probs[0]);
	for(int i=1; i<length; ++i) s += "," + toString(probs[i]);
	return  s +"}";
};

//------------------------------------------------------
//TPMD
//------------------------------------------------------
void TPMD::initializeFunction(std::string & pmdString, PMDType type){
	//parse string to get model.  options are
	// none
	// Skoglund[lambda,c]
	// Veeramah[a,b,c]
	std::string example = "Use either Skoglund[p,c], Exponential[a,b,c] or Empiric[0.2,0.3,...]";

	//check if function was initialized abefore
	if(functionsInitialized[type]) throw "PMD function has been initialized previously!";

	//check if it is none
	if(pmdString == "none"){
		myFunctions[type] = new TPMDFunction();
	} else {
		std::string::size_type pos = pmdString.find_first_of('[');
		if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
		std::string name = pmdString.substr(0,pos);

		//switch between functions
		if(name == "Empiric"){
			std::string::size_type endPos = pmdString.find_first_of(']');
			if(endPos == std::string::npos || endPos != pmdString.length()-1) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
			std::string list = pmdString.substr(pos+1, endPos-pos-1);
			myFunctions[type] = new TPMDEmpiric(list, example);
		} else {
			//prepare first value
			std::string tmp = pmdString.substr(pos+1, pmdString.length() - pos - 1);
			pos = tmp.find_first_of(',');
			if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
			double first = atof(tmp.substr(0, pos).c_str());

			if(name == "Skoglund"){
				//get lambda and
				if(first < 0.0) throw "Can not initialize Skoglund function with lambda < 0!";
				if(first > 1.0) throw "Can not initialize Skoglund function with lambda > 1!";
				double c = atof(tmp.substr(pos+1).c_str());
				if(c < 0.0) throw "Can not initialize Skoglund function with c < 0!";
				myFunctions[type] = new TPMDSkoglund(first, c);
			} else if(name == "Exponential"){
				//get a, b and c
				if(first < 0.0) throw "Can not initialize Exponential function with a < 0!";

				//get b
				tmp = tmp.substr(pos+1);
				pos = tmp.find_first_of(',');
				if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
				double b = atof(tmp.substr(0, pos).c_str());
				if(b < 0.0) throw "Can not initialize Exponential function with b < 0!";

				//get c
				double c = atof(tmp.substr(pos+1).c_str());
				if(c < 0.0) throw "Can not initialize Exponential function with c < 0!";

				//test if it can be > 1
				//if(first + c > 1) throw "Can not initialize Exponential function with a + c > 1!";

				//initialze
				myFunctions[type] = new TPMDExponential(first, b, c);
			} else throw "Can not initialize post mortem damage function '" + pmdString + "': wrong name!\n" + example;
		}
	}
	functionsInitialized[type] = true;
};



