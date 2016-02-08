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
void TPMDTable::fillFAndJacobian(arma::vec & F, arma::mat & J, double* oldParams){
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
		weight = counts[p][C][T] / tmp;
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
		weight = (sums[p][C] - counts[p][C][T]) / tmp;
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

void TPMDTable::fillF(arma::vec & F, double* oldParams){
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
		weight = counts[p][C][T] / (oldParams[0] + dExpMinusAlphaP);
		weight -= (sums[p][C] - counts[p][C][T]) / (1.0 - oldParams[0] - dExpMinusAlphaP);

		//add to F
		F(0) += weight;
		F(1) += weight * expMinusAlphaP;
		F(2) -= weight * (double) p * dExpMinusAlphaP;
	}
}

double TPMDTable::calcLL(double* oldParams){
	double LL = 0.0;
	double dExpMinusAlphaP;
	for(int p=0; p<maxLength; ++p){
		dExpMinusAlphaP = oldParams[1] * exp(-oldParams[2] * p);
		LL += counts[p][C][T] * log(oldParams[0] + dExpMinusAlphaP)
			+ (sums[p][C] - counts[p][C][T]) * log(1.0 - oldParams[0] - dExpMinusAlphaP);
	}
	return LL;
}

void TPMDTable::fitExponentialModel(int & numNRIterations, double & eps, TLog* logfile){
	//variables
	arma::mat J(3,3);
	arma::vec F(3);
	arma::vec oldF(3);
	arma::mat JxF;

	//set starting values
	double oldParams[3]; //mu, delta, alpha
	double newParams[3];
	oldParams[0] = 0.1;
	oldParams[1] = 0.3;
	oldParams[2] = 0.01;
	double LL, oldLL;
	oldLL = calcLL(oldParams);
	logfile->write("Initial LL = " + toString(oldLL));

	std::cout << "INITIAL Params:";
	for(int x=0; x<3; ++x)
		std::cout << "\t" << oldParams[x];
	std::cout << std::endl;

	//find last entry with counts
	int lastPositionToConsiderPlusOne = -1;
	for(int p=0; p<maxLength; ++p){
		if(sums[p][C] <= 0){
			lastPositionToConsiderPlusOne = p;
			break;
		}
	}

/*
	//use OLS to find starting point
	//------------------------------
	//some variables
	double alphaStep = 0.01;
	double alphaTmp = -alphaStep;
	int SSDold = lastPositionToConsiderPlusOne;
	int SSRdiff = -1.0;
	arma::vec vecOFOnes(lastPositionToConsiderPlusOne);
	vecOFOnes.ones(lastPositionToConsiderPlusOne);
	arma::mat X(lastPositionToConsiderPlusOne, 2);
	X.ones();

	//fill vector f to fit using OLS
	double* f = new double[lastPositionToConsiderPlusOne];
	for(int p=0; p<lastPositionToConsiderPlusOne; ++p){
		f[p] = counts[p][C][T] / sums[p][C];
	}

	//do until we get a small alpha
	while(fabs(alphaStep) > 0.00000001){
		while(SSRdiff < 0.0){
			//update alpha
			alphaTmp += alphaStep;

			//fill x
			for(int p=0; p<lastPositionToConsiderPlusOne; ++p){
				X(p,1) = exp(-alphaTmp * p);
			}

			//run OLS



		}
	}


	y <- f<-as.matrix(f, ncol=1);
	ones <-matrix(1,nrow=P,ncol=1);

	SSR.old <- P;
	SSR.diff <- -1;
	delta.alpha <- 0.01;

	while (abs(delta.alpha) > 10^{-6}){
	  while (SSR.diff < 0){

	    alpha <- alpha + delta.alpha;
	    x <- matrix(exp(-alpha*(1:P)),ncol=1);
	    X <- cbind(ones,x);

	    beta.hat <- solve(t(X)%*%X)%*%t(X)%*%y;
	    SSR.new <- sum(y^2) - t(beta.hat)%*%t(X)%*%y;
	    SSR.diff <- SSR.new - SSR.old;
	    SSR.old <- SSR.new;
	  }
	delta.alpha <- -0.1*delta.alpha;
	SSR.diff <- -1;
	}
	*/

	//do grid search to find starting values
	logfile->listFlush("Running grid search ...");
	int gridSize = 11;
	for(int i=0; i<gridSize; ++i){
		newParams[0] = i * 1.0/(gridSize - 1);
		for(int j=0; j<gridSize; ++j){
			newParams[1] = j * 1.0/(gridSize - 1);
			for(int k=0; k<gridSize; ++k){
				newParams[2] = k * 1.0/(gridSize - 1);
				//calc LL and decide if it is better than old
				LL = calcLL(newParams);
				if(LL > oldLL){
					for(int x=0; x<3; ++x)
						oldParams[x] = newParams[x];
					oldLL = LL;
				}
			}
		}
	}

	logfile->write(" done!");
	logfile->conclude("LL increased to " + toString(LL));

	std::cout << "GRID Params:";
	for(int x=0; x<3; ++x)
		std::cout << "\t" << oldParams[x];
	std::cout << std::endl;

	//then conduct a few steepst ascent runs
	logfile->listFlush("Running steepest ascent runs ...");
	double normF2;
	for(int i = 0; i<100000; ++i){
		fillF(F, oldParams);
		normF2 = 0;
		for(int x=0; x<3; ++x) normF2 += F(x) * F(x);
		normF2 = sqrt(normF2);
		F = F / normF2;
		for(int x=0; x<3; ++x)
			oldParams[x] = oldParams[x] +  0.001 * F(x);

		//test if we break
		if(i>100 && dot(F, oldF) < 0.0) break;
		oldF = F;
	}

	LL = calcLL(oldParams);
	logfile->write(" done!");
	logfile->conclude("LL increased to " + toString(LL));

	std::cout << "ASCENT Params:";
	for(int x=0; x<3; ++x)
		std::cout << "\t" << oldParams[x];
	std::cout << std::endl;


	//and finally conduct Newton-Raphson
	logfile->listFlush("refining using Newton-Raphson ...");
	double lambda; //used in backtracking
	bool acceptMove;
	bool NRconverged = false;

	for(int i = 0; i<numNRIterations; ++i){
		fillFAndJacobian(F, J, oldParams);
		if(solve(JxF, J, F)){
			//estimate new params
			for(int x=0; x<3; ++x)
				newParams[x] = oldParams[x] - lambda * JxF(x);

				//calculate LL at new location
				LL = calcLL(newParams);

				//check if we accept or backtrack
				if(LL > oldLL){
					//store new params
					for(int x=0; x<3; ++x)
						oldParams[x] = newParams[x];

					//check if we stop NR
					if(LL - oldLL < eps){
						oldLL = LL;
						logfile->conclude("Stopping Newton-Raphson: increase in LL was < " + toString(eps));
						break;
					}

					oldLL = LL;
				}
		} else {
			std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << J << std::endl << std::endl;
			throw "Issue solving JxF in TPMDTable::fitExponentialModel!";
		}

		std::cout << "TMP Params:";
		for(int x=0; x<3; ++x)
			std::cout << "\t" << oldParams[x];
		std::cout << std::endl;
	}

	//print params
	std::cout << "Params:";
	for(int x=0; x<3; ++x)
		std::cout << "\t" << oldParams[x];
	std::cout << std::endl;
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

void TPMDTables::fitExponentialModel(int numNRIterations, double eps, TLog* logfile){
	//loop over all read groups
	for(int i=0; i<readGroups->numGroups; ++i){
		forward[i]->fitExponentialModel(numNRIterations, eps, logfile);
		//reverse[i]->fitExponentialModel(numNRIterations, eps, logfile);
	}
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
TPMDVeeramah::TPMDVeeramah(double & A, double & B, double & C){
	a = A; b = B; c = C;
}
double TPMDVeeramah::getProb(int & pos){
	//Note: distance is zero based!
	return a * exp(- (double) pos * b) + c;
};
std::string TPMDVeeramah::getString(){
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
	std::string example = "Use either Skoglund[p,c], Veeramah[a,b,c] or Empiric[0.2,0.3,...]";

	//check if function was initialized abefore
	if(functionsInitialized[type]) throw "PMD function has been initialized previously!";

	std::cout << "PMD STRING = '" << pmdString << "'" << std::endl;

	//check if it is none
	if(pmdString == "none"){
		myFunctions[type] = new TPMDFunction();
	} else {
		std::string::size_type pos = pmdString.find_first_of('[');
		if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
		std::string name = pmdString.substr(0,pos);



		//switch between functions
		if(name == "Empiric"){

			std::cout << "EMPIRIC!!!!!" << std::endl;

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
			} else if(name == "Veeramah"){
				//get a, b and c
				if(first < 0.0) throw "Can not initialize Veeramah function with a < 0!";

				//get b
				tmp = tmp.substr(pos+1);
				pos = tmp.find_first_of(',');
				if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
				double b = atof(tmp.substr(0, pos).c_str());
				if(b < 0.0) throw "Can not initialize Veeramah function with b < 0!";

				//get c
				double c = atof(tmp.substr(pos+1).c_str());
				if(c < 0.0) throw "Can not initialize Veeramah function with c < 0!";

				//test if it can be > 1
				//if(first + c > 1) throw "Can not initialize Veeramah function with a + c > 1!";

				//initialze
				myFunctions[type] = new TPMDVeeramah(first, b, c);
			} else throw "Can not initialize post mortem damage function '" + pmdString + "': wrong name!\n" + example;
		}
	}
	functionsInitialized[type] = true;
};



