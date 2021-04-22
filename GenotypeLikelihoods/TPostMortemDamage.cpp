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

namespace GenotypeLikelihoods{



std::string TPMDTable::getPMDString(Base first, Base second){
	calculateSums();
	std::string s = "Empiric[";
	double tmpFirstToSecond, tmpSecondToFirst;  //tmpRefRead
	for(int p=0; p<_maxLength; ++p){
		if(p>0) s += ",";
		if(_sums[p][second] < 1 || _sums[p][first] < 1) s += "0.0";
		else {
			tmpFirstToSecond = (double) _counts[p][first][second] / (double) _sums[p][first];
			tmpSecondToFirst = (double) _counts[p][second][first] / (double) _sums[p][second];
			s += toString(std::max(0.0, (tmpFirstToSecond - tmpSecondToFirst)/(1.0 - tmpSecondToFirst)));
		}
	}
	return s + "]";
};
/*
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
*/
//function to fit an exponential model
void TPMDTable::fillFAndJacobian(arma::vec & F, arma::mat & J, Base & from, Base & to, double* oldParams){
	F.zeros();
	J.zeros();

	double weight, weightJ, tmp;
	double expMinusAlphaP;
	double dExpMinusAlphaP;
	for(int p=0; p<_maxLength; ++p){
		//exp
		expMinusAlphaP = exp(-oldParams[2] * p);
		dExpMinusAlphaP = oldParams[1] * expMinusAlphaP;

		//first term
		//----------
		tmp = oldParams[0] + dExpMinusAlphaP;
		weight = _counts[p][from][to] / tmp;
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
		weight = (_sums[p][from] - _counts[p][from][to]) / tmp;
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
	for(int p=0; p<_maxLength; ++p){
		//exp
		expMinusAlphaP = exp(-oldParams[2] * (double) p);
		dExpMinusAlphaP = oldParams[1] * expMinusAlphaP;

		//first term
		//----------
		weight = _counts[p][from][to] / (oldParams[0] + dExpMinusAlphaP);
		weight -= (_sums[p][from] - _counts[p][from][to]) / (1.0 - oldParams[0] - dExpMinusAlphaP);

		//add to F
		F(0) += weight;
		F(1) += weight * expMinusAlphaP;
		F(2) -= weight * (double) p * dExpMinusAlphaP;
	}
}

double TPMDTable::calcLL(Base & from, Base & to, double* oldParams){
	double LL = 0.0;
	double dExpMinusAlphaP;
	for(int p=0; p<_maxLength; ++p){
		dExpMinusAlphaP = oldParams[1] * exp(-oldParams[2] * p);
		LL += _counts[p][from][to] * log(oldParams[0] + dExpMinusAlphaP)
			+ (_sums[p][from] - _counts[p][from][to]) * log(1.0 - oldParams[0] - dExpMinusAlphaP);
	}
	return LL;
}

std::string TPMDTable::fitExponentialModel(Base from, Base to, int & numNRIterations, double & eps, std::string readGroupName, int maxReadLength, TLog* logfile){
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
	for(int p=_maxLength-1; p >= 0; --p){
		if(_sums[p][from] > 100){
			lastPositionToConsiderPlusOne = p;
			break;
		}
	}
	++lastPositionToConsiderPlusOne;

	if(lastPositionToConsiderPlusOne < 10) throw "Not sufficient data for PMD estimation in read group '" + readGroupName + "': less than the ten first positions have > 100 data points!";
	for(int p=0; p<lastPositionToConsiderPlusOne; ++p){
		if(_sums[p][from] == 0){
			throw "Not sufficient data for PMD estimation in read group '" + readGroupName + "': no observations for some reference alleles! Consider reducing the considered length.";
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
		y(p) = (double) _counts[p][from][to] / (double) _sums[p][from];
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
				for(int x=0; x<3; ++x){
					oldParams[x] = newParams[x];
				}

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
		mu += _counts[p][to][from];
		sum += _sums[p][to];
	}
	mu = mu / sum;

	double a = oldParams[1] / (1.0 - mu);
	double b = oldParams[2];
	double c = (oldParams[0] - mu) / (1.0 - mu);

	//calculate variance of estimator
	arma::mat Fisher = inv(J);

	//now return string
	if((mu + (1-mu) * ( a*exp(-b*maxReadLength) + c )) < 0)
		logfile->warning("Read group " + readGroupName + " was estimated to have a negative exponential damage pattern at the maxReadLength " + toString(maxReadLength) + ". This may be due to a lack of data");

	//if(a < 0) logfile->warning("Read group " + readGroupName + " was estimated to have an exponential damage pattern with a < 0! This may be due to a lack of data");
	return "Exponential[" + toString(a) + "," + toString(b) + "," + toString(c) + "]";
}

void TPMDTables::fitExponentialModel(int numNRIterations, double eps, std::string & filename, TLog* logfile){
	//TODO: do not have to fit exponential for all RG that are pooled
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups, fit and write exponential model
	for(size_t i=0; i<readGroups->size(); ++i){
		if(readGroups->readGroupInUse(i)){
			uint16_t index = readGroupMap->getIndex(i);
			out << readGroups->getName(i) << "\t" << forward[index]->fitExponentialModel(C, T, numNRIterations, eps, readGroups->getName(i), maxReadLength, logfile);
			out <<                           "\t" << reverse[index]->fitExponentialModel(G, A, numNRIterations, eps, readGroups->getName(i), maxReadLength, logfile) << "\n";
		}
	}
	out.close();
}

//---------------------------------------------------------------
//TPMDFunction
//---------------------------------------------------------------
void TPMDFunction::_parseParameters(const std::string & string){
	//expect string of the form NAME[P1,P2,...]
	//extract P1, P2, ... as a vector of doubles
	std::string tmp = readAfter(string, '[');
	fillVectorFromString(extractBefore(tmp, '['), _parameters, ',', true, true);
};

//---------------------------------------------------------------
// TPMDFunctionNoPMD
//---------------------------------------------------------------
TPMDFunctionNoPMD::TPMDFunctionNoPMD(const std::string & string){
	_parseParameters(string);

	if(_parameters.size() != 0){
		throw "Cannot initialize PMD function '" + PMDFunctionName_none + "': expected 0 but found " + toString(_parameters.size()) + " parameters!";
	}
};

//---------------------------------------------------------------
// TPMDFunctionSkoglund
//---------------------------------------------------------------
TPMDFunctionSkoglund::TPMDFunctionSkoglund(const std::string & string){
	_parseParameters(string);

	if(_parameters.empty()){
		//parameters missing: set to no PMD
		_parameters = {0.0, 0.0};
	} else {
		//parameters are provided
		if(_parameters.size() != 2){
			throw "Cannot initialize PMD function '" + PMDFunctionName_Skoglund + "': expected 2 (" + example() + ") but found " + toString(_parameters.size()) + " parameters!";
		}

		if(_parameters[0] < 0.0 || _parameters[0] > 1.0){
			throw "Cannot initialize PMD function '" + PMDFunctionName_Skoglund + "': p is outside [0,1]!";
		}
		if(_parameters[1] < 0.0){
			throw "Cannot initialize PMD function '" + PMDFunctionName_Skoglund + "': c must be >0!";
		}
	}
};

double TPMDFunctionSkoglund::prob(const uint16_t & pos) const {
	//Note: distance is zero based!
	return pow(1.0 - _parameters[0], (double) pos) * _parameters[0] + _parameters[1];
};

//---------------------------------------------------------------
// TPMDFunctionExponential
//--------------------------------------------------------------
TPMDFunctionExponential::TPMDFunctionExponential(const std::string & string){
	_parseParameters(string);

	if(_parameters.empty()){
		//parameters missing: set to no PMD
		_parameters = {0.0, 0.0, 0.0};
	} else {
		//parameters are provided
		if(_parameters.size() != 2){
			throw "Cannot initialize PMD function '" + PMDFunctionName_exponential + "': expected 3 (" + example() + ") but found " + toString(_parameters.size()) + " parameters!";
		}

		if(_parameters[1] < 0.0){
			throw "Cannot initialize PMD function '" + PMDFunctionName_exponential + "': b must be >0!";
		}

		if(_parameters[2] < 0.0){
			throw "Cannot initialize PMD function '" + PMDFunctionName_exponential + "': c must be >0!";
		}
	}
};

void TPMDFunctionExponential::learn(const TPMDTable & table, const Base & from, const Base & to){
	//resize parameters
	_parameters.resize(table.size() + 1); //include extra bin for sites beyond size (available in PMDTables)

	//extract counts in PMD direction and the inverse direction
	const countVec& pmdCounts = table[from][to];
	const countVec& pmdSums = table[from].sums();
	const countVec& invCounts = table[to][from];
	const countVec& invSums = table[to].sums();



};

double TPMDFunctionExponential::prob(const uint16_t & pos) const{
	//Note: distance is zero based!
	return _parameters[0] * exp(-_parameters[1] * (double) pos) + _parameters[2];
};

//---------------------------------------------------------------
// TPMDFunctionEmpiric
//---------------------------------------------------------------
TPMDFunctionEmpiric::TPMDFunctionEmpiric(const std::string & string){
	_parseParameters(string);

	if(_parameters.empty()){
		//parameters missing: set to no PMD
		_parameters = {0.0};
	} else {
		//parameters are provided
		for(auto& d : _parameters){
			if(d < 0.0 || d > 1.0){
				throw "Cannot initialize post mortem damage function '" + PMDFunctionName_empiric + "': some probabilities are outside [0,1]!";
			}
		}
	}
};

void TPMDFunctionEmpiric::learn(const TPMDTable & table, const Base & from, const Base & to){
	//resize parameters
	_parameters.resize(table.size() + 1); //include extra bin for sites beyond size (available in PMDTables)

	//extract counts in PMD direction and the inverse direction
	const countVec& pmdCounts = table[from][to];
	const countVec& pmdSums = table[from].sums();
	const countVec& invCounts = table[to][from];
	const countVec& invSums = table[to].sums();

	for(int p=0; p <= _parameters.size(); ++p){
		if(pmdSums[p] == 0 || invSums[p] == 0){
			_parameters[p] = 0.0;
		} else {
			double pmd = (double) pmdCounts[p] / pmdSums[p];
			double inv = (double) invCounts[p] / invSums[p];

			_parameters[p] = std::max(0.0, (pmd - inv)/(1.0 - inv)));
		}
	}
};

double TPMDFunctionEmpiric::prob(const uint16_t & pos) const{
	if(pos < _parameters.size()) return _parameters[pos];
	else return _parameters.back();
};

//------------------------------------------------------
//TPMDDoubleStrand
//------------------------------------------------------
void TPMDType::_initializeFunction(const std::string & pmdString, std::unique_ptr<TPMDFunction> & ptr){
	//parse string to get model. Options are
	// none
	// Empiric[0.5,0.3,...]
	// Skoglund[p,c]
	// Exponential[a,b,c]

	//extract function name
	std::string name = readBefore(pmdString , '[');

	if(name == PMDFunctionName_none){
		ptr = std::make_unique<TPMDFunctionNoPMD>(pmdString);
	} else if(name == PMDFunctionName_Skoglund){
		ptr = std::make_unique<TPMDFunctionSkoglund>(pmdString);
	} else if(name == PMDFunctionName_exponential){
		ptr = std::make_unique<TPMDFunctionExponential>(pmdString);
	} else if(name == PMDFunctionName_empiric){
		ptr = std::make_unique<TPMDFunctionEmpiric>(pmdString);
	} else {
		throw "Cannot initialize PMD function: unknown function '" + name + "'!. Use either " + PMDFunctionName_none + ", " + PMDFunctionName_Skoglund + ", " + PMDFunctionName_exponential + " or " + PMDFunctionName_empiric + ".";
	}
};

TPMDTypeDoubleStrand::TPMDTypeDoubleStrand(const std::string & string){
	//expect string pmdCT;pmdGA
	std::vector<std::string> tmp;
	fillVectorFromString(string, tmp, ';');
	if(tmp.size() != 2){
		throw "Cannot initialize PMD type " + PMDTypeName_doubleStrand + ": expect 2 functions but found " + tmp.size() + "!\nFunctions should be separated by ';'.\nProvided string: '" + string + "'.";
	}

	_initializeFunction(tmp[0], pmdCT);
	_initializeFunction(tmp[1], pmdGA);
};

std::string TPMDTypeDoubleStrand::functionString() const{
	return pmdCT->string() + ";" + pmdGA->string();
};

void TPMDTypeDoubleStrand::fillBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const{
	//Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	//no PMD for A and C
	baseLikelihoods[A] = baseLikelihoodsNoPMD[A];
	baseLikelihoods[T] = baseLikelihoodsNoPMD[T];

	//get relevant PMD probabilities
	double pmdProb_CT, pmdProb_GA;
	if(!base.isReverseStrand()){
		pmdProb_CT = pmdCT->prob(base.distFrom5Prime);
		pmdProb_GA = pmdGA->prob(base.distFrom3Prime);
	} else {
		pmdProb_CT = pmdCT->prob(base.distFrom3Prime);
		pmdProb_GA = pmdGA->prob(base.distFrom5Prime);
	}

	//add PMD
	baseLikelihoods[C] = (1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[C] + pmdProb_CT * baseLikelihoodsNoPMD[T];
	baseLikelihoods[G] = (1.0 - pmdProb_GA) * baseLikelihoodsNoPMD[G] + pmdProb_GA * baseLikelihoodsNoPMD[A];
};

//------------------------------------------------------
//TPostMortemDamage
//------------------------------------------------------
TPostMortemDamage::TPostMortemDamage(){
	_hasPMD = false;
};

void TPostMortemDamage::writeToFile(const BAM::TReadGroups & ReadGroups, const std::string filename){
	std::vector<std::string> header = {"ReadGroup", "Type", "Functions"};
	TOutputFile out(filename, header);

	//write for each read group
	for(auto& r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r){
		out << r->name_ID << _pmdObjects[r->id]->type() << _pmdObjects[r->id]->functionString() << std::endl;
	}
};

void TPostMortemDamage::_createPMDType(const std::string & type, const std::string & functions, std::unique_ptr<TPMDType> & ptr){
	//switch type
	if(type == PMDTypeName_none){
		ptr = std::make_unique<TPMDTypeNone>();
	} else if(type == PMDTypeName_doubleStrand){
		ptr = std::make_unique<TPMDTypeDoubleStrand>(functions);
	} else {
		throw "Cannot initialize PMD: unknown PMD type '" + type + "'!"
			 +"\nUse " + PMDTypeName_none + " or " + PMDTypeName_doubleStrand + ".";
	}
};

void TPostMortemDamage::_initializeFromFile(BAM::TReadGroups & ReadGroups, const std::string filename, TLog* logfile){
	//create an array of TPMD objects for each read group
	//also works if no parameters are provided (e.g. for estimation)
	_pmdObjects.resize(ReadGroups.size());

	//read from file for each read group
	TInputFile in(filename, header);
	if(in.numCols() != 3){
		throw "PMD file has the wrong format: expect three columns (read group, PMD type, PMD Functions) but found " + toString(in.numCols()) + "!";
	}

	//parse file that has structure: ReadGroup, Type, Functions
	std::vector<std::string> vec;
	while(in.read(vec)){
		if(ReadGroups.readGroupExists(vec[0])){ //ignore if it does not exist
			//get read group
			uint16_t readGroupId = ReadGroups.getId(vec[0]);

			//create type
			_createPMDType(vec[1], vec[2], _pmdObjects[readGroupId]);
		}
	}
	logfile->done();

	//test if we have a function for all read groups
	for(uint16_t i=0; i<ReadGroups.size(); ++i){
		if(!_pmdObjects[i]) throw "PMD not defined for read group '" + ReadGroups.getName(i) + "' in file '" + filename + "'!";
	}

	//check if there is PMD for at least one read group
	_hasPMD = false;
	for(auto& p : _pmdObjects){
		if(p->hasDamage()){
			_hasPMD = true;
			break;
		}
	}
};

void TPostMortemDamage::initialize(TParameters & params, BAM::TReadGroups & ReadGroups, TLog* logfile){
	if(params.parameterExists("pmd")){
		std::string pmdString = params.getParameterString("pmd");

		//expect either a file name or a type of the form "type:functions"
		//check if it is a file
		if(stringContains(pmdString, ":")){
			//not a file: initialize all read groups have the same pmd
			logfile->startIndent("PMD function used for all read groups (parameter pmd):");

			std::vector<std::string> tmp;
			fillVectorFromString(pmdString, tmp, ':');
			if(tmp.size() != 2){
				throw "Unable to understand PMD string '" + pmdString + "': expect format type:functions.";
			}

			_pmdObjects.resize(ReadGroups.size());
			for(size_t i=0; i<ReadGroups.size(); ++i){
				_createPMDType(tmp[0], tmp[1], _pmdObjects[i]);
			}
			_hasPMD = true;
			logfile->list(_pmdObjects[0]->type() + ":" + _pmdObjects[0]->functionString());
			logfile->endIndent();
		} else {
			//probably a file
			logfile->listFlush("Initializing PMD from file '" + pmdString + "' (parameter pmd) ...");
			_initializeFromFile(ReadGroups, pmdString, logfile);
			logfile->done();
		}
	} else {
		logfile->list("Assuming there is no PMD in the data. (use 'pmd' to add PMD definitions)");
	}
};

void TPostMortemDamage::calculateBaseLikelihoods(const BAM::TBase & base, const TBaseData & baseLikelihoodsNoPMD, TBaseData & baseLikelihoods) const{
	if(_hasPMD){
		_pmdObjects[ base.readGroupID ]->fillBaseLikelihoods(base, baseLikelihoodsNoPMD, baseLikelihoods);
	} else {
		//just copy
		baseLikelihoods = baseLikelihoodsNoPMD;
	}
};


}; //end namespace





