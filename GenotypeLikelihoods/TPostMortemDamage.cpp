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

using namespace coretools::str;

//Existing Functions
const std::string PMDFunctionName_none = "none";
const std::string PMDFunctionName_empiric = "Empiric";
const std::string PMDFunctionName_exponential = "Exponential";
//const std::string PMDFunctionName_Skoglund = "Skoglund";

//Existing types
const std::string PMDTypeName_none = "none";
const std::string PMDTypeName_singleStrand = "singleStrand";
const std::string PMDTypeName_doubleStrand = "doubleStrand";

//Estimation Parameters
const std::string PMDEstimationExponential_epsilon = "PMDExponentialEpsilon";
const std::string PMDEstimationExponential_numNR = "PMDExponentialNumNR";

//---------------------------------------------------------------
//TPMDFunction
//---------------------------------------------------------------
void TPMDFunction::_parseParameters(const std::string & string){
	//expect string of the form NAME[P1,P2,...]
	//extract P1, P2, ... as a vector of doubles
	std::string tmp = readAfter(string, '[');
	fillContainerFromString(extractBefore(tmp, '['), _parameters, ',', true, true, true);
};

//---------------------------------------------------------------
// TPMDFunctionNoPMD
//---------------------------------------------------------------
TPMDFunctionNoPMD::TPMDFunctionNoPMD(const std::string & string){
	_parseParameters(string);

	if(_parameters.size() != 0){
		throw "Cannot initialize PMD function '" + (std::string) PMDFunctionName_none + "': expected 0 but found " + toString(_parameters.size()) + " parameters!";
	}
};

//---------------------------------------------------------------
// TPMDFunctionSkoglund
//---------------------------------------------------------------
/*
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
*/

//---------------------------------------------------------------
// TPMDFunctionExponential
//--------------------------------------------------------------
TPMDFunctionExponential::TPMDFunctionExponential(const std::string & string){
	_parseParameters(string);

	if(_parameters.empty()){
		//parameters missing: set to no PMD
		_parameters = {1.0, 0.0, 0.0, 0.0};
	} else {
		//parameters are provided
		if(_parameters.size() != 2){
			throw "Cannot initialize PMD function '" + (std::string) PMDFunctionName_exponential + "': expected 4 (" + example() + ") but found " + toString(_parameters.size()) + " parameters!";
		}

		if(_parameters[0] < 1.0){
			throw "Cannot initialize PMD function '" + (std::string) PMDFunctionName_exponential + "': last position must be > 0!";
		}

		if(_parameters[2] < 0.0){
			throw "Cannot initialize PMD function '" + (std::string) PMDFunctionName_exponential + "': b must be > 0!";
		}

		if(_parameters[3] < 0.0){
			throw "Cannot initialize PMD function '" + (std::string) PMDFunctionName_exponential + "': c must be > 0!";
		}
	}

	_fillPMDProbabilities();
};

void TPMDFunctionExponential::_fillPMDProbabilities(){
	_lastPosition = round(_parameters[0]);
	_probs.resize(_lastPosition + 1);
	for(uint16_t p = 0; p <= _lastPosition; ++p){
		_probs[p] = _parameters[0] * exp(-_parameters[1] * (double) p) + _parameters[2];
	}
};

void TPMDFunctionExponential::parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile){
	if(!EstimationParameters.exists(PMDEstimationExponential_epsilon)){
		double eps = Params.getParameterWithDefault<double>(PMDEstimationExponential_epsilon, 0.001);
		EstimationParameters.add(PMDEstimationExponential_epsilon, eps);
		Logfile->list("Will consider the Newton-Raphson algorithm to have converged if the likelihood difference < " + toString(eps) + ". (parameter '" + PMDEstimationExponential_epsilon + "')");
	}

	if(!EstimationParameters.exists(PMDEstimationExponential_numNR)){
		double numNRIterations = Params.getParameterWithDefault<int>(PMDEstimationExponential_numNR, 100);
		EstimationParameters.add(PMDEstimationExponential_numNR, numNRIterations);
		Logfile->list("Will run up to " + toString(numNRIterations) + " Newton-Raphson iterations. (parameter '" + PMDEstimationExponential_numNR + ")");
	}
};

void TPMDFunctionExponential::_initialEstimatesOLS(const countVec & pmdCounts, const countVec & pmdSums, std::vector<double> & Parameters){
	//fill vector y to fit using OLS
	arma::vec y(_lastPosition + 1);
	double sumYSquared = 0.0;
	for(int p=0; p <= _lastPosition; ++p){
		y(p) = (double) pmdCounts[p] / (double) pmdSums[p];
		sumYSquared += y(p) * y(p);
	}

	//some variables
	double gammaStep = 0.01;
	double gammaTmp = -gammaStep + 0.00000001;
	double SSRold = pmdCounts.size() + 1;
	double SSRdiff = -1.0;
	arma::mat betaHat;
	arma::mat X(pmdCounts.size() + 1, 2);
	X.ones();

	//do until we get a small alpha
	while(fabs(gammaStep) > 0.00000001){
		while(SSRdiff < 0.0){
			//update gamma
			gammaTmp += gammaStep;

			//fill x
			for(int p = 0; p < pmdCounts.size(); ++p){
				X(p,1) = exp(-gammaTmp * (double) p);
			}

			betaHat = inv(X.t()*X) * X.t() * y;

			//calc sum of squares
			arma::mat tmp = (betaHat.t() * X.t() * y);
			double SSRnew = sumYSquared - tmp(0,0);
			SSRdiff = SSRnew - SSRold;
			SSRold = SSRnew;
		}

		//update alpha step
		gammaStep = -0.1 * gammaStep;
		SSRdiff = -1.0;
	}

	//set parameters
	Parameters = { betaHat(0), betaHat(1), gammaTmp };
};

void TPMDFunctionExponential::_fillFAndJacobian(arma::vec & F, arma::mat & J, const countVec & pmdCounts, const countVec& pmdSums, const std::vector<double> & Parameters){
	F.zeros();
	J.zeros();

	double weight, weightJ, tmp;
	double expMinusAlphaP;
	double dExpMinusAlphaP;
	for(int p = 0; p <= _lastPosition; ++p){
		//exp
		expMinusAlphaP = exp(-Parameters[2] * p);
		dExpMinusAlphaP = Parameters[1] * expMinusAlphaP;

		//first term
		//----------
		tmp = Parameters[0] + dExpMinusAlphaP;
		weight = pmdCounts[p] / tmp;
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
		J(1,2) -= weightJ * p * Parameters[0] * expMinusAlphaP;
		J(2,2) += weightJ * p * p * Parameters[0] * dExpMinusAlphaP;

		//second term
		//-----------
		tmp = (1.0 - Parameters[0] - dExpMinusAlphaP);
		weight = (pmdSums[p] - pmdCounts[p]) / tmp;
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
		J(1,2) += weightJ * p * (1.0 - Parameters[0]) * expMinusAlphaP;
		J(2,2) -= weightJ * p * p * (1.0 - Parameters[0]) * dExpMinusAlphaP;
	}

	//now fill in lower triangle of J
	J(1,0) = J(0,1);
	J(2,0) = J(0,2);
	J(2,1) = J(1,2);
};

double TPMDFunctionExponential::_calcLL(const countVec & pmdCounts, const countVec& pmdSums, const std::vector<double> & Parameters){
	double LL = 0.0;
	for(int p=0; p <= _lastPosition; ++p){
		double dExpMinusAlphaP = Parameters[1] * exp(-Parameters[2] * p);
		LL += pmdCounts[p] * log(Parameters[0] + dExpMinusAlphaP)
			+ (pmdSums[p] - pmdCounts[p]) * log(1.0 - Parameters[0] - dExpMinusAlphaP);
	}
	return LL;
};

void TPMDFunctionExponential::_estimateWithNewtonRaphson(const countVec & pmdCounts, const countVec& pmdSums, std::vector<double> & Parameters, const uint32_t & numNRIterations, const double & epsilon){
	//variables
	arma::mat J(3,3);
	arma::vec F(3);
	arma::mat JxF;

	//set starting values
	std::vector<double> newParams(3);

	//Conduct Newton-Raphson to refine
	//----------------------------------
	double oldLL = _calcLL(pmdCounts, pmdSums, Parameters);

	for(int i = 0; i<numNRIterations; ++i){
		_fillFAndJacobian(F, J, pmdCounts, pmdSums, Parameters);
		if(solve(JxF, J, F)){
			//estimate new params
			for(int x=0; x<3; ++x){
				newParams[x] = Parameters[x] - JxF(x);
			}

			//calculate LL at new location
			double LL = _calcLL(pmdCounts, pmdSums, newParams);

			//check if we accept or backtrack
			if(LL > oldLL){
				//store new params
				for(int x=0; x<3; ++x){
					Parameters = newParams;
				}

				//check if we stop NR
				if(LL - oldLL < epsilon){
					oldLL = LL;
					break;
				}
				oldLL = LL;
			}
		} else {
			std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << J << std::endl << std::endl;
			throw std::runtime_error("Issue solving JxF in TPMDTable::fitExponentialModel!");
		}
	}
};

void TPMDFunctionExponential::learn(const TPMDTable & Table, const genometools::Base & from, const genometools::Base & to, const TPMDEstimationParameters & EstimationParameters){
	//extract counts in PMD direction and the inverse direction
	const countVec& pmdCounts = Table[from][to];
	const countVec& pmdSums = Table[from].sums();
	const countVec& invCounts = Table[to][from];
	const countVec& invSums = Table[to].sums();

	//Check if we have sufficient data
	//find last entry with counts
	_lastPosition = -1;
	for(int p = pmdCounts.size()-1; p >= 0; --p){
		if(pmdSums[p] > 100){
			_lastPosition = p;
			break;
		}
	}
	_lastPosition;

	if(_lastPosition < 10) throw "Not sufficient data to fit exponential PMD model: less than the ten first positions have > 100 data points!\nConsider pooling read groups (parameter poolReadGroups).";
	for(int p=0; p<_lastPosition; ++p){
		if(pmdSums[p] == 0){
			throw "Not sufficient data to fit exponential PMD model: no observations for some reference alleles!<nConsider reducing the relevant length (parameter length).";
		}
	}

	//get initial estimates via OLS
	std::vector<double> Parameters;
	_initialEstimatesOLS(pmdCounts, pmdSums, Parameters);

	//run Newton-Raphson
	_estimateWithNewtonRaphson(pmdCounts, pmdSums, Parameters, EstimationParameters[PMDEstimationExponential_epsilon], EstimationParameters[PMDEstimationExponential_numNR]);

	//transform parameters
	//the exponential PMD model is f(C->T) = mu + (1-mu) *[ a*exp(-b * position) + c ]
	//but we fitted f(C->T) = alpha + beta * exp(-gamma * position).
	//Hence we have:
	// mu is estimated from T->C transitions
	// a =  beta / (1 - mu)
	// b = gamma
	// c = (alpha - mu) / (1 - mu)
	double mu = 0.0; double sum = 0.0;
	for(int p=0; p <= _lastPosition; ++p){
		mu += invCounts[p];
		sum += invSums[p];
	}
	mu = mu / sum;

	//store parameters, including lastPosition
	_parameters = { (double) _lastPosition,
					Parameters[1] / (1.0 - mu),
					Parameters[2],
				   (Parameters[0] - mu) / (1.0 - mu) };

	_fillPMDProbabilities();

	//check if pattern is negativ
	if(_probs[_lastPosition] < 0){
		throw "Estimation resulted in negative PMD at high positions!\nThis is likely be due to limited data. Consider pooling read groups (parameter poolReadGroups).";
	}
	if(Parameters[1] < 0){
		throw "Estimation resulted in a < 0!\nThis is likely due to limited data. Consider pooling read groups (parameter poolReadGroups).";
	}
};

double TPMDFunctionExponential::prob(const uint16_t & pos) const{
	//Note: distance is zero based!
	//model is fit up to _lastPosition. We assume constant PMD after that
	if(pos < _lastPosition){
		return _probs[pos];
	} else {
		return _probs[_lastPosition];
	}
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
				throw "Cannot initialize post mortem damage function '" + (std::string) PMDFunctionName_empiric + "': some probabilities are outside [0,1]!";
			}
		}
	}
};

void TPMDFunctionEmpiric::learn(const TPMDTable & Table, const genometools::Base & from, const genometools::Base & to, const TPMDEstimationParameters & EstimationParameters){
	//resize parameters
	_parameters.resize(Table.size() + 1); //include extra bin for sites beyond size (available in PMDTables)

	//extract counts in PMD direction and the inverse direction
	const countVec& pmdCounts = Table[from][to];
	const countVec& pmdSums = Table[from].sums();
	const countVec& invCounts = Table[to][from];
	const countVec& invSums = Table[to].sums();

	for(int p=0; p <= _parameters.size(); ++p){
		if(pmdSums[p] == 0 || invSums[p] == 0){
			_parameters[p] = 0.0;
		} else {
			double pmd = (double) pmdCounts[p] / pmdSums[p];
			double inv = (double) invCounts[p] / invSums[p];

			_parameters[p] = std::max(0.0, (pmd - inv)/(1.0 - inv));
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
	/*
	} else if(name == PMDFunctionName_Skoglund){
		ptr = std::make_unique<TPMDFunctionSkoglund>(pmdString);
	*/
	} else if(name == PMDFunctionName_exponential){
		ptr = std::make_unique<TPMDFunctionExponential>(pmdString);
	} else if(name == PMDFunctionName_empiric){
		ptr = std::make_unique<TPMDFunctionEmpiric>(pmdString);
	} else {
		throw "Cannot initialize PMD function: unknown function '" + name + "'!. Use either " + PMDFunctionName_none + ", " + PMDFunctionName_exponential + " or " + PMDFunctionName_empiric + ".";
	}
};

TPMDTypeDoubleStrand::TPMDTypeDoubleStrand(const std::vector<std::string> & Details){
	//expect three elements: type, pmdCT, pmdGA
	if(Details.size() != 2){
		throw "Cannot initialize PMD type " + (std::string)  PMDTypeName_doubleStrand + ": expect 3 entries but found " + toString(Details.size()) + "!"
			+ "\nProvided string: '" + concatenateString(Details, ':') + "'.";
			+ "\nExpect string of the form '" + PMDTypeName_doubleStrand + "':functionCT:functionGA'.";
	}

	_initializeFunction(Details[1], _pmdCT);
	_initializeFunction(Details[2], _pmdGA);
};

std::string TPMDTypeDoubleStrand::functionString() const{
	return PMDTypeName_doubleStrand + ":" + _pmdCT->string() + ":" + _pmdGA->string();
};

void TPMDTypeDoubleStrand::parseEstimationParameters(TPMDEstimationParameters & EstimationParameters, TParameters & Params, TLog* Logfile){
	_pmdCT->parseEstimationParameters(EstimationParameters, Params, Logfile);
	_pmdGA->parseEstimationParameters(EstimationParameters, Params, Logfile);
};

void TPMDTypeDoubleStrand::estimate(const TPMDTableReadGroup & PMDTable, const TPMDEstimationParameters & EstimationParameters){
	//Note: TPMDTables stores bases as during sequencing (not as after mapping)
	//Assumption: C->T pattern is the same for forward and reverse reads from their respective 5-prime ends.
	TPMDTable from5(PMDTable[forward5]);
	from5.add(PMDTable[reverse5]);
	_pmdCT->learn(from5, BAM::C, BAM::T, EstimationParameters);

	//Assumption: G->A pattern is the same for forward and reverse reads from their respective 3-prime ends.
	TPMDTable from3(PMDTable[forward3]);
	from3.add(PMDTable[reverse3]);
	_pmdGA->learn(from3, BAM::G, BAM::A, EstimationParameters);
};

void TPMDTypeDoubleStrand::fillBaseLikelihoods(const BAM::TSequencedBase & base, const TBaseProbabilities & baseLikelihoodsNoPMD, TBaseProbabilities & baseLikelihoods) const {
	//Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	//no PMD for A and C
	baseLikelihoods[BAM::A] = baseLikelihoodsNoPMD[BAM::A];
	baseLikelihoods[BAM::T] = baseLikelihoodsNoPMD[BAM::T];

	//get relevant PMD probabilities
	double pmdProb_CT, pmdProb_GA;
	if(!base.isReverseStrand()){
		pmdProb_CT = _pmdCT->prob(base.distFrom5Prime);
		pmdProb_GA = _pmdGA->prob(base.distFrom3Prime);
	} else {
		pmdProb_CT = _pmdCT->prob(base.distFrom3Prime);
		pmdProb_GA = _pmdGA->prob(base.distFrom5Prime);
	}

	//add PMD
	baseLikelihoods[BAM::C] = (1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[BAM::C].get() + pmdProb_CT * baseLikelihoodsNoPMD[BAM::T].get();
	baseLikelihoods[BAM::G] = (1.0 - pmdProb_GA) * baseLikelihoodsNoPMD[BAM::G].get() + pmdProb_GA * baseLikelihoodsNoPMD[BAM::A].get();
};

void TPMDTypeDoubleStrand::simulatePMD(BAM::TSequencedBase & base, TRandomGenerator & RandomGenerator) const {
	simulatePMD(base.base, base.distFrom5Prime, base.distFrom3Prime, base.isReverseStrand(), RandomGenerator);
};

void TPMDTypeDoubleStrand::simulatePMD(genometools::Base & base, const uint16_t & DistFrom5Prime, const uint16_t & DistFrom3Prime, const bool & IsReverseStrand, TRandomGenerator & RandomGenerator) const{
	//simulate PMD
	if(!IsReverseStrand){
		//forward strand
		if(base == BAM::C){
			if(RandomGenerator.getRand() < _pmdCT->prob(DistFrom5Prime)){
				base = BAM::T;
			}
		} else if(base == BAM::G){
			if(RandomGenerator.getRand() < _pmdGA->prob(DistFrom3Prime)){
				base = BAM::A;
			}
		}
	} else {
		//reverse strand
		if(base == BAM::C){
			if(RandomGenerator.getRand() < _pmdCT->prob(DistFrom3Prime)){
				base = BAM::T;
			}
		} else if(base == BAM::G){
			if(RandomGenerator.getRand() < _pmdGA->prob(DistFrom5Prime)){
				base = BAM::A;
			}
		}
	}
};

//------------------------------------------------------
//TPostMortemDamage
//------------------------------------------------------
TPostMortemDamage::TPostMortemDamage(){
	_hasPMD = false;
};

void TPostMortemDamage::writeToFile(const BAM::TReadGroups & ReadGroups, const std::string filename){
	std::vector<std::string> header = {"ReadGroup", "Type", "Functions"};
	coretools::TOutputFile out(filename, header);

	//write for each read group
	for(auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r){
		out << r->name_ID << _pmdObjects[r->id]->functionString() << std::endl;
	}
};

void TPostMortemDamage::writeToFile(const BAM::TReadGroups & ReadGroups, const BAM::TReadGroupMap & ReadGroupMap, const std::string filename){
	std::vector<std::string> header = {"ReadGroup", "Type", "Functions"};
	coretools::TOutputFile out(filename, header);

	//write for each read group
	for(auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r){
		out << r->name_ID << _pmdObjects[ ReadGroupMap.pooledIndex(r->id) ]->functionString() << std::endl;
	}
};

void TPostMortemDamage::_createPMDType(const std::string & pmdString, std::shared_ptr<TPMDType> & ptr){
	//split by ':'
	std::vector<std::string> details;
	fillContainerFromString(pmdString, details, ":");

	//switch type
	if(details[0] == PMDTypeName_none){
		ptr = std::make_shared<TPMDTypeNone>();
	} else if(details[0] == PMDTypeName_doubleStrand){
		ptr = std::make_shared<TPMDTypeDoubleStrand>(details);
	} else {
		throw "Cannot initialize PMD: unknown PMD type '" + details[0] + "'!"
			 +"\nUse " + PMDTypeName_none + " or " + PMDTypeName_doubleStrand + ".";
	}
};

void TPostMortemDamage::_initializeFromString(const std::string & pmdString, TLog* logfile){
	//not a file: initialize all read groups have the same pmd
	logfile->startIndent("PMD function used for all read groups:");

	for(auto& p : _pmdObjects){
		_createPMDType(pmdString, p);
	}

	//report
	logfile->list(_pmdObjects[0]->type() + ":" + _pmdObjects[0]->functionString());
	logfile->endIndent();
};

void TPostMortemDamage::_initializeFromFile(const BAM::TReadGroups & ReadGroups, const std::string & filename, TLog* logfile, std::vector<uint16_t> & ReadGroupsWithoutPMD){
	//create an array of TPMD objects for each read group
	//also works if no parameters are provided (e.g. for estimation)
	//read from file for each read group

	logfile->listFlush("Initializing PMD from file '" + filename + "' ...");
	coretools::TInputFile in(filename, {"readGroup", "pmd"}, "/t", "//");

	//parse file that has structure: ReadGroup, Type, Functions
	std::vector<std::string> vec;
	while(in.read(vec)){
		if(ReadGroups.readGroupExists(vec[0])){ //ignore if it does not exist
			//get read group
			uint16_t readGroupId = ReadGroups.getId(vec[0]);

			//create type
			_createPMDType(vec[1], _pmdObjects[readGroupId]);
		}
	}
	logfile->done();

	//check if we have a function for all read groups
	//create no-PMD types for all remaining ones and return their indexes
	for(uint16_t i=0; i<ReadGroups.size(); ++i){
		if(!_pmdObjects[i]){
			_createPMDType("non", _pmdObjects[i]);
			ReadGroupsWithoutPMD.push_back(i);
		}
	}
};

void TPostMortemDamage::_setHasDamage(){
	//check if there is PMD for at least one read group
	_hasPMD = false;
	for(auto& p : _pmdObjects){
		if(p->hasDamage()){
			_hasPMD = true;
			break;
		}
	}
};

void TPostMortemDamage::initialize(const std::string & pmdString, const BAM::TReadGroups & ReadGroups, TLog* Logfile, std::vector<uint16_t> & ReadGroupsWithoutPMD){
	if(_hasPMD){
		throw std::runtime_error("void TPostMortemDamage::initialize(const std::string & pmdString, const BAM::TReadGroups & ReadGroups, TLog* Logfile, std::vector<uint16_t> & ReadGroupsWithoutPMD): Models already initialized!");
	}

	//prepare objects
	ReadGroupsWithoutPMD.clear();
	_pmdObjects.resize(ReadGroups.size());

	//expect either a file name or a type of the form "type:functions"
	//check if it is a file
	if(stringContains(pmdString, ":")){
		_initializeFromString(pmdString, Logfile);
	} else {
		//probably a file
		_initializeFromFile(ReadGroups, pmdString, Logfile, ReadGroupsWithoutPMD);
	}

	//check if there is PMD for at least one read group
	_setHasDamage();
};

void TPostMortemDamage::fillBaseLikelihoods(const BAM::TSequencedBase & base, const TBaseProbabilities & baseLikelihoodsNoPMD, TBaseProbabilities & baseLikelihoods) const{
	if(_hasPMD){
		_pmdObjects[base.readGroupID]->fillBaseLikelihoods(base, baseLikelihoodsNoPMD, baseLikelihoods);
	} else {
		//just copy
		baseLikelihoods = baseLikelihoodsNoPMD;
	}
};


}; //end namespace





