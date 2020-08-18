/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TSequencingErrorModel.h"

namespace GenotypeLikelihoods{

//--------------------------------------------------------------------
// TRecalibrationEMModelCovariateDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
bool TSequencingErrorCovariateDefinition::parse(const std::string & modelString, std::string & error){
	std::vector<std::string> tmp;
	fillVectorFromString(modelString, tmp, ';', true);
	for(std::string s : tmp){
		size_t pos = s.find(':');
		if(pos == std::string::npos){
			//intercept?
			size_t pos_1 = s.find('[');
			size_t pos_2 = s.find(']');
			if(s.substr(0, pos_1) != "intercept"){
				error = "Unable to understand model string '" + modelString + "': expecting an '=' or the function 'intercept'!";
				return false;
			} else {
				if(pos_1 == std::string::npos){
					if(pos_2 != std::string::npos){
						error = "Unable to understand model string '" + modelString + "': missing '['";
						return false;
					}
					intercept = "0.0";
				} else {
					if(pos_2 == std::string::npos){
						error = "Unable to understand model string '" + modelString + "': missing ']'";
						return false;
					}
					intercept = s.substr(pos_1+1, pos_2-pos_1-1);
				}
			}
		} else {
			covariateFunctions.emplace_back(s.substr(0, pos), s.substr(pos+1));
		}
	}
	return true;
};

void TSequencingErrorCovariateDefinition::setIntercept(const double Intercept){
	intercept = toString(Intercept);
};

void TSequencingErrorCovariateDefinition::addCovariate(const std::string covariate, const std::string function){
	covariateFunctions.emplace_back(covariate, function);
};

std::string TSequencingErrorCovariateDefinition::getModelString(){
	std::string modelString = "intercept[" + intercept + "]";
	for(auto& it : covariateFunctions){

		modelString += ";" + it.covariate + ":" + it.function;
	}
	return modelString;
};

//--------------------------------------------------------------------
// TRecalibrationEMModelCovariateList
//--------------------------------------------------------------------
TSequencingErrorCovariateList::TSequencingErrorCovariateList(){
	numParameters = intercept.numParameters();
};

TSequencingErrorCovariateList::~TSequencingErrorCovariateList(){
	_clear();
};

TSequencingErrorCovariateList::TSequencingErrorCovariateList(TSequencingErrorCovariateList&& other){
	//copy from other
	covariates = std::move(other.covariates);
	pointerToCovariateFunctions = std::move(other.pointerToCovariateFunctions);
	intercept = std::move(other.intercept);
	numParameters = other.numParameters;

	//set intercept of other to zero
	other.intercept.setIntercept(0.0);
};

TSequencingErrorCovariateList& TSequencingErrorCovariateList::operator=(TSequencingErrorCovariateList&& other){
	if(&other != this){ //don't copy yourself
		//clear
		_clear();

		//copy from other
		covariates = std::move(other.covariates);
		pointerToCovariateFunctions = std::move(other.pointerToCovariateFunctions);
		intercept = std::move(other.intercept);
		numParameters = other.numParameters;

		//set intercept of other to zero
		other.intercept.setIntercept(0.0);
	}
	//return
	return *this;
};

void TSequencingErrorCovariateList::_clear(){
	for(auto* it : covariates){
		delete it;
	}
	covariates.clear();
	pointerToCovariateFunctions.clear();
	intercept.setIntercept(0.0);
};

void TSequencingErrorCovariateList::createCovariatesAndIntercept(TSequencingErrorCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable){
	//include intercept
	if(!covariateMap.intercept.empty()){
		std::vector<std::string> vec = {covariateMap.intercept};
		intercept.initialize(0, vec);
	}

	//create covariates
	numParameters = intercept.numParameters();
	for(TRecalibrationEMModelCovariateDefinitionIterator it = covariateMap.begin(); it != covariateMap.end(); ++it){
		//create function for each covariate
		if(it->covariate == SequencingErrorCovariateName_none){
			continue;
		} else if(it->covariate == SequencingErrorCovariateName_quality){
			covariates.emplace_back(new TSequencingErrorCovariate_quality(numParameters, it->function, dataTable));
		} else if(it->covariate == SequencingErrorCovariateName_position){
			covariates.emplace_back(new TSequencingErrorCovariate_position(numParameters, it->function, dataTable));
		} else if(it->covariate == SequencingErrorCovariateName_context){
			covariates.emplace_back(new TSequencingErrorCovariate_context(numParameters, it->function, dataTable));
		} else if(it->covariate == SequencingErrorCovariateName_fragmentLength){
			covariates.emplace_back(new TSequencingErrorCovariate_fragmentLength(numParameters, it->function, dataTable));
		} else if(it->covariate == SequencingErrorCovariateName_mappingQuality){
			covariates.emplace_back(new TSequencingErrorCovariate_mappingQuality(numParameters, it->function, dataTable));
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}

		//add new parameters
		numParameters += covariates.back()->numParameters();
	}

	//summarize
	_storePointersToCovariateFunctions();
};

void TSequencingErrorCovariateList::createCovariatesAndIntercept(TSequencingErrorCovariateDefinition & covariateMap){
	//include intercept
	std::vector<std::string> vec = {covariateMap.intercept};
	intercept.initialize(0, vec);

	//create covariates
	numParameters = intercept.numParameters();
	for(TRecalibrationEMModelCovariateDefinitionIterator it = covariateMap.begin(); it != covariateMap.end(); ++it){
		//create function for each covariate
		if(it->covariate == SequencingErrorCovariateName_none){
			continue;
		} else if(it->covariate == SequencingErrorCovariateName_quality){
			covariates.emplace_back(new TSequencingErrorCovariate_quality(numParameters, it->function));
		} else if(it->covariate == SequencingErrorCovariateName_position){
			covariates.emplace_back(new TSequencingErrorCovariate_position(numParameters, it->function));
		} else if(it->covariate == SequencingErrorCovariateName_context){
			covariates.emplace_back(new TSequencingErrorCovariate_context(numParameters, it->function));
		} else if(it->covariate == SequencingErrorCovariateName_fragmentLength){
			covariates.emplace_back(new TSequencingErrorCovariate_fragmentLength(numParameters, it->function));
		} else if(it->covariate == SequencingErrorCovariateName_mappingQuality){
			covariates.emplace_back(new TSequencingErrorCovariate_mappingQuality(numParameters, it->function));
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}

		//add new parameters
		numParameters += covariates.back()->numParameters();
	}

	//summarize
	_storePointersToCovariateFunctions();
};

void TSequencingErrorCovariateList::_storePointersToCovariateFunctions(){
	//add intercept
	pointerToCovariateFunctions.emplace_back(&intercept);

	//add covariates
	for(auto & cov : covariates){
		//store function pointer
		pointerToCovariateFunctions.push_back(cov->getPointerToFunction());
	}
};

TSequencingErrorCovariateDefinition TSequencingErrorCovariateList::getCovariateDefinition() const{
	TSequencingErrorCovariateDefinition def;
	def.setIntercept(intercept.getIntercept());
	for(const auto & cov : covariates){
		def.addCovariate(cov->name(), cov->functionString());
	}
	return def;
}

//--------------------------------------------------------------------
// TSequencingErrorRhoStorage
//--------------------------------------------------------------------
TSequencingErrorRhoStorage::TSequencingErrorRhoStorage(){
	reset();
};

TSequencingErrorRhoStorage::TSequencingErrorRhoStorage(const std::string & def, std::string & error){
	set(def, error);
};

TSequencingErrorRhoStorage::TSequencingErrorRhoStorage(const TSequencingErrorRhoStorage & other){
	rho = other.rho;
};

void TSequencingErrorRhoStorage::operator=(const TSequencingErrorRhoStorage & other){
	rho = other.rho;
};

double TSequencingErrorRhoStorage::operator()(const uint8_t & from, const uint8_t & to){
	return rho[from][to];
};

void TSequencingErrorRhoStorage::reset(){
	for(int a=0; a<4; ++a){
		for(int b=0; b<4; ++b){
			if(a==b){
				rho[a][b] = 0.0;
			} else {
				rho[a][b] = 1.0 / 3.0;
			}
		}
	}
};

bool TSequencingErrorRhoStorage::set(const std::string & def, std::string & error){
	std::vector<std::string> vec;
	std::string s = def;
	fillVectorFromString(def, vec, ';');
	if(vec.size() != 4){
		error = "Rho matrix has " + toString(vec.size()) + " instead of 4 rows!";
		return false;
	}

	//parse rows
	std::vector<double> r;
	for(size_t a = 0; a<vec.size(); ++a){
		std::string& row = vec[a];
		trimString(row, "()");
		fillVectorFromString(row, r, ',');
		if(r.size() != 4){
			error = "Rho matrix has " + toString(r.size()) + " instead of 4 columns for row " + toString(a+1) + "!";
			return false;
		}

		//fill
		for(size_t b=0; b<4; ++b){
			if(a == b){
				rho[a][b] = 0.0;
			} else {
				rho[a][b] = r[b];
			}
		}
		++a;
	}

	return true;
};

std::string TSequencingErrorRhoStorage::getDefinition() const{
	std::string def;
	for(int a=0; a<4; ++a){
		if(a>0){
			def += ";";
		}
		for(int b=0; b<4; ++b){
			if(b>0){
				def += ",";
			}
			if(a!=b){
				def += toString(rho[a][b]);
			} else {
				def += "-";
			}
		}
	}
	return def;
};

//--------------------------------------------------------------------
// TSequencingErrorRho
//--------------------------------------------------------------------
void TSequencingErrorRho::operator=(const TSequencingErrorRhoStorage & other){
	TSequencingErrorRhoStorage::operator =(other);
};

void TSequencingErrorRho::fillBaseLikelihoods(const Base base, const double epsilon, TBaseData & baseLikelihoods) const{
	baseLikelihoods[base] = 1.0 - epsilon;
	if(base == A){
		baseLikelihoods[C] = epsilon * rho[C][A];
		baseLikelihoods[G] = epsilon * rho[G][A];
		baseLikelihoods[T] = epsilon * rho[T][A];
	} else if(base == C){
		baseLikelihoods[A] = epsilon * rho[A][C];
		baseLikelihoods[G] = epsilon * rho[G][C];
		baseLikelihoods[T] = epsilon * rho[T][C];
	} else if(base == G){
		baseLikelihoods[A] = epsilon * rho[A][G];
		baseLikelihoods[C] = epsilon * rho[C][G];
		baseLikelihoods[T] = epsilon * rho[T][G];
	} else {
		baseLikelihoods[A] = epsilon * rho[A][T];
		baseLikelihoods[C] = epsilon * rho[C][T];
		baseLikelihoods[G] = epsilon * rho[G][T];
	}
};

void TSequencingErrorRho::prepareEstimationFromEMWeights(){
	for(int b=0; b<4; ++b){
		for(int a=0; a<4; ++a){
			rho[a][b] = 0.0;
		}
	}
};

void TSequencingErrorRho::addBaseForEstimation(const Base & base, const TBaseData & EMWeights){
	if(base == A){
		rho[C][A] += EMWeights.at(C);
		rho[G][A] += EMWeights.at(G);
		rho[T][A] += EMWeights.at(T);
	} else if(base == C){
		rho[A][C] += EMWeights.at(A);
		rho[G][C] += EMWeights.at(G);
		rho[T][C] += EMWeights.at(T);
	} else if(base == G){
		rho[A][G] += EMWeights.at(A);
		rho[C][G] += EMWeights.at(C);
		rho[T][G] += EMWeights.at(T);
	} else {
		rho[A][T] += EMWeights.at(A);
		rho[C][T] += EMWeights.at(C);
		rho[G][T] += EMWeights.at(G);
	}
};

void TSequencingErrorRho::estimate(){
	//calculate denominators
	std::vector<double> denom(4, 0.0);
	for(int a=0; a<4; ++a){
		for(int b=0; b<4; ++b){
			if(a!=b){
				denom[a] += rho[a][b];
			}
		}
	}

	//scale rho
	for(int a=0; a<4; ++a){
		for(int b=0; b<4; ++b){
			if(a!=b){
				rho[a][b] /= denom[a];
			}
		}
	}

	//set denominators = 0.0
	for(int a=0; a<4; ++a){
		rho[a][a] = 0.0;
	}
};

//--------------------------------------------------------------------
// TSequencingErrorModel
//--------------------------------------------------------------------
TSequencingErrorModel::TSequencingErrorModel(TSequencingErrorModelDefinition & modelDef, TLog* Logfile){
	logfile = Logfile;

	//create covariates
	_covariates.createCovariatesAndIntercept(modelDef.covariates);
	_rho = modelDef.rho;

	//prepare Newton-Raphson variables
	setNewtonRaphosnParamsToZero();
};

TSequencingErrorModel::TSequencingErrorModel(TSequencingErrorModelDefinition & modelDef, TRecalibrationEMDataTable* dataTable, TLog* Logfile){
	logfile = Logfile;

	//create covariates
	_covariates.createCovariatesAndIntercept(modelDef.covariates, dataTable);
	_rho = modelDef.rho;

	//prepare Newton-Raphson variables
	setNewtonRaphosnParamsToZero();
};

bool TSequencingErrorModel::checkParameterRange(TRecalibrationEMDataTable* dataTable, std::string & error){
	for(auto & cov : _covariates.covariates){
		if(!cov->checkParameterRange(dataTable)){
			error = "Function for covariate " + cov->name() + " does not cover full range of data";
			return false;
		}
	}
	return true;
};

void TSequencingErrorModel::_initializeDerivatives(){
	//intercept
	size_t numNonZeroFirstDeriv = _covariates.intercept.numNonZeroFirstDerivatives();
	size_t numNonZeroSecondDeriv = _covariates.intercept.numNonZeroSecondDerivatives();

	//covariates
	for(const auto & cov : _covariates.covariates){
		numNonZeroFirstDeriv += cov->numNonZeroFirstDerivatives();
		numNonZeroSecondDeriv += cov->numNonZeroSecondDerivatives();
	}
	_firstDerivatives.resize(numNonZeroFirstDeriv);
	_secondDerivatives.resize(numNonZeroSecondDeriv);
};

double TSequencingErrorModel::_calcEpsilon(const double eta) const{
	if(eta > 23.03){
		return 0.9999999999;
	}
	if(eta < -23.03){
		return 0.0000000001;
	}

	return 1.0 / (1.0 + exp(-eta));
};

double TSequencingErrorModel::getErrorRate(const BAM::TBase & base) const{
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = _covariates.intercept.getEtaTerm();

	for(const auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(base);
	}

	return _calcEpsilon(eta);
};

void TSequencingErrorModel::fillBaseLikelihoods(const BAM::TBase & base, TBaseData & baseLikelihoods) const{
	//first calculate epsilon
	double eta = _covariates.intercept.getEtaTerm();
	for(const auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(base);
	}

	//then calculate base likelihoods
	_rho.fillBaseLikelihoods(base.base, _calcEpsilon(eta), baseLikelihoods);
};

std::string TSequencingErrorModel::getCovariateDefinition() const{
	return _covariates.getCovariateDefinition().getModelString();
};

std::string TSequencingErrorModel::getRhoDefinition() const{
	return _rho.getDefinition();
};

//-------------------------------------------------
//functions to estimate rho
//-------------------------------------------------
void TSequencingErrorModel::prepareRhoEstimationFromEMWeights(){
	_rho.prepareEstimationFromEMWeights();
};

void TSequencingErrorModel::addBaseForRhoEstimation(BAM::TBase & base, const TBaseData & EMWeights){
	_rho.addBaseForEstimation(base.base, EMWeights);
};

void TSequencingErrorModel::estimateRho(){
	_rho.estimate();
};

//-------------------------------------------------
//functions for estimation
//-------------------------------------------------
void TSequencingErrorModel::setNewtonRaphosnParamsToZero(){
	_Jacobian.resize(_covariates.numParameters, _covariates.numParameters);
	_F.resize(_covariates.numParameters);
	_JxF.resize(_covariates.numParameters, 1);

	_Jacobian.zeros();
	_F.zeros();

	_initializeDerivatives();

	_numSitesAdded = 0;
	_NRconverged = false;
	_NRStepAccepted = false;
};

void TSequencingErrorModel::setQToZero(){
	if(!_NRconverged){
		_oldQ = _Q;
		_Q = 0.0;
	}
};

void TSequencingErrorModel::addToQ(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d){
	if(!_NRconverged){
		//get error rate
		double eps = getErrorRate(base);
		//calculate sum_bbar [ Ind(bbar=d)log(1-eps) + Ind(bbar!=d)log(eps) ]
		_Q += EM_weights_bbar_given_d.at( base.base ) * log(1.0 - eps) + (1.0 - EM_weights_bbar_given_d.at( base.base )) * log(eps);
	}
};

void TSequencingErrorModel::addToFandJacobian(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d){
	// 1) Calculate epsilon
	//--------------------
	double eps = getErrorRate(base);

	// 2 ) fill derivatives
	//--------------------
	_firstDerivatives.restart();
	_secondDerivatives.restart();

	//fill derivatives of intercept
	_covariates.intercept.fillDerivatives(0.0, _firstDerivatives, _secondDerivatives);

	//fill derivatives of covariates
	for(const auto & cov : _covariates.covariates){
		cov->fillDerivatives(base, _firstDerivatives, _secondDerivatives);
	}

	// 3) add derivatives to F and Jacobian
	//calculate weights
	double weight1 = (1.0 - EM_weights_bbar_given_d.at(base.base))*(1.0 - eps) - EM_weights_bbar_given_d.at(base.base) * eps;

	std::cout << "weight1 = " << weight1 << ", EM_weight = " << EM_weights_bbar_given_d.at(base.base) << ", eps = " << eps << std::endl;

	double weight2 = (1.0 - eps) * eps;

	//add first derivatives
	for(TRecalibrationEMFirstDerivativesIterator it = _firstDerivatives.begin(); it != _firstDerivatives.end(); ++it){
		//add to F
		_F(it->index) += weight1 * it->derivative;

		//add to J
		for(TRecalibrationEMFirstDerivativesIterator it2 = it; it2 != _firstDerivatives.end(); ++it2){
			_Jacobian(it->index, it2->index) += weight2 * it->derivative * it2->derivative;
		}
	}

	//add second derivatives to Jacobian (happens to have the same weigth as F!)
	for(auto& it : _secondDerivatives){
		_Jacobian(it.index1, it.index2) += weight1 * it.derivative;
	}

	++_numSitesAdded;
};

bool TSequencingErrorModel::solveJxF(){
	if(_NRconverged){
		return true;
	} else {
		//Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
		for(int i=0; i<(_covariates.numParameters-1); ++i){
			for(unsigned int j=i+1; j<_covariates.numParameters; ++j){
				//copy from upper triangle to lower triangle
				_Jacobian(j,i) = _Jacobian(i,j);
			}
		}

		//scale F and J by 1/#sites
		_Jacobian = _Jacobian / (double) _numSitesAdded;
		_F = _F / (double) _numSitesAdded;

		std::cout << "DEBUG!!!" << std::endl;
		printJacobianToStdOut();
		printFToStdOut();

		//now solve J^-1 x F
		return solve(_JxF, _Jacobian, _F);
	}
};

void TSequencingErrorModel::proposeNewParameters(double & lambda){
	if(!_NRStepAccepted){
		uint16_t index = 0;
		for(const auto it : _covariates.pointerToCovariateFunctions){
			it->proposeNewParameters(_JxF, index, lambda);
		}
	}
};

bool TSequencingErrorModel::acceptProposedParametersBasedOnQ(){
	if(_NRStepAccepted) return true;
	if(_Q > _oldQ){
		_NRStepAccepted = true;
	} else {
		_NRStepAccepted = false;
		_Q = _oldQ;

		for(const auto it : _covariates.pointerToCovariateFunctions){
			it->rejectProposedParameters();
		}
	}
	return _NRStepAccepted;
};

void TSequencingErrorModel::adjustParametersPostEstimation(){
	for(const auto it : _covariates.pointerToCovariateFunctions){
		_covariates.intercept.addToIntercept(it->adjustParametersPostEstimation());
	}
};

double TSequencingErrorModel::getSteepestGradient(){
	if(_NRStepAccepted) return 0.0;
	double maxF = 0.0;
	for(unsigned int i=0; i<_covariates.numParameters; ++i){
		if(fabs(_F(i)) > maxF) maxF = fabs(_F(i));
	}
	return maxF;
};

void TSequencingErrorModel::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << _Jacobian << std::endl << std::endl;
};

void TSequencingErrorModel::printFToStdOut(){
	std::cout << std::endl << std::endl << "F:" << std::endl << _F << std::endl << std::endl;
};

void TSequencingErrorModel::printJxFToStdOut(){
	std::cout << std::endl << std::endl << "JxF:" << std::endl << _JxF << std::endl << std::endl;
};

}; //end namespace


