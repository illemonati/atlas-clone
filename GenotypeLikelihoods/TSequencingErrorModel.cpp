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
	fillVectorFromStringAnySkipEmpty(modelString, tmp, ";");
	for(std::string s : tmp){
		size_t pos = s.find('=');
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

		modelString += ";" + it.covariate + "=" + it.function;
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
	std::vector<std::string> vec = {covariateMap.intercept};
	intercept.initialize(0, vec);


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
// TSequencingErrorRho
//--------------------------------------------------------------------
TSequencingErrorRho::TSequencingErrorRho(){
	for(int b=0; b<4; ++b){
		for(int a=0; a<4; ++a){
			if(a==b){
				rho[b][a] = 0.0;
			} else {
				rho[b][a] = 1.0 / 3.0;
			}
		}
	}
};

void TSequencingErrorRho::fillBaseLikelihoods(const Base base, const double epsilon, TBaseData & baseLikelihoods) const{
	baseLikelihoods[base] = 1.0 - epsilon;
	if(base == A){
		baseLikelihoods[C] = epsilon * rho[A][C];
		baseLikelihoods[G] = epsilon * rho[A][G];
		baseLikelihoods[T] = epsilon * rho[A][T];
	} else if(base == C){
		baseLikelihoods[A] = epsilon * rho[C][A];
		baseLikelihoods[G] = epsilon * rho[C][G];
		baseLikelihoods[T] = epsilon * rho[C][T];
	} else if(base == G){
		baseLikelihoods[A] = epsilon * rho[G][A];
		baseLikelihoods[C] = epsilon * rho[G][C];
		baseLikelihoods[T] = epsilon * rho[G][T];
	} else {
		baseLikelihoods[A] = epsilon * rho[T][A];
		baseLikelihoods[C] = epsilon * rho[T][C];
		baseLikelihoods[G] = epsilon * rho[T][G];
	}
};

//--------------------------------------------------------------------
// TSequencingErrorModel
//--------------------------------------------------------------------
TSequencingErrorModel::TSequencingErrorModel(TSequencingErrorCovariateDefinition & covariateMap, TLog* Logfile){
	logfile = Logfile;
	setEMParamsToZero();

	//create covariates
	_covariates.createCovariatesAndIntercept(covariateMap);
};

TSequencingErrorModel::TSequencingErrorModel(TSequencingErrorCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable, TLog* Logfile){
	logfile = Logfile;
	setEMParamsToZero();

	//create covariates
	_covariates.createCovariatesAndIntercept(covariateMap, dataTable);
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

double TSequencingErrorModel::getErrorRate(const TBase & base) const{
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = _covariates.intercept.getEtaTerm();

	for(const auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(base);
	}

	return _calcEpsilon(eta);
};

double TSequencingErrorModel::getErrorRate(const TRecalibrationEMReadData & data) const{
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = _covariates.intercept.getEtaTerm();

	for(auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(data);
	}

	return _calcEpsilon(eta);
};

void TSequencingErrorModel::fillBaseLikelihoods(const TBase & base, TBaseData & baseLikelihoods) const{
	//first calculate epsilon
	double eta = _covariates.intercept.getEtaTerm();
	for(const auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(base);
	}

	std::cout << "error = " <<  _calcEpsilon(eta) << std::endl;

	//then calculate base likelihoods
	rho.fillBaseLikelihoods(base.base, _calcEpsilon(eta), baseLikelihoods);
};

TSequencingErrorCovariateDefinition TSequencingErrorModel::getCovariateDefinition() const{
	return _covariates.getCovariateDefinition();
};

//-------------------------------------------------
//functions for estimation
void TSequencingErrorModel::setEMParamsToZero(){
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

void TSequencingErrorModel::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	if(!_NRconverged){
		double eps = getErrorRate(data);
		_Q += _calcQ(eps, knownGenotype, data);
	}
};

void TSequencingErrorModel::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
	if(!_NRconverged){
		//get error rate
		double eps = getErrorRate(data);

		//add this data for all genotypes
		_Q += P_g_given_d_oldBeta[0] * _calcQ(eps, A, data);
		_Q += P_g_given_d_oldBeta[1] * _calcQ(eps, C, data);
		_Q += P_g_given_d_oldBeta[2] * _calcQ(eps, G, data);
		_Q += P_g_given_d_oldBeta[3] * _calcQ(eps, T, data);
	}
};

void TSequencingErrorModel::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	//fill derivatives
	_firstDerivatives.restart();
	_secondDerivatives.restart();

	//fill derivatives of intercept
	_covariates.intercept.fillDerivatives(0.0, _firstDerivatives, _secondDerivatives);

	//fill derivatives of covariates
	for(const auto & cov : _covariates.covariates){
		cov->fillDerivatives(data, _firstDerivatives, _secondDerivatives);
	}

	//add first derivatives to F and Jacobian
	for(TRecalibrationEMFirstDerivativesIterator it = _firstDerivatives.begin(); it != _firstDerivatives.end(); ++it){
		//add to F
		_F(it->index) += weightF * it->derivative;

		//add to J
		for(TRecalibrationEMFirstDerivativesIterator it2 = it; it2 != _firstDerivatives.end(); ++it2){
			_Jacobian(it->index, it2->index) += weightJacobian * it->derivative * it2->derivative;
		}
	}

	//add second derivatives to Jacobian (happens to have the same weigth as F!)
	for(auto& it : _secondDerivatives){
		_Jacobian(it.index1, it.index2) += weightF * it.derivative;
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


