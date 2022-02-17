/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TSequencingErrorModel.h"
#include "probability.h"

namespace GenotypeLikelihoods{
namespace SequencingError {
using coretools::Probability;

//*********************************************************
// TRecalibrationEMModelCovariateDefinition
// class to store model definition. Used when parsing files
//*********************************************************
void TCovariateDefinition::clear(){
	_intercept = "";
	_covariateFunctions.clear();
};

bool TCovariateDefinition::parse(const std::string & modelString, std::string & error){
	//make sure it is empty
	clear();

	//split string
	std::vector<std::string> tmp;
	coretools::str::fillContainerFromString(modelString, tmp, ';', true);

	//loop over entries
	for(std::string s : tmp){
		size_t pos = s.find(':');
		if(pos == std::string::npos){
			//intercept?
			size_t pos_1 = s.find('[');
			size_t pos_2 = s.find(']');
			if(s.substr(0, pos_1) != "intercept"){
				error = "Unable to understand model string '" + modelString + "': expecting an ':' or the function 'intercept'!";
				return false;
			} else {
				if(pos_1 == std::string::npos){
					if(pos_2 != std::string::npos){
						error = "Unable to understand model string '" + modelString + "': missing '['";
						return false;
					}
					_intercept = "0.0";
				} else {
					if(pos_2 == std::string::npos){
						error = "Unable to understand model string '" + modelString + "': missing ']'";
						return false;
					}
					_intercept = s.substr(pos_1+1, pos_2-pos_1-1);
				}
			}
		} else {
			_covariateFunctions.emplace_back(s.substr(0, pos), s.substr(pos+1));
		}
	}
	return true;
};

void TCovariateDefinition::setIntercept(const double Intercept){
	_intercept = coretools::str::toString(Intercept);
};

void TCovariateDefinition::addCovariate(const std::string covariate, const std::string function){
	_covariateFunctions.emplace_back(covariate, function);
};

std::string TCovariateDefinition::getModelString() const{
	std::string modelString = "intercept[" + _intercept + "]";
	for(auto& it : _covariateFunctions){

		modelString += ";" + it.covariate + ":" + it.function;
	}
	return modelString;
};

//*********************************************************
// TRecalibrationEMModelCovariateList
//*********************************************************
TCovariateList::TCovariateList(){
	numParameters = intercept.numParameters();
};

TCovariateList::~TCovariateList(){
	_clear();
};

TCovariateList::TCovariateList(TCovariateList&& other){
	//copy from other
	covariates = std::move(other.covariates);
	pointerToCovariateFunctions = std::move(other.pointerToCovariateFunctions);
	intercept = std::move(other.intercept);
	numParameters = other.numParameters;

	//set intercept of other to zero
	other.intercept.setIntercept(0.0);
};

TCovariateList& TCovariateList::operator=(TCovariateList&& other){
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

void TCovariateList::_clear(){
	for(auto* it : covariates){
		delete it;
	}
	covariates.clear();
	pointerToCovariateFunctions.clear();
	intercept.setIntercept(0.0);
};

void TCovariateList::createCovariatesAndIntercept(const TCovariateDefinition & covariateMap, const RecalEstimatorTools::TRecalDataTable & DataTable){
	//include intercept
	if(!covariateMap.intercept().empty()){
		std::vector<std::string> vec = {covariateMap.intercept()};
		intercept.initialize(0, vec);
	}

	//create covariates
	numParameters = intercept.numParameters();
	for(auto it = covariateMap.cbegin(); it != covariateMap.cend(); ++it){
		//create function for each covariate
		if(it->covariate == TCovariate::name){
			continue;
		} else if(it->covariate == TCovariate_quality::name){
			covariates.emplace_back(new TCovariate_quality(numParameters, it->function, DataTable));
		} else if(it->covariate == TCovariate_position::name){
			covariates.emplace_back(new TCovariate_position(numParameters, it->function, DataTable));
		} else if(it->covariate == TCovariate_context::name){
			covariates.emplace_back(new TCovariate_context(numParameters, it->function, DataTable));
		} else if(it->covariate == TCovariate_fragmentLength::name){
			covariates.emplace_back(new TCovariate_fragmentLength(numParameters, it->function, DataTable));
		} else if(it->covariate == TCovariate_mappingQuality::name){
			covariates.emplace_back(new TCovariate_mappingQuality(numParameters, it->function, DataTable));
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}

		//add new parameters
		numParameters += covariates.back()->numParameters();
	}

	//summarize
	_storePointersToCovariateFunctions();
};

void TCovariateList::createCovariatesAndIntercept(const TCovariateDefinition & covariateMap){
	//include intercept
	std::vector<std::string> vec = {covariateMap.intercept()};
	intercept.initialize(0, vec);

	//create covariates
	numParameters = intercept.numParameters();
	for(auto it = covariateMap.cbegin(); it != covariateMap.cend(); ++it){
		//create function for each covariate
		if(it->covariate == TCovariate::name){
			continue;
		} else if(it->covariate == TCovariate_quality::name){
			covariates.emplace_back(new TCovariate_quality(numParameters, it->function));
		} else if(it->covariate == TCovariate_position::name){
			covariates.emplace_back(new TCovariate_position(numParameters, it->function));
		} else if(it->covariate == TCovariate_context::name){
			covariates.emplace_back(new TCovariate_context(numParameters, it->function));
		} else if(it->covariate == TCovariate_fragmentLength::name){
			covariates.emplace_back(new TCovariate_fragmentLength(numParameters, it->function));
		} else if(it->covariate == TCovariate_mappingQuality::name){
			covariates.emplace_back(new TCovariate_mappingQuality(numParameters, it->function));
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}

		//add new parameters
		numParameters += covariates.back()->numParameters();
	}

	//summarize
	_storePointersToCovariateFunctions();
};

void TCovariateList::_storePointersToCovariateFunctions(){
	//add intercept
	pointerToCovariateFunctions.emplace_back(&intercept);

	//add covariates
	for(auto & cov : covariates){
		//store function pointer
		pointerToCovariateFunctions.push_back(cov->getPointerToFunction());
	}
};

TCovariateDefinition TCovariateList::getCovariateDefinition() const{
	TCovariateDefinition def;
	def.setIntercept(intercept.getIntercept());
	for(const auto & cov : covariates){
		def.addCovariate(cov->typeString(), cov->functionString());
	}
	return def;
}

//*********************************************************
// TRho
//*********************************************************

void TRho::reset() noexcept {
	for (int a = 0; a < 4; ++a) {
		for (int b = 0; b < 4; ++b) {
			if(a==b){
				rho[a][b] = 0.0;
			} else {
				rho[a][b] = 1.0 / 3.0;
			}
		}
	}
}

bool TRho::set(const std::string & def, std::string & error){
	using coretools::str::toString;
	//"default" implies default rho
	if(def == "default"){
		reset();
		return true;
	}

	//otherwise: full matrix is provided
	std::vector<std::string> vec;
	std::string s = def;
	coretools::str::fillContainerFromString(def, vec, ';');
	if(vec.size() != 4){
		error = "Rho matrix has " + toString(vec.size()) + " instead of 4 rows!";
		return false;
	}

	//parse rows
	std::vector<double> r;
	for(size_t a = 0; a<vec.size(); ++a){
		std::string& row = vec[a];
		coretools::str::trimString(row, "()");
		coretools::str::fillContainerFromString(row, r, ',');
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
	}
	return true;
}

std::string TRho::getDefinition() const noexcept{
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
				def += coretools::str::toString(rho[a][b]);
			} else {
				def += "-";
			}
		}
	}
	return def;
}

void TRho::fillBaseLikelihoods(const genometools::Base base, const Probability & epsilon, TBaseLikelihoods & baseLikelihoods) const noexcept {
	using genometools::Base;
	using genometools::index;
	if(base == Base::N){
		baseLikelihoods.reset();
	} else {
		baseLikelihoods[base] = epsilon.complement();
		if(base == Base::A){
			baseLikelihoods[Base::C] = epsilon * rho[index(Base::C)][index(Base::A)];
			baseLikelihoods[Base::G] = epsilon * rho[index(Base::G)][index(Base::A)];
			baseLikelihoods[Base::T] = epsilon * rho[index(Base::T)][index(Base::A)];
		} else if(base == Base::C){
			baseLikelihoods[Base::A] = epsilon * rho[index(Base::A)][index(Base::C)];
			baseLikelihoods[Base::G] = epsilon * rho[index(Base::G)][index(Base::C)];
			baseLikelihoods[Base::T] = epsilon * rho[index(Base::T)][index(Base::C)];
		} else if(base == Base::G){
			baseLikelihoods[Base::A] = epsilon * rho[index(Base::A)][index(Base::G)];
			baseLikelihoods[Base::C] = epsilon * rho[index(Base::C)][index(Base::G)];
			baseLikelihoods[Base::T] = epsilon * rho[index(Base::T)][index(Base::G)];
		} else {
			baseLikelihoods[Base::A] = epsilon * rho[index(Base::A)][index(Base::T)];
			baseLikelihoods[Base::C] = epsilon * rho[index(Base::C)][index(Base::T)];
			baseLikelihoods[Base::G] = epsilon * rho[index(Base::G)][index(Base::T)];
		}
	}
}

void TRho::prepareEstimationFromEMWeights() noexcept {
	rho.fill({0., 0., 0., 0.});
}

void TRho::addBaseForEstimation(const genometools::Base & base, const TBaseLikelihoods & EMWeights) noexcept {
	using genometools::Base;
	using genometools::index;
	if(base == Base::A){
		rho[index(Base::C)][index(Base::A)] += EMWeights[Base::C].get();
		rho[index(Base::G)][index(Base::A)] += EMWeights[Base::G].get();
		rho[index(Base::T)][index(Base::A)] += EMWeights[Base::T].get();
	} else if(base == Base::C){
		rho[index(Base::A)][index(Base::C)] += EMWeights[Base::A].get();
		rho[index(Base::G)][index(Base::C)] += EMWeights[Base::G].get();
		rho[index(Base::T)][index(Base::C)] += EMWeights[Base::T].get();
	} else if(base == Base::G){
		rho[index(Base::A)][index(Base::G)] += EMWeights[Base::A].get();
		rho[index(Base::C)][index(Base::G)] += EMWeights[Base::C].get();
		rho[index(Base::T)][index(Base::G)] += EMWeights[Base::T].get();
	} else if(base == Base::N){
		rho[index(Base::A)][index(Base::T)] += EMWeights[Base::A].get();
		rho[index(Base::C)][index(Base::T)] += EMWeights[Base::C].get();
		rho[index(Base::G)][index(Base::T)] += EMWeights[Base::G].get();
	}
};

void TRho::estimate() noexcept{
	//calculate denominators
	std::array<double, 4> denom;
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
}

//*********************************************************
// TModelNoRecal
//*********************************************************
Probability TModelNoRecal::getErrorRate(const BAM::TSequencedBase & base) const{
	if(base == genometools::Base::N){
		return Probability(1.0);
	} else {
		return (Probability) base.originalQuality_phredInt;
	}
}

genometools::PhredIntProbability TModelNoRecal::getPhredInt(const BAM::TSequencedBase & base) const{
	if(base == genometools::Base::N){
		return genometools::PhredIntProbability(0); //Todo: change to maxProb() one available.
	} else {
		return base.originalQuality_phredInt;
	}
}

void TModelNoRecal::fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const{
	using genometools::Base;
	if(base == Base::N){
		baseLikelihoods.reset();
	} else {
		const auto eps = static_cast<Probability>(base.originalQuality_phredInt);
		for (auto other = Base::min; other < Base::max; ++other) {
			if (other == base.base) baseLikelihoods[other] = eps.complement();
			else baseLikelihoods[other] = (1./3)*eps;
		}
	}
}

//*********************************************************
// TModelRecal
//*********************************************************

TModelRecal::TModelRecal(const TModelDefinition & modelDef){
	//create covariates
	_covariates.createCovariatesAndIntercept(modelDef.covariates);
	_rho = modelDef.rho;

	//prepare Newton-Raphson variables
	setNewtonRaphsonParamsToZero();
};

TModelRecal::TModelRecal(const TModelDefinition & ModelDef, const RecalEstimatorTools::TRecalDataTable & DataTable){
	//create covariates
	_covariates.createCovariatesAndIntercept(ModelDef.covariates, DataTable);
	_rho = ModelDef.rho;

	//prepare Newton-Raphson variables
	setNewtonRaphsonParamsToZero();
};

TModelDefinition TModelRecal::getModelDefinition() const{
	std::string error;
	TModelDefinition modelDef(getCovariateDefinition(), getRhoDefinition(), error);
	return modelDef;
};

//-------------------------------------------------
//functions to calculate error rates
//-------------------------------------------------

Probability TModelRecal::_calcEpsilon(double eta) const{
	if(eta > 23.03){
		return Probability(0.9999999999);
	}
	if(eta < -23.03){
		return Probability(0.0000000001);
	}

	return coretools::logistic(eta);
};

Probability TModelRecal::_calcErrorRate(const BAM::TSequencedBase & base) const{
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = _covariates.intercept.getEtaTerm();

	for(const auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(base);
	}

	return _calcEpsilon(eta);
};

Probability TModelRecal::getErrorRate(const BAM::TSequencedBase & base) const{
	if(base == genometools::Base::N){
		return Probability::highest();
	} else {
		return _calcErrorRate(base);
	}
};

genometools::PhredIntProbability TModelRecal::getPhredInt(const BAM::TSequencedBase & base) const{
	if(base == genometools::Base::N){
		return genometools::PhredIntProbability(0); //Todo: change to maxProb() one available.
	} else {
		return genometools::PhredIntProbability(_calcErrorRate(base));
	}
};

void TModelRecal::fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const{
	if(base == genometools::Base::N){
		baseLikelihoods.reset();
	} else {
		_rho.fillBaseLikelihoods(base.base, _calcErrorRate(base), baseLikelihoods);
	}
};

//-------------------------------------------------
//functions to estimate rho
//-------------------------------------------------
void TModelRecal::prepareRhoEstimationFromEMWeights(){
	_rho.prepareEstimationFromEMWeights();
};

void TModelRecal::addBaseForRhoEstimation(BAM::TSequencedBase & base, const TBaseLikelihoods & EMWeights){
	_rho.addBaseForEstimation(base.base, EMWeights);
};

void TModelRecal::estimateRho(){
	_rho.estimate();
};

//-------------------------------------------------
//functions for estimation
//-------------------------------------------------
bool TModelRecal::checkParameterRange(RecalEstimatorTools::TRecalDataTable & DataTable, std::string & error){
	for(auto & cov : _covariates.covariates){
		if(!cov->checkParameterRange(DataTable)){
			error = "Function for covariate " + cov->typeString() + " does not cover full range of data";
			return false;
		}
	}
	return true;
};

void TModelRecal::_initializeDerivatives(){
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

void TModelRecal::setNewtonRaphsonParamsToZero(){
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

void TModelRecal::setQToZero(){
	if(!_NRconverged){
		_oldQ = _Q;
		_Q = 0.0;
	}
};

void TModelRecal::addToQ(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d){
	if(!_NRconverged){
		//get error rate
		Probability eps = getErrorRate(base);
		//calculate sum_bbar [ Ind(bbar=d)log(1-eps) + Ind(bbar!=d)log(eps) ]
		_Q += EM_weights_bbar_given_d[ base.base ].get() * log(eps.complement().get())
			+ EM_weights_bbar_given_d[ base.base ].complement().get() * log(eps.get());
	}
};

void TModelRecal::addToFandJacobian(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d){
	// 1) Calculate epsilon
	//--------------------
	double eps = getErrorRate(base).get();

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
	double weight1 = 1.0 - eps - EM_weights_bbar_given_d[base.base].get();

	std::cout << "weight1 = " << weight1 << ", EM_weight = " << EM_weights_bbar_given_d[base.base] << ", eps = " << eps << std::endl;

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

bool TModelRecal::solveJxF(){
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

void TModelRecal::proposeNewParameters(double & lambda){
	if(!_NRStepAccepted){
		uint16_t index = 0;
		for(const auto it : _covariates.pointerToCovariateFunctions){
			it->proposeNewParameters(_JxF, index, lambda);
		}
	}
};

bool TModelRecal::acceptProposedParametersBasedOnQ(){
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

void TModelRecal::adjustParametersPostEstimation(){
	for(const auto it : _covariates.pointerToCovariateFunctions){
		_covariates.intercept.addToIntercept(it->adjustParametersPostEstimation());
	}
};

double TModelRecal::getSteepestGradient(){
	if(_NRStepAccepted) return 0.0;
	double maxF = 0.0;
	for(unsigned int i=0; i<_covariates.numParameters; ++i){
		if(fabs(_F(i)) > maxF) maxF = fabs(_F(i));
	}
	return maxF;
};

void TModelRecal::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << _Jacobian << std::endl << std::endl;
};

void TModelRecal::printFToStdOut(){
	std::cout << std::endl << std::endl << "F:" << std::endl << _F << std::endl << std::endl;
};

void TModelRecal::printJxFToStdOut(){
	std::cout << std::endl << std::endl << "JxF:" << std::endl << _JxF << std::endl << std::endl;
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods
