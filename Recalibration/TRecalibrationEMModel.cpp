/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include <TRecalibrationEMModel.h>

//--------------------------------------------------------------------
// TRecalibrationEMModelCovariateDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
bool TRecalibrationEMModelCovariateDefinition::parse(const std::string & modelString, std::string & error){
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

void TRecalibrationEMModelCovariateDefinition::setIntercept(const double Intercept){
	intercept = toString(Intercept);
};

void TRecalibrationEMModelCovariateDefinition::addCovariate(const std::string covariate, const std::string function){
	covariateFunctions.emplace_back(covariate, function);
};

std::string TRecalibrationEMModelCovariateDefinition::getModelString(){
	std::string modelString = "intercept[" + intercept + "]";
	for(auto& it : covariateFunctions){

		modelString += ";" + it.covariate + "=" + it.function;
	}
	return modelString;
};

//--------------------------------------------------------------------
// TRecalibrationEMModelCovariateList
//--------------------------------------------------------------------
TRecalibrationEMModelCovariateList::TRecalibrationEMModelCovariateList(){
	numParameters = intercept.numParameters();
};

TRecalibrationEMModelCovariateList::~TRecalibrationEMModelCovariateList(){
	_clear();
};

TRecalibrationEMModelCovariateList::TRecalibrationEMModelCovariateList(TRecalibrationEMModelCovariateList&& other){
	//copy from other
	covariates = std::move(other.covariates);
	pointerToCovariateFunctions = std::move(other.pointerToCovariateFunctions);
	intercept = std::move(other.intercept);
	numParameters = other.numParameters;

	//set intercept of other to zero
	other.intercept.setIntercept(0.0);
};

TRecalibrationEMModelCovariateList& TRecalibrationEMModelCovariateList::operator=(TRecalibrationEMModelCovariateList&& other){
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

void TRecalibrationEMModelCovariateList::_clear(){
	for(TRecalibrationEMCovariate* it : covariates){
		delete it;
	}
	covariates.clear();
	pointerToCovariateFunctions.clear();
	intercept.setIntercept(0.0);
};

void TRecalibrationEMModelCovariateList::_createCovariatesAndIntercept(TRecalibrationEMModelCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable){
	//include intercept
	std::vector<std::string> vec = {covariateMap.intercept};
	intercept.initialize(0, vec);

	//create covariates
	numParameters = intercept.numParameters();
	for(TRecalibrationEMModelCovariateDefinitionIterator it = covariateMap.begin(); it != covariateMap.end(); ++it){
		//create function for each covariate
		if(it->covariate == RecalCovariateName_none){
			continue;
		} else if(it->covariate == RecalCovariateName_quality){
			covariates.emplace_back(new TRecalibrationEMCovariate_quality(numParameters, it->function, dataTable));
		} else if(it->covariate == RecalCovariateName_position){
			covariates.emplace_back(new TRecalibrationEMCovariate_position(numParameters, it->function, dataTable));
		} else if(it->covariate == RecalCovariateName_context){
			covariates.emplace_back(new TRecalibrationEMCovariate_context(numParameters, it->function, dataTable));
		} else if(it->covariate == RecalCovariateName_fragmentLength){
			throw "Covariate " + (std::string) RecalCovariateName_fragmentLength + " is currently not implemented!";
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}

		//add new parameters
		numParameters += covariates.back()->numParameters();
	}

	//summarize
	_storePointersToCovariateFunctions();
};

void TRecalibrationEMModelCovariateList::_createCovariatesAndIntercept(TRecalibrationEMModelCovariateDefinition & covariateMap){
	//include intercept
	std::vector<std::string> vec = {covariateMap.intercept};
	intercept.initialize(0, vec);

	//create covariates
	numParameters = intercept.numParameters();
	for(TRecalibrationEMModelCovariateDefinitionIterator it = covariateMap.begin(); it != covariateMap.end(); ++it){
		//create function for each covariate
		if(it->covariate == RecalCovariateName_none){
			continue;
		} else if(it->covariate == RecalCovariateName_quality){
			covariates.emplace_back(new TRecalibrationEMCovariate_quality(numParameters, it->function));
		} else if(it->covariate == RecalCovariateName_position){
			covariates.emplace_back(new TRecalibrationEMCovariate_position(numParameters, it->function));
		} else if(it->covariate == RecalCovariateName_context){
			covariates.emplace_back(new TRecalibrationEMCovariate_context(numParameters, it->function));
		} else if(it->covariate == RecalCovariateName_fragmentLength){
			throw "Covariate " + (std::string) RecalCovariateName_fragmentLength + " is currently not implemented!";
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}

		//add new parameters
		numParameters += covariates.back()->numParameters();
	}

	//summarize
	_storePointersToCovariateFunctions();
};

void TRecalibrationEMModelCovariateList::_storePointersToCovariateFunctions(){
	//add intercept
	pointerToCovariateFunctions.emplace_back(&intercept);

	//add covariates
	for(auto & cov : covariates){
		//store function pointer
		pointerToCovariateFunctions.push_back(cov->getPointerToFunction());
	}
};

TRecalibrationEMModelCovariateDefinition TRecalibrationEMModelCovariateList::getCovariateDefinition(){
	TRecalibrationEMModelCovariateDefinition def;
	def.setIntercept(intercept.getIntercept());
	for(const auto & cov : covariates){
		def.addCovariate(cov->name(), cov->functionString());
	}
	return def;
}

//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
TRecalibrationEMModel::TRecalibrationEMModel(TRecalibrationEMModelCovariateDefinition & covariateMap, TLog* Logfile){
	logfile = Logfile;
	setEMParamsToZero();

	//create covariates
	_covariates._createCovariatesAndIntercept(covariateMap);
};

TRecalibrationEMModel::TRecalibrationEMModel(TRecalibrationEMModelCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable, TLog* Logfile){
	logfile = Logfile;
	setEMParamsToZero();

	//create covariates
	_covariates._createCovariatesAndIntercept(covariateMap, dataTable);
};

bool TRecalibrationEMModel::checkParameterRange(TRecalibrationEMDataTable* dataTable, std::string & error){
	for(auto & cov : _covariates.covariates){
		if(!cov->checkParameterRange(dataTable)){
			error = "Function for covariate " + cov->name() + " does not cover full range of data";
			return false;
		}
	}
	return true;
};

bool TRecalibrationEMModel::checkParameterRange(std::vector<uint16_t> & Qualities, uint16_t maxPos, std::string & error){
	for(auto & cov : _covariates.covariates){
		if(!cov->checkParameterRange(Qualities, maxPos)){
			error = "Function for covariate " + cov->name() + " does not cover full range of data";
			return false;
		}
	}
	return true;
};

void TRecalibrationEMModel::_initializeDerivatives(){
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

double TRecalibrationEMModel::_calcEpsilon(const double eta){
	if(eta > 23.03){
		return 0.9999999999;
	}
	if(eta < -23.03){
		return 0.0000000001;
	}

	return 1.0 / (1.0 + exp(-eta));
};

double TRecalibrationEMModel::getErrorRate(const TBase & base){
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = _covariates.intercept.getEtaTerm();

	for(const auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(base);
	}

	return _calcEpsilon(eta);
};

double TRecalibrationEMModel::getErrorRate(const TRecalibrationEMReadData & data){
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = _covariates.intercept.getEtaTerm();

	for(auto & cov : _covariates.covariates){
		eta += cov->getEtaTerm(data);
	}

	return _calcEpsilon(eta);
};

TRecalibrationEMModelCovariateDefinition TRecalibrationEMModel::getCovariateDefinition(){
	return _covariates.getCovariateDefinition();
};

//-------------------------------------------------
//functions for estimation
void TRecalibrationEMModel::setEMParamsToZero(){
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

void TRecalibrationEMModel::setQToZero(){
	if(!_NRconverged){
		_oldQ = _Q;
		_Q = 0.0;
	}
};

void TRecalibrationEMModel::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	if(!_NRconverged){
		double eps = getErrorRate(data);
		_Q += _calcQ(eps, knownGenotype, data);
	}
};

void TRecalibrationEMModel::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
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

void TRecalibrationEMModel::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
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

bool TRecalibrationEMModel::solveJxF(){
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

void TRecalibrationEMModel::proposeNewParameters(double & lambda){
	if(!_NRStepAccepted){
		uint16_t index = 0;
		for(const auto it : _covariates.pointerToCovariateFunctions){
			it->proposeNewParameters(_JxF, index, lambda);
		}
	}
};

bool TRecalibrationEMModel::acceptProposedParametersBasedOnQ(){
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

void TRecalibrationEMModel::adjustParametersPostEstimation(){
	for(const auto it : _covariates.pointerToCovariateFunctions){
		_covariates.intercept.addToIntercept(it->adjustParametersPostEstimation());
	}
};

double TRecalibrationEMModel::getSteepestGradient(){
	if(_NRStepAccepted) return 0.0;
	double maxF = 0.0;
	for(unsigned int i=0; i<_covariates.numParameters; ++i){
		if(fabs(_F(i)) > maxF) maxF = fabs(_F(i));
	}
	return maxF;
};

void TRecalibrationEMModel::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << _Jacobian << std::endl << std::endl;
};

void TRecalibrationEMModel::printFToStdOut(){
	std::cout << std::endl << std::endl << "F:" << std::endl << _F << std::endl << std::endl;
};

void TRecalibrationEMModel::printJxFToStdOut(){
	std::cout << std::endl << std::endl << "JxF:" << std::endl << _JxF << std::endl << std::endl;
};

//--------------------------------------------------------------------
// TRecalibrationEMModels
//--------------------------------------------------------------------
//Has model definitions, which are read from file or command line. For each model definition a model is created (taking into account RG pooling)
TRecalibrationEMModels::TRecalibrationEMModels(TReadGroups* ReadGroups,  TReadGroupMap* ReadGroupMap, TLog* Logfile){
	readGroups = ReadGroups;
	readGroupMap = ReadGroupMap;
	totNumParameters = 0;
	readGroupIndex.initialize(readGroups, readGroupMap);
	logfile = Logfile;
};

//general function to add and remove models
//-----------------------------------------
void TRecalibrationEMModels::_readRecalFile(const std::string filename, std::vector<TRecalibrationEMModelDefinition> & modelDefs){
	//read parameters from file
	logfile->listFlush("Reading recalibration models from '" + filename + "' ...");

	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "' for reading!";

	//skip header
	std::string line;
	std::getline(file, line);

	//tmp variables for reading
	int lineNum = 0;
	std::vector<std::string> vec;

	//parse file to read details for each read group
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 3)
				throw "Wrong number of entries in file '" + filename + "' on line " + toString(lineNum) + ": expected 3 (read group, mate, recal model) but found " + toString(vec.size()) + "!";

			//check if read group exists
			if(!readGroups->readGroupExists(vec[0])){
				logfile->warning("Read group '" + vec[0] + "' does not exist in the BAM header! Are you using the correct recal file?");
			} else {
				//get read group index (handles pooling)
				uint16_t rg = readGroupMap->getIndex(vec[0]);

				if(readGroups->readGroupInUse(rg)){
					bool isSecondMate;
					if(vec[1] == "second")
						isSecondMate = true;
					else if(vec[1] == "first")
						isSecondMate = false;
					else
						throw "Unknown mate '" + vec[1] + "' in file '" + filename + "' on line " + toString(lineNum) + "! Must be 'first' or 'second'.";

					//save
					modelDefs.emplace_back(rg, isSecondMate);
					std::string error;
					if(!modelDefs.back().parseModel(vec[2], error)){
						throw error + " in file '" + filename + "' on line " + toString(lineNum) + "!";
					}
				}
			}
		}
	}
	logfile->done();
};

void TRecalibrationEMModels::removeModel(int readGroupId, bool isSecondMate){
	//get index of model
	int index = readGroupIndex.index(readGroupId, isSecondMate);

	//adjust total number of parameters
	totNumParameters -= models[index].numParameters();

	//erase
	models.erase(models.begin() + index);

	//update read group index
	readGroupIndex.setAsNotUsed(readGroupId, isSecondMate);
};

// Functions to add models for estimation: dataTable is provided
//-------------------------------------------------------------
void TRecalibrationEMModels::addModel(const uint16_t readGroupId, const bool isSecondMate, TRecalibrationEMModelCovariateDefinition & covariates, TRecalibrationEMDataTable* dataTable){
	if(readGroupIndex.inUse(readGroupId,isSecondMate)){
		//check model
		std::string error;
		if(!models[readGroupIndex.index(readGroupId, isSecondMate)].checkParameterRange(dataTable, error)){
			//get names of matching read groups
			std::vector<std::string> rgNames;
			readGroupMap->fillNamesOfReadgroups(readGroupId, rgNames);

			//compile error
			error += " for ";
			if(isSecondMate){
				error += "second";
			} else {
				error += "first";
			}
			error += " mate of read group";
			if(rgNames.size() > 1){
				error += "s " + concatenateString(rgNames, ", ") + "!";
			} else {
				error += " " + rgNames[0] + "!";
			}
			throw error;
		}
	} else {
		readGroupIndex.setAsUsed(readGroupId, isSecondMate);

		models.emplace_back(covariates, dataTable, logfile);

		totNumParameters += models.back().numParameters();
	}
};

void TRecalibrationEMModels::addModelsFromFile(const std::string filename, TRecalibrationEMDataTables* dataTables){
	//read recal file
	std::vector<TRecalibrationEMModelDefinition> modelDefs;
	_readRecalFile(filename, modelDefs);

	//now create models
	logfile->listFlush("Creating recalibration models ...");
	for(auto& def: modelDefs){
		//if read group is pooled, only create model using the values of the first read group of the pool
		if(!modelExists(def)){
			addModel(def.readGroupId, def.isSecondMate, def.covariates, dataTables->getTable(def.readGroupId, def.isSecondMate));
		}
	}
	logfile->done();
	logfile->conclude("Created a total of " + toString(modelDefs.size()) + " models.");
};

// Functions to add models for recalibration: no dataTable is provided
//-------------------------------------------------------------------
void TRecalibrationEMModels::createModels(std::string s){
	//initialize from string or file
	size_t pos = s.find_first_of('=');
	if(pos == std::string::npos){
		_createModelsFromFile(s);
	} else {
		_createModelsFromString(s);
	}
};

void TRecalibrationEMModels::addModel(const uint16_t readGroupId, const bool isSecondMate, TRecalibrationEMModelCovariateDefinition & covariateFunctions){
	if(!readGroupIndex.inUse(readGroupId,isSecondMate)){
		readGroupIndex.setAsUsed(readGroupId, isSecondMate);
		models.emplace_back(covariateFunctions, logfile);
		totNumParameters += models.back().numParameters();
	}
};

void TRecalibrationEMModels::_createModelsFromString(const std::string & s){
	//s has format model[param1, param2, param3, ...]
	logfile->startIndent("Initializing recal with string '" + s + "' for all read groups:");

	std::string error;
	TRecalibrationEMModelCovariateDefinition modelDef(s, error);
	if(!error.empty()){
		throw error + "!";
	}

	//initialize same model for all read groups
	std::pair<int, bool> missingReadGroupInfo;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo))
		addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, modelDef);

	logfile->endIndent();
};

void TRecalibrationEMModels::_createModelsFromFile(std::string filename){
	//read recal file
	std::vector<TRecalibrationEMModelDefinition> modelDefs;
	_readRecalFile(filename, modelDefs);

	//now create models
	logfile->listFlush("Creating recalibration models ...");
	for(auto& def: modelDefs){
		//if read group is pooled, only create model using the values of the first read group of the pool
		if(!modelExists(def)){
			addModel(def.readGroupId, def.isSecondMate, def.covariates);
		}
	}
	logfile->done();
	logfile->conclude("Created a total of " + toString(modelDefs.size()) + " models.");

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	warningForMissingReadGroups();
	reportReadGroupsConsideredSingleEnd();
	addNoRecalModelIfMissing();
};


// functions for reporting / writing
//-------------------------------------------------------------------
bool TRecalibrationEMModels::hasReadGroupsWithoutModel(){
	return readGroupIndex.hasCasesWithoutIndex();
};

void TRecalibrationEMModels::addNoRecalModelIfMissing(){
	//report read groups for which no recal model was given and initialize them as "no_recal" model
	std::pair<int, bool> missingReadGroupInfo;
	TRecalibrationEMModelCovariateDefinition empty;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo)){
		addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, empty);
	}
};

void TRecalibrationEMModels::reportReadGroupsNotUsed(){
	readGroupIndex.reportReadGroupsNotUsed(logfile);
};

void TRecalibrationEMModels::reportReadGroupsConsideredSingleEnd(){
	if(readGroupIndex.hasCasesWithoutSecondMate()){
		logfile->startIndent("Will assume the following read groups to be single end (no recalibration provided for second mate):");
		readGroupIndex.reportReadGroupsConsideredSingleEnd(logfile);
		logfile->endIndent();
	}
};

void TRecalibrationEMModels::warningForMissingReadGroups(){
	readGroupIndex.warningForMissingReadGroups(logfile);
};

void TRecalibrationEMModels::writeRecalFile(const std::string filename){
	//open file and write header
	TOutputFilePlain out(filename);
	out.writeHeader({"ReadGroup", "Mate", "Model"});

	//add models
	for(size_t r=0; r<readGroups->size(); ++r){
		int index = readGroupMap->getIndex(r);
		if(readGroupIndex.inUse(index, false)){
			_writeParameters(out, readGroups->getName(r), index, false);
		}
		if(readGroupIndex.inUse(index, true)){
			_writeParameters(out, readGroups->getName(r), index, true);
		}
	}
};

void TRecalibrationEMModels::_writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate){
	if(readGroupIndex.inUse(readGroup, isSecondMate)){
		out << readGroupName;
		if(isSecondMate) out << "second";
		else out << "first";
		TRecalibrationEMModelCovariateDefinition def = models[ readGroupIndex.index(readGroup, isSecondMate) ].getCovariateDefinition();
		out << def.getModelString() << std::endl;
	}
};

// functions for estimation
//-------------------------------------------------------------------
void TRecalibrationEMModels::setEMParamsToZero(){
	for(auto& model : models){
		model.setEMParamsToZero();
	}
};

void TRecalibrationEMModels::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	models[ readGroupIndex.index(data) ].addToFandJacobian(data, weightF, weightJacobian);
};

void TRecalibrationEMModels::setQToZero(){
	for(auto& model : models){
		model.setQToZero();
	}
};

void TRecalibrationEMModels::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
	models[ readGroupIndex.index(data) ].addToQ(data, P_g_given_d_oldBeta);
};

void TRecalibrationEMModels::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	models[ readGroupIndex.index(data) ].addToQ(data, knownGenotype);
};

double TRecalibrationEMModels::curQ(){
	double Q = 0.0;
	for(auto& model : models){
		Q += model.curQ();
	}
	return Q;
};

bool TRecalibrationEMModels::solveJxF(){
	bool couldSolve = true;
	for(auto& model : models){
		if(!model.solveJxF()){
			model.printJacobianToStdOut();
			throw "Issue solving JxF! This may be due to a lack of data. Consider adding more sites.";
			couldSolve = false;
		}
	}
	return couldSolve;
};

void TRecalibrationEMModels::proposeNewParameters(double lambda){
	for(auto& model : models){
		model.proposeNewParameters(lambda);
	}
};

unsigned int TRecalibrationEMModels::acceptProposedParametersBasedOnQ(){
	unsigned int numAccepted = 0;
	for(auto& model : models){
		numAccepted += (unsigned int) model.acceptProposedParametersBasedOnQ();
	}
	return numAccepted;
};

void TRecalibrationEMModels::adjustParametersPostEstimation(){
	for(auto& model : models){
		model.adjustParametersPostEstimation();
	}
};

double TRecalibrationEMModels::getSteepestGradient(){
	double maxF = 0.0;
	for(auto& model : models){
		double tmp = model.getSteepestGradient();
		if(fabs(tmp) > maxF) maxF = fabs(tmp);
	}
	return maxF;
};




