/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include <TRecalibrationEMModel.h>

//--------------------------------------------------------------------
// TRecalibrationEMModelNEW
//--------------------------------------------------------------------

TRecalibrationEMModelNEW::TRecalibrationEMModelNEW(std::map<std::string, std::string> & modules, TLog* Logfile, bool verbose){
	//set parameters
	logfile = Logfile;

	//create modules







	//set number of parameters
	_numParameters = 0;

};

void TRecalibrationEMModelNEW::_createModule(TRecalibrationEMModule* module, std::string str, bool verbose){
	std::string orig = str;
	std::string format = "Expected format is NAME(PARAMETERS)[VALUES], where (PARAMETERS) and [VALUES] are optional.";

	//split string into parameters and values
	size_t pos = str.find('(');
	std::string name;
	std::vector<std::string> parameters;
	std::vector<std::string> values;

	if(pos != std::string::npos){
		//extract name
		name = str.substr(0, pos);
		str.erase(0, pos + 1);
		if(!stringContainsOnlyLetters(name)){
			throw "Wrong format for recal module '" + orig + "'! " + format;
		}

		//extract parameters
		size_t pos = str.find(')');
		if(pos == std::string::npos){
			throw "Wrong format for recal module '" + orig + "': missing ')'! " + format;
		}
		fillVectorFromStringAnySkipEmpty(str.substr(0, pos), parameters, ",");
		str.erase(0, pos + 1);
	}



		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
	std::string tmp =


}

double TRecalibrationEMModelNEW::createModules(std::map<std::string, std::string> & modules, bool verbose){
	for(std::map<std::string, std::string>::iterator it = modules.begin(); it != modules.end(); ++it){
		//which module?
		if(it->first == "quality"){

		}
	}





}

double TRecalibrationEMModelNEW::_calcEpsilon(const double eta){
	if(eta > 23.03){
		return 0.9999999999;
	}
	if(eta < -23.03){
		return 0.0000000001;
	}

	double tmp = exp(eta);
	return tmp / (1.0 + tmp);
};

double TRecalibrationEMModelNEW::getErrorRate(const TBase & base){
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as modules
	double eta = intercept.getEtaTerm();

	eta += quality->getEtaTerm(base.qualityAsPhredInt);
	eta += position->getEtaTerm(base.distFrom5Prime);
	eta += context->getEtaTerm(base.context);
	eta += fragmentLength->getEtaTerm(base.fragmentLength);

	return _calcEpsilon(eta);
};

double TRecalibrationEMModelNEW::getErrorRate(const TRecalibrationEMReadData & data){
	//eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as modules
	double eta = intercept.getEtaTerm();

	eta += quality->getEtaTerm(data.quality);
	eta += position->getEtaTerm(data.position);
	eta += context->getEtaTerm(data.context);
	eta += fragmentLength->getEtaTerm(data.fragmentLength);

	return _calcEpsilon(eta);
};

//-------------------------------------------------
//functions for estimation
void TRecalibrationEMModel_Base::setEMParamsToZero(){
	Jacobian.resize(_numParameters, _numParameters);
	F.resize(_numParameters);
	JxF.resize(_numParameters, 1);

	Jacobian.zeros();
	F.zeros();

	_numSitesAdded = 0;
	_NRconverged = false;
	_NRStepAccepted = false;
};

void TRecalibrationEMModelNEW::setQToZero(){
	if(!_NRconverged){
		_oldQ = _Q;
		_Q = 0.0;
	}
};

double TRecalibrationEMModelNEW::_calcQ(const int & genotype, TRecalibrationEMReadData & data){
	double eps = getErrorRate(data);
	double B = 1.33333333333333333333 * data.D[genotype] - 1.0;
	double P_d_given_g_beta = B * eps - data.D[genotype] + 1.0;
	return log(P_d_given_g_beta);
};

void TRecalibrationEMModelNEW::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	if(!_NRconverged){
		_Q += _calcQ(knownGenotype, data);
	}
};

void TRecalibrationEMModelNEW::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
	if(!_NRconverged){
		//add this data for all genotypes
		for(int g=0; g<4; ++g){
			_Q += P_g_given_d_oldBeta[g] * _calcQ(g, data);
		}
	}
};

bool TRecalibrationEMModelNEW::solveJxF(){
	if(_NRconverged){
		return true;
	} else {
		//Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
		for(unsigned int i=0; i<(_numParameters-1); ++i){
			for(unsigned int j=i+1; j<_numParameters; ++j){
				//copy from upper triangle to lower triangle
				Jacobian(j,i) = Jacobian(i,j);
			}
		}

		//scale F and J by 1/#sites
		Jacobian = Jacobian / (double) _numSitesAdded;
		F = F / (double) _numSitesAdded;

		//now solve J^-1 x F
		return solve(JxF, Jacobian, F);
	}
};

void TRecalibrationEMModelNEW::proposeNewParameters(double & lambda){
	if(!_NRStepAccepted){
		size_t index = 0;

		for(TRecalibrationEMModule* it : activeModules){
			it->proposeNewParameters(JxF, index, lambda);
		}
	}
};

bool TRecalibrationEMModelNEW::acceptProposedParametersBasedOnQ(){
	if(_NRStepAccepted) return true;
	if(_Q > _oldQ){
		_NRStepAccepted = true;
	} else {
		_NRStepAccepted = false;
		_Q = _oldQ;

		for(TRecalibrationEMModule* it : activeModules){
			it->rejectProposedParameters();
		}
	}
	return _NRStepAccepted;
};

double TRecalibrationEMModelNEW::getSteepestGradient(){
	if(_NRStepAccepted) return 0.0;
	double maxF = 0.0;
	for(unsigned int i=0; i<_numParameters; ++i){
		if(fabs(F(i)) > maxF) maxF = fabs(F(i));
	}
	return maxF;
};

void TRecalibrationEMModelNEW::printJacobianToStdOut(){
	std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << Jacobian << std::endl << std::endl;
};

void TRecalibrationEMModelNEW::printFToStdOut(){
	std::cout << std::endl << std::endl << "F:" << std::endl << F << std::endl << std::endl;
};

void TRecalibrationEMModelNEW::printJxFToStdOut(){
	std::cout << std::endl << std::endl << "JxF:" << std::endl << JxF << std::endl << std::endl;
};


//--------------------------------------------------------------------
// TRecalibrationEMModels
//--------------------------------------------------------------------
TRecalibrationEMModels::TRecalibrationEMModels(int numReadGroups, TLog* Logfile){
	totNumParameters = 0;
	readGroupIndex.initialize(numReadGroups);
	logfile = Logfile;
};

TRecalibrationEMModels::~TRecalibrationEMModels(){
	for(TRecalibrationEMModel_Base* it : models)
		delete it;
};

void TRecalibrationEMModels::_addModel(std::map<std::string, std::string> & modules, bool verbose){

	models.push_back(TRecalibrationEMModel(std::map<std::string, std::string> & modules, logfile, verbose));

	totNumParameters += models.back()->numParameters();
};

void TRecalibrationEMModels::addSameModelForAllReadGroups(std::string modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	std::pair<int, bool> missingReadGroupInfo;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo))
		addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, modelTag, values, false);
};

void TRecalibrationEMModels::addModel(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	//add to read group index
	readGroupIndex.setAsUsed(readGroupId, isSecondMate);

	//add model according to tag
	_addModel(modelTag, values, verbose);
};

void TRecalibrationEMModels::addModel(int readGroupId, bool isSecondMate, std::string modelTag, int maxPos){
	//add to read group index
	readGroupIndex.setAsUsed(readGroupId, isSecondMate);

	//create model
	models.push_back(createTRecalibrationEMModel(modelTag, maxPos, false, logfile));
	totNumParameters += models.back()->numParameters();
};

void TRecalibrationEMModels::addModelIfItDoesNotExist(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<int> & Qualities, int maxPos){
	if(readGroupIndex.inUse(readGroupId,isSecondMate)){
		//check model
		models[readGroupIndex.index(readGroupId, isSecondMate)]->checkParameterRange(Qualities, maxPos);
	} else {
		//add to read group index
		readGroupIndex.setAsUsed(readGroupId, isSecondMate);

		//create model
		models.push_back(createTRecalibrationEMModel(modelTag, maxPos, false, logfile));
		totNumParameters += models.back()->numParameters();
	}
};

void TRecalibrationEMModels::removeModel(int readGroupId, bool isSecondMate){
	//get index of model
	int index = readGroupIndex.index(readGroupId, isSecondMate);

	//adjust total number of parameters
	totNumParameters -= models[index]->numParameters();

	//erase
	models.erase(models.begin() + index);

	//update read group index
	readGroupIndex.setAsNotUsed(readGroupId, isSecondMate);
};

void TRecalibrationEMModels::createModels(std::string s, TReadGroups & readGroups, TReadGroupMap & readGroupMap){
	//initialize from string or file
	size_t pos = s.find_first_of('[');
	if(pos == std::string::npos){
		_createModelsFromFile(s, readGroups, readGroupMap);
	} else {
		_createModelsFromString(s, readGroups);
	}
};

void TRecalibrationEMModels::_createModelsFromString(std::string & s, TReadGroups & readGroups){
	//s has format model[param1, param2, param3, ...]
	logfile->startIndent("Initializing recal with string '" + s + "' for all read groups:");

	//read model tag
	size_t pos = s.find_first_of('[');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing '['!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::string modelTag = s.substr(0,pos);
	s.erase(0, pos+1);

	//read parameters: quality, position and context separted by semicolon (;)
	pos = s.find_first_of(']');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing ']'!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::vector<std::string> tmpVec;
	fillVectorFromString(s.substr(0, pos), tmpVec, ";");

	//if(tmpVec.size() != 3)
	//	throw "Failed to understand recal string: wrong number of parameter sets (" + toString(tmpVec.size()) + " instead of 3)!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";

	//initialize model
	addSameModelForAllReadGroups(modelTag, tmpVec, true);

	logfile->endIndent();
};

void TRecalibrationEMModels::_createModelsFromFile(std::string filename, TReadGroups & readGroups, TReadGroupMap & readGroupMap){
	//read parameters from file
	logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");

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
			if(vec.size() != 6)
				throw "Wrong number of entries in file '" + filename + "' on line " + toString(lineNum) + ": expected 6 but found " + toString(vec.size()) + "!";

			//check if read group exists
			if(readGroups.readGroupExists(vec[0])){
				//read read group, mate and model
				int rg = readGroupMap.getIndex(readGroups.find(vec[0]));

				bool isSecondMate;
				if(vec[1] == "second")
					isSecondMate = true;
				else if(vec[1] == "first")
					isSecondMate = false;
				else
					throw "Unknown mate '" + vec[1] + "' in file '" + filename + "' on line " + toString(lineNum) + "!";

				//if read group is pooled. And if so, only create model using the values of the first read group of the pool
				if(!modelExists(rg, isSecondMate)){
					//clean up vec to only contain parameters (remove read group, mate, model and LL)
					vec.erase(vec.begin(), vec.begin() + 2);

					//create model
					addModel(rg, isSecondMate, modelTag, vec, false);
				}
			} else {
				logfile->warning("Read group '" + vec[0] + "' does not exist in the BAM header! Are you using the correct recal file?");
			}
		}
	}
	logfile->done();

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	if(hasReadGroupsWithoutModel()){
		warningForMissingReadGroups(readGroups, readGroupMap);
		logfile->startIndent("Will assume the following read groups to be single end (no recalibration provided for second mate):");
		reportReadGroupsConsideredSingleEnd(readGroups, readGroupMap);
		addNoRecalModelIfMissing();
		logfile->endIndent();
	}
};

bool TRecalibrationEMModels::hasReadGroupsWithoutModel(){
	return readGroupIndex.hasCasesWithoutIndex();
};

void TRecalibrationEMModels::addNoRecalModelIfMissing(){
	//report read groups for which no recal model was given and initialize them as "no_recal" model
	std::pair<int, bool> missingReadGroupInfo;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo))
		addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, noRecal_name, 0);
};

void TRecalibrationEMModels::reportReadGroupsNotUsed(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	readGroupIndex.reportReadGroupsNotUsed(logfile, readGroups, ReadGroupMap);
};

void TRecalibrationEMModels::reportReadGroupsNotUsed(TReadGroups & readGroups){
	readGroupIndex.reportReadGroupsNotUsed(logfile, readGroups);
};

void TRecalibrationEMModels::reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	readGroupIndex.reportReadGroupsConsideredSingleEnd(logfile, readGroups, ReadGroupMap);
};

void TRecalibrationEMModels::reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups){
	readGroupIndex.reportReadGroupsConsideredSingleEnd(logfile, readGroups);
};

void TRecalibrationEMModels::warningForMissingReadGroups(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	readGroupIndex.warningForMissingReadGroups(logfile, readGroups, ReadGroupMap);
};

void TRecalibrationEMModels::warningForMissingReadGroups(TReadGroups & readGroups){
	readGroupIndex.warningForMissingReadGroups(logfile, readGroups);
};

void TRecalibrationEMModels::setEMParamsToZero(){
	for(TRecalibrationEMModel_Base* model : models){
		model->setEMParamsToZero();
	}
};

void TRecalibrationEMModels::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	models[ readGroupIndex.index(data) ]->addToFandJacobian(data, weightF, weightJacobian);
};

void TRecalibrationEMModels::setQToZero(){
	for(TRecalibrationEMModel_Base* model : models){
		model->setQToZero();
	}
};

void TRecalibrationEMModels::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
	models[ readGroupIndex.index(data) ]->addToQ(data, P_g_given_d_oldBeta);
};

void TRecalibrationEMModels::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	models[ readGroupIndex.index(data) ]->addToQ(data, knownGenotype);
};

double TRecalibrationEMModels::curQ(){
	double Q = 0.0;
	for(TRecalibrationEMModel_Base* model : models){
		Q += model->curQ();
	}
	return Q;
};

bool TRecalibrationEMModels::solveJxF(){
	bool couldSolve = true;
	for(TRecalibrationEMModel_Base* model : models){
		if(!model->solveJxF()){
			model->printJacobianToStdOut();
			throw "Issue solving JxF! This may be due to a lack of data. Consider adding more sites.";
			couldSolve = false;
		}
	}
	return couldSolve;
};

void TRecalibrationEMModels::proposeNewParameters(double lambda){
	for(TRecalibrationEMModel_Base* it : models){
		it->proposeNewParameters(lambda);
	}
};

unsigned int TRecalibrationEMModels::acceptProposedParametersBasedOnQ(){
	unsigned int numAccepted = 0;
	for(TRecalibrationEMModel_Base* model : models){
		numAccepted += model->acceptProposedParametersBasedOnQ();
	}
	return numAccepted;
};

void TRecalibrationEMModels::rejectProposedParameters(){
	for(TRecalibrationEMModel_Base* it : models)
		it->rejectProposedParameters();
};

double TRecalibrationEMModels::getSteepestGradient(){
	double maxF = 0.0;
	for(TRecalibrationEMModel_Base* model : models){
		double tmp = model->getSteepestGradient();
		if(fabs(tmp) > maxF) maxF = fabs(tmp);
	}
	return maxF;
};

void TRecalibrationEMModels::writeHeader(TOutputFilePlain & out){
	out.writeHeader({"readGroup", "mate", "model", "quality", "position", "context"});
};

void TRecalibrationEMModels::writeParameters(TOutputFilePlain & out, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap){
	for(size_t r=0; r<readGroups.size(); ++r){
		int index = ReadGroupMap[r];
		if(readGroupIndex.inUse(index, false)){
			_writeParameters(out, readGroups.getName(r), index, false);
		}
		if(readGroupIndex.inUse(index, true)){
			_writeParameters(out, readGroups.getName(r), index, true);
		}
	}
};

void TRecalibrationEMModels::_writeParameters(TOutputFilePlain & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate){
	if(readGroupIndex.inUse(readGroup, isSecondMate)){
		out << readGroupName;
		if(isSecondMate) out << "second";
		else out << "first";
		models[ readGroupIndex.index(readGroup, isSecondMate) ]->writeParametersToFile(out);
		out << std::endl;
	}
};

