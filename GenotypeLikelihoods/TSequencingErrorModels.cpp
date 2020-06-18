/*
 * TSequencingErrorModels.cpp
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#include "TSequencingErrorModels.h"

namespace GenotypeLikelihoods{


//--------------------------------------------------------------------
// TErrorModels
//--------------------------------------------------------------------
//Has model definitions, which are read from file or command line. For each model definition a model is created (taking into account RG pooling)
TSequencingErrorModels::TSequencingErrorModels(){
	readGroups = nullptr;
	readGroupMap = nullptr;
	totNumParameters = 0;
	logfile = nullptr;
	doRecalibration = false;
};

void TSequencingErrorModels::_init(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile){
	if(doRecalibration)
		throw "TSequencingErrorModels already initialized!";
	readGroups = ReadGroups;
	readGroupMap = ReadGroupMap;
	totNumParameters = 0;
	readGroupIndex.initialize(readGroups, readGroupMap);
	logfile = Logfile;
	doRecalibration = true;
};

//general function to add and remove models
//-----------------------------------------
void TSequencingErrorModels::_readRecalFile(const std::string filename, std::vector<TSequencingErrorModelDefinition> & modelDefs){
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

void TSequencingErrorModels::removeModel(int readGroupId, bool isSecondMate){
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
void TSequencingErrorModels::prepareModelsForEstimation(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile){
	_init(ReadGroups, ReadGroupMap, Logfile);
};

void TSequencingErrorModels::addModel(const uint16_t readGroupId, const bool isSecondMate, TSequencingErrorCovariateDefinition & covariates, TRecalibrationEMDataTable* dataTable){
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

void TSequencingErrorModels::addModelsFromFile(const std::string filename, TRecalibrationEMDataTables* dataTables){
	//read recal file
	std::vector<TSequencingErrorModelDefinition> modelDefs;
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
void TSequencingErrorModels::createModels(std::string s, BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile){
	_init(ReadGroups, ReadGroupMap, Logfile);

	//initialize from string or file
	size_t pos = s.find_first_of('=');
	if(pos == std::string::npos){
		_createModelsFromFile(s);
	} else {
		_createModelsFromString(s);
	}
};

void TSequencingErrorModels::createEmptyModels(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile){
	_init(ReadGroups, ReadGroupMap, Logfile);
	_addNoRecalModelIfMissing();
};

void TSequencingErrorModels::_addModel(const uint16_t readGroupId, const bool isSecondMate, TSequencingErrorCovariateDefinition & covariateFunctions){
	if(!readGroupIndex.inUse(readGroupId,isSecondMate)){
		readGroupIndex.setAsUsed(readGroupId, isSecondMate);
		models.emplace_back(covariateFunctions, logfile);
		totNumParameters += models.back().numParameters();
	}
};

void TSequencingErrorModels::_createModelsFromString(const std::string & s){
	//s has format model[param1, param2, param3, ...]
	logfile->startIndent("Initializing recal with string '" + s + "' for all read groups:");

	std::string error;
	TSequencingErrorCovariateDefinition modelDef(s, error);
	if(!error.empty()){
		throw error + "!";
	}

	//initialize same model for all read groups
	std::pair<int, bool> missingReadGroupInfo;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo))
		_addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, modelDef);

	logfile->endIndent();
};

void TSequencingErrorModels::_createModelsFromFile(std::string filename){
	//read recal file
	std::vector<TSequencingErrorModelDefinition> modelDefs;
	_readRecalFile(filename, modelDefs);

	//now create models
	logfile->listFlush("Creating recalibration models ...");
	for(auto& def: modelDefs){
		//if read group is pooled, only create model using the values of the first read group of the pool
		if(!modelExists(def)){
			_addModel(def.readGroupId, def.isSecondMate, def.covariates);
		}
	}
	logfile->done();
	logfile->conclude("Created a total of " + toString(modelDefs.size()) + " models.");

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	warningForMissingReadGroups();
	reportReadGroupsConsideredSingleEnd();
	_addNoRecalModelIfMissing();
};


// functions for reporting / writing
//-------------------------------------------------------------------
bool TSequencingErrorModels::hasReadGroupsWithoutModel() const{
	return readGroupIndex.hasCasesWithoutIndex();
};

void TSequencingErrorModels::_addNoRecalModelIfMissing(){
	//create no-recal model: only covariate is quality and beta is 1
	std::string error; //needed by TSequencingErrorCovariateDefinition to write errors to
	TSequencingErrorCovariateDefinition noRecal("quality=polynomial[1]", error);

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	std::pair<int, bool> missingReadGroupInfo;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo)){
		_addModel(missingReadGroupInfo.first, missingReadGroupInfo.second, noRecal);
	}
};

void TSequencingErrorModels::reportReadGroupsNotUsed() const{
	readGroupIndex.reportReadGroupsNotUsed(logfile);
};

void TSequencingErrorModels::reportReadGroupsConsideredSingleEnd() const{
	if(readGroupIndex.hasCasesWithoutSecondMate()){
		logfile->startIndent("Will assume the following read groups to be single end (no recalibration provided for second mate):");
		readGroupIndex.reportReadGroupsConsideredSingleEnd(logfile);
		logfile->endIndent();
	}
};

void TSequencingErrorModels::warningForMissingReadGroups() const{
	readGroupIndex.warningForMissingReadGroups(logfile);
};


//functions to get error rates
//-------------------------------------------------------

double TSequencingErrorModels::getErrorRate(const TRecalibrationEMReadData & data) const{
	if(doRecalibration){
		return models[ readGroupIndex.index(data) ].getErrorRate(data);
	} else {
		return qualMap.phredIntToError(data.qualityAsPhredInt);
	}
};

double TSequencingErrorModels::getErrorRate(const BAM::TBase & base) const{
	if(base.base == N){
		return 1.0;
	} else if(doRecalibration){
		return models[ readGroupIndex.index(base) ].getErrorRate(base);
	} else {
		return qualMap.phredIntToError(base.originalQuality_phredInt);
	}
};

uint8_t TSequencingErrorModels::getPhredInt(const BAM::TBase & base) const{
	if(base.base == N){
		return 0;
	} else if(doRecalibration){
		return qualMap.errorToPhredInt( models[ readGroupIndex.index(base) ].getErrorRate(base) );
	} else {
		return base.originalQuality_phredInt;
	}
};

void TSequencingErrorModels::recalibrate(BAM::TBase & base) const{
	if(base.base == N){
		base.recalibratedQualityAsPhredInt = 0;
	} else if(doRecalibration){
		base.recalibratedQualityAsPhredInt = qualMap.errorToPhredInt(models[ readGroupIndex.index(base) ].getErrorRate(base));
	} else {
		base.recalibratedQualityAsPhredInt = base.originalQuality_phredInt;
	}
};

void TSequencingErrorModels::recalibrate(std::vector<BAM::TBase> & bases, const uint16_t length) const{
	if(doRecalibration){
		TSequencingErrorModel& model = models[ readGroupIndex.index(bases[0]) ];
		for(uint16_t i=0; i<length; ++i){
			if(bases[i].base == N){
				bases[i].recalibratedQualityAsPhredInt = 0;
			} else {
				bases[i].recalibratedQualityAsPhredInt = qualMap.errorToPhredInt(model.getErrorRate(bases[i]));
			}
		}
	} else {
		for(uint16_t i=0; i<length; ++i){
			bases[i].recalibratedQualityAsPhredInt = bases[i].originalQuality_phredInt;
		}
	}
};

void TSequencingErrorModels::recalibrate(std::vector<BAM::TBase> & bases) const{
	if(doRecalibration){
		TSequencingErrorModel& model = models[ readGroupIndex.index(bases[0]) ];
		for(auto& b : bases){
			if(b.base == N){
				b.recalibratedQualityAsPhredInt = 0;
			} else {
				b.recalibratedQualityAsPhredInt = qualMap.errorToPhredInt(model.getErrorRate(b));
			}
		}
	} else {
		for(auto& b : bases){
			b.recalibratedQualityAsPhredInt = b.originalQuality_phredInt;
		}
	}
};

void TSequencingErrorModels::calculateBaseLikelihoods(const BAM::TBase & base, TBaseData & baseLikelihoods) const{
	if(base.base == N){
		baseLikelihoods.reset();
	} else if(doRecalibration){
		models[ readGroupIndex.index(base) ].fillBaseLikelihoods(base, baseLikelihoods);
	} else {
		defaultRho.fillBaseLikelihoods(base.base, qualMap.phredIntToError(base.originalQuality_phredInt), baseLikelihoods);
	}
};


// functions for estimation
//-------------------------------------------------------------------
void TSequencingErrorModels::setEMParamsToZero(){
	for(auto& model : models){
		model.setEMParamsToZero();
	}
};

void TSequencingErrorModels::addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian){
	models[ readGroupIndex.index(data) ].addToFandJacobian(data, weightF, weightJacobian);
};

void TSequencingErrorModels::setQToZero(){
	for(auto& model : models){
		model.setQToZero();
	}
};

void TSequencingErrorModels::addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta){
	models[ readGroupIndex.index(data) ].addToQ(data, P_g_given_d_oldBeta);
};

void TSequencingErrorModels::addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype){
	models[ readGroupIndex.index(data) ].addToQ(data, knownGenotype);
};

double TSequencingErrorModels::curQ(){
	double Q = 0.0;
	for(auto& model : models){
		Q += model.curQ();
	}
	return Q;
};

bool TSequencingErrorModels::solveJxF(){
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

void TSequencingErrorModels::proposeNewParameters(double lambda){
	for(auto& model : models){
		model.proposeNewParameters(lambda);
	}
};

unsigned int TSequencingErrorModels::acceptProposedParametersBasedOnQ(){
	unsigned int numAccepted = 0;
	for(auto& model : models){
		numAccepted += (unsigned int) model.acceptProposedParametersBasedOnQ();
	}
	return numAccepted;
};

void TSequencingErrorModels::adjustParametersPostEstimation(){
	for(auto& model : models){
		model.adjustParametersPostEstimation();
	}
};

double TSequencingErrorModels::getSteepestGradient(){
	double maxF = 0.0;
	for(auto& model : models){
		double tmp = model.getSteepestGradient();
		if(fabs(tmp) > maxF) maxF = fabs(tmp);
	}
	return maxF;
};


// functions to write file
//-------------------------------------------------------------------
void TSequencingErrorModels::writeRecalFile(const std::string filename) const{
	//open file and write header
	TOutputFile out(filename);
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

void TSequencingErrorModels::_writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate) const{
	if(readGroupIndex.inUse(readGroup, isSecondMate)){
		out << readGroupName;
		if(isSecondMate) out << "second";
		else out << "first";
		TSequencingErrorCovariateDefinition def = models[ readGroupIndex.index(readGroup, isSecondMate) ].getCovariateDefinition();
		out << def.getModelString() << std::endl;
	}
};

}; //end namespace
