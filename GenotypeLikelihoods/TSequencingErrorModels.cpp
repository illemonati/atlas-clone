/*
 * TSequencingErrorModels.cpp
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#include "TSequencingErrorModels.h"

namespace GenotypeLikelihoods{


//-------------------------------------
// TSequencingErrorModelsOneReadGroup
//-------------------------------------
TSequencingErrorModelsOneReadGroup::TSequencingErrorModelsOneReadGroup(std::shared_ptr<TSequencingErrorModelNoRecal> & NoRecalModel){
	_modelFirstMate(NoRecalModel);
	_modelSecondMate(NoRecalModel);
};

void TSequencingErrorModelsOneReadGroup::addRecalModel(const TSequencingErrorModelDefinition & ModelDef, const bool & IsSecondMate){
	if(!IsSecondMate){
		_recalModelFirstMate = std::make_shared<TSequencingErrorModelRecal>(ModelDef);
		_modelFirstMate(_recalModelFirstMate);
	} else {
		_recalModelSecondMate = std::make_shared<TSequencingErrorModelRecal>(ModelDef);
		_modelSecondMate(_recalModelSecondMate);
	}
};

const TSequencingErrorModel& TSequencingErrorModelsOneReadGroup::operator[](const bool & IsSecondMate) const{
	if(!IsSecondMate){
		return *_modelFirstMate;
	} else {
		return *_modelSecondMate;
	}
};

const std::shared_ptr<TSequencingErrorModelRecal>& TSequencingErrorModelsOneReadGroup::getPointerToRecalModel(bool IsSecondMate) const{
	if(!IsSecondMate){
		return _recalModelFirstMate;
	} else {
		return _recalModelSecondMate;
	}
};

double TSequencingErrorModelsOneReadGroup::getErrorRate(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const {
	if(!base.isSecondMate()){
		return _modelFirstMate->getErrorRate(base, qualMap);
	} else {
		return _modelSecondMate->getErrorRate(base, qualMap);
	}
};

uint8_t TSequencingErrorModelsOneReadGroup::getPhredInt(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const{
	if(!base.isSecondMate()){
		return _modelFirstMate->getPhredInt(base, qualMap);
	} else {
		return _modelSecondMate->getPhredInt(base, qualMap);
	}
};

void TSequencingErrorModelsOneReadGroup::fillBaseLikelihoods(const BAM::TBase & base, const BAM::TQualityMap & qualMap, TBaseData & baseLikelihoods) const{
	if(!base.isSecondMate()){
		_modelFirstMate->fillBaseLikelihoods(base, qualMap, baseLikelihoods);
	} else {
		_modelSecondMate->fillBaseLikelihoods(base, qualMap, baseLikelihoods);
	}
};

//--------------------------------------------------------------------
// TSequencingErrorModels
//--------------------------------------------------------------------
TSequencingErrorModels::TSequencingErrorModels(){
	doRecalibration = false;
	_noRecalModel = std::make_shared<TSequencingErrorModelNoRecal>();
};

void TSequencingErrorModels::initialize(const std::string & RecalString, const std::string & RhoString, const BAM::TReadGroups & ReadGroups, TLog* Logfile){
	if(doRecalibration)
		throw std::runtime_error("void TSequencingErrorModels::initialize(const std::string & RecalString, const std::string & RhoString, const BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	//prepare objects
	models.resize(ReadGroups.size());

	//create model definition
	TSequencingErrorModelDefinition modelDef;
	std::string error;
	if(!modelDef.parse(RecalString, RhoString, error)){
		throw error + "!";
	}

	//initialize models
	for(auto& m : models){
		m.addRecalModel(modelDef, false);
		m.addRecalModel(modelDef, true);
	}
};

void TSequencingErrorModels::initializeFromFile(const std::string & Filename, const BAM::TReadGroups & ReadGroups, TLog* Logfile,
		                                        std::vector<uint16_t> & ReadGroupsWithoutRecal,
												std::vector<uint16_t> & ReadGroupsLikelySingleEnd){
	if(doRecalibration)
		throw std::runtime_error("void TSequencingErrorModels::initializeFromFile(const std::string & filename, const BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	//read parameters from file
	Logfile->listFlush("Initializing recalibration models from '" + Filename + "' ...");
	TInputFile in(Filename, {"readGroup", "mate", "covariates", "rho"}, "/t", "//");

	//prepare objects
	models.resize(ReadGroups.size(), TSequencingErrorModelsOneReadGroup(_noRecalModel));

	//tmp variables for reading
	std::vector<std::string> vec;
	TSequencingErrorModelDefinition modelDef;
	std::string error;

	//parse file to read details for each read group
	while(in.read(vec)){
		if(ReadGroups.readGroupExists(vec[0])){ //ignore if it does not exist
			//get read group
			uint16_t readGroupId = ReadGroups.getId(vec[0]);

			//parse model definition (allows us to get errors right here)
			if(!modelDef.parse(vec[2], vec[3], error)){
				throw error + " in file '" + Filename + "' on line " + toString(in.lineNumber()) + "!";
			}

			//add model
			if(vec[1] == "first")
				models[readGroupId].addRecalModel(modelDef, false);
			else if(vec[1] == "second")
				models[readGroupId].addRecalModel(modelDef, true);
			else
				throw "Unknown mate '" + vec[1] + "' in file '" + Filename + "' on line " + toString(in.lineNumber()) + "! Must be 'first' or 'second'.";
		}
	}
	Logfile->done();

	//check if some read groups do not have recal
	ReadGroupsWithoutRecal.clear();
	ReadGroupsLikelySingleEnd.clear();
	for(uint16_t r = 0; r < ReadGroups.size(); ++r){
		if(!models[r][0].recalibrates()){
			ReadGroupsWithoutRecal.push_back(r);
		} else {
			if(!models[r][1].recalibrates()){
				ReadGroupsLikelySingleEnd.push_back(r);
			}
		}
	}
};

//functions to get error rates
//-------------------------------------------------------
double TSequencingErrorModels::getErrorRate(const BAM::TBase & base) const{
	if(doRecalibration){
		return models[base.readGroupID].getErrorRate(base, qualMap);
	} else {
		return qualMap.phredIntToError(base.originalQuality_phredInt);
	}
};

uint8_t TSequencingErrorModels::getPhredInt(const BAM::TBase & base) const{
	if(doRecalibration){
		return models[base.readGroupID].getPhredInt(base, qualMap);
	} else {
		return base.originalQuality_phredInt;
	}
};

void TSequencingErrorModels::recalibrate(BAM::TBase & base) const{
	base.recalibratedQualityAsPhredInt = getPhredInt(base);
};

void TSequencingErrorModels::recalibrate(std::vector<BAM::TBase> & bases) const{
	for(auto& b : bases){
		recalibrate(b);
	}
};

void TSequencingErrorModels::fillBaseLikelihoods(const BAM::TBase & base, TBaseData & baseLikelihoods) const{
	if(doRecalibration){
		models[base.readGroupID].fillBaseLikelihoods(base, qualMap, baseLikelihoods);
	} else {
		_noRecalModel->fillBaseLikelihoods(base, qualMap, baseLikelihoods);
	}
};


// functions to write file
//-------------------------------------------------------------------
void TSequencingErrorModels::writeRecalFile(const std::string Filename, const BAM::TReadGroups & ReadGroups) const{
	//open file and write header
	TOutputFile out(Filename);
	out.writeHeader({"readGroup", "mate", "covariates", "rho"});

	//add models
	for(uint16_t r=0; r<ReadGroups.size(); ++r){
		_writeParameters(out, ReadGroups.getName(r), r, false);
		_writeParameters(out, ReadGroups.getName(r), r, true);
	}
};

DO I NEED THSI FUNCTION???

void TSequencingErrorModels::_writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate) const{
	out << readGroupName;
	if(isSecondMate) out << "second";
	else out << "first";
	out << models[readGroup][isSecondMate].getCovariateDefinition() << models[readGroup][isSecondMate].getRhoDefinition() << std::endl;
};

void TSequencingErrorModels::print() const{
	for(size_t r=0; r<readGroups->size(); ++r){
		int index = readGroupMap->getIndex(r);
		if(readGroupIndex.inUse(index, false)){
			auto& model = models[ readGroupIndex.index(index, false) ];
			std::cout << "Model rg = " << index << ", first: " << model.getCovariateDefinition() << "\t" << model.getRhoDefinition() << std::endl;

		}
		if(readGroupIndex.inUse(index, true)){
			auto& model = models[ readGroupIndex.index(index, true) ];
			std::cout << "Model rg = " << index << ", second: " << model.getCovariateDefinition() << "\t" << model.getRhoDefinition() << std::endl;
		}
	}
};



}; //end namespace
