/*
 * TSequencingErrorModels.cpp
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#include "TSequencingErrorModels.h"

namespace GenotypeLikelihoods{


//-------------------------------------
// TSequencingErrorModelEntry / TSequencingErrorModelsOneReadGroup
//-------------------------------------
void TSequencingErrorModelEntry::addModel(const TSequencingErrorModelDefinition & ModelDef){
	_recalModel = std::make_shared<TSequencingErrorModelRecal>(ModelDef);
};

const TSequencingErrorModel& TSequencingErrorModelEntry::model() const{
	if(_recalModel){
		return *_recalModel;
	} else {
		return _noRecalModel;
	}
};

TSequencingErrorModelRecal* TSequencingErrorModelEntry::getPointerToRecalModel(){
	return _recalModel.get();
};

std::shared_ptr<TSequencingErrorModelRecal>& TSequencingErrorModelEntry::getSharedPointerToRecalModel(){
	return _recalModel;
};

Probability TSequencingErrorModelEntry::getErrorRate(const BAM::TSequencedBase & base) const {
	if(_recalModel){
		return _recalModel->getErrorRate(base);
	} else {
		return _noRecalModel.getErrorRate(base);
	}
};

BAM::PhredIntErrorRate TSequencingErrorModelEntry::getPhredInt(const BAM::TSequencedBase & base) const{
	if(_recalModel){
		return _recalModel->getPhredInt(base);
	} else {
		return _noRecalModel.getPhredInt(base);
	}
};

void TSequencingErrorModelEntry::fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const{
	if(_recalModel){
		_recalModel->fillBaseLikelihoods(base, baseLikelihoods);
	} else {
		_noRecalModel.fillBaseLikelihoods(base, baseLikelihoods);
	}
};

//TSequencingErrorModelsOneReadGroup
void TSequencingErrorModelsOneReadGroup::addRecalModel(const TSequencingErrorModelDefinition & ModelDef, const bool & IsSecondMate){
	_models[IsSecondMate].addModel(ModelDef);
};

const TSequencingErrorModel& TSequencingErrorModelsOneReadGroup::operator[](const bool & IsSecondMate) const{
	return _models[IsSecondMate].model();
};

TSequencingErrorModelRecal* TSequencingErrorModelsOneReadGroup::getPointerToRecalModel(const bool & IsSecondMate){
	return _models[IsSecondMate].getPointerToRecalModel();
};

std::shared_ptr<TSequencingErrorModelRecal>& TSequencingErrorModelsOneReadGroup::getSharedPointerToRecalModel(const bool & IsSecondMate){
	return _models[IsSecondMate].getSharedPointerToRecalModel();
};

Probability TSequencingErrorModelsOneReadGroup::getErrorRate(const BAM::TSequencedBase & base) const {
	return _models[base.isSecondMate()].getErrorRate(base);
};

BAM::PhredIntErrorRate TSequencingErrorModelsOneReadGroup::getPhredInt(const BAM::TSequencedBase & base) const{
	return _models[base.isSecondMate()].getPhredInt(base);
};

void TSequencingErrorModelsOneReadGroup::fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const{
	_models[base.isSecondMate()].fillBaseLikelihoods(base, baseLikelihoods);
};

//--------------------------------------------------------------------
// TSequencingErrorModels
//--------------------------------------------------------------------
TSequencingErrorModels::TSequencingErrorModels(){
	_doRecalibration = false;
};

void TSequencingErrorModels::initialize(const std::string & RecalString, const std::string & RhoString, const BAM::TReadGroups & ReadGroups, TLog* Logfile){
	if(_doRecalibration)
		throw std::runtime_error("void TSequencingErrorModels::initialize(const std::string & RecalString, const std::string & RhoString, const BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	//prepare objects
	_models.resize(ReadGroups.size());

	//create model definition
	TSequencingErrorModelDefinition modelDef;
	std::string error;
	if(!modelDef.parse(RecalString, RhoString, error)){
		throw error + "!";
	}

	//initialize models
	for(auto& m : _models){
		m.addRecalModel(modelDef, false);
		m.addRecalModel(modelDef, true);
	}
};

void TSequencingErrorModels::initializeFromFile(const std::string & Filename, const BAM::TReadGroups & ReadGroups, TLog* Logfile){
	if(_doRecalibration)
		throw std::runtime_error("void TSequencingErrorModels::initializeFromFile(const std::string & filename, const BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	//read parameters from file
	Logfile->listFlush("Initializing recalibration models from '" + Filename + "' ...");
	TInputFile in(Filename, {"readGroup", "mate", "covariates", "rho"}, "/t", "//");

	//prepare objects
	_models.resize(ReadGroups.size());

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
				_models[readGroupId].addRecalModel(modelDef, false);
			else if(vec[1] == "second")
				_models[readGroupId].addRecalModel(modelDef, true);
			else
				throw "Unknown mate '" + vec[1] + "' in file '" + Filename + "' on line " + toString(in.lineNumber()) + "! Must be 'first' or 'second'.";
		}
	}
	Logfile->done();
};

void TSequencingErrorModels::checkReadGroups(const BAM::TReadGroups & ReadGroups,
											 std::vector<uint16_t> & ReadGroupsWithoutRecal,
											 std::vector<uint16_t> & ReadGroupsLikelySingleEnd){
	ReadGroupsWithoutRecal.clear();
	ReadGroupsLikelySingleEnd.clear();
	for(uint16_t r = 0; r < ReadGroups.size(); ++r){
		if(!_models[r][0].recalibrates()){
			ReadGroupsWithoutRecal.push_back(r);
		} else {
			if(!_models[r][1].recalibrates()){
				ReadGroupsLikelySingleEnd.push_back(r);
			}
		}
	}
};

//functions to get error rates
//-------------------------------------------------------
Probability TSequencingErrorModels::getErrorRate(const BAM::TSequencedBase & base) const{
	if(_doRecalibration){
		return _models[base.readGroupID].getErrorRate(base);
	} else {
		return base.originalQuality_phredInt;
	}
};

BAM::PhredIntErrorRate TSequencingErrorModels::getPhredInt(const BAM::TSequencedBase & base) const{
	if(_doRecalibration){
		return _models[base.readGroupID].getPhredInt(base);
	} else {
		return base.originalQuality_phredInt;
	}
};

void TSequencingErrorModels::recalibrate(BAM::TSequencedBase & base) const{
	base.recalibratedQualityAsPhredInt = getPhredInt(base);
};

void TSequencingErrorModels::recalibrate(std::vector<BAM::TSequencedBase> & bases) const{
	for(auto& b : bases){
		recalibrate(b);
	}
};

void TSequencingErrorModels::fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const{
	if(_doRecalibration){
		_models[base.readGroupID].fillBaseLikelihoods(base, baseLikelihoods);
	} else {
		_noRecalModel.fillBaseLikelihoods(base, baseLikelihoods);
	}
};


// functions to write file
//-------------------------------------------------------------------
void TSequencingErrorModels::writeRecalFile(const BAM::TReadGroups & ReadGroups, const std::string Filename) const{
	//open file and write header
	TOutputFile out(Filename);
	out.writeHeader({"readGroup", "mate", "covariates", "rho"});

	//add models
	for(uint16_t r=0; r<ReadGroups.size(); ++r){
		for(uint8_t mate=0; mate < 2; ++mate){
			out << ReadGroups.getName(r);

			if(mate == 0){
				out << "first";
			} else {
				out << "second";
			}

			out << _models[r][mate].getCovariateDefinition() << _models[r][mate].getRhoDefinition() << std::endl;
		}
	}
};


}; //end namespace
