/*
 * TModels.cpp
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#include "TSequencingErrorModels.h"
#include "mathFunctions.h"
#include "probability.h"
#include "stringFunctions.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using coretools::Probability;
using coretools::str::toString;

//-------------------------------------
// TModelEntry / TModelsOneReadGroup
//-------------------------------------
TModelNoRecal TModelEntry::_noRecalModel;

void TModelEntry::addModel(const TModelDefinition &ModelDef) {
	_recalModel = std::make_shared<TModelRecal>(ModelDef);
};

const TModel &TModelEntry::model() const {
	if (_recalModel) return *_recalModel;
	return _noRecalModel;
};

TModelRecal *TModelEntry::getPointerToRecalModel() { return _recalModel.get(); }

std::shared_ptr<TModelRecal> &TModelEntry::getSharedPointerToRecalModel() {
	return _recalModel;
}

Probability TModelEntry::getErrorRate(const BAM::TSequencedBase &base) const {
	if (_recalModel) return _recalModel->getErrorRate(base);
	return _noRecalModel.getErrorRate(base);
}

genometools::PhredIntProbability TModelEntry::getPhredInt(const BAM::TSequencedBase &base) const {
	if (_recalModel) return _recalModel->getPhredInt(base);
	return _noRecalModel.getPhredInt(base);
}

void TModelEntry::fillBaseLikelihoods(const BAM::TSequencedBase &base,
						     TBaseLikelihoods &baseLikelihoods) const {
	if (_recalModel) {
		_recalModel->fillBaseLikelihoods(base, baseLikelihoods);
	} else {
		_noRecalModel.fillBaseLikelihoods(base, baseLikelihoods);
	}
}

// TModelsOneReadGroup
void TModelsOneReadGroup::addRecalModel(const TModelDefinition &ModelDef,
						       const bool &IsSecondMate) {
	_models[IsSecondMate].addModel(ModelDef);
}

const TModel &TModelsOneReadGroup::operator[](const bool &IsSecondMate) const {
	return _models[IsSecondMate].model();
}

TModelRecal *TModelsOneReadGroup::getPointerToRecalModel(const bool &IsSecondMate) {
	return _models[IsSecondMate].getPointerToRecalModel();
}

std::shared_ptr<TModelRecal> &
TModelsOneReadGroup::getSharedPointerToRecalModel(const bool &IsSecondMate) {
	return _models[IsSecondMate].getSharedPointerToRecalModel();
}

bool TModelsOneReadGroup::recalibrates() const {
	return _models[0].recalibrates() || _models[1].recalibrates();
}

Probability TModelsOneReadGroup::getErrorRate(const BAM::TSequencedBase &base) const {
	return _models[base.isSecondMate()].getErrorRate(base);
}

genometools::PhredIntProbability
TModelsOneReadGroup::getPhredInt(const BAM::TSequencedBase &base) const {
	return _models[base.isSecondMate()].getPhredInt(base);
}

void TModelsOneReadGroup::fillBaseLikelihoods(const BAM::TSequencedBase &base,
							     TBaseLikelihoods &baseLikelihoods) const {
	_models[base.isSecondMate()].fillBaseLikelihoods(base, baseLikelihoods);
}

//--------------------------------------------------------------------
// TModels
//--------------------------------------------------------------------
bool TModels::recalStringIsLikelyAModel(const std::string &RecalString) {
	// check if it contains a ';', ':', '[' or ']'
	return coretools::str::stringContainsAny(RecalString, ";:[]");
}

void TModels::initialize(const std::string &RecalString, const std::string &RhoString,
					const BAM::TReadGroups &ReadGroups, coretools::TLog *) {
	if (!_models.empty())
		throw std::runtime_error(
			"void TModels::initialize(const std::string & RecalString, const std::string & RhoString, "
		    "const BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	// prepare objects
	_models.resize(ReadGroups.size());

	// create model definition
	TModelDefinition modelDef;
	std::string error;
	if (!modelDef.parse(RecalString, RhoString, error)) throw error + "!";

	// initialize models
	for (auto &m : _models) {
		m.addRecalModel(modelDef, false);
		m.addRecalModel(modelDef, true);
	}
}

void TModels::initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups,
						coretools::TLog *Logfile) {
	if (!_models.empty())
		throw std::runtime_error("void TModels::initializeFromFile(const std::string & filename, const "
					 "BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	// read parameters from file
	Logfile->listFlush("Initializing recalibration models from '" + Filename + "' ...");
	coretools::TInputFile in(Filename, {"readGroup", "mate", "covariates", "rho"}, "/t", "//");

	// prepare objects
	_models.resize(ReadGroups.size());

	// tmp variables for reading
	std::vector<std::string> vec;

	// parse file to read details for each read group
	while (in.read(vec)) {
		if (ReadGroups.readGroupExists(vec[0])) { // ignore if it does not exist
			// get read group
			const uint16_t readGroupId = ReadGroups.getId(vec[0]);
			std::string error;
			TModelDefinition modelDef;

			// parse model definition (allows us to get errors right here)
			if (!modelDef.parse(vec[2], vec[3], error)) 
				throw error + " in file '" + Filename + "' on line " + toString(in.lineNumber()) + "!";

			// add model
			if (vec[1] == "first")
				_models[readGroupId].addRecalModel(modelDef, false);
			else if (vec[1] == "second")
				_models[readGroupId].addRecalModel(modelDef, true);
			else
				throw "Unknown mate '" + vec[1] + "' in file '" + Filename + "' on line " + toString(in.lineNumber()) +
					"! Must be 'first' or 'second'.";
		}
	}
	Logfile->done();
}

void TModels::checkReadGroups(const BAM::TReadGroups &ReadGroups,
					     std::vector<uint16_t> &ReadGroupsWithoutRecal,
					     std::vector<uint16_t> &ReadGroupsLikelySingleEnd) {
	ReadGroupsWithoutRecal.clear();
	ReadGroupsLikelySingleEnd.clear();
	for (uint16_t r = 0; r < ReadGroups.size(); ++r) {
		if (!_models[r][0].recalibrates()) {
			ReadGroupsWithoutRecal.push_back(r);
		} else if (!_models[r][1].recalibrates()) {
			ReadGroupsLikelySingleEnd.push_back(r);
		}
	}
}

// functions to get error rates
//-------------------------------------------------------
Probability TModels::getErrorRate(const BAM::TSequencedBase &base) const {
	if (!_models.empty()) return _models[base.readGroupID].getErrorRate(base);
	return (Probability)base.originalQuality_phredInt;
}

genometools::PhredIntProbability TModels::getPhredInt(const BAM::TSequencedBase &base) const {
	if (!_models.empty()) return _models[base.readGroupID].getPhredInt(base);
	return base.originalQuality_phredInt;
}

void TModels::recalibrate(BAM::TSequencedBase &base) const {
	base.recalibratedQualityAsPhredInt = getPhredInt(base);
}

void TModels::recalibrate(std::vector<BAM::TSequencedBase> &bases) const {
	for (auto &b : bases) recalibrate(b);
}

void TModels::fillBaseLikelihoods(const BAM::TSequencedBase &base,
						 TBaseLikelihoods &baseLikelihoods) const {
	if (!_models.empty()) _models[base.readGroupID].fillBaseLikelihoods(base, baseLikelihoods);
}

// functions to write file
//-------------------------------------------------------------------
void TModels::writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const {
	// open file and write header
	coretools::TOutputFile out(Filename);
	out.writeHeader({"readGroup", "mate", "covariates", "rho"});

	// add models
	for (uint16_t r = 0; r < ReadGroups.size(); ++r) {
		for (uint8_t mate = 0; mate < 2; ++mate) {
			const auto s = mate == 0 ? "first" : "second";
			out << ReadGroups.getName(r) << s << _models[r][mate].getCovariateDefinition()
				<< _models[r][mate].getRhoDefinition() << std::endl;
		}
	}
}

} // namespace SequencingError
}; // namespace GenotypeLikelihoods
