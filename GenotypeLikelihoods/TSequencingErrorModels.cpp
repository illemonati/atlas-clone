/*
 * TModels.cpp
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#include "TSequencingErrorModels.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <ostream>
#include <stdexcept>

#include "TFile.h"
#include "TLog.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "TSequencingErrorModel.h"
#include "probability.h"
#include "stringFunctions.h"

namespace GenotypeLikelihoods { class TBaseLikelihoods; }

namespace GenotypeLikelihoods {
namespace SequencingError {

using coretools::Probability;
using coretools::str::toString;
using coretools::instances::logfile;

//--------------------------------------------------------------------
// TModels
//--------------------------------------------------------------------
bool TModels::recalStringIsLikelyAModel(const std::string &RecalString) const noexcept {
	// check if it contains a ';', ':', '[' or ']'
	return coretools::str::stringContainsAny(RecalString, ";:[]");
}

void TModels::initialize(const std::string &RecalString, const std::string &RhoString,
					const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		throw std::runtime_error(
			"void TModels::initialize(const std::string & RecalString, const std::string & RhoString, "
		    "const BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	// prepare objects
	_models.resize(ReadGroups.size());

	// create model definition
	TModelDefinition modelDef(RecalString, RhoString);

	// initialize models
	for (auto &m : _models) {
		m[0] = std::make_shared<TModelRecal>(modelDef);
		m[1] = std::make_shared<TModelRecal>(modelDef);
	}
}

void TModels::initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		throw std::runtime_error("void TModels::initializeFromFile(const std::string & filename, const "
					 "BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	// read parameters from file
	logfile().listFlush("Initializing recalibration models from '" + Filename + "' ...");
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
			try {
				const TModelDefinition modelDef(vec[2], vec[3]);

				// add model
				if (vec[1] == "first")
					_models[readGroupId][0] = std::make_shared<TModelRecal>(modelDef);
				else if (vec[1] == "second")
					_models[readGroupId][1] = std::make_shared<TModelRecal>(modelDef);
				else
					throw "Unknown mate '" + vec[1] + "! Must be 'first' or 'second'.";
			} catch (const char *error) {
				throw std::string(error) + " in file '" + Filename + "' on line " + toString(in.lineNumber()) + "!";
			}
		}
	}
	logfile().done();
}

void TModels::checkReadGroups(const BAM::TReadGroups &ReadGroups, std::vector<uint16_t> &ReadGroupsWithoutRecal,
			      std::vector<uint16_t> &ReadGroupsLikelySingleEnd) const noexcept {
	ReadGroupsWithoutRecal.clear();
	ReadGroupsLikelySingleEnd.clear();
	for (uint16_t r = 0; r < ReadGroups.size(); ++r) {
		if (!_models[r][0]->recalibrates()) {
			ReadGroupsWithoutRecal.push_back(r);
		} else if (!_models[r][1]->recalibrates()) {
			ReadGroupsLikelySingleEnd.push_back(r);
		}
	}
}

// functions to get error rates
//-------------------------------------------------------
Probability TModels::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	if (!_models.empty()) return _models[base.readGroupID][base.isSecondMate()]->getErrorRate(base);
	return _noRecal.getErrorRate(base);
}

genometools::PhredIntProbability TModels::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	if (!_models.empty()) return _models[base.readGroupID][base.isSecondMate()]->getPhredInt(base);
	return _noRecal.getPhredInt(base);
}

void TModels::recalibrate(BAM::TSequencedBase &base) const noexcept {
	base.recalibratedQualityAsPhredInt = getPhredInt(base);
}

void TModels::recalibrate(std::vector<BAM::TSequencedBase> &bases) const noexcept {
	for (auto &b : bases) recalibrate(b);
}

void TModels::fillBaseLikelihoods(const BAM::TSequencedBase &base,
						 TBaseLikelihoods &baseLikelihoods) const noexcept {
	if (!_models.empty()) _models[base.readGroupID][base.isSecondMate()]->fillBaseLikelihoods(base, baseLikelihoods);
	else _noRecal.fillBaseLikelihoods(base, baseLikelihoods);
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
			out << ReadGroups.getName(r) << s << _models[r][mate]->getCovariateDefinition()
				<< _models[r][mate]->getRhoDefinition() << std::endl;
		}
	}
}

} // namespace SequencingError
}; // namespace GenotypeLikelihoods
