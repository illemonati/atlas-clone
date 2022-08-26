/*
 * TModels.cpp
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#include "SequencingError/TModels.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <ostream>
#include <stdexcept>

#include "TError.h"
#include "TFile.h"
#include "TLog.h"
#include "TReadGroupInfo.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "SequencingError/TModel.h"
#include "probability.h"
#include "stringFunctions.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using coretools::Probability;
using coretools::str::toString;
using coretools::instances::logfile;

//--------------------------------------------------------------------
// TModels
//--------------------------------------------------------------------

namespace impl {

auto epsRho(const std::string &s) {
	// Format: intercept[];cov1:function1[];cov2:function2[];...;rho[[]]
	const auto rBegin = s.find("rho");
	if (rBegin == std::string::npos) UERROR("Recal string ", s, " does not contain 'rho'");
	return std::make_pair(s.substr(0, rBegin-1), s.substr(rBegin + 3, s.size()));
}
} // namespace impl

void TModels::initializeNoRecal(const BAM::TReadGroups &ReadGroups) {
	_models.resize(ReadGroups.size());
	for (auto &m : _models) {
		m[0] = std::make_unique<TModelNoRecal>();
		m[1] = std::make_unique<TModelNoRecal>();
	}
}

std::vector<std::array<std::unique_ptr<TModel>, 2>> TModels::forget() {
	std::vector<std::array<std::unique_ptr<TModel>, 2>> forgottenModels(_models.size());
	for (auto &m : forgottenModels) {
		m[0] = std::make_unique<TModelNoRecal>();
		m[1] = std::make_unique<TModelNoRecal>();
	}
	std::swap(_models, forgottenModels);
	return forgottenModels;
}
void TModels::remember(std::vector<std::array<std::unique_ptr<TModel>, 2>>& forgottenModels) {
	if (forgottenModels.size() != _models.size()) DEVERROR("Forgotten models are not correct size!");
	std::swap(_models, forgottenModels);
}

void TModels::initialize(const std::string &RecalString, const std::string &RhoString,
					const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		throw std::runtime_error(
			"void TModels::initialize(const std::string & RecalString, const std::string & RhoString, "
		    "const BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	// prepare objects
	_models.resize(ReadGroups.size());

	// initialize models
	for (auto &m : _models) {
		m[0] = std::make_unique<TModelRecal>(RhoString, RecalString);
		m[1] = std::make_unique<TModelRecal>(RhoString, RecalString);
	}
}

void TModels::initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		throw std::runtime_error("void TModels::initializeFromFile(const std::string & filename, const "
					 "BAM::TReadGroups & ReadGroups, TLog* Logfile): Models already initialized!");

	// read parameters from file
	logfile().listFlush("Initializing recalibration models from '" + Filename + "' ...");
	coretools::TInputFile in(Filename, {"readGroup", "mate", "covariates", "rho"}, "\t", "//");

	// prepare objects
	_models.resize(ReadGroups.size());

	// tmp variables for reading
	std::vector<std::string> vec;

	// parse file to read details for each read group
	while (in.read(vec)) {
		if (ReadGroups.readGroupExists(vec[0])) { // ignore if it does not exist
			// get read group
			const uint16_t readGroupId = ReadGroups.getId(vec[0]);
			const auto& rhoDef = vec[3];
			const auto& epsilonDef = vec[2];
			try {
				// add model
				if (vec[1] == "first")
					_models[readGroupId][0] = std::make_unique<TModelRecal>(rhoDef, epsilonDef);
				else if (vec[1] == "second")
					_models[readGroupId][1] = std::make_unique<TModelRecal>(rhoDef, epsilonDef);
				else
					throw "Unknown mate '" + vec[1] + "! Must be 'first' or 'second'.";
			} catch (const char *error) {
				throw std::string(error) + " in file '" + Filename + "' on line " + toString(in.lineNumber()) + "!";
			}
		}
	}

	// Fill remaining with NoRecal
	for (auto &m : _models) {
		if (!m[0]) m[0] = std::make_unique<TModelNoRecal>();
		if (!m[1]) m[1] = std::make_unique<TModelNoRecal>();
	}

	logfile().done();
}

void TModels::initialize(BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;

	_models.resize(RgInfo.size());
	RgInfo.parse(InfoType::recal1, InfoType::recal2);

	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		const auto s = RgInfo.get(rg, InfoType::recal1);
		if (s == "-") { _models[rg][0] = std::make_unique<TModelNoRecal>(); }
		const auto [e, r] = impl::epsRho(s);
		_models[rg][0]    = std::make_unique<TModelRecal>(e, r);
	}
	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		const auto s = RgInfo.get(rg, InfoType::recal2);
		if (s == "-") { _models[rg][1] = std::make_unique<TModelNoRecal>(); }
		const auto [e, r] = impl::epsRho(s);
		_models[rg][1]    = std::make_unique<TModelRecal>(e, r);
	}
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
	return _models[base.readGroupID][base.isSecondMate()]->getErrorRate(base);
}

genometools::PhredIntProbability TModels::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()]->getPhredInt(base);
}

void TModels::recalibrate(BAM::TSequencedBase &base) const noexcept {
	base.recalibratedQualityAsPhredInt = getPhredInt(base);
}

void TModels::recalibrate(std::vector<BAM::TSequencedBase> &bases) const noexcept {
	for (auto &b : bases) recalibrate(b);
}

TBaseLikelihoods TModels::getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()]->getBaseLikelihoods(base);
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
