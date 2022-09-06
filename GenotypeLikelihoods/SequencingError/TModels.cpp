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
using namespace std::string_literals;

//-------------------------------------
// TReadGroupModels
//-------------------------------------

namespace impl {

std::pair<std::string, std::string> rhoEps(const std::string &s) {
	// Format: intercept[];cov1:function1[];cov2:function2[];...;rho[[]]
	const auto rBegin = s.find("rho");
	if (rBegin == std::string::npos) {
		// no rho definition
		return std::make_pair("default"s, s);
	}
	return std::make_pair(s.substr(rBegin + 3, s.size()), s.substr(0, rBegin-1));
}

void initModel(std::unique_ptr<TModel> & model, const BAM::RGInfo::TReadGroupInfoEntry & Info, const BAM::RGInfo::InfoType Type){
	if(Info.has(BAM::RGInfo::InfoType::recal1)){
			const auto [r, e] = rhoEps(Info.get(Type));
			model = std::make_unique<TModelRecal>(r, e);
		} else {
			model = std::make_unique<TModelNoRecal>();
		}
}

} // namespace impl

TReadGroupModels::TReadGroupModels(){
	_models[0] = std::make_unique<TModelNoRecal>();
	_models[1] = std::make_unique<TModelNoRecal>();
}

TReadGroupModels::TReadGroupModels(const std::string &RecalString, const std::string &RhoString){
	if(RecalString.empty() || RecalString == "-" || RecalString == "default"){
		_models[0] = std::make_unique<TModelNoRecal>();
		_models[1] = std::make_unique<TModelNoRecal>();
	} else {
		_models[0] = std::make_unique<TModelRecal>(RhoString, RecalString);
		_models[1] = std::make_unique<TModelRecal>(RhoString, RecalString);
	}
}

TReadGroupModels::TReadGroupModels(const std::string &RecalString1, const std::string &RhoString1, const std::string &RecalString2, const std::string &RhoString2){
	if(RecalString1.empty() || RecalString1 == "-" || RecalString1 == "default"){
		_models[0] = std::make_unique<TModelNoRecal>();
	} else {
		_models[0] = std::make_unique<TModelRecal>(RhoString1, RecalString1);
	}
	if(RecalString2.empty() || RecalString2 == "-" || RecalString2 == "default"){
		_models[1] = std::make_unique<TModelNoRecal>();
	} else {
		_models[1] = std::make_unique<TModelRecal>(RhoString2, RecalString2);
	}
}

TReadGroupModels::TReadGroupModels(const BAM::RGInfo::TReadGroupInfoEntry & Info){
	using BAM::RGInfo::InfoType;
	impl::initModel(_models[0], Info, InfoType::recal1);
	impl::initModel(_models[1], Info, InfoType::recal2);
}

void TReadGroupModels::simulate(BAM::TAlignment & Alignment) const {
	const TModel& mod = _models[Alignment.isSecondMate()];
	for (auto & b : Alignment) {
		mod.simulate(b);
	}
}

//--------------------------------------------------------------------
// TModels
//--------------------------------------------------------------------

namespace impl {

void fillNoRecal(std::vector<TReadGroupModels> & vec, const size_t size){
	vec.reserve(size);
	for(uint16_t i = 0; i < size; ++i){
		vec.emplace_back(); //constructor without arguments
	}
}

} // namespace impl

void TModels::initializeNoRecal(const BAM::TReadGroups &ReadGroups) {
	impl::fillNoRecal(_models, ReadGroups.size());
}

std::vector<TReadGroupModels> TModels::forget() {
	std::vector<TReadGroupModels> forgottenModels;
	impl::fillNoRecal(forgottenModels, _models.size());
	std::swap(_models, forgottenModels);
	return forgottenModels;
}

void TModels::remember(std::vector<TReadGroupModels>& forgottenModels) {
	if (forgottenModels.size() != _models.size()) DEVERROR("Forgotten models are not correct size!");
	std::swap(_models, forgottenModels);
}

void TModels::initialize(const std::string &RecalString, const std::string &RhoString,
					const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		DEVERROR("Models already initialized!");

	// prepare objects
	_models.reserve(ReadGroups.size());

	// initialize models
	for(uint16_t i = 0; i < ReadGroups.size(); ++i){
		_models.emplace_back(RhoString, RecalString);
	}
}

void TModels::initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		DEVERROR("Models already initialized!");

	// read parameters from file
	logfile().listFlush("Initializing recalibration models from '" + Filename + "' ...");
	coretools::TInputFile in(Filename, {"readGroup", "covariates1", "rho1", "covariates2", "rho2"}, "\t", "//");

	// prepare objects
	_models.resize(ReadGroups.size());

	// tmp variables for reading
	std::vector<std::string> vec;

	//store per RG info
	std::vector< std::vector<std::string> > info(ReadGroups.size());

	// parse file to read details for each read group
	while (in.read(vec)) {
		if (ReadGroups.readGroupExists(vec[0])) { // ignore if it does not exist
			// store by read group id
			info[ReadGroups.getId(vec[0])] = vec;
		}
	}

	//create models
	_models.reserve(ReadGroups.size());
	for(uint16_t r = 0; r < ReadGroups.size(); ++r){
		try {
			_models.emplace_back(info[r][1], info[r][2], info[r][3], info[r][4]);
		} catch (const char *error) {
			throw std::string(error) + " for read group " + ReadGroups[r].name_ID + " in file '" + Filename + "!";
		}
	}
	logfile().done();
}

void TModels::initialize(BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;

	RgInfo.parse(InfoType::recal1, InfoType::recal2);
	_models.reserve(RgInfo.size());

	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		_models.emplace_back(RgInfo[rg]);
	}
}

void TModels::checkReadGroups(const BAM::TReadGroups &ReadGroups, std::vector<uint16_t> &ReadGroupsWithoutRecal,
			      std::vector<uint16_t> &ReadGroupsLikelySingleEnd) const noexcept {
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
Probability TModels::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()].getErrorRate(base);
}

genometools::PhredIntProbability TModels::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()].getPhredInt(base);
}

void TModels::recalibrate(BAM::TSequencedBase &base) const noexcept {
	base.recalibratedQualityAsPhredInt = getPhredInt(base);
}

void TModels::recalibrate(std::vector<BAM::TSequencedBase> &bases) const noexcept {
	for (auto &b : bases) recalibrate(b);
}

TBaseLikelihoods TModels::getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()].getBaseLikelihoods(base);
}

// functions to write file
//-------------------------------------------------------------------
void TModels::writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const {
	// open file and write header
	coretools::TOutputFile out(Filename);
	out.writeHeader({"readGroup", "covariates1", "rho1", "covariates2", "rho2"});

	// add models
	for (uint16_t r = 0; r < ReadGroups.size(); ++r) {
		out << ReadGroups.getName(r);
		for (uint8_t mate = 0; mate < 2; ++mate) {
			out << _models[r][mate].getCovariateDefinition()
				<< _models[r][mate].getRhoDefinition();
		}
		out.endLine();
	}
}

} // namespace SequencingError
}; // namespace GenotypeLikelihoods
