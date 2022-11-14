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

#include "coretools/Main/TError.h"
#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "TReadGroupInfo.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "SequencingError/TModel.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using coretools::Probability;
using coretools::instances::logfile;
	using namespace std::string_literals;

//--------------------------------------------------------------------
// TModels
//--------------------------------------------------------------------

namespace impl {

std::pair<std::string, std::string> epsRho(const std::string &s) {
	// Format: intercept[];cov1:function1[];cov2:function2[];...;rho[[]]
	const auto rBegin = s.find("rho");
	if (rBegin == std::string::npos) {
		// no rho definition
		return std::make_pair(s, "default"s);
	}
	return std::make_pair(s.substr(0, rBegin-1), s.substr(rBegin + 3, s.size()));
}

void initModel(std::unique_ptr<TModel> & model, const BAM::RGInfo::TReadGroupInfoEntry & Info, const BAM::RGInfo::InfoType Type){
	if(Info.has(Type)){
		const auto [e, r] = epsRho(Info.get(Type));
		model = std::make_unique<TModelRecal>(e, r);
	} else {
		model = std::make_unique<TModelNoRecal>();
	}
}

void initModel(std::unique_ptr<TModel> & model, const BAM::RGInfo::TInfo & info){
	if(info.empty() || (info.is_string() && info.get<std::string_view>() == "default")){
		model = std::make_unique<TModelNoRecal>();
	} else {
		// TODO
		model = std::make_unique<TModelRecal>(info.get<std::string>());
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
		_models[0] = std::make_unique<TModelRecal>(RecalString, RhoString);
		_models[1] = std::make_unique<TModelRecal>(RecalString, RhoString);
	}
}

TReadGroupModels::TReadGroupModels(const std::string &RecalString1, const std::string &RhoString1, const std::string &RecalString2, const std::string &RhoString2){
	if(RecalString1.empty() || RecalString1 == "-" || RecalString1 == "default"){
		_models[0] = std::make_unique<TModelNoRecal>();
	} else {
		_models[0] = std::make_unique<TModelRecal>(RecalString1, RhoString1);
	}
	if(RecalString2.empty() || RecalString2 == "-" || RecalString2 == "default"){
		_models[1] = std::make_unique<TModelNoRecal>();
	} else {
		_models[1] = std::make_unique<TModelRecal>(RecalString2, RhoString2);
	}
}

void TReadGroupModels::initialize(size_t mate, const std::string &RecalString, const std::string &RhoString) {
	if(RecalString.empty() || RecalString == "-" || RecalString == "default"){
		_models[mate] = std::make_unique<TModelNoRecal>();
	} else {
		_models[mate] = std::make_unique<TModelRecal>(RecalString, RhoString);
	}
}

TReadGroupModels::TReadGroupModels(const BAM::RGInfo::TReadGroupInfoEntry & Info){
	using BAM::RGInfo::InfoType;

	//check if recal is provided
	if(Info.has(InfoType::recal)){
		auto& json = Info[InfoType::recal];

		//is this a single-end read group?
		bool single = false;
		if(Info.has(InfoType::seqType) && Info[InfoType::seqType] == BAM::RGInfo::seqType::single){
			single = true;
			_models[0] = std::make_unique<TModelNoRecal>();
		}

		//check if two mates are provided
		if(json.contains("first")){
			impl::initModel(_models[0], json["first"]);
			if(json.contains("second")){
				if(single){
					UERROR("Recal provided for second mate of single-end read group '", Info.name(), "'!");
				} else {
					impl::initModel(_models[0], json["second"]);
				}
			} else {
				_models[1] = std::make_unique<TModelNoRecal>();
			}
		} else {
			//assume a single recal model is provided
			impl::initModel(_models[0], json);
			if(!single){
				impl::initModel(_models[1], json);
			}
		}

	} else {
		//initialize no recal model
		_models[0] = std::make_unique<TModelNoRecal>();
		_models[1] = std::make_unique<TModelNoRecal>();
	}
}

void TReadGroupModels::simulate(BAM::TAlignment & Alignment) const {
	const TModel& mod = *_models[Alignment.isSecondMate()];
	for (auto & b : Alignment) {
		mod.simulate(b);
	}
}

BAM::RGInfo::TInfo TReadGroupModels::info() const{
	BAM::RGInfo::TInfo info;
	info["Mate1"]  = _models[0]->info();
	info["Mate2"] = _models[1]->info();
	return info;
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
		_models.emplace_back(RecalString, RhoString);
	}
}

void TModels::initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		DEVERROR("Models already initialized!");

	// read parameters from file
	logfile().listFlush("Initializing recalibration models from '" + Filename + "' ...");
	coretools::TInputFile in(Filename, {"readGroup", "mate", "covariates", "rho"}, "\t", "//");

	// prepare objects

	// tmp variables for reading
	std::vector<std::string> vec;

	//store per RG info
	std::vector< std::vector<std::string> > info(ReadGroups.size());

	_models.reserve(ReadGroups.size());
	// parse file to read details for each read group
	while (in.read(vec)) {
		OUT(vec);
		try {
			const auto readGroup = ReadGroups.getId(vec[0]);
			const auto mate      = vec[1] == "first" ? 0 : 1;
			WINK();
			_models[readGroup].initialize(mate, vec[2], vec[3]);
			WINK();
		} catch (const char *error) {
			throw std::string(error) + " for read group " + vec[0] + " in file '" + Filename + "!";
		}
	}
	WINK();

	logfile().done();
}

void TModels::initialize(BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;

	RgInfo.parse(InfoType::recal);
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
Probability TModels::errorRate(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()].errorRate(base);
}

genometools::PhredIntProbability TModels::phredInt(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()].phredInt(base);
}

void TModels::recalibrate(BAM::TSequencedBase &base) const noexcept {
	base.recalibratedQualityAsPhredInt = phredInt(base);
}

void TModels::recalibrate(std::vector<BAM::TSequencedBase> &bases) const noexcept {
	for (auto &b : bases) recalibrate(b);
}

TBaseLikelihoods TModels::baseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	return _models[base.readGroupID][base.isSecondMate()].baseLikelihoods(base);
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
			out << _models[r][mate].epsilonDefinition()
				<< _models[r][mate].rhoDefinition();
		}
		out.endln();
	}
}

void TModels::addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const {
	for(size_t r = 0; r < _models.size(); ++r){
		RgInfo.set(r, BAM::RGInfo::InfoType::recal, _models[r].info());
	}
}

} // namespace SequencingError
}; // namespace GenotypeLikelihoods
