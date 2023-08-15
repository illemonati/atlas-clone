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

using coretools::instances::logfile;

//--------------------------------------------------------------------
// TModels
//--------------------------------------------------------------------

namespace impl {

std::pair<std::string_view, std::string_view> epsRho(std::string_view s) {
	// Format: intercept[];cov1:function1[];cov2:function2[];...;rho[[]]
	const auto rBegin = s.find("rho");
	if (rBegin == std::string::npos) {
		// no rho definition
		return std::make_pair(s, "default");
	}
	return std::make_pair(s.substr(0, rBegin-1), s.substr(rBegin + 3, s.size()));
}

void initModel(std::unique_ptr<TModel> & model, const BAM::RGInfo::TInfo & info){
	if(info.empty() || (info.is_string() && info.get<std::string_view>() == "default")){
		model = std::make_unique<TModelNoRecal>();
	} else if (info.is_string()) {
		auto [recal, rho] = epsRho(info.get<std::string_view>());
		model = std::make_unique<TModelRecal>(recal, rho);
	}
	else {
		model = std::make_unique<TModelRecal>(info);
	}
}

std::vector<TReadGroupModels *> fillPointers(std::vector<TReadGroupModels> &Models) {
	std::vector<TReadGroupModels *> pModels;
	for (auto &m : Models) pModels.push_back(&m);
	return pModels;
}

} // namespace impl

TReadGroupModels::TReadGroupModels(){
	_models[0] = std::make_unique<TModelNoRecal>();
	_models[1] = std::make_unique<TModelNoRecal>();
}

TReadGroupModels::TReadGroupModels(std::string_view RecalString, std::string_view RhoString){
	if(RecalString.empty() || RecalString == "-" || RecalString == "default"){
		_models[0] = std::make_unique<TModelNoRecal>();
		_models[1] = std::make_unique<TModelNoRecal>();
	} else {
		_models[0] = std::make_unique<TModelRecal>(RecalString, RhoString);
		_models[1] = std::make_unique<TModelRecal>(RecalString, RhoString);
	}
}

TReadGroupModels::TReadGroupModels(std::string_view RecalString1, std::string_view RhoString1, std::string_view RecalString2, std::string_view RhoString2){
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

void TReadGroupModels::initialize(size_t mate, std::string_view RecalString, std::string_view RhoString) {
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
		}

		//check if two mates are provided
		if(json.contains("Mate1")){
			impl::initModel(_models[0], json["Mate1"]);
			if(json.contains("Mate2")){
				if(single){
					// what should we do?
				}
				impl::initModel(_models[1], json["Mate2"]);
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

void TModels::initializeNoRecal(const BAM::TReadGroups &ReadGroups) {
	_models.clear();
	_models.emplace_back(); // default-constructor -> NoRecal
	_pModels.clear();

	for(size_t i = 0; i < ReadGroups.size(); ++i){
		_pModels.push_back(&_models.front());
	}
	_recalibrates=false;
}

std::vector<TReadGroupModels> TModels::forget() {
	std::vector<TReadGroupModels> forgottenModels(_models.size());
	//impl::fillNoRecal(forgottenModels, _models.size());
	std::swap(_models, forgottenModels);
	_pModels = impl::fillPointers(_models);
	_recalibrates=false;
	return forgottenModels;
}

void TModels::remember(std::vector<TReadGroupModels>& forgottenModels) {
	if (forgottenModels.size() != _models.size()) DEVERROR("Forgotten models are not correct size!");
	std::swap(_models, forgottenModels);
	_pModels = impl::fillPointers(_models);
}

void TModels::initialize(std::string_view RecalString, std::string_view RhoString,
					const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		DEVERROR("Models already initialized!");

	// prepare objects
	_models.reserve(ReadGroups.size());
	_recalibrates=false;

	// initialize models
	for(size_t i = 0; i < ReadGroups.size(); ++i){
		_models.emplace_back(RecalString, RhoString);
		if (_models.back().recalibrates()) _recalibrates = true;
	}
	_pModels = impl::fillPointers(_models);
}

void TModels::initializeFromFile(std::string_view Filename, const BAM::TReadGroups &ReadGroups) {
	if (!_models.empty())
		DEVERROR("Models already initialized!");

	// read parameters from file
	logfile().listFlush("Initializing recalibration models from '", Filename, "' ...");
	coretools::TInputFile in(Filename, {"readGroup", "mate", "covariates", "rho"}, "\t", "//");

	// prepare objects

	// tmp variables for reading
	std::vector<std::string> vec;

	_models.resize(ReadGroups.size());
	_recalibrates=false;
	// parse file to read details for each read group
	while (in.read(vec)) {
		try {
			const auto readGroup = ReadGroups.getId(vec[0]);
			const auto mate      = vec[1] == "first" ? 0 : 1;
			_models[readGroup].initialize(mate, vec[2], vec[3]);
			if (_models[readGroup][mate].recalibrates()) _recalibrates = true;
		} catch (const char *error) {
			UERROR(error, " for read group ", vec[0], " in file '", Filename, "!");
		}
	}
	_pModels = impl::fillPointers(_models);

	logfile().done();
}

void TModels::initialize(BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;

	_models.reserve(RgInfo.size());
	_recalibrates=false;

	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		_models.emplace_back(RgInfo[rg]);
		if (_models.back().recalibrates()) _recalibrates = true;
	}
	_pModels = impl::fillPointers(_models);
}

void TModels::checkReadGroups(const BAM::TReadGroups &ReadGroups, std::vector<size_t> &ReadGroupsWithoutRecal,
			      std::vector<size_t> &ReadGroupsLikelySingleEnd) const noexcept {
	ReadGroupsWithoutRecal.clear();
	ReadGroupsLikelySingleEnd.clear();
	for (size_t r = 0; r < ReadGroups.size(); ++r) {
		if (!_model(r)[0].recalibrates()) {
			ReadGroupsWithoutRecal.push_back(r);
		} else if (!_model(r)[1].recalibrates()) {
			ReadGroupsLikelySingleEnd.push_back(r);
		}
	}
}

// functions to get error rates
//-------------------------------------------------------

void TModels::recalibrate(std::vector<BAM::TSequencedBase> &datas) const noexcept {
	const auto & front = datas.front();
	const auto & model = _model(front); // assuming all datas are in same readgroup and same mate
	for (auto &b : datas) {
		if (model.recalibrates()) {
			b.recalibratedQualityAsPhredInt = model.phredInt(b);
		} else {
			b.recalibratedQualityAsPhredInt = b.originalQuality_phredInt;
		}
	}
}

// functions to write file
//-------------------------------------------------------------------
void TModels::writeRecalFile(const BAM::TReadGroups &ReadGroups, std::string_view  Filename) const {
	// open file and write header
	coretools::TOutputFile out(Filename);
	out.writeHeader({"readGroup", "covariates1", "rho1", "covariates2", "rho2"});

	// add models
	for (size_t r = 0; r < ReadGroups.size(); ++r) {
		out << ReadGroups.getName(r);
		for (size_t mate = 0; mate < 2; ++mate) {
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
