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

void TReadGroupModels::initialize(bool isSecondMate, std::string_view RecalString, std::string_view RhoString) {
	if(RecalString.empty() || RecalString == "-" || RecalString == "default"){
		_models[isSecondMate] = std::make_unique<TModelNoRecal>();
	} else {
		_models[isSecondMate] = std::make_unique<TModelRecal>(RecalString, RhoString);
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

void TModels::pool(const BAM::TReadGroupMap& rgMap) {
	if (!recalibrates()) DEVERROR("Cannot pool models that do not recalibrate!");
	_pModels.resize(_models.size());
	for (size_t i = 0; i < _models.size(); ++i) {
		_pModels[i] = &_models[rgMap.pooledIndex(i)];
	}
}

void TModels::unPool() {
	_pModels.clear();
	_recalibrates = false;
	for (auto &m : _models) {
		_pModels.push_back(&m);
		if (m.recalibrates()) _recalibrates = true;
	}
}

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
	std::vector<TReadGroupModels> forgottenModels;
	//impl::fillNoRecal(forgottenModels, _models.size());
	const auto N = _pModels.size();
	std::swap(_models, forgottenModels);

	_models.clear();
	_models.emplace_back(); // default-constructor -> NoRecal
	_pModels.clear();

	for(size_t i = 0; i < N; ++i){
		_pModels.push_back(&_models.front());
	}
	_recalibrates=false;
	return forgottenModels;
}

void TModels::remember(std::vector<TReadGroupModels>& forgottenModels) {
	if (forgottenModels.size() != _pModels.size()) DEVERROR("Forgotten models are not correct size!");
	std::swap(_models, forgottenModels);
	unPool();
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
	}
	unPool();
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
		} catch (const char *error) {
			UERROR(error, " for read group ", vec[0], " in file '", Filename, "!");
		}
	}
	unPool();
	logfile().done();
}

void TModels::initialize(BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;

	_models.reserve(RgInfo.size());
	_recalibrates=false;

	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		_models.emplace_back(RgInfo[rg]);
	}
	unPool();
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
	const auto & m = model(front); // assuming all datas are in same readgroup and same mate
	for (auto &b : datas) {
		if (m.recalibrates()) {
			b.recalibratedQualityAsPhredInt = m.phredInt(b);
		} else {
			b.recalibratedQualityAsPhredInt = b.originalQuality_phredInt;
		}
	}
}
void TModels::writeRecalFile(const BAM::TReadGroups &ReadGroups, std::string_view Filename) const {
	// open file and write header
	coretools::TOutputFile out(Filename, {"readGroup", "mate", "covariates", "rho"});

	// add models
	for (size_t r = 0; r < ReadGroups.size(); ++r) {
		for (size_t mate = 0; mate < 2; ++mate) {
			out.write(ReadGroups.getName(r), std::array{"first", "second"}[mate]);
			out.writeln(_model(r)[mate].epsilonDefinition(), _model(r)[mate].rhoDefinition());
		}
	}
}

void TModels::addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const {
	for(size_t r = 0; r < _pModels.size(); ++r){
		RgInfo.set(r, BAM::RGInfo::InfoType::recal, _model(r).info());
	}
}

} // namespace SequencingError
}; // namespace GenotypeLikelihoods
