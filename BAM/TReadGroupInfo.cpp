/*
 * TSimulatorReadGroupInfo.cpp
 *
 *  Created on: Jul 14, 2022
 *      Author: phaentu
 */

#include "TReadGroupInfo.h"
#include "coretools/Files/TStdWriter.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Types/commonWeakTypes.h"
#include "coretools/enum.h"
#include <vector>

using coretools::instances::logfile;
using coretools::instances::parameters;

namespace BAM::RGInfo {

//------------------------------------------------
// Functions to initialize: only visible in cpp file
//------------------------------------------------
namespace implJson {
InfoType argument2InfoType(std::string_view Argument) {
	for (auto i = InfoType::min; i < InfoType::max; ++i) {
		if (infos[i].argument == Argument) return i;
	}
	return InfoType::max;
}
} // namespace implJson

//------------------------------------------------
// TReadGroupInfo
//------------------------------------------------
void TReadGroupInfo::_setAllReadGroups(InfoType Info, std::string_view Val) {
	for (auto &i : _info) { i[Info] = Val; }
}

void TReadGroupInfo::_setDefault(InfoType Info) {
	// use default values
	logfile().write("default value '", infos[Info].defaults, "' for all read groups. (set with '",
	                TReadGroupInfo::RGInfoArgument, "' or '", infos[Info].argument, "')");
	_setAllReadGroups(Info, infos[Info].defaults);
}

void TReadGroupInfo::_setFromCommandLine(InfoType Info) {
	// read from command line
	const std::string &key = infos[Info].argument;
	const auto arg         = parameters().get(key, "");

	if (arg.empty()) {
		logfile().write("using default for all read groups. (parameter '", key, "')");
		_setAllReadGroups(Info, "");
	} else {
		logfile().write("using '", arg, "' for all read groups. (parameter '", key, "')");
		_setAllReadGroups(Info, arg);
	}
}

void TReadGroupInfo::_setFromRGInfoFile(InfoType Info) {
	// present in file -> read for each read group
	logfile().write("reading read group specific settings from read group info file '", _filename,
	                "'. (overwrite with '", infos[Info].argument, "')");
	for (auto &r : _info) {
		if (_json.contains(r.name()) && _json[r.name()].contains(infos[Info].argument)) {
			r[Info] = _json[r.name()][infos[Info].argument];
		} else {
			r[Info] = infos[Info].defaults;
		}
	}
}

bool TReadGroupInfo::_readGroupExists(std::string_view Name) {
	for (auto &r : _info) {
		if (r[InfoType::RGName] == Name) { return true; }
	}
	return false;
}

void TReadGroupInfo::_readFile(std::string_view Filename) {
	try {
		_json = nlohmann::ordered_json::parse(std::ifstream(std::string(Filename)));
	} catch (nlohmann::json::parse_error &ex) {
		UERROR("Failed to parse read group info file '", Filename, "': JSON error '", ex.what(), " at byte ", ex.byte,
		       "!");
	}

	// warn if file contains inexisting read groups
	if (!_info.empty()) {
		std::vector<std::string> ignoredRGs;
		for (auto it = _json.begin(); it != _json.end(); ++it) {
			if (!_readGroupExists(it.key())) { ignoredRGs.push_back(it.key()); }
		}

		if (ignoredRGs.size() > 0) {
			logfile().warning("The following read groups are given in the read group info file '", Filename,
			                  "' but are not present in the BAM file!");
		}
	}

	_filename = Filename;
}

void TReadGroupInfo::_createReadGroupInfoEntries(const BAM::TReadGroups &ReadGroups) {
	if (!_info.empty()) { DEVERROR("Read group info already read!"); }

	// create read group info entries
	_info.reserve(ReadGroups.size());
	for (size_t i = 0; i < ReadGroups.size(); ++i) { _info.emplace_back(this, ReadGroups[i].name_ID); }
	_parsed.reset();
	_parsed.set<InfoType::RGName>();
}

bool TReadGroupInfo::isParsed() const {
	// Not counting RGName
	for (auto info_t = coretools::next(InfoType::min); info_t < InfoType::max; ++info_t) {
		if (_parsed[info_t]) return true;
	}
	return false;
}

void TReadGroupInfo::_parse(const InfoType Info) {
	if (!_parsed[Info]) {
		logfile().listFlush(coretools::str::capitalizeFirst(infos[Info].description), ": ");
		std::string arg = infos[Info].argument;

		// check if info is provided on the command line -> overwrites file
		if (parameters().exists(arg)) {
			_setFromCommandLine(Info);
		} else {
			// check if provided in file
			if (fileHasInfo(Info)) {
				_setFromRGInfoFile(Info);
			} else {
				_setDefault(Info);
			}
		}
		_parsed[Info] = true;
	}
}

TReadGroupInfo::TReadGroupInfo(const BAM::TReadGroups &ReadGroups) {
	_createReadGroupInfoEntries(ReadGroups);
	if (parameters().exists(RGInfoArgument)) { _readFile(parameters().get(RGInfoArgument)); }
}

// either: read info from file and match with TReadGroups (used for analyzes)
TReadGroupInfo::TReadGroupInfo(const BAM::TReadGroups &ReadGroups, std::string_view Filename) {
	_createReadGroupInfoEntries(ReadGroups);
	_readFile(Filename);
}

// or: read info and fill TReadGroups (used for simulations)
BAM::TReadGroups TReadGroupInfo::createReadGroups(std::string_view RgInfoFileName) {
	if (!_info.empty()) { DEVERROR("Read group info already read!"); }

	// create empty read groups
	BAM::TReadGroups readGroups;

	// Info is provided as a) a RG info file OR b) as the number of read groups and default arguments
	if (!RgInfoFileName.empty()) {
		_readFile(RgInfoFileName);

		// create RGs from RG info file
		for (auto it = _json.begin(); it != _json.end(); ++it) { readGroups.add(it.key()); }
	} else {
		// create identical read groups from command line
		const auto numRG = parameters().get<coretools::StrictlyPositiveInt>(numRGArgument, 1);
		if (numRG == 1) {
			logfile().list("Initializing one read group from arguments. (parameter '", numRGArgument, "')");
		} else {
			logfile().list("Initializing ", numRG, " identical read groups from arguments (parameter '", numRGArgument,
			               "').");
		}

		// create read groups
		for (int i = 0; i < numRG; ++i) { readGroups.add("SimReadGroup" + coretools::str::toString(i + 1)); }
	}
	_createReadGroupInfoEntries(readGroups);

	return readGroups;
}

bool TReadGroupInfo::fileHasInfo(const InfoType Info) const {
	// return true if at least one RG has this info in file
	if (hasFile()) {
		const std::string &arg = infos[Info].argument;
		for (auto &j : _json) {
			if (j.contains(arg)) { return true; }
		}
	}
	return false;
}

std::vector<std::string> TReadGroupInfo::getUnusedAttributesInFile() {
	std::vector<std::string> ret;
	if (hasFile()) {
		for (auto &j : _json) {
			for (auto it = j.begin(); it != j.end(); ++it) {
				const InfoType &arg = implJson::argument2InfoType(it.key());
				if (arg == InfoType::max || !_parsed[arg]) { ret.push_back(it.key()); }
			}
		}
	}
	return ret;
}

void TReadGroupInfo::warnAboutUnusedColumnsInFile() {
	if (hasFile()) {
		auto up = getUnusedAttributesInFile();
		if (up.size() == 1) {
			logfile().warning("The following attribute in read group info file '", _filename,
			                  "' was never used: ", coretools::str::concatenateString(up, ", "), "!");
		} else if (up.size() > 1) {
			logfile().warning("The following attributes in read group info file '", _filename,
			                  "' were never used: ", coretools::str::concatenateString(up, ", "), "!");
		}
	}
};

void TReadGroupInfo::set(const uint16_t RGIndex, const InfoType Type, const TInfo &Value) {
	// check if info was already parsed. Else, add
	_parsed[Type] = true;

	// now add to specific
	assert(RGIndex < _info.size());
	_info[RGIndex][Type] = Value;
}

void TReadGroupInfo::write(std::string_view Filename) {
	// write RG info file
	logfile().list("Writing read group info to file '", Filename, "'.");
	for (auto &r : _info) {
		// make sure json has entry for that read group
		if (!_json.contains(r.name())) { _json[r.name()]; }

		// add (or overwrite) all parsed attributes
		TInfo &x = _json[r.name()];

		for (auto i = InfoType::min; i < InfoType::max; ++i) {
			if (_parsed[i]) { x[infos[i].argument] = r[i]; }
		}
	}

	// open file
	coretools::TStdWriter writer(Filename);
	writer.write(_json.dump(2));
}

//-------------------------------------------
// TTask_testReadGroupInfo
//-------------------------------------------
void TReadGroupInfoTest::run() {

	std::string filename = parameters().get("json");

	TReadGroupInfo r;
	r.createReadGroups(filename);

	r[0][InfoType::cycles] = "200,200";
	r.write("new.json");
}

} // namespace BAM::RGInfo
