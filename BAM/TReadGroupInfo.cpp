/*
 * TSimulatorReadGroupInfo.cpp
 *
 *  Created on: Jul 14, 2022
 *      Author: phaentu
 */


#include "TFile.h"
#include "TError.h"
#include <vector>
#include "TReadGroupInfo.h"
#include "TLog.h"
#include "TParameters.h"
#include "commonWeakTypes.h"
#include "devtools.h"

using coretools::instances::parameters;
using coretools::instances::logfile;
using nlohmann::json;

namespace BAM {
namespace RGInfo{

//------------------------------------------------
// Functions to initialize: only visible in cpp file
//------------------------------------------------
namespace implJson{
	InfoType argument2InfoType(const std::string & Argument){
		for(auto i = InfoType::min; i < InfoType::max; ++i){
			if(infos[i].argument == Argument)
				return i;
		}
		return InfoType::max;
	}
}

//------------------------------------------------
// TReadGroupInfo
//------------------------------------------------
void TReadGroupInfo::_setAllReadGroups(InfoType Info, const std::string & Val){
	for(auto& i : _info){
		i.set(Info, Val);
	}
}

void TReadGroupInfo::_setDefault(InfoType Info){
	//use default values
	logfile().write("default value '",
				   infos[Info].defaults,
				   "' for all read groups. (set with '",
				   TReadGroupInfo::_RGInfoArgument,
				   "' or '",
				   infos[Info].argument,
				   "')");
	_setAllReadGroups(Info, infos[Info].defaults);
}

void TReadGroupInfo::_setFromCommandLine(InfoType Info){
	//read from command line
	const std::string& arg = infos[Info].argument;

	std::vector<std::string> tmp, argVec;
	parameters().fillParameterIntoContainer(arg, tmp, true);
	coretools::str::repeatIndexes(tmp, argVec);
	if (argVec.size() == 1){
		logfile().write("using '", argVec[0], "' for all read groups. (argument '", arg, "')");
		_setAllReadGroups(Info, argVec[0]);
	} else if (argVec.size() == _info.size()){
		logfile().write("using read group specific settings provided on the command line. (argument '", arg, "')");
		for(size_t i = 0; i < _info.size(); ++i){
			_info[i].set(Info, argVec[i]);
		}
	} else {
		UERROR("Number of provided values does not match number of read groups!");
	}
}

void TReadGroupInfo::_setFromRGInfoFile(InfoType Info){
 	 //present in file -> read for each read group
	logfile().write("reading read group specific settings from read group info file '", _filename, "'. (overwrite with '", infos[Info].argument, "')");
	for(auto& r : _info){
		if(_json.contains(r.name()) && _json[r.name()].contains(infos[Info].argument)){
			r.set(Info, _json[r.name()][infos[Info].argument]);
		} else {
			r.set(Info, infos[Info].defaults);
		}
	}
}

bool TReadGroupInfo::_readGroupExists(const std::string & Name){
	for(auto& r : _info){
		if(r[InfoType::RGName] == Name){
			return true;
		}
	}
	return false;
}

void TReadGroupInfo::_readFile(const std::string & Filename){
	std::ifstream in(Filename);
	if(!in){
		UERROR("Failed to open read group info file '", Filename, "' for reading!");
	}

	try {
		_json = nlohmann::ordered_json::parse(in);
	}
	catch (json::parse_error& ex)
	{
		UERROR("Failed to parse read group info file '", Filename, "': JSON error '", ex.what(), " at byte ", ex.byte, "!");
	}

	// warn if file contains inexisting read groups
	if(!_info.empty()){
		std::vector<std::string> ignoredRGs;
		for (auto it = _json.begin(); it != _json.end(); ++it){
			if(!_readGroupExists(it.key())){
				ignoredRGs.push_back(it.key());
			}
		}

		if(ignoredRGs.size() > 0){
			logfile().warning("The following read groups are given in the read group info file '", Filename, "' but are not present in the BAM file!");
		}
	}

	_filename = Filename;
}

void TReadGroupInfo::_createReadGroupInfoEntries(const BAM::TReadGroups & ReadGroups){
	if(!_info.empty()){
		DEVERROR("Read group info already read!");
	}

	//create read group info entries
	_info.reserve(ReadGroups.size());
	for(auto i = 0; i < ReadGroups.size(); ++i){
		_info.push_back(ReadGroups[i].name_ID);
	}
	std::fill(_parsed.begin(), _parsed.end(), false);
	_parsed[InfoType::RGName] = true;
}

// either: read info from file and match with TReadGroups (used for analyzes)
void TReadGroupInfo::readInfoAndMatchReadGroups(const BAM::TReadGroups & ReadGroups, const std::string & Filename){
	if(!_info.empty()){
		DEVERROR("Read group info already read!");
	}

	_createReadGroupInfoEntries(ReadGroups);
	if(Filename.empty()){
		_readFile(parameters().getParameter<std::string>(_RGInfoArgument, false));
	} else {
		_readFile(Filename);
	}
}

// or: read info and fill TReadGroups (used for simulations)
BAM::TReadGroups TReadGroupInfo::readInfoAndCreateReadGroups(const std::string & RgInfoFileName){
	if(!_info.empty()){
		DEVERROR("Read group info already read!");
	}

	//create empty read groups
	BAM::TReadGroups readGroups;

	// Info is provided as a) a RG info file OR b) as the number of read groups and default arguments
	if(!RgInfoFileName.empty()){
		_readFile(RgInfoFileName);

		// create RGs from RG info file
		for (auto it = _json.begin(); it != _json.end(); ++it){
			readGroups.add(it.key());
		}
	} else {
		// create identical read groups from command line
		uint16_t numRG;
		if(parameters().parameterExists(_numRGArgument)){
			numRG = parameters().getParameter<coretools::StrictlyPositive<int>>(_numRGArgument);
			if (numRG == 1) {
				logfile().list("Initializing one read group from arguments. (parameter '", _numRGArgument, "')");
			} else {
				logfile().list("Initializing ", numRG, " identical read groups from arguments (parameter '", _numRGArgument, "').");
			}
		} else {
			numRG = 1;
			logfile().list("Initializing one read group from arguments. (set with '", _numRGArgument, "')");
		}

		//create read groups
		for (int i = 0; i < numRG; ++i) {
			readGroups.add("SimReadGroup" + coretools::str::toString(i + 1));
		}
	}
	_createReadGroupInfoEntries(readGroups);

	return readGroups;
}

void TReadGroupInfo::parse(const InfoType Info){
	if(!_parsed[Info]){
		logfile().listFlush(coretools::str::capitalizeFirst(infos[Info].description), ": ");
		std::string arg = infos[Info].argument;

		//check if info is provided on the command line -> overwrites file
		if(parameters().parameterExists(arg)){
			_setFromCommandLine(Info);
		} else {
			//check if provided in file
			if(fileHasInfo(Info)){
				_setFromRGInfoFile(Info);
			} else {
				_setDefault(Info);
			}
		}
		_parsed[Info] = true;
	}
}

bool TReadGroupInfo::fileHasInfo(const InfoType Info) const {
	//return true if at least one RG has thsi info in file
	if(hasFile()){
		const std::string& arg = infos[Info].argument;
		for(auto& j : _json){
			if(j.contains(arg)){
				return true;
			}
		}
	}
	return false;
}

std::vector<std::string> TReadGroupInfo::getUnusedAttributesInFile(){
	std::vector<std::string> ret;
	if(hasFile()){
		for(auto& j : _json){
			for (auto it = j.begin(); it != j.end(); ++it){
				const InfoType& arg = implJson::argument2InfoType(it.key());
				if(arg == InfoType::max || !_parsed[arg]){
					 ret.push_back(it.key());
				}
			}
		}
	}
	return ret;
}

void TReadGroupInfo::warnAboutUnusedColumnsInFile(){
	if(hasFile()){
		auto up = getUnusedAttributesInFile();
		if(up.size() == 1){
			logfile().warning("The following attribute in read group info file '", _filename, "' was never used: ", coretools::str::concatenateString(up, ", "), "!");
		} else if(up.size() > 1){
			logfile().warning("The following attributes in read group info file '", _filename, "' were never used: ", coretools::str::concatenateString(up, ", "), "!");
		}
	}
};


void TReadGroupInfo::set(const uint16_t RGIndex, const InfoType Info, const ordered_json & Value){
	//check if info was already parsed. Else, add
	_parsed[Info] = true;

	//now add to specific
	_info[RGIndex].set(Info, Value);
}

void TReadGroupInfo::write(const std::string & Filename){
	//write RG info file
	logfile().listFlush("Writing read group info to file '", Filename, "' ...");
	for(auto& r : _info){
		//make sure json has entry for that read group
		if(!_json.contains(r.name())){
			_json[r.name()] = json::object();
		}

		// add (or overwrite) all parsed attributes
		ordered_json& x = _json[r.name()];

		for(auto i = InfoType::min; i < InfoType::max; ++i){
			if(_parsed[i]){
				x[infos[i].argument] = r[i];
			}
		}

	}

	//open file
	std::ofstream out(Filename);
	if(!out){
		UERROR("Failed to open file '", Filename, "' for writing!");
	}
	out << std::setw(4) << _json << std::endl;
	out.close();
	logfile().done();
}

//-------------------------------------------
// TTask_testReadGroupInfo
//-------------------------------------------
void TTask_testReadGroupInfo::run(){

	std::string filename = parameters().getParameter("json");

	TReadGroupInfo r;
	r.readInfoAndCreateReadGroups(filename);

	r.parse(InfoType::mappingQuality, InfoType::cycles);

	r.set(0, InfoType::cycles, "200,200");

	r.write("new.json");

}


} //end namespace RGInfo
} //end namespace BAM




