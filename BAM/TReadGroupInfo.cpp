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

using coretools::instances::parameters;
using coretools::instances::logfile;


namespace BAM {
namespace RGInfo{

//------------------------------------------------
// TFileData
//------------------------------------------------

TFileData::TFileData(const std::string & Filename){
	//open RG file
	_filename = Filename;
	logfile().listFlush("Reading read group info from file '" + _filename + "' ... ");
	coretools::TInputFile in(Filename, coretools::TFile_Filetype::header, "\t", "//");
	_header = in.header();

	//extract RG column
	//check that file has a column named "readGroup"
	if(!in.hasColname(infos[InfoType::RGName].argument)){
		UERROR("Column '", infos[InfoType::RGName].argument, "' missing in file '", Filename, "'!");
	}
	auto rgCol = in.getIndexOfColname(infos[InfoType::RGName].argument);

	//read file and create read group entries
	std::vector<std::string> tmp;
	while (in.read(tmp)) {
		_fileData.push_back(tmp);
		_rgNames.push_back(tmp[rgCol]);
	}
	logfile().done();
	logfile().conclude("Found ", _rgNames.size(), " read groups.");

	// ensure names are unique
	std::vector<std::string> rgNames = _rgNames;
	sort(rgNames.begin(), rgNames.end());
	if(std::adjacent_find(rgNames.cbegin(), rgNames.cend()) != rgNames.cend()){
		DEVERROR("Duplicate read group names in file '", Filename, "'!");
	}
}

bool TFileData::hasInfo(InfoType Info) const noexcept {
	auto it = find(_header.cbegin(), _header.cend(), infos[Info].argument);
	return it != _header.cend();
}

size_t TFileData::getInfoCol(InfoType Info) const {
	auto it = find(_header.cbegin(), _header.cend(), infos[Info].argument);
	if(it != _header.cend()){
		return it - _header.cbegin();
	}
	DEVERROR("Info '", infos[Info].argument, "' not present in read group info file!");
}

size_t TFileData::getRow(const std::string & ReadGroupName) const {
	auto it = std::find(_rgNames.cbegin(), _rgNames.cend(), ReadGroupName);
	if(it == _rgNames.cend()){
		UERROR("Read group '", ReadGroupName, "' missing in read group info file '", _filename, "'!");
	}
	return it - _rgNames.cbegin();
}

//------------------------------------------------
// Functions to initialize: only visible in cpp file
//------------------------------------------------
namespace impl{

	using InfoVec = std::vector<TReadGroupInfoEntry>;

	void setAllReadGroups(InfoVec & Vec, InfoType Info, const std::string & Val){
		for(auto& i : Vec){
			i.set(Info, Val);
		}
	}

	void setDefault(InfoVec & Vec, InfoType Info){
		//use default values
		logfile().list("Initializing ",
					   infos[Info].description,
					   " with default value '",
					   infos[Info].defaults,
					   "' for all read groups. (set with '",
					   TReadGroupInfo::_RGInfoArgument,
					   "' or '",
					   infos[Info].argument,
					   "')");
		setAllReadGroups(Vec, Info, infos[Info].defaults);
	}

	void setFromCommandLine(InfoVec & Vec, InfoType Info){
		//read from command line
		const std::string& arg = infos[Info].argument;
		std::string val = parameters().getParameter<std::string>(arg);
		logfile().list("Initializing ", infos[Info].description, " with '", val, "' for all read groups. (argument '", arg, "')");
		setAllReadGroups(Vec, Info, val);
	}

	void setFromRGInfoFile(InfoVec & Vec, InfoType Info, const TFileData & FileData){
	 	 //present in file -> read for each read group
		logfile().list("Initializing ", infos[Info].description, " from read group info file. (overwrite with '", infos[Info].argument, "')");
		auto col = FileData.getInfoCol(Info);
		for(auto& r : Vec){
			auto row = FileData.getRow(r[InfoType::RGName]);
			r.set(Info, FileData[row][col]);
		}
	}

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

// either: read info from file and match with TReadGroups (used for analyzes)
void TReadGroupInfo::readInfoAndMatchReadGroups(const BAM::TReadGroups & ReadGroups){
	if(!_info.empty()){
		DEVERROR("Read group info already read!");
	}

	//create read group info entries
	_info.reserve(ReadGroups.size());
	for(auto i = 0; i < ReadGroups.size(); ++i){
		_info.push_back(ReadGroups[i].name_ID);
	}
	_parsed[InfoType::RGName];
	_readFileIfProvided();
}

// or: read info and fill TReadGroups (used for simulations)
BAM::TReadGroups TReadGroupInfo::readInfoAndCreateReadGroups(){
	if(!_info.empty()){
		DEVERROR("Read group info already read!");
	}

	//create empty read groups
	BAM::TReadGroups readGroups;

	// Info is provided as a) a RG info file OR b) as the number of read groups and default arguments
	_readFileIfProvided();
	if(_fileData){
		// create RGs from RG info file
		//create read groups
		auto col = _fileData->getInfoCol(InfoType::RGName);
		for(size_t i = 0; i < _fileData->size(); i++){
			readGroups.add((*_fileData)[i][col]);
		}
	} else {
		// create identical read groups from command line
		const auto numRG = parameters().getParameterWithDefault<coretools::StrictlyPositive<int>>(_numRGArgument, 1);
		if (numRG == 1) {
			logfile().list("Initializing one read group from arguments. (parameter '", _numRGArgument, "')");
		} else {
			logfile().list("Initializing ", numRG, " identical read groups from arguments (parameter '", _numRGArgument, "').");
		}

		//create read groups
		for (int i = 0; i < numRG; ++i) {
			readGroups.add("SimReadGroup" + coretools::str::toString(i + 1));
		}
	}
	_parsed[InfoType::RGName];
	return readGroups;
}

void TReadGroupInfo::_readFileIfProvided(){
	if(parameters().parameterExists(_RGInfoArgument)){
		_fileData = std::make_unique<TFileData>(parameters().getParameter<std::string>(_RGInfoArgument));
	}
}

void TReadGroupInfo::parse(const InfoType Info){
	//check if info is provided on the command line -> overwrites file
	std::string arg = infos[Info].argument;
	if(parameters().parameterExists(arg)){
		impl::setFromCommandLine(_info, Info);
	} else {
		//check if provided in file
		if(_fileData && _fileData->hasInfo(Info)){
			impl::setFromRGInfoFile(_info, Info, *_fileData);
		} else {
			impl::setDefault(_info, Info);
		}
	}
	_parsed[Info];
}

std::vector<std::string> TReadGroupInfo::getUnusedColumnsInFile(){
	std::vector<std::string> ret;
	if(_fileData){
		for(auto& s : _fileData->header()){
			InfoType arg = impl::argument2InfoType(s);
			if(arg == InfoType::max || !_parsed[arg]){
				 ret.push_back(s);
			}
		}
	}
	return ret;
}

void TReadGroupInfo::warnAboutUnusedColumnsInFile(){
	if(_fileData){
		auto up = getUnusedColumnsInFile();
		if(!up.empty()){
			logfile().warning("The following columns in read group info file '", _fileData->filename(), "' were never used: ", coretools::str::concatenateString(up, ", "), "!");
		}
	}
};

void TReadGroupInfo::set(const uint16_t RGIndex, const InfoType Info, const std::string & Value){
	//check if info was already parsed. Else, add
	_parsed[Info] = true;

	//now add to specific
	_info[RGIndex].set(Info, Value);
}

void TReadGroupInfo::write(const std::string & Filename){
	//write RG info file
	if(_fileData){
		// keep order of original file and add new columns at end

		// 1) compile header and open file
		// keep original header without RGName columns (we will put that first)
		std::vector<std::string> header;
		std::vector<std::string> fileHeader = _fileData->header();
		header.push_back(infos[InfoType::RGName].argument);
		for(auto& s : fileHeader){
			if(s != infos[InfoType::RGName].argument){
				 header.push_back(s);
			 }
		}

		//add novel columns
		for(auto i = InfoType::min; i < InfoType::max; ++i){
			if(_parsed[i] && (!_fileData || !_fileData->hasInfo(i))){
				header.push_back(infos[i].argument);
			}
		}

		//open file
		coretools::TOutputFile out(Filename, header);

		// 2) figure out what to write from where
		std::vector<bool> colFromFile(header.size());
		std::vector<InfoType> colInfoType(header.size());
		std::vector<size_t> colInFile(header.size());

		for(size_t c = 0; c < header.size(); ++c){
			auto i = impl::argument2InfoType(header[c]);
			if(i != InfoType::max && _parsed[i]){
				colFromFile[c] = false;
				colInfoType[c] = i;
			} else {
				colFromFile[c] = true;
				colInFile[c] = std::find(fileHeader.begin(), fileHeader.end(), header[c]) - fileHeader.begin();
			}
		}

		// 3) write file
		for(auto rg : _info){
			//find RG info in file
			size_t row = _fileData->getRow(rg[InfoType::RGName]);
			for(size_t c = 0; c < header.size(); ++c){
				if(colFromFile[c]){
					out << (*_fileData)[row][colInFile[c]];
				} else {
					out << rg[colInfoType[c]];
				}
			}
			out << std::endl;
		}
	} else {
		//no file was written: use order of enum
		std::vector<std::string> header;
		std::vector<InfoType> colInfoType;

		for(auto i = InfoType::min; i < InfoType::max; ++i){
			if(_parsed[i]){
				header.push_back(infos[i].argument);
				colInfoType.push_back(i);
			}
		}

		//open file
		coretools::TOutputFile out(Filename, header);

		//and write
		for(auto rg : _info){
			for(size_t c = 0; c < header.size(); ++c){
				out << rg[colInfoType[c]];
			}
			out << std::endl;
		}
	}
}

} //end namespace RGInfo
} //end namespace BAM




