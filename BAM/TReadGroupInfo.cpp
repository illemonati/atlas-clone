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

namespace impl{
	//------------------------------------------------
	// class / functions to initialize: only visible in cpp file
	//------------------------------------------------

	class TFileData{
	private:
		std::string _filename;
		std::vector< std::vector<std::string> > _fileData;
		std::vector<std::string> _header;
		std::vector<std::string> _rgNames;

	public:

		TFileData(const std::string & Filename){
			//open RG file
			_filename = Filename;
			logfile().listFlush("Reading read group info from file '" + _filename + "' ... ");
			coretools::TInputFile in(Filename, coretools::TFile_Filetype::header, "\t", "//");
			_header = in.header();

			//extract RG column
			//check that file has a column named "readGroup"
			if(!in.hasColname(infoType2ArgumentString(InfoType::RGName))){
				UERROR("Column '", infoType2ArgumentString(InfoType::RGName), "' missing in file '", Filename, "'!");
			}
			auto rgCol = in.getIndexOfColname(infoType2ArgumentString(InfoType::RGName));

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

		const std::vector<std::string> header() const { return _header; }

		bool hasInfo(InfoType Info) const noexcept {
			auto it = find(_header.cbegin(), _header.cend(), infoType2ArgumentString(Info));
			return it != _header.cend();
		}

		size_t getInfoCol(InfoType Info) const {
			auto it = find(_header.cbegin(), _header.cend(), infoType2ArgumentString(Info));
			if(it != _header.cend()){
				return it - _header.cbegin();
			}
			DEVERROR("Info '", infoType2ArgumentString(Info), "' not present in read group info file!");
		}

		size_t size() const noexcept {
			return _fileData.size();
		}

		const std::vector<std::string> operator[](size_t row) const noexcept {
			return _fileData[row];
		}

		size_t getRow(const std::string & ReadGroupName) const {
			auto it = std::find(_rgNames.cbegin(), _rgNames.cend(), ReadGroupName);
			if(it == _rgNames.cend()){
				UERROR("Read group '", ReadGroupName, "' missing in rea dgroup info file '", _filename, "'!");
			}
			return it - _rgNames.cbegin();
		}
	};

	//------------------------------------------------
	// Functions to set info
	//------------------------------------------------

	using InfoVec = std::vector<TReadGroupInfoEntry>;

	void setAllReadGroups(InfoVec & Vec, InfoType Info, const std::string & Val){
		for(auto& i : Vec){
			i.set(Info, Val);
		}
	}

	void setDefault(InfoVec & Vec, InfoType Info){
		//use default values
		logfile().list("Initializing ",
				       infoType2Description(Info),
					   " with default value '",
					   infoType2DefaultValue(Info),
					   "' for all read groups. (set with '",
					   TReadGroupInfo::_RGInfoArgument,
					   "' or '",
					   infoType2ArgumentString(Info),
					   "')");
		setAllReadGroups(Vec, Info, infoType2DefaultValue(Info));
	}

	void setFromCommandLine(InfoVec & Vec, InfoType Info){
		//read from command line
		std::string arg = infoType2ArgumentString(Info);
		std::string val = parameters().getParameter<std::string>(arg);
		logfile().list("Initializing ", infoType2Description(Info), " with '", val, "' for all read groups. (argument '", arg, "')");
		setAllReadGroups(Vec, Info, val);
	}

	void setFromRGInfoFile(InfoVec & Vec, InfoType Info, const TFileData & FileData){
	 	 //present in file -> read for each read group
		logfile().list("Initializing ", infoType2Description(Info), " from read group info file. (overwrite with '", infoType2ArgumentString(Info), "')");
		auto col = FileData.getInfoCol(Info);
		for(auto& r : Vec){
			auto row = FileData.getRow(r.get(InfoType::RGName));
			r.set(Info, FileData[row][col]);
		}
	}

	void setPerReadGroup(InfoVec Vec, InfoType Info, const TFileData & FileData){
		//check if info is provided on the command line -> overwrites file
		std::string arg = infoType2ArgumentString(Info);
		if(parameters().parameterExists(arg)){
			setFromCommandLine(Vec, Info);
		} else {
			//check if provided in file
			if(FileData.hasInfo(Info)){
				setFromRGInfoFile(Vec, Info, FileData);
			} else {
				setDefault(Vec, Info);
			}
		}
	}

	void setPerReadGroup(InfoVec Vec, InfoType Info){
		//check if info is provided on the command line -> overwrites file
		std::string arg = infoType2ArgumentString(Info);
		if(parameters().parameterExists(arg)){
			setFromCommandLine(Vec, Info);
		} else {
			setDefault(Vec, Info);
		}
	}

	void setAllPerReadGroup(InfoVec Vec, const TFileData & FileData){
		for(int i = static_cast<int>(InfoType::seqType); i < static_cast<int>(InfoType::COUNT); i++){
			setPerReadGroup(Vec, static_cast<InfoType>(i), FileData);
		}
	}

	void setAllPerReadGroup(InfoVec Vec){
		for(int i = static_cast<int>(InfoType::seqType); i < static_cast<int>(InfoType::COUNT); i++){
			setPerReadGroup(Vec, static_cast<InfoType>(i));
		}
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
	for(auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r){
		_info.push_back(r->name_ID);
	}

	//read file, if provided
	if(parameters().parameterExists(_RGInfoArgument)){
		impl::TFileData data(parameters().getParameter<std::string>(_RGInfoArgument));

		//initialize from RG file or command line
		impl::setAllPerReadGroup(_info, data);
	} else {
		impl::setAllPerReadGroup(_info);
	}
}

// or: read info and fill TReadGroups (used for simulations)
void TReadGroupInfo::readInfoAndCreateReadGroups(BAM::TReadGroups & ReadGroups){
	if(!_info.empty()){
		DEVERROR("Read group info already read!");
	}

	//empty read groups
	ReadGroups.clear();

	// Info is provided as a) a RG info file OR b) as the number of read groups and default arguments
	if(parameters().parameterExists(_RGInfoArgument)){
		// create RGs from RG info file
		// read RG info file
		impl::TFileData data(parameters().getParameter<std::string>(_RGInfoArgument));

		//create read groups
		auto col = data.getInfoCol(InfoType::RGName);
		for(auto i = 0; i < data.size(); i++){
			ReadGroups.add(data[i][col]);
		}

		//initialize from RG file or command line
		impl::setAllPerReadGroup(_info, data);
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
			ReadGroups.add("SimReadGroup" + coretools::str::toString(i + 1));
		}

		//initialize from command line or with default
		impl::setAllPerReadGroup(_info);
	}
}

} //end namespace RGInfo
} //end namespace BAM




