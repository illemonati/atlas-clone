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
		std::vector< std::vector<std::string> > _fileData;
		std::vector<std::string> _header;

	public:

		TFileData(const std::string & Filename){
			//open RG file
			logfile().listFlush("Reading read group info from file '" + Filename + "' ... ");
			coretools::TInputFile in(Filename, coretools::header, "\t", "//");
			header = in.header();

			//extract RG column
			//check that file has a column named "readGroup"
			if(!in.hasColname(infoType2ArgumentString(InfoType::RGName))){
				UERROR("Column '", infoType2ArgumentString(InfoType::RGName), "' missing in file '", Filename, "'!");
			}
			auto rgCol = in.getIndexOfColname(infoType2ArgumentString(InfoType::RGName));

			//read file and create read group entries
			std::vector<std::string> tmp;
			std::vector< std::vector<std::string> > fileData;
			std::vector<std::string> rgNames;

			while (in.read(tmp)) {
				fileData.push_back(tmp);
				rgNames.push_back(tmp[rgCol]);
			}
			logfile().done();
			logfile().conclude("Found ", _readGroups.size(), " read groups.");

			// ensure names are unique
			sort(rgNames.begin(), rgNames.end());
			if(std::adjacent_find(rgNames.cbegin(), rgNames.cend()) != rgNames.cend()){
				DEVERROR("Duplicate read group names in file '", Filename, "'!");
			}
		}

		const std::vector<std::string> header() const { return _header; }

		bool hasInfo(InfoType Info){
			auto it = find(_header.cbegin(), _header.cend(), infoType2ArgumentString(Info));
			return it != header.cend();
		}

		size_t getInfoCol(InfoType Info){
			auto it = find(_header.cbegin(), _header.cend(), infoType2ArgumentString(Info));
			if(it != header.cend()){
				return it - _header.cbegin();
			}
			DEVERROR("Info '", infoType2ArgumentString(Info), "' not present in read group info file!");
		}



		const std::vector<std::string> operator[](size_t row) const {
			return _fileData[row];
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
		logfile().list("Initializing ", infoType2Description(Info), " with default value '", infoType2DefaultValue(Info), "' for all read groups. (set with '", infoType2ArgumentString(Info), "' or '", _RGInfoArgument, "')");
		_setAllReadGroups(Vec, Info, infoType2DefaultValue(Info));
	}

	void setFromCommandLine(InfoVec & Vec, InfoType Info){
		//read from command line
		std::string arg = infoType2ArgumentString(Info);
		std::string val = parameters().getParameter<std::string>(arg);
		logfile().list("Initializing ", infoType2Description(Info), " with '", val, "' for all read groups. (argument '", arg, "')");
		_setAllReadGroups(Vec, Info, val);
	}

	void setFromFile(InfoVec & Vec, InfoType Info, const TFileData & FileData){
	 	 //present in file -> read for each read group
		logfile().list("Initializing ", infoType2Description(Info), " from read group info file. (overwrite with '", infoType2ArgumentString(Info), "')");
		size_t col = FileData.getInfoCol(Info);


			Vec[r].set(Info, FileData[r][col]);
			++rg;
		}
	}


	void setPerReadGroup(InfoVec Vec, InfoType Info, const TFileData & FileData){
		//check if it is provided on the command line -> overwrites file
		std::string arg = infoType2ArgumentString(Info);
		if(parameters().parameterExists(arg)){
			setFromCommandLine(Info);
		} else {
			//check if provided in file
			if(FileData.hasInfo(Info)){
				_readInfoFromFile(Info, fileData, it-header.cbegin());
			} else {
				setDefault(Vec, Info);
			}
		}
	}

}


void TReadGroupInfo::_readInfoFromRGInfoFile(InfoType Info, const std::vector< std::vector<std::string> > & fileData, size_t col){
 	 //present in RG info file -> read for each read group
	logfile().list("Initializing ", infoType2Description(Info), " from read group info file. (overwrite with '", infoType2ArgumentString(Info), "')");
	size_t rg = 0;
	for(auto& row : fileData){
		_info[rg].set(Info, row[col]);
		++rg;
	}
}

void TReadGroupInfo::_readInfoPerReadGroup(InfoType Info){
	//check if it is provded on the command line
	if(parameters().parameterExists(infoType2ArgumentString(Info))){
		_readInfoFromCommandLine(Info);
	} else {
		_setDefault(Info, infoType2DefaultValue(Info));
	}
}

void TReadGroupInfo::_readInfoPerReadGroup(InfoType Info, const std::vector<std::string> & header, const std::vector< std::vector<std::string> > & fileData){
	//check if it is provided on the command line -> overwrites RG info file
	std::string arg = infoType2ArgumentString(Info);
	if(parameters().parameterExists(arg)){
		_readInfoFromCommandLine(Info);
	} else {
		//check if provided in file
		auto it = find(header.cbegin(), header.cend(), arg);
		if(it != header.cend()){
			_readInfoFromFile(Info, fileData, it-header.cbegin());
		} else {
			_setDefault(Info, Default);
		}
	}
}
/*
void TReadGroupInfo::_readInfoFile(){
	//open RG file
	impl::TD
	std::string Filename = parameters().getParameter<std::string>(_RGInfoArgument);
	logfile().listFlush("Reading read group info from file '" + Filename + "' ... ");
	coretools::TInputFile in(Filename, coretools::header, "\t", "//");

	//extract RG column
	//check that file has a column names "readGroup"
	if(!in.hasColname(infoType2ArgumentString(InfoType::RGName))){
		UERROR("Column '", infoType2ArgumentString(InfoType::RGName), "' missing in file '", Filename, "'!");
	}
	auto rgCol = in.getIndexOfColname(infoType2ArgumentString(InfoType::RGName));

	//read file and create read group entries
	std::vector<std::string> tmp;
	std::vector< std::vector<std::string> > fileData;
	std::vector<std::string> rgNames;

	while (in.read(tmp)) {
		fileData.push_back(tmp);
		rgNames.push_back(tmp[rgCol]);
	}
	logfile().done();
	logfile().conclude("Found ", _readGroups.size(), " read groups.");

	// ensure names are unique
	sort(rgNames.begin(), rgNames.end());
	if(std::adjacent_find(rgNames.cbegin(), rgNames.cend()) != rgNames.cend()){
		DEVERROR("Duplicate read group names in file '", Filename, "'!");
	}

	// create read groups
	_info.clear();
	_info.reserve(fileData.size());
	for(auto& f : fileData){
		_info.push_back(f[rgCol]);
	}

	// parse other data, either from command line or file
	for(int i = InfoType::seqType; i < InfoType::COUNT; i++){
		_readInfoPerReadGroup(static_cast<InfoType>(i));
	}
}

void TReadGroupInfo::_initializeFromRGInfoFile(){
	//open RG file
	std::string Filename = parameters().getParameter<std::string>(_RGInfoArgument);
	logfile().listFlush("Reading read group info from file '" + Filename + "' ... ");
	coretools::TInputFile in(Filename, coretools::header, "\t", "//");

	//extract RG column
	//check that file has a column names "readGroup"
	if(!in.hasColname(infoType2ArgumentString(InfoType::RGName))){
		UERROR("Column '", infoType2ArgumentString(InfoType::RGName), "' missing in file '", Filename, "'!");
	}
	auto rgCol = in.getIndexOfColname(infoType2ArgumentString(InfoType::RGName));

	//read file and create read group entries
	std::vector<std::string> tmp;
	std::vector< std::vector<std::string> > fileData;
	std::vector<std::string> rgNames;

	while (in.read(tmp)) {
		fileData.push_back(tmp);
		rgNames.push_back(tmp[rgCol]);
	}
	logfile().done();
	logfile().conclude("Found ", _readGroups.size(), " read groups.");

	// ensure names are unique
	sort(rgNames.begin(), rgNames.end());
	if(std::adjacent_find(rgNames.cbegin(), rgNames.cend()) != rgNames.cend()){
		DEVERROR("Duplicate read group names in file '", Filename, "'!");
	}

	// create read groups
	_info.clear();
	_info.reserve(fileData.size());
	for(auto& f : fileData){
		_info.push_back(f[rgCol]);
	}

	// parse other data, either from command line or file
	_readInfoPerReadGroup(InfoType::seqType);
	_readInfoPerReadGroup(InfoType::numCycles);
	_readInfoPerReadGroup(InfoType::fragmentLengthDistr);
	_readInfoPerReadGroup(InfoType::baseQualityDistr);
	_readInfoPerReadGroup(InfoType::mappingQualityDistr);
	_readInfoPerReadGroup(InfoType::softClipDistr);
}
*/

void TReadGroupInfo::_initializeFromCommandLine(){
	//parse number of read groups
	const auto numRG = parameters().getParameterWithDefault<StrictlyPositive<int>>(_numRGArgument, 1);
	if (numRG == 1) {
		logfile().list("Initializing one read group from arguments. (parameter '", _numRGArgument, "')");
	} else {
		logfile().list("Initializing ", numRG, " identical read groups from arguments (parameter '", _numRGArgument, "').");
	}

	//create read groups
	for (int i = 0; i < numRG; ++i) {
		_readGroups.add("SimReadGroup" + coretools::str::toString(i + 1));
	}

	//initialize info
	_readInfoPerReadGroup(InfoType::seqType);
	_readInfoPerReadGroup(InfoType::numCycles);
	_readInfoPerReadGroup(InfoType::fragmentLengthDistr);
	_readInfoPerReadGroup(InfoType::baseQualityDistr);
	_readInfoPerReadGroup(InfoType::mappingQualityDistr);
	_readInfoPerReadGroup(InfoType::softClipDistr);
}


void TReadGroupInfo::_matchReadGroups(BAM::TReadGroups & ReadGroups){
	// Ensure info is in the same oder as in ReadGroups so it can be looked up by read group index
	auto v = _info;
	std::vector<bool> infoUsed(v.size(), false);
	_info.clear();
	_info.resize(ReadGroups.size());
	for(size_t i = 0; i < ReadGroups.size(); ++i){
		//find name
		size_t k = 0;
		for(; k < v.size(); ++k){
			if(v[k].get(InfoType::RGName) == ReadGroups[i].name_ID){
				_info[i] = v[k];
				infoUsed[k] = true;
				break;
			}
		}
		if(k == v.size()){
			UERROR("No info provided for read group '", ReadGroups[i].name_ID, "'!");
		}
	}

	//warn if some info was not used
	if(v.size() > ReadGroups.size()){
		logfile().warning("Additional read groups present in read group info file!");
	}
}

TReadGroupInfo::TSimulatorReadGroupInfo(){
	//initialize from command line parameters
	if(parameters().parameterExists(_RGInfoargument)){
		_initializeFromRGInfoFile();
	} else {
		_initializeFromCommandLine();
	}
};


//**************************************************************************************************33

// either: read info from file and match with TReadGroups (used for analyzes)
void TReadGroupInfo::readInfoAndMatchReadGroups(const BAM::TReadGroups & ReadGroups){
	if(!_info.empty()){
		DEVERROR("Read group info already read!")
	}

	//create read group info entries
	_info.reserve(ReadGroups.size());
	for(auto& r : ReadGroups){
		_info.push_back(r.name_ID);
	}

	//read file, if provided



}

// or: read info and fill TReadGroups (used for simulations)
void TReadGroupInfo::readInfoAndCreateReadGroups(BAM::TReadGroups & ReadGroups){
	if(!_info.empty()){
		DEVERROR("Read group info already read!")
	}

	// Info is provided as a) a RG info file OR b) as the number of read groups and default arguments
	if(parameters().parameterExists(_RGInfoArgument)){
		// create RGs from RG info file
		// read RG info file
		impl::TFileData(parameters().getParameter<std::string>(_RGInfoArgument));

		//create read groups


	} else {
		//
	}





}



void TReadGroupInfo::createReadGroups(BAM::TReadGroups & ReadGroups) const {
	//create read groups from info entries
	_readGroups.clear();
	for(auto& i : _info){
		_readGroups.add(i.get(InfoType::RGName));
	}
}

void TReadGroupInfo::matchReadGroups(BAM::TReadGroups & ReadGroups){
	// Ensure info is in the same oder as in ReadGroups so it can be looked up by read group index
	auto v = _info;
	std::vector<bool> infoUsed(v.size(), false);
	_info.clear();
	_info.resize(ReadGroups.size());
	for(size_t i = 0; i < ReadGroups.size(); ++i){
		//find name
		size_t k = 0;
		for(; k < v.size(); ++k){
			if(v[k].get(InfoType::RGName) == ReadGroups[i].name_ID){
				_info[i] = v[k];
				infoUsed[k] = true;
				break;
			}
		}
		if(k == v.size()){
			UERROR("No info provided for read group '", ReadGroups[i].name_ID, "'!");
		}
	}

	//warn if some info was not used
	if(v.size() > ReadGroups.size()){
		logfile().warning("Additional read groups present in read group info file!");
	}
}

} //end namespace RGInfo
} //end namespace BAM




