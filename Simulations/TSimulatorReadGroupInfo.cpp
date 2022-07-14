/*
 * TSimulatorReadGroupInfo.cpp
 *
 *  Created on: Jul 14, 2022
 *      Author: phaentu
 */


#include "TSimulatorReadGroupInfo.h"
#include "TFile.h"
#include "TError.h"
#include <vector>

using coretools::instances::parameters;
using coretools::instances::logfile;


namespace Simulations {

namespace RGInfo{

void TSimulatorReadGroupInfo::_setAllReadGroups(InfoType Info, const std::string & Val){
	for(auto& i : _info){
		i.set(Info, Val);
	}
}

void TSimulatorReadGroupInfo::_readInfoFromCommandLine(InfoType Info){
	//read from command line
	std::string arg = infoType2ArgumentString(Info);
	std::string val = parameters().getParameter<std::string>(arg);
	logfile().list("Initializing ", infoType2Description(Info), " with '", val, "' for all read groups. (argument '", arg, "')");
	_setAllReadGroups(Info, val);
}

void TSimulatorReadGroupInfo::_readInfoFromFile(InfoType Info, const std::vector< std::vector<std::string> > & fileData, size_t col){
 	 //present in file -> read for each read group
	logfile().list("Initializing ", infoType2Description(Info), " from read group info file. (overwrite with '", infoType2ArgumentString(Info), "')");
	size_t rg = 0;
	for(auto& row : fileData){
		_info[rg].set(Info, row[col]);
		++rg;
	}
}

void TSimulatorReadGroupInfo::_setDefault(InfoType Info){
	//use default values
	logfile().list("Initializing ", infoType2Description(Info), " with default value '", infoType2DefaultValue(Info), "' for all read groups. (set with '", infoType2ArgumentString(Info), "' or '", _RGInfoArgument, "')");
	_setAllReadGroups(Info, infoType2DefaultValue(Info));
}

void TSimulatorReadGroupInfo::_readInfoPerReadGroup(InfoType Info){
	//check if it is provded on the command line
	if(parameters().parameterExists(infoType2ArgumentString(Info))){
		_readInfoFromCommandLine(Info);
	} else {
		_setDefault(Info, infoType2DefaultValue(Info));
	}
}

void TSimulatorReadGroupInfo::_readInfoPerReadGroup(InfoType Info, const std::vector<std::string> & header, const std::vector< std::vector<std::string> > & fileData){
	//check if it is provded on the command line -> overwrites file
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

void TSimulatorReadGroupInfo::_initializeFromRGInfoFile(){
	//open RG file
	std::string Filename = parameters().getParameter<std::string>("readGroupInfo");
	logfile().listFlush("Initializing read groups from file '" + Filename + "' ... ");
	coretools::TInputFile in(Filename, coretools::header, "\t", "//");

	//extract RG column
	//check that file has a column names "readGroup"
	if(!in.hasColname(infoType2ArgumentString(InfoType::RGName))){
		UERROR("Column '", infoType2ArgumentString(InfoType::RGName), "' missing in file '", Filename, "'!");
	}
	auto rgCol = in.getIndexOfColname(infoType2ArgumentString(InfoType::RGName));

	//read file and create read groups
	std::vector<std::string> tmp;
	std::vector< std::vector<std::string> > fileData;
	_readGroups.clear();
	while (in.read(tmp)) {
		_readGroups.add(tmp[rgCol]);
		fileData.push_back(tmp);
	}
	logfile().done();
	logfile().conclude("Found ", _readGroups.size(), " read groups.");

	//initialize info entries
	_info.resize(_readGroups.size());
	for(size_t i = 0; i < _readGroups.size(); ++i){
		_info[i].set(InfoType::RGName, _readGroups[i].name_ID);
	}

	//now parse other data, either from command line or file
	_readInfoPerReadGroup(InfoType::seqType);
	_readInfoPerReadGroup(InfoType::numCycles);
	_readInfoPerReadGroup(InfoType::fragmentLengthDistr);
	_readInfoPerReadGroup(InfoType::baseQualityDistr);
	_readInfoPerReadGroup(InfoType::mappingQualityDistr);
	_readInfoPerReadGroup(InfoType::softClipDistr);
}

void TSimulatorReadGroupInfo::_initializeFromCommandLine(){
	//parse number of read groups
	const auto numRG = parameters().getParameterWithDefault<StrictlyPositive<int>>(_numRGArgument, 1);
	if (numRG == 1) {
		logfile().list("Initializing one read group. (parameter '", _numRGArgument, "')");
	} else {
		logfile().list("Initializing ", numRG, " identical read groups (parameter '", _numRGArgument, "').");
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

TSimulatorReadGroupInfo::TSimulatorReadGroupInfo(){
	//initialize from command line parameters
	if(parameters().parameterExists("readGroupInfo")){
		_initializeFromRGInfoFile();
	} else {
		_initializeFromCommandLine();
	}
};


} //end namespace RGInfo
} //end namespace Simulations




