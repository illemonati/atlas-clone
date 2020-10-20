/*
 * TSamHeader.cpp
 *
 *  Created on: Jun 9, 2020
 *      Author: wegmannd
 */

#include "TSamHeader.h"

namespace BAM{

//---------------------------------
// TSamHeader_HD
// Stores the HD line
//---------------------------------
TSamHeader_HD::TSamHeader_HD(){
	_version_VN = "1.6";
	_sortOrder_SO = "unknown";
	_grouping_GO = "none";
	_subSorting_SS = "";
};

void TSamHeader_HD::setSortOrder(const std::string SortOrder){
	//check if valid
	if(SortOrder.empty()){
		_sortOrder_SO = "unknown";
	} else if(SortOrder == "unknown" || SortOrder == "unsorted" || SortOrder == "queryname" || SortOrder == "coordinate"){
		_sortOrder_SO = SortOrder;
	} else {
		throw "Unknow BAM sort order '" + SortOrder + "'! Must be either 'unknown', 'unsorted', 'queryname' or 'coordinate'.";
	}
};

void TSamHeader_HD::setGrouping(const std::string Grouping){
	//check if valid
	if(Grouping.empty()){
		_grouping_GO = "none";
	} else if(Grouping == "none" || Grouping == "query" || Grouping == "reference"){
		_grouping_GO = Grouping;
	} else {
		throw "Unknow BAM grouping '" + Grouping + "'! Must be either 'none', 'query', or 'reference'.";
	}
};

std::string TSamHeader_HD::compileSamHeader() const{
	std::string header = "@HD\tVN:" + _version_VN;
	if(!_sortOrder_SO.empty()){
		header += "\tSO:" + _sortOrder_SO;
	}
	if(!_grouping_GO.empty()){
		header += "\tGO:" + _grouping_GO;
	}
	if(!_subSorting_SS.empty()){
		header += "\tSS:" + _subSorting_SS;
	}
	return header + "\n";
};

//---------------------------------
// TSamProgram
// Stores programs used
//---------------------------------
TSamProgram::TSamProgram(const std::string ID, const std::string Name){
	_ID = ID;
	_name_PN = Name;
	_hasPrevious = false;
	_previous_PP = nullptr;
	_hasNext = false;
	_next_PP = nullptr;
};

TSamProgram::TSamProgram(const std::string ID, const std::string Name, const std::string CommandLine, const std::string Description, const std::string Version){
	_ID = ID;
	_name_PN = Name;
	_commandLine_CL = CommandLine;
	_description_DS = Description;
	_version_VN = Version;
	_hasPrevious = false;
	_previous_PP = nullptr;
	_hasNext = false;
			_next_PP = nullptr;
};

void TSamProgram::addPrevious(const TSamProgram & Previous) const{
	_previous_PP = &Previous;
	_hasPrevious = true;
};

void TSamProgram::addNext(const TSamProgram & Next) const{
	_next_PP = &Next;
	_hasNext = true;
};

std::string TSamProgram::compileSamHeader() const{
	std::string header = "@PG\tID:" + _ID;

	if(!_name_PN.empty()){
		header += "\tPN:" + _name_PN;
	}

	if(!_commandLine_CL.empty()){
		header += "\tCL:" + _commandLine_CL;
	}

	if(_hasPrevious){
		header += "\tPP:" + _previous_PP->_ID;
	}

	if(!_description_DS.empty()){
		header += "\tDS:" + _description_DS;
	}

	if(!_version_VN.empty()){
		header += "\tVN:" + _version_VN;
	}

	return header + '\n';
};

bool operator<(const std::string & left, const TSamProgram & right){
	return left < right.id();
};

//---------------------------------
// TSamHeader
// main class to interact with
//---------------------------------
void TSamHeader::set(const std::string Version,
			 	 	 const std::string SortOrder,
					 const std::string Grouping,
					 const std::string SubSorting){
	_HD.setVersion(Version);
	_HD.setSortOrder(SortOrder);
	_HD.setGrouping(Grouping);
	_HD.setSubSorting(SubSorting);
};

void TSamHeader::addProgram(const std::string ID, const std::string Name){
	auto it = _programs_PG.emplace(ID, Name);
	if(!it.second){
		throw "Failed to add program to BAM header: duplicate ID '" + ID + "'!";
	}
};

void TSamHeader::addProgram(const std::string ID, const std::string Name, const std::string CommandLine, const std::string Description, const std::string Version){
	auto it = _programs_PG.emplace(ID, Name, CommandLine, Description, Version);
	if(!it.second){
		throw "Failed to add program to BAM header: duplicate ID '" + ID + "'!";
	}
};

void TSamHeader::addPreviousProgramInChain(const std::string ID, const std::string previousID){
	//search if previousID exists
	auto prev = _programs_PG.find(previousID);
	if(prev == _programs_PG.end()){
		throw "Failed ot add PP tag to BAM header: previous ID '" + previousID + "' does not exists!";
	}

	//add
	auto nex = _programs_PG.find(ID);
	nex->addPrevious(*prev);
	prev->addNext(*nex);
};

std::string TSamHeader::compileSamHeader(const TReadGroups & ReadGroups) const{
	//HD line
	std::string header = _HD.compileSamHeader();

	//read groups
	header += ReadGroups.compileSamHeader();

	//programs
	for(auto& p : _programs_PG){
		header += p.compileSamHeader();
	}

	//comments
	for(auto& c : _comments_CO){
		header += "@CO\t" + c + "\n";
	}

	return header;
};

std::string TSamHeader::compileSamHeader(const TReadGroups & ReadGroups, const TChromosomes & Chromosomes) const{
	std::string header = compileSamHeader(ReadGroups);

	return header + Chromosomes.compileSamHeader();
};

}; //end namespace


