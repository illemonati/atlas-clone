/*
 * TSamHeader.cpp
 *
 *  Created on: Jun 9, 2020
 *      Author: wegmannd
 */

#include "TSamHeader.h"

#include <utility>

#include "genometools/GenomePositions/TChromosomes.h"
#include "TReadGroups.h"

namespace BAM{

//---------------------------------
// TSamHeader_HD
// Stores the HD line
//---------------------------------
void TSamHeader_HD::setSortOrder(std::string_view SortOrder){
	//check if valid
	if(SortOrder.empty()){
		_sortOrder_SO = "unknown";
	} else if(SortOrder == "unknown" || SortOrder == "unsorted" || SortOrder == "queryname" || SortOrder == "coordinate"){
		_sortOrder_SO = SortOrder;
	} else {
		UERROR("Unknow BAM sort order '", SortOrder, "'! Must be either 'unknown', 'unsorted', 'queryname' or 'coordinate'.");
	}
}

void TSamHeader_HD::setGrouping(std::string_view Grouping){
	//check if valid
	if(Grouping.empty()){
		_grouping_GO = "none";
	} else if(Grouping == "none" || Grouping == "query" || Grouping == "reference"){
		_grouping_GO = Grouping;
	} else {
		UERROR("Unknow BAM grouping '", Grouping, "'! Must be either 'none', 'query', or 'reference'.");
	}
}

std::string TSamHeader_HD::compileSamHeader() const {
	std::string header = "@HD\tVN:" + _version_VN;
	if (!_sortOrder_SO.empty()) { header += "\tSO:" + _sortOrder_SO; }
	if(!_grouping_GO.empty()){
		header += "\tGO:" + _grouping_GO;
	}
	if(!_subSorting_SS.empty()){
		header += "\tSS:" + _subSorting_SS;
	}
	return header + "\n";
}

//---------------------------------
// TSamProgram
// Stores programs used
//---------------------------------
TSamProgram::TSamProgram(std::string_view ID, std::string_view Name, std::string_view CommandLine, std::string_view Description, std::string_view Version){
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
}

bool operator<(std::string_view left, const TSamProgram & right){
	return left < right.id();
}

//---------------------------------
// TSamHeader
// main class to interact with
//---------------------------------
void TSamHeader::set(std::string_view Version,
			 	 	 std::string_view SortOrder,
					 std::string_view Grouping,
					 std::string_view SubSorting){
	_HD.setVersion(Version);
	_HD.setSortOrder(SortOrder);
	_HD.setGrouping(Grouping);
	_HD.setSubSorting(SubSorting);
};

void TSamHeader::addProgram(std::string_view ID, std::string_view Name, std::string_view CommandLine, std::string_view Description, std::string_view Version){
	auto it = _programs_PG.emplace(ID, Name, CommandLine, Description, Version);
	if(!it.second){
		UERROR("Failed to add program to BAM header: duplicate ID '", ID, "'!");
	}
};

void TSamHeader::addPreviousProgramInChain(std::string_view ID, std::string_view previousID){
	//search if previousID exists
	auto prev = _programs_PG.find(previousID);
	if(prev == _programs_PG.end()){
		UERROR("Failed ot add PP tag to BAM header: previous ID '", previousID, "' does not exists!");
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

std::string TSamHeader::compileSamHeader(const TReadGroups & ReadGroups, const genometools::TChromosomes & Chromosomes) const{
	std::string header = compileSamHeader(ReadGroups);

	return header + Chromosomes.compileSamHeader();
};

}; //end namespace


