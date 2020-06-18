/*
 * TSamHeader.h
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#ifndef BAM_TSAMHEADER_H_
#define BAM_TSAMHEADER_H_

#include "stringFunctions.h"
#include "TReadGroups.h"
#include "TChromosomes.h"

namespace BAM{

//-----------------------------------------------------
//TSamHeader
//A class to store, access and manipulate the SAM header
// does currently not store SG (chromosomes) and RG (ReadGroups) tags. These are store in their own class
//-----------------------------------------------------

//---------------------------------
// TSamHeader_HD
// Stores the HD line
//---------------------------------
class TSamHeader_HD{
private:
	std::string _version_VN;
	std::string _sortOrder_SO;
	std::string _grouping_GO;
	std::string _subSorting_SS;

public:
	TSamHeader_HD(){
		_version_VN = "1.6";
		_sortOrder_SO = "unknown";
		_grouping_GO = "none";
		_subSorting_SS = "";
	};

	void setVersion(const std::string Version){ _version_VN = Version; };
	void setSortOrder(const std::string SortOrder);
	void setGrouping(const std::string Grouping);
	void setSubSorting(const std::string SubSorting){ _subSorting_SS = SubSorting; };

	//getters
	std::string version() const{ return _version_VN; };
	std::string sortOrder() const{ return _sortOrder_SO; };
	std::string grouping() const{ return _grouping_GO; };
	std::string subSorting() const{ return _subSorting_SS; };

	std::string compileSamHeader() const;
};

//---------------------------------
// TSamProgram
// Stores programs used
//---------------------------------
class TSamProgram{
private:
	std::string _ID;
	mutable std::string _name_PN;
	mutable std::string _commandLine_CL;

	mutable bool _hasPrevious;
	mutable const TSamProgram* _previous_PP; //nullpointer indicates
	mutable bool _hasNext;
	mutable const TSamProgram* _next_PP; //nullpointer indicates
	mutable std::string _description_DS;
	mutable std::string _version_VN;

public:
	TSamProgram(const std::string ID, const std::string Name);
	TSamProgram(const std::string ID, const std::string Name, const std::string CommandLine, const std::string Description, const std::string Version);
	void addPrevious(const TSamProgram & Previous) const;
	void addNext(const TSamProgram & Next) const;
	std::string compileSamHeader() const;

	bool operator<(const TSamProgram & other) const{
		return this->_ID < other._ID;
	};
	bool operator<(const std::string & other) const{
		return this->_ID < other;
	};
};

//---------------------------------
// TSamHeader
// main class to interact with
//---------------------------------
class TSamHeader{
private:
	TSamHeader_HD _HD;
	std::set<TSamProgram, std::less<>> _programs_PG;
	std::vector<std::string> _comments_CO;

public:
	TSamHeader(){};
	TSamHeader(const std::string Version, const std::string SortOrder, const std::string Grouping, const std::string SubSorting=""){
		set(Version, SortOrder, Grouping, SubSorting);
	};

	//add info
	void set(const std::string Version, const std::string SortOrder, const std::string Grouping, const std::string SubSorting="");
	void addProgram(const std::string ID, const std::string Name);
	void addProgram(const std::string ID, const std::string Name, const std::string CommandLine, const std::string Description, const std::string Version);
	void addPreviousProgramInChain(const std::string ID, const std::string previousID);
	void addComment(const std::string Comment){ _comments_CO.push_back(Comment); };
	std::string compileSamHeader(const TReadGroups & ReadGroups) const;
	std::string compileSamHeader(const TReadGroups & ReadGroups, const TChromosomes & Chromosomes) const;
};

}; //end namespace


#endif /* BAM_TSAMHEADER_H_ */
