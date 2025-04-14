/*
 * TSamHeader.h
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#ifndef BAM_TSAMHEADER_H_
#define BAM_TSAMHEADER_H_

#include <set>
#include <string>
#include <vector>

namespace BAM { class TReadGroups; }
namespace genometools { class TChromosomes; }

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
struct TSamHeader_HD{
private:
	std::string _version_VN    = "1.6";
	std::string _sortOrder_SO  = "unknown";
	std::string _grouping_GO   = "none";
	std::string _subSorting_SS = "";

public:
	void setVersion(std::string_view Version){ _version_VN = Version; };
	void setSortOrder(std::string_view SortOrder);
	void setGrouping(std::string_view Grouping);
	void setSubSorting(std::string_view SubSorting){ _subSorting_SS = SubSorting; };

	//getters
	const std::string& version() const noexcept { return _version_VN; };
	const std::string& sortOrder() const noexcept { return _sortOrder_SO; };
	const std::string& grouping() const noexcept { return _grouping_GO; };
	const std::string& subSorting() const noexcept { return _subSorting_SS; };

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
	TSamProgram(std::string_view ID, std::string_view Name, std::string_view CommandLine, std::string_view Description,
				std::string_view Version);
	void addPrevious(const TSamProgram &Previous) const;
	void addNext(const TSamProgram &Next) const;
	std::string id() const{ return _ID; };
	std::string compileSamHeader() const;

	bool operator<(const TSamProgram & other) const{
		return this->_ID < other._ID;
	};
	bool operator<(std::string_view other) const{
		return this->_ID < other;
	};
};

bool operator<(std::string_view left, const TSamProgram & right);

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
	TSamHeader() = default;
	TSamHeader(std::string_view Version, std::string_view SortOrder, std::string_view Grouping,
			   std::string_view SubSorting = "") {
		set(Version, SortOrder, Grouping, SubSorting);
	};

	//add info
	void set(std::string_view Version, std::string_view SortOrder, std::string_view Grouping, std::string_view SubSorting="");
	void addProgram(std::string_view ID, std::string_view Name, std::string_view CommandLine, std::string_view Description, std::string_view Version);
	void addPreviousProgramInChain(std::string_view ID, std::string_view previousID);
	void addComment(std::string_view Comment){ _comments_CO.emplace_back(Comment); };
	std::string compileSamHeader(const TReadGroups & ReadGroups) const;
	std::string compileSamHeader(const TReadGroups & ReadGroups, const genometools::TChromosomes & Chromosomes) const;

    //getters
    const std::string& version() const{ return _HD.version(); };
    const std::string& sortOrder() const{ return _HD.sortOrder(); };
    const std::string& grouping() const{ return _HD.grouping(); };
    const std::string& subSorting() const{ return _HD.subSorting(); };
};

}; //end namespace


#endif /* BAM_TSAMHEADER_H_ */
