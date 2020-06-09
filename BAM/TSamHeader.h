/*
 * TSamHeader.h
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#ifndef BAM_TSAMHEADER_H_
#define BAM_TSAMHEADER_H_

#include "TReadGroups.h"
#include "TChromosomes.h"

namespace BAM{

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
	void setSortOrder(const std::string SortOrder){
		//check if valid
		if(SortOrder == "unknown" || SortOrder == "unsorted" || SortOrder == "queryname" || SortOrder == "coordinate"){
			_sortOrder_SO = SortOrder;
		} else {
			throw "Unknow BAM sort order '" + SortOrder + "'! Must be either 'unknown', 'unsorted', 'queryname' or 'coordinate'.";
		}
	};
	void setGrouping(const std::string Grouping){
		//check if valid
		if(Grouping == "none" || Grouping == "query" || Grouping == "reference"){
			_grouping_GO = Grouping;
		} else {
			throw "Unknow BAM grouping '" + Grouping + "'! Must be either 'none', 'query', or 'reference'.";
		}
	};
	void setSubSorting(const std::string SubSorting){
		_subSorting_SS = SubSorting;
	};

	//getters
	std::string version() const{ return _version_VN; };
	std::string sortOrder() const{ return _sortOrder_SO; };
	std::string grouping() const{ return _grouping_GO; };
	std::string subSorting() const{ return _subSorting_SS; };
};

class TSamHeader{
private:
	TSamHeader_HD _HD;

	TReadGroups& _readGroups_RG;
	TChromosomes& _chromosomes_SQ;

	//programs and comments are currently nor parsed
	//TODO: add parsing of these
	std::vector<std::string> _programs_PG;
	std::vector<std::string> _comments_CO;

public:
	TSamHeader(TChromosomes & Chromosomes, TReadGroups & ReadGroups);

	//add info
	void set(const std::string Version,
			 const std::string SortOrder,
			 const std::string Grouping,
			 const std::string SubSorting){
		_HD.setVersion(Version);
		_HD.setSortOrder(SortOrder);
		_HD.setGrouping(Grouping);
		_HD.setSubSorting(SubSorting);
	};

	void addProgram(const std::string Program){
		_programs_PG.push_back(Program);
	};

	void addComment(const std::string Comment){
		_comments_CO.push_back(Comment);
	};

};





#endif /* BAM_TSAMHEADER_H_ */
