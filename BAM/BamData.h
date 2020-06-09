/*
 * TSamFlags.h
 *
 *  Created on: May 20, 2020
 *      Author: phaentu
 */

#ifndef BAM_BAMDATA_H_
#define BAM_BAMDATA_H_

#include <bitset>
#include <vector>
#include <set>
#include "TReadGroups.h"
#include "TChromosomes.h"

//-----------------------------------------------------
//TSamFlags
//A class to store, access and manipulate SAM flags.
//SAM flags are bitwise flags and this class has the only purpose to provided named access
//-----------------------------------------------------

class TSamFlags{
private:
	std::bitset<16> flags;

public:
	TSamFlags(const uint16_t Flags){
		flags = std::bitset<16>(Flags);
	};

	//getters
	uint16_t asInt() const{ return flags.to_ulong(); };
	bool operator[](const uint16_t index){
		return flags[index];
	};
	bool isPaired() const{ return flags[0]; };
	bool isProperPair() const{ return flags[1]; };
	bool isUnmapped() const{ return flags[2]; };
	bool mateIsUnmapped() const{ return flags[3]; };

	bool isReverseStrand() const{ return flags[4]; };
	bool mateIsReverseStrand() const{ return flags[5]; };
	bool isFirstMate() const{ return flags[6]; }; //first segment in template
	bool isSecondMate() const{ return flags[7]; }; //last segment in template

	bool isSecondary() const{ return flags[8]; };
	bool isQCFailed() const{ return flags[9]; };
	bool isDuplicate() const{ return flags[10]; };
	bool isSupplementary() const{ return flags[11]; };

	//setters
	void set(const uint16_t & Flags){ flags = std::bitset<16>(Flags); };
	void operator=(const uint16_t & Flags){ flags = std::bitset<16>(Flags); };
	void operator=(const TSamFlags & other){ flags = other.flags; };
	void setIsPaired(const bool ok){flags[0] = ok; };
	void setIsProperPair(const bool ok){flags[1] = ok; };
	void setIsUnmapped(const bool ok){flags[2] = ok; };
	void setMateIsUnmapped(const bool ok){flags[3] = ok; };

	void setIsReverseStrand(const bool ok){flags[4] = ok; };
	void setMateIsReverseStrand(const bool ok){flags[5] = ok; };
	void setIsRead1(const bool ok){flags[6] = ok; };
	void setIsRead2(const bool ok){flags[7] = ok; };

	void setIsSecondary(const bool ok){flags[8] = ok; };
	void setIsQCFailed(const bool ok){flags[9] = ok; };
	void setIsDuplicate(const bool ok){flags[10] = ok; };
	void setIsSupplementary(const bool ok){flags[11] = ok; };
};

//-----------------------------------------------------
//TCigar
//A class to store, access and manipulate CIGAR operators
//-----------------------------------------------------
struct CigarOperator {
    char     type;
    uint32_t length;

    CigarOperator(){
    	type = '\0';
    	length = 0;
    };
    CigarOperator(const char & Type, const uint32_t & Length){
    	type = Type;
    	length = Length;
    };
};

class TCigar{
private:
	std::vector<CigarOperator> _cigar;
	uint16_t _lengthAligned;
	uint16_t _lengthInserted;
	uint16_t _lengthDeleted;
	uint16_t _lengthSoftClippedLeft;
	uint16_t _lengthSoftClippedRight;
	bool _addSoftClippedLeft;

public:
	TCigar(){
		clear();
	};

	void clear();
	std::vector<CigarOperator>::const_iterator begin() const { return _cigar.begin(); }
	std::vector<CigarOperator>::const_iterator end() const { return _cigar.end(); }
    void add(const char & Type, const uint32_t & Length);
    void removeSoftClips();

    uint16_t lengthAligned() const{ return _lengthAligned; };
    uint16_t lengthInserted() const{ return _lengthInserted; };
    uint16_t lengthDeleted() const{ return _lengthDeleted; };
    uint16_t lengthSoftClippedLeft() const{ return _lengthSoftClippedLeft; };
    uint16_t lengthSoftClippedRight() const{ return _lengthSoftClippedRight; };
    uint16_t lengthSoftClipped() const{ return _lengthSoftClippedLeft + _lengthSoftClippedRight; };
    uint16_t lengthSequenced() const{ return _lengthAligned + _lengthInserted; };
    uint16_t lengthRead() const{ return _lengthAligned + _lengthInserted + _lengthSoftClippedLeft + _lengthSoftClippedRight; };
};

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

	friend bool operator<(const TSamProgram & other) const{
		return this->_ID < other->_ID;
	};
	friend bool operator<(const std::string & other) const{
		return this->_ID < other;
	};
	friend bool operator<(const std::string & left, const TSamProgram & right) const{
	    return left < right._ID;
	};
	friend bool operator<(const TSamProgram & left, const std::string & right) const{
	  	return left._ID < right;
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

#endif /* BAM_BAMDATA_H_ */
