/*
 * TReadGroups.h
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

#ifndef TREADGROUPS_H_
#define TREADGROUPS_H_

#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace BAM{

//---------------------------------------------------------------
//TReadGroup
//---------------------------------------------------------------
struct TReadGroup {
public:
	size_t id; //internal ID
	std::string name_ID;
	std::string barcodeSequence_BC;
	std::string sequencingCenter_CN;
	std::string description_DS;
	std::string productionDate_DT;
	std::string flowOrder_FO;
	std::string keySequence_KS;
	std::string library_LB;
	std::string program_PG;
	std::string predictedInsertSize_PI;
	std::string sequencingTechnology_PL;
	std::string platformModel_PM;
	std::string platformUnit_PU;
	std::string sample_SM;

    //flags
    bool inUse; 						//read groups not in use are ignored when reading
    bool writeToHeader;                 //is false if read group is not in use or replaced by new one

    TReadGroup();
    TReadGroup(size_t ID, std::string_view Name);

    std::string compileSamHeader() const;

    bool operator<(const TReadGroup & right) const;
    bool operator<(std::string_view right) const;
    bool operator==(std::string_view name) const;
};

bool operator<(std::string_view left, const TReadGroup & right);

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------
class TReadGroups{
private:
	const TReadGroup _noReadGroup;
	std::vector<TReadGroup> _readGroups;
	std::vector<size_t> _readGroupsById;
	bool _limitReadGroups;

	std::vector<TReadGroup>::iterator _getReadGroup(std::string_view Name);
	std::vector<TReadGroup>::const_iterator _getReadGroup(std::string_view Name) const;
	void _fillLookupFromId();

public:
	TReadGroups();
	~TReadGroups(){};

	TReadGroups(const TReadGroups && other);
	TReadGroups(const TReadGroups & other);
	TReadGroups& operator=(const TReadGroups & other);
	TReadGroups& operator=(const TReadGroups && other);

	static constexpr size_t noReadGroupId = -1;

	// add and remove read groups
	void clear();
	TReadGroup& add(std::string_view name);
	TReadGroup& addAlternativeRG(std::string_view Name, std::string_view original);
	size_t size() const;
	bool empty() const;

	// access read groups
	size_t getId(std::string_view name) const; // returns noReadGroupId if read group does not exist.
	const TReadGroup& getReadGroup(std::string_view name) const;
	TReadGroup& getReadGroup(std::string_view name);
	const TReadGroup& getReadGroup(size_t ReadGroupId) const;
	TReadGroup& getReadGroup(size_t ReadGroupId);
	const TReadGroup& operator[](size_t readGroupId) const; //no checking

	bool readGroupExists(std::string_view name) const;
	bool readGroupExists(size_t readGroupId) const;

	//looping over
	std::vector<TReadGroup>::iterator begin(){ return _readGroups.begin(); };
	std::vector<TReadGroup>::iterator end(){ return _readGroups.end(); };
	std::vector<TReadGroup>::const_iterator cbegin() const{ return _readGroups.cbegin(); };
	std::vector<TReadGroup>::const_iterator cend() const{ return _readGroups.cend(); };

	//getters of specific entries
	template <typename T> bool readGroupInUse(T Identifier) const {
		return getReadGroup(Identifier).inUse;
	}
	const std::string& getName (size_t readGroupId) const;
	std::vector<std::string> getNames(std::vector<size_t> & readGroupIds) const;

	//some additional tasks
	void filterReadGroups(std::string_view readGroupList);
	template <typename T> void removeFromHeader(T Identifier){
		auto rg = getReadGroup(Identifier);
		rg.writeToHeader = false;
	}

	void printReadgroupsInUse() const;
	void fillVectorWithNames(std::vector<std::string> & vec) const;
	std::string compileSamHeader() const;
};

//--------------------------------------------------------------------------------------
//TReadGroupMap
//Maps bam file read group index to internal index, which may differ in case of pooling
//--------------------------------------------------------------------------------------
class TReadGroupMap{
private:
	static const size_t ReadGroupMapNotInitializedIndex; //largest possible values
	std::vector<size_t> _readGroupMap; //maps read group index to pooled index
	std::vector< std::vector<size_t> > _reverseReadGroupMap; //IDs of all Rgs pooled with that index. Includes itself.
	std::vector<size_t> _readGroupsInUse;

	void _resize(const TReadGroups & ReadGroups);
	void _markAsInUse(size_t index);
	void _fillWithoutPooling(const TReadGroups & ReadGroups);
	void _fillFromFile(const TReadGroups & ReadGroups, std::string_view filename);

public:
	TReadGroupMap(const TReadGroups & ReadGroups, std::string_view filename = "");

	size_t size() const { return _readGroupMap.size(); };
	size_t numReadGroupsInUse() const { return _readGroupsInUse.size(); };

	size_t operator[](size_t rg) const { return _readGroupMap[rg]; };
	size_t pooledIndex(size_t rg) const { return _readGroupMap[rg]; };
	bool inUse(size_t rg) const { return _readGroupMap[rg] == rg; };

	const std::vector<size_t>& readGroupsInUse() const { return _readGroupsInUse; };
	const std::vector<size_t>& readGroupsPooledWith(size_t rg) const { return _reverseReadGroupMap[rg]; };
};

}; //end namespace

#endif /* TREADGROUPS_H_ */
