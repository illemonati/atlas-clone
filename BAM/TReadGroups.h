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
//TReadGroupMaxLength
// ToDo: remove?
//---------------------------------------------------------------
struct TReadGroupMaxLength{
public:
	int maxLen;
	uint16_t truncatedOrMergedReadGroupID;
	std::string truncatedOrMergedReadGroup;
	int sequencingType; //0 = single, 1 = mixed, 2 = paired

	TReadGroupMaxLength(int MaxLen, int TruncatedOrMergedReadGroupID, std::string_view TruncatedOrMergedReadGroup, int Type){
		maxLen = MaxLen;
		truncatedOrMergedReadGroupID = TruncatedOrMergedReadGroupID;
		truncatedOrMergedReadGroup = TruncatedOrMergedReadGroup;
		sequencingType = Type;
	};
};

//---------------------------------------------------------------
//TReadGroup
//---------------------------------------------------------------
class TReadGroup{
private:
	uint16_t _id; //internal ID

public:
	std::string name_ID;
	mutable std::string barcodeSequence_BC;
	mutable std::string sequencingCenter_CN;
	mutable std::string description_DS;
	mutable std::string productionDate_DT;
	mutable std::string flowOrder_FO;
	mutable std::string keySequence_KS;
	mutable std::string library_LB;
	mutable std::string program_PG;
	mutable std::string predictedInsertSize_PI;
	mutable std::string sequencingTechnology_PL;
	mutable std::string platformModel_PM;
	mutable std::string platformUnit_PU;
	mutable std::string sample_SM;

    //flags
    mutable bool inUse; 						//read groups not in use are ignored when reading
    mutable bool writeToHeader;                 //is false if read group is not in use or replaced by new one

    TReadGroup(const uint16_t ID, std::string_view Name);
    TReadGroup(const TReadGroup & other) = default;
    TReadGroup* getPointer(){ return this; };
    std::string compileSamHeader() const;

    //getters
    uint16_t id() const { return _id; };
    void setId(const uint16_t id){ _id = id; };

    bool operator<(const TReadGroup & right) const;
    bool operator<(std::string_view right) const;
};

bool operator<(std::string_view left, const TReadGroup & right);

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------
class TReadGroups{
private:
	std::set<TReadGroup, std::less<>> _readGroups;
	std::vector<const TReadGroup*> _readGroupsById;
	bool _limitReadGroups;

	void _fillLookupFromId();

public:
	TReadGroups();
	~TReadGroups(){};

	TReadGroups(const TReadGroups && other);
	TReadGroups(const TReadGroups & other);
	TReadGroups& operator=(const TReadGroups & other);
	TReadGroups& operator=(const TReadGroups && other);

	void clear();
	const TReadGroup& add(std::string_view name);
	const TReadGroup& addAlternativeRG(std::string_view Name, std::string_view original);
	uint16_t size() const;
	bool empty() const;

	uint16_t getId(std::string_view name) const;
	const std::string& getName (uint16_t readGroupId) const;
	std::vector<std::string> getNames(std::vector<uint16_t> & readGroupIds) const;
	const TReadGroup& getReadGroup(std::string_view name);
	const TReadGroup& operator[](uint16_t readGroupId) const;
	bool readGroupExists(std::string_view name) const;
	bool readGroupInUse(uint16_t readGroupId) const;
	bool readGroupInUse(std::string_view name) const;

	//looping over
	std::set<TReadGroup, std::less<>>::iterator begin(){ return _readGroups.begin(); };
	std::set<TReadGroup, std::less<>>::iterator end(){ return _readGroups.end(); };
	std::set<TReadGroup, std::less<>>::iterator cbegin() const{ return _readGroups.cbegin(); };
	std::set<TReadGroup, std::less<>>::iterator cend() const{ return _readGroups.cend(); };

	void filterReadGroups(std::string_view readGroupList);
	void removeFromHeader(std::string_view name);
	void removeFromHeader(const uint16_t readGroupId);
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
	static const uint16_t ReadGroupMapNotInitializedIndex; //largest possible values
	std::vector<uint16_t> _readGroupMap; //maps read group index to pooled index
	std::vector< std::vector<uint16_t> > _reverseReadGroupMap; //IDs of all Rgs pooled with that index. Includes itself.
	std::vector<uint16_t> _readGroupsInUse;

	void _resize(const TReadGroups & ReadGroups);
	void _markAsInUse(uint16_t index);
	void _fillWithoutPooling(const TReadGroups & ReadGroups);
	void _fillFromFile(const TReadGroups & ReadGroups, std::string_view filename);

public:
	TReadGroupMap(const TReadGroups & ReadGroups);
	TReadGroupMap(const TReadGroups & ReadGroups, std::string_view filename);

	~TReadGroupMap() = default;

	uint16_t size() const { return _readGroupMap.size(); };
	uint16_t numReadGroupsInUse() const { return _readGroupsInUse.size(); };

	uint16_t operator[](uint16_t rg) const { return _readGroupMap[rg]; };
	uint16_t pooledIndex(uint16_t rg) const { return _readGroupMap[rg]; };
	bool inUse(uint16_t rg) const { return _readGroupMap[rg] == rg; };

	const std::vector<uint16_t>& readGroupsInUse() const { return _readGroupsInUse; };
	const std::vector<uint16_t>& readGroupsPooledWith(uint16_t rg) const { return _reverseReadGroupMap[rg]; };
};

}; //end namespace

#endif /* TREADGROUPS_H_ */
