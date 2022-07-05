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

namespace coretools { class TLog; }

namespace BAM{

//------------------------------------------------
// ReadGroupType
//------------------------------------------------
enum ReadGroupType : uint8_t { unchanged=0, single, mixed, paired};

constexpr char readGroupType2String(ReadGroupType type) noexcept {
	static std::array<std::string, 4> name = {"unchanged", "single-end", "mixed", "paired-end"};
	return name[type];
};

//---------------------------------------------------------------
//TReadGroupMaxLength
//---------------------------------------------------------------
struct TReadGroupMaxLength{
public:
	int maxLen;
	uint16_t truncatedOrMergedReadGroupID;
	std::string truncatedOrMergedReadGroup;
	int sequencingType; //0 = single, 1 = mixed, 2 = paired

	TReadGroupMaxLength(int MaxLen, int TruncatedOrMergedReadGroupID, std::string & TruncatedOrMergedReadGroup, int Type){
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
public:
	uint16_t id; //internal ID

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

    TReadGroup(const uint16_t ID, const std::string Name);
    TReadGroup(const TReadGroup & other) = default;
    TReadGroup* getPointer(){ return this; };
    std::string compileSamHeader() const;

    bool operator<(const TReadGroup & right) const;
    bool operator<(const std::string & right) const;
};

bool operator<(const std::string & left, const TReadGroup & right);

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

	void clear();
	const TReadGroup& add(const std::string name);
	const TReadGroup& addAlternativeRG(const std::string Name, const std::string original);
	const TReadGroup& addAlternativeRG(const std::string Name, const uint16_t original);
	const TReadGroup& addAlternativeRG(const std::string Name, const TReadGroup original);

	uint16_t size() const;
	bool empty() const;

	uint16_t getId(const std::string & name) const;
	const std::string& getName (uint16_t readGroupId) const;
	std::vector<std::string> getNames(std::vector<uint16_t> & readGroupIds) const;
	const TReadGroup& getReadGroup(const std::string & name);
	const TReadGroup& operator[](uint16_t readGroupId) const;
	bool readGroupExists(const std::string & name) const;
	bool readGroupInUse(uint16_t readGroupId) const;
	bool readGroupInUse(const std::string name) const;

	//looping over
	std::set<TReadGroup, std::less<>>::iterator begin(){ return _readGroups.begin(); };
	std::set<TReadGroup, std::less<>>::iterator end(){ return _readGroups.end(); };
	std::set<TReadGroup, std::less<>>::iterator cbegin() const{ return _readGroups.cbegin(); };
	std::set<TReadGroup, std::less<>>::iterator cend() const{ return _readGroups.cend(); };

	void filterReadGroups(std::string readGroupList);
	void removeFromHeader(const std::string name);
	void removeFromHeader(const uint16_t readGroupId);
	void printReadgroupsInUse(coretools::TLog* logfile) const;
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
	void _fillFromFile(const TReadGroups & ReadGroups, const std::string & filename, coretools::TLog* logfile);

public:
	TReadGroupMap(const TReadGroups & ReadGroups);
	TReadGroupMap(const TReadGroups & ReadGroups, const std::string filename, coretools::TLog* logfile);

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
