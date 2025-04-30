/*
 * TReadGroups.h
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

#ifndef TREADGROUPS_H_
#define TREADGROUPS_H_

#include <string>
#include <vector>

#include "TReadGroup.h"

namespace BAM{

class TReadGroups {
private:
	TReadGroup _noReadGroup;
	std::vector<TReadGroup> _readGroups;
	std::vector<size_t> _readGroupsById;
	bool _limitReadGroups = false;
	int _readGroupIdForReadsWithoutReadGroup; // is noReadGroupId if it does not exist (default)

	std::vector<TReadGroup>::iterator _getReadGroup(std::string_view Name);
	std::vector<TReadGroup>::const_iterator _getReadGroup(std::string_view Name) const;
	void _fillLookupFromId();

public:
	static constexpr size_t noReadGroupId = -1;

	TReadGroups();

	// add and remove read groups
	void clear();	
	TReadGroup& add(std::string_view name);
	TReadGroup& addAlternativeRG(std::string_view Name, std::string_view original);
	void addReadGroupForReadsWithoutReadGroup(std::string_view Name);
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

}; //end namespace

#endif /* TREADGROUPS_H_ */
