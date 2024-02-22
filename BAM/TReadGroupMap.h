#ifndef TREADGROUPMAP_H_
#define TREADGROUPMAP_H_

#include <cstddef>
#include <vector>

#include "TReadGroups.h"

namespace BAM {

class TReadGroups;

//--------------------------------------------------------------------------------------
// TReadGroupMap
// Maps bam file read group index to internal index, which may differ in case of pooling
//--------------------------------------------------------------------------------------
class TReadGroupMap {
private:
	static constexpr size_t ReadGroupMapNotInitializedIndex = -1; // largest possible values

	std::string _name;
	std::vector<size_t> _readGroupMap;                            // maps read group index to pooled index
	std::vector<std::vector<size_t>> _reverseReadGroupMap; // IDs of all Rgs pooled with that index. Includes itself.
	std::vector<size_t> _readGroupsInUse;

	void _resize(const TReadGroups &ReadGroups);
	void _markAsInUse(size_t index);
	void _noPooling(const TReadGroups &ReadGroups);
	void _poolAll(const TReadGroups &ReadGroups);
	void _fromFile(const TReadGroups &ReadGroups, std::string_view filename);

public:
	TReadGroupMap(std::string_view Name, const TReadGroups &ReadGroups, std::string_view Type = "");

	size_t size() const { return _readGroupMap.size(); };
	size_t pooledIndex(size_t rg) const { return _readGroupMap[rg]; };

	const std::vector<size_t> &readGroupsInUse() const { return _readGroupsInUse; };
	const std::vector<size_t> &readGroupsPooledWith(size_t rg) const { return _reverseReadGroupMap[rg]; };
};
} // namespace BAM

#endif
