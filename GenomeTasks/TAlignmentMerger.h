/*
 * TAlignmentMerger.h
 *
 *  Created on: Oct 20, 2022
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TALIGNMENTMERGER_H_
#define GENOMETASKS_TALIGNMENTMERGER_H_

#include <memory>
#include <set>

#include "TBamFilter.h"
#include "TWaitingListBamTraverser.h"

namespace BAM { class TAlignment; }
namespace BAM { class TReadGroups; }

namespace GenomeTasks{

namespace AlignmentMerger {


//-----------------------------------------
// TAlignmentMergerReadGroupSettings
//-----------------------------------------
	enum class ReadGroupType : uint8_t { min=0, unchanged=min, single, mixed, paired, max};

struct TAlignmentMergerReadGroupSetting{
	size_t readGroupId;
	size_t altReadGroupId;
	ReadGroupType type;
	size_t maxCycles;

	constexpr TAlignmentMergerReadGroupSetting(size_t ReadGroupId, const ReadGroupType Type, size_t MaxCycles)
		: readGroupId(ReadGroupId), altReadGroupId(ReadGroupId), type(Type), maxCycles(MaxCycles){};

	constexpr TAlignmentMergerReadGroupSetting(size_t ReadGroupId, size_t AltReadGroupId, const ReadGroupType Type, size_t MaxCycles)
		: readGroupId(ReadGroupId), altReadGroupId(AltReadGroupId), type(Type), maxCycles(MaxCycles){};

	constexpr bool operator<(const TAlignmentMergerReadGroupSetting & right) const noexcept { return readGroupId < right.readGroupId; };
	constexpr bool operator<(size_t right) const noexcept { return readGroupId < right; };
};

constexpr bool operator<(size_t left, const TAlignmentMergerReadGroupSetting & right) noexcept {
	return left < right.readGroupId;
};

class TAlignmentMergerReadGroupSettings{
private:
	std::set<TAlignmentMergerReadGroupSetting, std::less<> > _settings;

	void _printSummary();
public:
	void initialize(BAM::TReadGroups & readGroups);
	void setAllAsUnchanged(const BAM::TReadGroups & readGroups);
	bool needTruncation() const;
	bool needsMerging() const;
	ReadGroupType getType(size_t readGroupId) const;
	size_t getMaxCycles(size_t readGroupId) const;
	const TAlignmentMergerReadGroupSetting& getSettings(size_t readGroupId) const;
};

//-----------------------------------------
// TAlignmentMergerType
// base class does not merge
//-----------------------------------------
class TAlignmentMerger{
public:
	TAlignmentMerger(){};
	virtual ~TAlignmentMerger(){};
	virtual size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
	size_t determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate);
	virtual size_t overlapLengthAndMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_randomRead final : public TAlignmentMerger {
public:
	TAlignmentMerger_randomRead();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_middle final :public TAlignmentMerger{
public:
	TAlignmentMerger_middle();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
	std::pair<size_t,bool> determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate);
	void sameDirectionMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate, std::pair<size_t,bool> overlapLength);
};

class TAlignmentMerger_firstMate final : public TAlignmentMerger {
public:
	TAlignmentMerger_firstMate();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};

class TAlignmentMerger_secondMate final : public TAlignmentMerger {
public:
	TAlignmentMerger_secondMate();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};

class TAlignmentMerger_highestQuality final : public TAlignmentMerger {
public:
	TAlignmentMerger_highestQuality();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
	size_t overlapLengthAndMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};

//-----------------------------------------
// TAlignmentSplitMerger
//-----------------------------------------
class TAlignmentSplitMerger final
	: public TWaitingListBamTraverser {
private:
	std::unique_ptr<TAlignmentMerger> _merger;
	TAlignmentMergerReadGroupSettings _rgSettings;
	bool _allowForLarger;

	void _initializeMerger();
	void _handleMates(BAM::TAlignment & alignment, iterator mate) override;
	void _handleSingle(BAM::TAlignment & alignment) override;
	bool _alignmentCanBeWrittenUnchanged() override;

public:
	TAlignmentSplitMerger();
	void run() {
		traverseBAM();
	}
};

//-----------------------------------------
// TOverlapQuantifier
//-----------------------------------------
class TOverlapQuantifier{
private:
	TGenome _genome;
	TAlignmentMerger _merger;
	std::vector<TWaitingAlignment> _alignmentStorage;

public:
	TOverlapQuantifier();
	void run();
};

} //end namespace AlignmentMerger
} //end namespace GenomeTasks


#endif /* GENOMETASKS_TALIGNMENTMERGER_H_ */
