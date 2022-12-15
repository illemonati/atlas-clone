/*
 * TAlignmentMerger.h
 *
 *  Created on: Oct 20, 2022
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TALIGNMENTMERGER_H_
#define GENOMETASKS_TALIGNMENTMERGER_H_

#include <functional>
#include <memory>
#include <set>
#include <stdint.h>
#include <string>

#include "TBamFilter.h"
#include "TBamFile.h"
#include "TGenome.h"
#include "coretools/Main/TTask.h"


namespace BAM { class TAlignment; }
namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }
namespace genometools { class TGenomePosition; }

namespace GenomeTasks{

namespace AlignmentMerger {

using namespace GenomeTasks::BamFilter;

//-----------------------------------------
// TAlignmentMergerReadGroupSettings
//-----------------------------------------
enum ReadGroupType : uint8_t { unchanged=0, single, mixed, paired};

struct TAlignmentMergerReadGroupSetting{
	uint16_t readGroupId;
	uint16_t altReadGroupId;
	ReadGroupType type;
	uint16_t maxCycles;

	constexpr TAlignmentMergerReadGroupSetting(const uint16_t ReadGroupId, const ReadGroupType Type, const uint16_t MaxCycles)
		: readGroupId(ReadGroupId), altReadGroupId(ReadGroupId), type(Type), maxCycles(MaxCycles){};

	constexpr TAlignmentMergerReadGroupSetting(const uint16_t ReadGroupId, uint16_t AltReadGroupId, const ReadGroupType Type, const uint16_t MaxCycles)
		: readGroupId(ReadGroupId), altReadGroupId(AltReadGroupId), type(Type), maxCycles(MaxCycles){};

	constexpr bool operator<(const TAlignmentMergerReadGroupSetting & right) const noexcept { return readGroupId < right.readGroupId; };
	constexpr bool operator<(const uint16_t right) const noexcept { return readGroupId < right; };
};

constexpr bool operator<(const uint16_t left, const TAlignmentMergerReadGroupSetting & right) noexcept {
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
	ReadGroupType getType(const uint16_t readGroupId) const;
	uint16_t getMaxCycles(const uint16_t readGroupId) const;
	const TAlignmentMergerReadGroupSetting& getSettings(const uint16_t readGroupId) const;
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

class TAlignmentMerger_randomRead:public TAlignmentMerger{
public:
	TAlignmentMerger_randomRead();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_middle:public TAlignmentMerger{
public:
	TAlignmentMerger_middle();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
	std::pair<size_t,bool> determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate);
	void sameDirectionMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate, std::pair<size_t,bool> overlapLength);
};

class TAlignmentMerger_firstMate:public TAlignmentMerger{
public:
	TAlignmentMerger_firstMate();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_secondMate:public TAlignmentMerger{
public:
	TAlignmentMerger_secondMate();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_highestQuality:public TAlignmentMerger{
public:
	TAlignmentMerger_highestQuality();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
	size_t overlapLengthAndMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};


//-----------------------------------------
// TAlignmentSplitMerger
//-----------------------------------------
class TAlignmentSplitMerger final:public TGenomeParsedWithAlignmentStorage<TAlignmentStorageSorted, TAlignmentStorageSortedIterator>{
private:
	std::unique_ptr<TAlignmentMerger> _merger;
	TAlignmentMergerReadGroupSettings _rgSettings;
	bool _allowForLarger;

	void _initializeMerger();
	void _openBamFileForWriting() override;
	void _handleMates(BAM::TAlignment & alignment, TAlignmentStorageSortedIterator mate) override;
	void _handleSingle(BAM::TAlignment & alignment) override;
	bool _alignmentCanBeWrittenUnchanged() override;

public:
	TAlignmentSplitMerger();
};

//-----------------------------------------
// TOverlapQuantifier
//-----------------------------------------
class TOverlapQuantifier:public TGenome_filtered{
private:
	TAlignmentMerger _merger;
	TAlignmentStorage _alignmentStorage;
	void _handleAlignment() override {};

public:
	void quantifyOverlap();
};

}; //end namespace AlignmentMerger


//-----------------------------------------
// Tasks
//-----------------------------------------

class TTask_splitMerge:public coretools::TTask{
public:
	TTask_splitMerge(){ _explanation = "Splitting single-end reads and merging paired-end reads in BAM file"; };

	void run(){
		AlignmentMerger::TAlignmentSplitMerger splitMerger;
		splitMerger.traverseBAM();
	};
};

class TTask_overlapQuantifier:public coretools::TTask{
public:
	TTask_overlapQuantifier(){ _explanation = "Estimating distribution of overlap of paired reads in BAM file"; };

	void run(){
		AlignmentMerger::TOverlapQuantifier overlapQuantifier;
		overlapQuantifier.quantifyOverlap();
	};
};


}; //end namespace GenomeTasks


#endif /* GENOMETASKS_TALIGNMENTMERGER_H_ */
