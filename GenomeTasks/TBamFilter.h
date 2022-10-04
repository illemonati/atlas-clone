/*
 * TAlignmentMerger.h
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#ifndef TBAMFILTER_H_
#define TBAMFILTER_H_

#include <functional>
#include <memory>
#include <set>
#include <stdint.h>
#include <string>

#include "TAlignmentStorage.h"
#include "TBamFile.h"
#include "TGenome.h"
#include "TTask.h"

namespace BAM { class TAlignment; }
namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }
namespace genometools { class TGenomePosition; }

namespace GenomeTasks{

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
// TBamFilter
//-----------------------------------------
class TBamFilter:public TGenome_parsed{
protected:
	TAlignmentMergerReadGroupSettings _rgSettings;
	BAM::TAlignmentList _blacklist; //used to keep track of filtered out mates
	TAlignmentStorage _alignmentStorage;

	BAM::TOutputBamFile _outBam;

	uint32_t _maxDistanceBetweenMates;
	bool _recalibrate;
	bool _incorporatePMD;
	bool _keepOrphans;

	virtual void _openBamFileForWriting();
	void _writeAlignment(TAlignmentInStorage & it);
	void _writeOrFilterAsOrphan(TAlignmentInStorage & it);
	void _writeAll();
	void _writeUpTo(const genometools::TGenomePosition & position);

	BAM::TAlignment* _parseIntoNewAlignment();
	virtual void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	virtual void _handleSingle(BAM::TAlignment* alignment);
	void _handleAlignment() override {}

public:
	TBamFilter();
	virtual void traverseBAM();
};

//-----------------------------------------
// TAlignmentMergerType
// base class does not merge
//-----------------------------------------
class TAlignmentMerger{
public: 
	TAlignmentMerger(){};
	virtual ~TAlignmentMerger(){};
	virtual uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
	virtual uint32_t calculateFirstReadPositionWithOverlap(BAM::TAlignment & alignment) const;
	virtual uint32_t calculateLastReadPositionWithOverlap(BAM::TAlignment & alignment) const;
	std::pair<genometools::TGenomePosition,genometools::TGenomePosition> findFirstAndLastReadPos(BAM::TAlignment & alignment) const;
	std::pair<uint16_t,bool> determineOverlapLength(std::pair<genometools::TGenomePosition,genometools::TGenomePosition> alignmentEndPositions, std::pair<genometools::TGenomePosition,genometools::TGenomePosition> mateEndPositions);
};
/*
class TAlignmentMerger_randomBase:public TAlignmentMerger{
protected:
	bool _adaptQuality;

	void _mergeBasesCore(BAM::TSequencedBase & bestBase, BAM::TSequencedBase & worstBase);
	virtual void _mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate);

public:
	TAlignmentMerger_randomBase(const bool AdaptQuality);
	virtual ~TAlignmentMerger_randomBase(){};
};
*/
class TAlignmentMerger_randomRead:public TAlignmentMerger{
private:
	bool _keepMate;
public:
	TAlignmentMerger_randomRead();
	uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_highestQuality:public TAlignmentMerger{
public:
	TAlignmentMerger_highestQuality();
	uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> getMinQuals(BAM::TAlignment & alignment, BAM::TAlignment & mate, std::pair<uint16_t,bool> overlapLength) const;
std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> minQual(BAM::TAlignment & firstRead, BAM::TAlignment & secondRead, uint16_t &overlapLength) const;
};

//-----------------------------------------
// TAlignmentSplitMerger
//-----------------------------------------
class TAlignmentSplitMerger:public TBamFilter{
private:
	std::unique_ptr<TAlignmentMerger> _merger;
	bool _allowForLarger;

	void _initializeMerger();
	void _openBamFileForWriting();
	void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	void _handleSingle(BAM::TAlignment* alignment);

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


//--------------------------------------
// Tasks
//--------------------------------------

	class TTask_filterBAM:public coretools::TTask{
public:
	TTask_filterBAM(){ _explanation = "Writing reads that pass filters to BAM file"; };

	void run(){
		TBamFilter filter; 
		filter.traverseBAM();
	};
};

class TTask_splitMerge:public coretools::TTask{
public:
	TTask_splitMerge(){ _explanation = "Splitting single-end reads and merging paired-end reads in BAM file"; };

	void run(){
		TAlignmentSplitMerger splitMerger;
		splitMerger.traverseBAM();
	};
};

class TTask_overlapQuantifier:public coretools::TTask{
public:
	TTask_overlapQuantifier(){ _explanation = "Estimating distribution of overlap of paired reads in BAM file"; };

	void run(){
		TOverlapQuantifier overlapQuantifier;
		overlapQuantifier.quantifyOverlap();
	};
};

}; //end namespace

#endif /* TBAMFILTER_H_ */
