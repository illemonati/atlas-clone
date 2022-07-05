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
#include "TReadGroups.h"

namespace BAM { class TAlignment; }
namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }
namespace genometools { class TGenomePosition; }

using BAM::ReadGroupType;
using BAM::readGroupType2String;

namespace GenomeTasks{

//-----------------------------------------
// TAlignmentMergerReadGroupSettings
//-----------------------------------------
class TAlignmentMergerReadGroupSetting{
private:
	uint16_t _readGroupId;
	ReadGroupType _type;
	uint16_t _maxCycles;

	uint16_t _splittingThreshold;
	std::vector<uint16_t> _altReadGroupIds; //indexed by length > threshold

public:
	constexpr TAlignmentMergerReadGroupSetting(const uint16_t ReadGroupId, const ReadGroupType Type, const uint16_t MaxCycles):
		_readGroupId(ReadGroupId),
		_type(Type),
		_maxCycles(MaxCycles),
		_splittingThreshold(MaxCycles+1){};

	constexpr TAlignmentMergerReadGroupSetting(const uint16_t ReadGroupId, const ReadGroupType Type, const uint16_t MaxCycles, uint16_t NumSplitCategories, BAM::TReadGroups & readGroups):
		_readGroupId(ReadGroupId),
		_type(Type),
		_maxCycles(MaxCycles){

		//create alternative read groups
		if(NumSplitCategories > _maxCycles){
			NumSplitCategories = _maxCycles;
		}
		_splittingThreshold = _maxCycles - NumSplitCategories + 2;

		_altReadGroupIds.resize(NumSplitCategories - 1);

		for(size_t i=0; i<(NumSplitCategories-1); ++i){
			_altReadGroupIds[i] = readGroups.addAlternativeRG(readGroups[ReadGroupId].name_ID + "_Length" + toString(_maxCycles-i), ReadGroupId);
		}
	};

	constexpr bool operator<(const TAlignmentMergerReadGroupSetting & right) const noexcept { return _readGroupId < right._readGroupId; };
	constexpr bool operator<(const uint16_t right) const noexcept { return _readGroupId < right; };

	//getters
	ReadGroupType type(){ return _type; };
	uint16_t maxCycles(){ return _maxCycles; };
	uint16_t readGroupId(){ return _readGroupId; };
	void splitReadGroup(BAM::TAlignment* alignment){
		if(alignment->length() >= _splittingThreshold){
			alignment->setReadGroup(_altReadGroupIds[alignment->length() - _splittingThreshold]);
		}
	};
};

constexpr bool operator<(const uint16_t left, const TAlignmentMergerReadGroupSetting & right) noexcept {
	return left < right._readGroupId;
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
protected:
	virtual void _mergeBases(BAM::TSequencedBase &, BAM::TSequencedBase &){};

public:
	TAlignmentMerger(){};
	virtual ~TAlignmentMerger(){};

	virtual uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_randomBase:public TAlignmentMerger{
protected:
	bool _adaptQuality;

	void _mergeBasesCore(BAM::TSequencedBase & bestBase, BAM::TSequencedBase & worstBase);
	virtual void _mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate);

public:
	TAlignmentMerger_randomBase(const bool AdaptQuality);
	virtual ~TAlignmentMerger_randomBase(){};
};

class TAlignmentMerger_randomRead:public TAlignmentMerger_randomBase{
private:
	bool _keepMate;

	void _mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate);
	bool _merge(BAM::TAlignment* alignment, BAM::TAlignment* mate);

public:
	TAlignmentMerger_randomRead(const bool AdaptQuality);

	uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_highestQuality:public TAlignmentMerger_randomBase{
private:
	void _mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate);
public:
	TAlignmentMerger_highestQuality(const bool AdaptQuality);
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
