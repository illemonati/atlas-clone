/*
 * TAlignmentMerger.h
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#ifndef TBAMFILTER_H_
#define TBAMFILTER_H_

#include "TGenome.h"
#include "TAlignmentStorage.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"

namespace GenomeTasks{

//-----------------------------------------
// TAlignmentMergerReadGroupSettings
//-----------------------------------------
enum ReadGroupType : uint8_t { unchanged=0, single, mixed, paired};

struct TAlignmentMergerReadGroupSetting{
	uint16_t readGroupId;
	mutable uint16_t altReadGroupId;
	mutable ReadGroupType type;
	mutable uint16_t maxCycles;

	TAlignmentMergerReadGroupSetting(const uint16_t ReadGroupId, const ReadGroupType Type, const uint16_t MaxCycles){
		readGroupId = ReadGroupId;
		altReadGroupId = readGroupId;
		type = Type;
		maxCycles = MaxCycles;
	};

	bool operator<(const TAlignmentMergerReadGroupSetting & right) const { return readGroupId < right.readGroupId; };
	bool operator<(const uint16_t right) const { return readGroupId < right; };
};

bool operator<(const uint16_t left, const TAlignmentMergerReadGroupSetting & right);

class TAlignmentMergerReadGroupSettings{
private:
	std::set<TAlignmentMergerReadGroupSetting, std::less<> > _settings;

	void _printSummary(coretools::TLog* logfile);

public:
	TAlignmentMergerReadGroupSettings(){};

	void initialize(coretools::TParameters & Params, coretools::TLog* logfile, BAM::TReadGroups & readGroups);
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
	void _writeUpTo(const BAM::TGenomePosition & position);

	BAM::TAlignment* _parseIntoNewAlignment();
	virtual void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	virtual void _handleSingle(BAM::TAlignment* alignment);

public:
	TBamFilter(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	virtual ~TBamFilter(){};
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
	coretools::TRandomGenerator* _randomGenerator;
	bool _adaptQuality;

	void _mergeBasesCore(BAM::TSequencedBase & bestBase, BAM::TSequencedBase & worstBase);
	virtual void _mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate);

public:
	TAlignmentMerger_randomBase(coretools::TRandomGenerator* RandomGenerator, const bool AdaptQuality);
	virtual ~TAlignmentMerger_randomBase(){};
};

class TAlignmentMerger_randomRead:public TAlignmentMerger_randomBase{
private:
	bool _keepMate;

	void _mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate);
	bool _merge(BAM::TAlignment* alignment, BAM::TAlignment* mate);

public:
	TAlignmentMerger_randomRead(coretools::TRandomGenerator* RandomGenerator, const bool AdaptQuality);

	uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
};

class TAlignmentMerger_highestQuality:public TAlignmentMerger_randomBase{
private:
	void _mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate);
public:
	TAlignmentMerger_highestQuality(coretools::TRandomGenerator* RandomGenerator, const bool AdaptQuality);
};

//-----------------------------------------
// TAlignmentSplitMerger
//-----------------------------------------
class TAlignmentSplitMerger:public TBamFilter{
private:
	std::unique_ptr<TAlignmentMerger> _merger;
	uint64_t _numReadsMerged;
	uint64_t _numBasesMerged;
	uint8_t _considerAtMaxUpToDist;
	bool _allowForLarger;

	void _initializeMerger(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void _openBamFileForWriting();
	void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	void _handleSingle(BAM::TAlignment* alignment);

public:
	TAlignmentSplitMerger(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);

};

//-----------------------------------------
// TOverlapQuantifier
//-----------------------------------------
class TOverlapQuantifier:public TGenome_filtered{
private:
	TAlignmentMerger _merger;
	TAlignmentStorage _alignmentStorage;

public:
	TOverlapQuantifier(coretools::TParameters & Params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void quantifyOverlap();
};


//--------------------------------------
// Tasks
//--------------------------------------

	class TTask_filterBAM:public coretools::TTask{
public:
	TTask_filterBAM(){ _explanation = "Writing reads that pass filters to BAM file"; };

	void run(){
		using namespace coretools::instances;
		TBamFilter filter(parameters(), &logfile(), &randomGenerator());
		filter.traverseBAM();
	};
};

class TTask_splitMerge:public coretools::TTask{
public:
	TTask_splitMerge(){ _explanation = "Splitting single-end reads and merging paired-end reads and in BAM file"; };

	void run(){
		using namespace coretools::instances;
		TAlignmentSplitMerger splitMerger(parameters(), &logfile(), &randomGenerator());
		splitMerger.traverseBAM();
	};
};

class TTask_overlapQuantifier:public coretools::TTask{
public:
	TTask_overlapQuantifier(){ _explanation = "Estimating distribution of overlap of paired reads in BAM file"; };

	void run(){
		using namespace coretools::instances;
		TOverlapQuantifier overlapQuantifier(parameters(), &logfile(), &randomGenerator());
		overlapQuantifier.quantifyOverlap();
	};
};

}; //end namespace

#endif /* TBAMFILTER_H_ */
