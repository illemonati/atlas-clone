/*
 * TAlignmentMerger.h
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#ifndef TBAMFILTER_H_
#define TBAMFILTER_H_

#include "TReadList.h"
#include "TGenome.h"

namespace AlignmentSplitMerge{


//-----------------------------------------
// TAlignmentMergerReadGroupSettings
//-----------------------------------------
enum ReadGroupType : uint8_t { unchanged=0, single, mixed, paired};

struct TAlignmentMergerReadGroupSetting{
	uint16_t readGroupId;
	uint16_t altReadGroupId;
	ReadGroupType type;
	uint16_t maxCycles;

	TAlignmentMergerReadGroupSetting(const uint16_t ReadGroupId, const ReadGroupType Type, const uint16_t MaxCycles){
		readGroupId = ReadGroupId;
		altReadGroupId = readGroupId;
		type = Type;
		maxCycles = MaxCycles;
	};

	bool operator<(const TAlignmentMergerReadGroupSetting & left, const TAlignmentMergerReadGroupSetting & right){ return left.readGroupId < right.readGroupId; };
	bool operator<(const TAlignmentMergerReadGroupSetting & left, const uint16_t & right){ return left.readGroupId < right; };
	bool operator<(const uint16_t & left, const TAlignmentMergerReadGroupSetting & right){ return left < right.readGroupId; };

};

class TAlignmentMergerReadGroupSettings{
private:
	std::set<TAlignmentMergerReadGroupSetting, std::less<> > _settings;

	void _printSummary(TLog* logfile);

public:
	TAlignmentMergerReadGroupSettings(){};

	void initialize(TParameters & Params, TLog* logfile, BAM::TReadGroups & readGroups);
	void setAllAsUnchanged(BAM::TReadGroups & readGroups);
	bool needTruncation() const;
	bool needsMerging() const;
	ReadGroupType getType(const uint16_t readGroupId) const;
	uint16_t getMaxCycles(const uint16_t readGroupId) const;
	const TAlignmentMergerReadGroupSetting& getSettings(const uint16_t readGroupId) const;
};

//-----------------------------------------
// TAlignmentStorage
//-----------------------------------------
class TAlignmentMergerEntry{
private:

public:
	BAM::TAlignment* alignment;
	bool ready;

	TAlignmentMergerEntry(BAM::TAlignment* Alignment, bool readyForWriting);
	TAlignmentMergerEntry(TAlignmentMergerEntry && other);
	TAlignmentMergerEntry& operator=(TAlignmentMergerEntry && other);
	~TAlignmentMergerEntry();
	void setAsNonProperPair();
};

typedef std::vector< TAlignmentMergerEntry >::iterator TAlignmentInStorage;

class TAlignmentStorage{
private:
	std::vector< TAlignmentMergerEntry > _storage;

public:
	TAlignmentStorage(){};
	~TAlignmentStorage();

	bool empty();
	TAlignmentInStorage& begin();
	TAlignmentInStorage& end();
	void emplace_back(BAM::TAlignment* alignment, const bool ready);
	TAlignmentInStorage& find(const std::string & name);
	TAlignmentInStorage& erase(TAlignmentInStorage it);
	void clear();
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
	void _writeUpTo(const int position);

	BAM::TAlignment* _parseIntoNewAlignment();
	virtual void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	virtual void _handleSingle(BAM::TAlignment* alignment);

public:
	TBamFilter(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TBamFilter(){};
	virtual void traverseBAM();
};

//-----------------------------------------
// TAlignmentMergerType
// base class does not merge
//-----------------------------------------
class TAlignmentMerger{
protected:
	virtual void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap){};

public:
	TAlignmentMerger(){};
	virtual ~TAlignmentMerger(){};

	virtual uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate, const TQualityMap & qualMap);
};

class TAlignmentMerger_randomBase:public TAlignmentMerger{
protected:
	TRandomGenerator* _randomGenerator;
	bool _adaptQuality;

	void _mergeBasesCore(TBase & bestBase, TBase & worstBase, const TQualityMap & qualMap);
	virtual void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap);

public:
	TAlignmentMerger_randomBase(TRandomGenerator* RandomGenerator, const bool AdaptQuality);
	virtual ~TAlignmentMerger_randomBase(){};
};

class TAlignmentMerger_randomRead:public TAlignmentMerger_randomBase{
private:
	bool _keepMate;

	void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap);
	bool _merge(BAM::TAlignment* alignment, BAM::TAlignment* mate);

public:
	TAlignmentMerger_randomRead(TRandomGenerator* RandomGenerator, const bool AdaptQuality);

	uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate, const TQualityMap & qualMap);
};

class TAlignmentMerger_highestQuality:public TAlignmentMerger_randomBase{
private:
	void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap);
public:
	TAlignmentMerger_highestQuality(TRandomGenerator* RandomGenerator, const bool AdaptQuality);
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

	void _initializeMerger(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void _openBamFileForWriting();
	void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	void _handleSingle(BAM::TAlignment* alignment);

public:
	TAlignmentSplitMerger(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);

};

//-----------------------------------------
// TOverlapQuantifier
//-----------------------------------------
class TOverlapQuantifier:public TGenome_filtered{
private:
	TAlignmentMerger _merger;
	TAlignmentStorage _alignmentStorage;

public:
	TOverlapQuantifier(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void quantifyOverlap();
};

}; //end namespace

#endif /* TBAMFILTER_H_ */
