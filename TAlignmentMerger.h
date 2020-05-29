/*
 * TAlignmentMerger.h
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#ifndef TALIGNMENTMERGER_H_
#define TALIGNMENTMERGER_H_

#include "TReadList.h"
#include "TBamFile.h"

//-----------------------------------------
// TAlignmentMergerEntry
//-----------------------------------------
class TAlignmentMergerEntry{
private:

public:
	BAM::TAlignment* alignment;
	bool ready;

	TAlignmentMergerEntry(BAM::TAlignment* Alignment, bool readyForWriting){
		alignment = Alignment;
		ready = readyForWriting;
	};

	TAlignmentMergerEntry(TAlignmentMergerEntry && other){
		//copy from other
		alignment = other.alignment;
		ready = other.ready;

		//set other to default
		other.alignment = nullptr;
		other.ready = false;
	};

	TAlignmentMergerEntry& operator=(TAlignmentMergerEntry && other){
		if(this != &other){
			//free object
			delete alignment;

			//copy from other
			alignment = other.alignment;
			ready = other.ready;

			//set other to default
			other.alignment = nullptr;
			other.ready = false;
		}

		return *this;
	};

	~TAlignmentMergerEntry(){
		delete alignment;
	};

	void setAsNonProperPair(){
		alignment->setIsProperPair(false);
		ready = true;
	};
};

typedef std::vector< TAlignmentMergerEntry > TAlignmentStorage;
typedef std::vector< TAlignmentMergerEntry >::iterator TAlignmentInStorage;

//-----------------------------------------
// TMateFinder
//-----------------------------------------
class TMateFilter{
protected:
	BAM::TBamFile& _bamFile;
	BAM::TAlignmentBlacklist _blacklist; //used to keep track of filtered out mates
	TAlignmentStorage _alignmentStorage;

	TGenotypeMap genoMap;
	TQualityMap qualMap;

	TLog* logfile;
	uint32_t _maxDistanceBetweenMates;
	bool _keepOrphans;

	TAlignmentInStorage _findMate(const std::string & name);
	void _writeAlignment(TAlignmentInStorage & it);
	void _writeOrFilterAsOrphan(TAlignmentInStorage & it);
	void _writeAll();
	void _writeUpTo(const int position);

	virtual void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);

public:
	TMateFilter(BAM::TBamFile& BamFile, TParameters & params, TLog* Logfile);
	virtual ~TMateFilter(){};
	void traverseBAM(const std::string outputName);
};

//-----------------------------------------
// TAlignmentMergerType
// base class does not merge
//-----------------------------------------
class TAlignmentMergerType{
protected:
	virtual void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap){};

public:
	TAlignmentMergerType(){};
	virtual ~TAlignmentMergerType(){};

	uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate, const TQualityMap & qualMap);
};

class TAlignmentMergerType_randomBase:public TAlignmentMergerType{
protected:
	TRandomGenerator* _randomGenerator;
	bool _adaptQuality;

	void _mergeBasesCore(TBase & bestBase, TBase & worstBase, const TQualityMap & qualMap);
	virtual void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap);

public:
	TAlignmentMergerType_randomBase(TRandomGenerator* RandomGenerator, const bool AdaptQuality);
	virtual ~TAlignmentMergerType_randomBase(){};
};

class TAlignmentMergerType_randomRead:public TAlignmentMergerType_randomBase{
private:
	bool _keepMate;

	void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap);
	bool _merge(BAM::TAlignment* alignment, BAM::TAlignment* mate);

public:
	TAlignmentMergerType_randomRead(TRandomGenerator* RandomGenerator, const bool AdaptQuality);

	uint16_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate, const TQualityMap & qualMap);
};

class TAlignmentMergerType_highestQuality:public TAlignmentMergerType_randomBase{
private:
	void _mergeBases(TBase & alignment, TBase & mate, const TQualityMap & qualMap);
public:
	TAlignmentMergerType_highestQuality(TRandomGenerator* RandomGenerator, const bool AdaptQuality);
};

//-----------------------------------------
// TAlignmentMerger
//-----------------------------------------
class TAlignmentMerger:public TMateFilter{
protected:
	std::unique_ptr<TAlignmentMergerType> _merger;
	uint64_t _numReadsMerged;
	uint64_t _numBasesMerged;

	void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	void _mergeBases(TBase & bestBase, TBase & worstBase);

public:
	TAlignmentMerger(BAM::TBamFile& BamFile, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TAlignmentMerger()

	void mergeBAM(const std::string outputName);
};

//-----------------------------------------
// TAlignmentSplitMerger
//-----------------------------------------
class TAlignmentSplitMerger:public TAlignmentMerger{
private:

	bool _allowForLarger;

public:


};


#endif /* TALIGNMENTMERGER_H_ */
