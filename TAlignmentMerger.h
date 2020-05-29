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

enum AlignmentMergerType : uint8_t {none=0, randomRead, randomBase, highestQuality};

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
	void traverseBAM(const std::string outputName);
};

//-----------------------------------------
// TAlignmentMergerType
//-----------------------------------------
class TAlignmentMergerType_none{
protected:
	uint16_t _numOverlap;
	bool _adaptQuality;

	virtual void _prepare();
	bool TAlignmentMergerType_none::_mergeBases(TBase & bestBase, TBase & worstBase);

public:
	TAlignmentMergerType_none();
	~TAlignmentMergerType_none(){};

	void merge(BAM::TAlignment* alignment, BAM::TAlignment* mate);
};

class TAlignmentMergerType_randomRead:public TAlignmentMergerType_none{
protected:
	TRandomGenerator* _randomGenerator;
	bool _keepMate;

	virtual void _prepare();
	virtual bool _merge(BAM::TAlignment* alignment, BAM::TAlignment* mate);

public:
	TAlignmentMergerType_randomRead(TRandomGenerator* RandomGenerator);
	~TAlignmentMergerType_randomRead(){};
	virtual bool merge(BAM::TAlignment* alignment, BAM::TAlignment* mate);
};

//-----------------------------------------
// TAlignmentMerger
//-----------------------------------------
class TAlignmentMerger:public TMateFilter{
protected:
	TRandomGenerator* _randomGenerator;
	AlignmentMergerType _type;

	bool _adaptQuality;

	void _handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate);
	void _mergeBases(TBase & bestBase, TBase & worstBase);

public:
	TAlignmentMerger(BAM::TBamFile& BamFile, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void mergeBAM(const std::string outputName);

	void addToBeMerged(TAlignment & alignment, TRandomGenerator* randomGenerator);
	void addToBeSplit(TAlignment & alignment, std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT);
	void checkForMateAndWriteUnmerged(TAlignment & alignment);
	void addAsImproperPair(TAlignment & alignment);
	void addReadyToBeWritten(TAlignment & alignment);
	void writeUpTo(const int position);
	void clear();

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
