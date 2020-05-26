/*
 * TAlignmentMerger.h
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#ifndef TALIGNMENTMERGER_H_
#define TALIGNMENTMERGER_H_

#include "TReadList.h"
#include "TAlignmentParser.h"

using namespace BAM;


//-----------------------------------------
// TAlignmentMergerEntry
//-----------------------------------------
class TAlignmentMergerEntry{
private:

public:
	TAlignment* alignment;
	bool ready;

	TAlignmentMergerEntry(TAlignment & Alignment, bool readyForWriting){
		alignment = new TAlignment(Alignment);
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
// TAlignmentMerger
//-----------------------------------------
class TAlignmentMerger{
private:
	TAlignmentStorage alignmentStorage; //bool indicates whether read is ready for writing
	BAM::TBamFile& bamFile;
	TQualityMap qualMap;
	TGenotypeMap genoMap;

	int _maxDistanceBetweenMates;
	bool _filterOrphans;
	bool _adaptQuality;
	bool _keepRandomBase;
	bool _keepRandomRead;
	bool _allowForLarger;
	bool _keepHigherQuality;

	void _writeAlignment(TAlignmentInStorage & it);
	void _addToBlacklist(TAlignmentInStorage & it, std::string error);
	void _writeAllThatAreReady();
	TAlignmentInStorage _findMate(TAlignment & alignment);

public:
	TAlignmentMerger(BamTools::BamWriter* Writer, TAlignmentParser* Parser, int MaxDistanceBetweenMates);

	void keepOrphans(){ _filterOrphans = false; };
	void keepOriginalQuality(){ _adaptQuality = false; };
	void keepRandomBase(){ _keepRandomBase = true; };
	void keepRandomRead(){ _keepRandomRead = true; };
	void allowForLarger(){ _allowForLarger = true; };
	void keepHigherQuality(){_keepHigherQuality = true; };

	void addToBeMerged(TAlignment & alignment, TRandomGenerator* randomGenerator);
	void addToBeSplit(TAlignment & alignment, std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT);
	void checkForMateAndWriteUnmerged(TAlignment & alignment);
	void addAsImproperPair(TAlignment & alignment);
	void addReadyToBeWritten(TAlignment & alignment);
	void writeUpTo(const int position);
	void clear();

};



#endif /* TALIGNMENTMERGER_H_ */
