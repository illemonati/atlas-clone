/*
 * TAlignmentMerger.h
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#ifndef TALIGNMENTMERGER_H_
#define TALIGNMENTMERGER_H_

#include "../IOTools/IOAbstractClasses/BamWriter.h"
#include "../TReadList.h"
#include "TAlignmentParser.h"


class TAlignmentMergerEntry{
private:

public:
	TAlignment* alignment;
	bool ready;

    TAlignmentMergerEntry(TAlignment & Alignment, bool readyForWriting);
    TAlignmentMergerEntry(TAlignmentMergerEntry && other);
    ~TAlignmentMergerEntry();

    TAlignmentMergerEntry& operator=(TAlignmentMergerEntry && other);

    void setAsNonProperPair();
};


class TAlignmentMerger{
private:
	std::vector< TAlignmentMergerEntry > alignmentStorage; //bool indicates wheter read is ready for writing
    BamWriter* writer;
	TAlignmentParser* parser;

	int _maxDistanceBetweenMates;
	bool _filterOrphans;
	bool _adaptQuality;

	void _writeAlignment(std::vector< TAlignmentMergerEntry >::iterator & it);
	void _addToBlacklist(std::vector< TAlignmentMergerEntry >::iterator & it, std::string error);
	void _writeAllThatAreReady();
	std::vector< TAlignmentMergerEntry >::iterator _findMate(TAlignment & alignment);

public:
    TAlignmentMerger(BamWriter* Writer, TAlignmentParser* Parser, int MaxDistanceBetweenMates);

	void keepOrphans(){ _filterOrphans = false; };
	void keepOriginalQuality(){ _adaptQuality = false; };

	void addToBeMerged(TAlignment & alignment);
	void checkForMateAndWriteUnmerged(TAlignment & alignment);
	void addAsImproperPair(TAlignment & alignment);
	void addReadyToBeWritten(TAlignment & alignment);
	void writeUpTo(const int position);
	void clear();

};



#endif /* TALIGNMENTMERGER_H_ */
