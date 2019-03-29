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

class TAlignmentMergerEntry{
private:

public:
	TAlignment* alignment;
	bool ready;

	TAlignmentMergerEntry(TAlignment & Alignment, bool readyForWriting){
		alignment = new TAlignment(Alignment);
		ready = readyForWriting;
	};

	constexpr TAlignmentMergerEntry(TAlignmentMergerEntry && other):alignment(nullptr),ready(false){
		//copy from other
		alignment = other.alignment;
		ready = other.ready;

		//set other to default
		other.alignment = nullptr;
		other.ready = false;
	};

	constexpr TAlignmentMergerEntry& operator=(TAlignmentMergerEntry && other){
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


class TAlignmentMerger{
private:
	std::vector< TAlignmentMergerEntry > alignmentStorage; //bool indicates wheter read is ready for writing
	BamTools::BamWriter* writer;
	TAlignmentParser* parser;

	int _maxDistanceBetweenMates;
	bool _filterOrphans;
	bool _adaptQuality;

	void _writeAlignment(std::vector< TAlignmentMergerEntry >::iterator & it);
	void _addToBlacklist(std::vector< TAlignmentMergerEntry >::iterator & it, std::string error);
	void _writeAllThatAreReady();
	std::vector< TAlignmentMergerEntry >::iterator _findMate(TAlignment & alignment);

public:
	TAlignmentMerger(BamTools::BamWriter* Writer, TAlignmentParser* Parser, int MaxDistanceBetweenMates);

	void keepOrphans(){ _filterOrphans = false; };
	void keepOriginalQuality(){ _adaptQuality = false; };

	void addToBeMerged(TAlignment & alignment);
	void addAsImproperPair(TAlignment & alignment);
	void addReadyToBeWritten(TAlignment & alignment);
	void writeUpTo(const int position);
	void clear();

};



#endif /* TALIGNMENTMERGER_H_ */
