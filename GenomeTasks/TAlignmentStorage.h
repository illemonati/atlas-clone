/*
 * TAlignmentStorage.h
 *
 *  Created on: Jul 12, 2021
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TALIGNMENTSTORAGE_H_
#define GENOMETASKS_TALIGNMENTSTORAGE_H_

#include "TAlignment.h"
#include <vector>

namespace GenomeTasks{

//-----------------------------------------
// TAlignmentMergerEntry
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

//-----------------------------------------
// TAlignmentStorage
//-----------------------------------------

class TAlignmentStorage{
private:
	std::vector< TAlignmentMergerEntry > _storage;

public:
	TAlignmentStorage(){};
	~TAlignmentStorage();

	bool empty();
	TAlignmentInStorage begin();
	TAlignmentInStorage end();
	void emplace_back(BAM::TAlignment* alignment, const bool ready);
	TAlignmentInStorage find(const std::string & name);
	TAlignmentInStorage erase(TAlignmentInStorage it);
	void clear();
};


}; //end namesapce


#endif /* GENOMETASKS_TALIGNMENTSTORAGE_H_ */
