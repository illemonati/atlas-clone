/*
 * TAlignmentStorage.cpp
 *
 *  Created on: Jul 12, 2021
 *      Author: wegmannd
 */

#include "TAlignmentStorage.h"
#include <algorithm>
#include "TAlignment.h"

namespace GenomeTasks{

//-----------------------------------------
// TAlignmentMergerEntry
//-----------------------------------------

TAlignmentMergerEntry::TAlignmentMergerEntry(BAM::TAlignment* Alignment, bool readyForWriting){
	alignment = Alignment;
	ready = readyForWriting;
};

TAlignmentMergerEntry::TAlignmentMergerEntry(TAlignmentMergerEntry && other){
	//copy from other
	alignment = other.alignment;
	ready = other.ready;

	//set other to default
	other.alignment = nullptr;
	other.ready = false;
};

TAlignmentMergerEntry& TAlignmentMergerEntry::operator=(TAlignmentMergerEntry && other){
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

TAlignmentMergerEntry::~TAlignmentMergerEntry(){
	delete alignment;
};

void TAlignmentMergerEntry::setAsNonProperPair(){
	alignment->setIsProperPair(false);
	ready = true;
};

//------------------------------------------------
// TAlignmentStorage
//------------------------------------------------

TAlignmentStorage::~TAlignmentStorage(){ _storage.clear(); };

bool TAlignmentStorage::empty(){
	return _storage.empty();
};

TAlignmentInStorage TAlignmentStorage::begin(){
	return _storage.begin();
};

TAlignmentInStorage TAlignmentStorage::end(){
	return _storage.end();
};

void TAlignmentStorage::emplace_back(BAM::TAlignment* alignment, const bool ready){
	_storage.emplace_back(alignment, ready);
};

TAlignmentInStorage TAlignmentStorage::find(const std::string & name){
	for(auto it=_storage.begin(); it!=_storage.end(); ++it){
		if(it->alignment->name() == name){
			return it;
		}
	}

	return _storage.end();
};

TAlignmentInStorage TAlignmentStorage::erase(TAlignmentInStorage it){
	return _storage.erase(it);
};

void TAlignmentStorage::clear(){
	_storage.clear();
};

}; //end namespace
