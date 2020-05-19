/*
 * TAlignmentBlacklist.cpp
 *
 *  Created on: May 19, 2020
 *      Author: phaentu
 */

#include "TAlignmentBlacklist.h"

void TAlignmentBlacklist::enableWriting(const std::string filename){
	_write = true;
	out.open(filename);
	out.writeHeader({"Alignment", "isReverse", "Error"});
};

void TAlignmentBlacklist::addToBlacklist(const std::string alignmentName, const bool isReverse, const std::string errorMessage){
	//TODO: should check if read already exists in list (could be case in paired-end data) -> remove
	blacklist.insert(alignmentName);
	if(_write){
		out << alignmentName << isReverse << errorMessage << std::endl;
	}
};

void TAlignmentBlacklist::removeFromBlacklist(const std::string alignmentName, const bool isReverse, const std::string errorMessage){
	blacklist.erase(alignmentName);
	if(_write){
		out << alignmentName << isReverse << errorMessage << std::endl;
	}
};

bool TAlignmentBlacklist::isInBlacklist(const std::string & alignmentName){
	return blacklist.find(alignmentName);
};


