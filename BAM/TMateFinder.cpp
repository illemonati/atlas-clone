/*
 * TAlignmentBlacklist.cpp
 *
 *  Created on: May 19, 2020
 *      Author: phaentu
 */

#include <TMateFinder.h>

namespace BAM{

void TMateFinder::enableWriting(const std::string filename){
	_write = true;
	out.open(filename);
	out.writeHeader({"Alignment", "isSecondMate", "Error"});
};

void TMateFinder::addToBlacklist(const std::string alignmentName, const bool isSecondMate, const std::string errorMessage){
	//TODO: should check if read already exists in list (could be case in paired-end data) -> remove
	blacklist.insert(alignmentName);
	if(_write){
		out << alignmentName << isSecondMate << errorMessage << std::endl;
	}
};

void TMateFinder::removeFromBlacklist(const std::string alignmentName, const bool isRisSecondMateverse, const std::string errorMessage){
	blacklist.erase(alignmentName);
	if(_write){
		out << alignmentName << isSecondMate << errorMessage << std::endl;
	}
};

bool TMateFinder::isInBlacklist(const std::string & alignmentName){
	return blacklist.find(alignmentName);
};

}; //end namesapce
