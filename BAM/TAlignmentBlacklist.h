/*
 * TAlignmentBlacklist.h
 *
 *  Created on: May 19, 2020
 *      Author: phaentu
 */

#ifndef BAM_TALIGNMENTBLACKLIST_H_
#define BAM_TALIGNMENTBLACKLIST_H_

#include <set>
#include "TFile.h"

class TAlignmentBlacklist{
private:
	std::set<std::string> blacklist;
	TOutputFileZipped out;
	bool _write;


public:
	TAlignmentBlacklist(){
		_write = false;
	};

	void enableWriting(const std::string filename);
	void removeFromBlacklist(const std::string alignmentName, const bool isReverse, const std::string errorMessage);
	bool isInBlacklist(const std::string & alignmentName);
};


#endif /* BAM_TALIGNMENTBLACKLIST_H_ */
