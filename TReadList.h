/*
 * TReadList.h
 *
 *  Created on: Mar 27, 2019
 *      Author: linkv
 */

#ifndef TREADLIST_H_
#define TREADLIST_H_

#include "TAlignment.h"

class TReadList{
private:
	bool _writeToFile;
	gz::ogzstream ignoredReads;
	std::string filename;
	std::map <std::string, int> readList;

	void _write(const std::string & Name, const std::string & Error, const bool & isSecondMate){
		if(_writeToFile){
			if(isSecondMate){
				ignoredReads << "Read " << Name << " (second mate): " << Error << "\n";
			} else {
				ignoredReads << "Read " << Name << " (first mate): " << Error << "\n";
			}
		}
	};

public:
	TReadList(){
		_writeToFile = false;
	}

	TReadList(const bool & writeToFile, std::string & filename){
		_writeToFile = writeToFile;
	}

	void setWriteReadListToFileToTrue(TLog* logfile, std::string & Filename){
		_writeToFile = true;
		filename = Filename + "_ignoredReads.txt.gz";
		logfile->list("Writing sequencing depth estimates to '" + filename + "'");
		ignoredReads.open(filename.c_str());
		if(!ignoredReads) throw "Failed to open output file '" + filename + "'!";
	}

	void addToReadList(TAlignment & alignment, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		readList.emplace(alignment.name(), 1);
		_write(alignment.name(), errorMessage, alignment.isSecondMate);
	};

	void addToReadList(std::string & alignmentName, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		readList.emplace(alignmentName, 1);
		if(_writeToFile){
			ignoredReads << "Read " << alignmentName << ": " << errorMessage << "\n";
		}
	};

	void removeFromReadList(TAlignment & alignment, const std::string & errorMessage){
		readList.erase(alignment.name());
		_write(alignment.name(), errorMessage, alignment.isSecondMate);
	};

	bool isInReadList(const std::string & alignmentName){
		if(readList.count(alignmentName) > 0)
			return true;
		return false;
	};

};

////TODO: this blacklist does nothing!!! write blacklist class that is then used in readAlignment
//if(params.parameterExists("blacklist")){
//	//open blacklist file
//	std::string blacklist = params.getParameterString("blacklist");
//	logfile->listFlush("Reading reads to be omitted from '" + blacklist + "...");
//	std::ifstream file(blacklist.c_str());
//	if(!file) throw "Failed to open file '" + blacklist + "!";
//
//	int lineNum = 0;
//	std::vector<std::string> vec;
//
//	//fill list of reads to omit
//	while(file.good() && !file.eof()){
//		++lineNum;
//		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
//		if(!vec.empty())
//			addToReadList(vec[0], "user-defined blacklist");
//	}
//	logfile->write("done! Read " + toString(lineNum) + " read names");
//}


#endif /* TREADLIST_H_ */
