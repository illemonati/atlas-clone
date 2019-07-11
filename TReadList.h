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

    void _write(const std::string & Name, const std::string & Error, const bool & isSecondMate);

public:
    TReadList();
    TReadList(const bool & writeToFile, std::string & filename);

    void setWriteReadListToFileToTrue(TLog* logfile, std::string & Filename);
    void addToReadList(TAlignment & alignment, const std::string & errorMessage);
    void addToReadList(std::string & alignmentName, const std::string & errorMessage);
    void removeFromReadList(TAlignment & alignment, const std::string & errorMessage);
    bool isInReadList(const std::string & alignmentName);

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
