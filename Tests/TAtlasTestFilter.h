/*
 * TAtlasTestFilter.h
 *
 *  Created on: May 8, 2019
 *      Author: vivian
 */

#ifndef TESTS_TATLASTESTFILTER_H_
#define TESTS_TATLASTESTFILTER_H_

#include "TAtlasTest.h"
#include "../TAlignmentParser.h"

class TAtlasTest_filter:public TAtlasTest{
private:
	std::string filenameTag;
	std::string bamFileName;
	int readLength;
	int chrLength;
	int phredError;
	std::string readGroupName;
	bool keepReadsLongerThanFragment;
	bool keepOrphanedReads;
	bool keepImproperPairs;
	bool keepUnmappedReads;
	bool keepFailedQC;
	bool keepSecondary;
	bool keepSupplementary;
	bool keepDuplicates;
	bool filterSoftClips;
	TQualityMap qualMap;
	std::vector<std::string> shouldKeep;
	std::vector<std::string> trueIgnoredReadMessages;

	void writeBAM();
	void setToSingleEnd(BamTools::BamAlignment & bamAlignment);
	void setToProperPairEtc(BamTools::BamAlignment & bamAlignment);
	void setToFwdMate(BamTools::BamAlignment & bamAlignment);
	void setToRevMate(BamTools::BamAlignment & bamAlignment);
	bool checkfilteredBAMFile();



public:
	TAtlasTest_filter(TParameters & params, TLog* logfile);
	bool run();
};


#endif /* TESTS_TATLASTESTFILTER_H_ */
