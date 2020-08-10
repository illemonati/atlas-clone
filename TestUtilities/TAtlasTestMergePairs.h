/*
 * TAtlasTestMergePairs.h
 *
 *  Created on: Jan 16, 2019
 *      Author: linkv
 */

#ifndef TATLASTESTMERGEPAIRS_H_
#define TATLASTESTMERGEPAIRS_H_

#include "../TestUtilities/TAtlasTest.h"

class TAtlasTest_mergePairs:public TAtlasTest{
protected:
	std::string filenameTag;
	std::string bamFileName;
	std::string readGroupName;
	int readLength;
	int chrLength;
	int phredError;
	bool filterOrphanedReads;
	BAM::TQualityMap qualMap;
	std::vector<std::string> trueQueryBases;
	std::vector<std::string> trueQualities;
	std::vector<bool> trueIsProper;
	std::vector<std::string> trueIgnoredReadMessages;

	void setToProperPairEtc(BamTools::BamAlignment & bamAlignment);
	void setToSingleEndEtc(BamTools::BamAlignment & bamAlignment);
	void setToFwdMate(BamTools::BamAlignment & bamAlignment);
	void setToRevMate(BamTools::BamAlignment & bamAlignment);
	void writeAlignments(BamTools::BamWriter & bamWriter);

	bool basicChecks(BamTools::BamAlignment & bamAlignment, const int pairNumber);
	bool checkMergedBAMFile();
	void writePairedEndReads(BamTools::BamWriter & bamWriter);
	void writeBAM();
	void setVariables(TParameters & params, TLog* Logfile, TTaskList * TaskList);

public:
	TAtlasTest_mergePairs();
	virtual ~TAtlasTest_mergePairs(){};
	virtual bool run(TParameters & params, TLog* logfile, TTaskList * TaskList);
};

//class TAtlasTest_mergePairsAdaptQual:public TAtlasTest_mergePairs{
//
//};

class TAtlasTest_mergeSplitPairs:public TAtlasTest_mergePairs{
private:
	void writeSingleEndReads(BamTools::BamWriter & bamWriter);
	void writeAlignments(BamTools::BamWriter & bamWriter);
	void writeRGSpecsFile();
	void setVariables(TParameters & params, TLog* Logfile, TTaskList* TaskList);

public:
	TAtlasTest_mergeSplitPairs();
	~TAtlasTest_mergeSplitPairs(){};
	bool run(TParameters & params, TLog* Logfile, TTaskList* TaskList);
};

#endif /* TATLASTESTMERGEPAIRS_H_ */
