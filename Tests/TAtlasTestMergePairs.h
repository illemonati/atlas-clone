/*
 * TAtlasTestMergePairs.h
 *
 *  Created on: Jan 16, 2019
 *      Author: linkv
 */

#ifndef TATLASTESTMERGEPAIRS_H_
#define TATLASTESTMERGEPAIRS_H_

#include "TAtlasTest.h"

class TAtlasTest_mergePairs:public TAtlasTest{
private:
	std::string _name;
	std::string filenameTag;
	std::string bamFileName;
	std::string readGroupName;
	int readLength;
	int chrLength;
	int phredError;
//	TGenotypeMap genoMap;
	TQualityMap qualMap;
	std::vector<std::string> trueQueryBases;
	std::vector<std::string> trueQualities;

	void setToProperPairEtc(BamTools::BamAlignment & bamAlignment);
	bool basicChecks(BamTools::BamAlignment & bamAlignment, const int pairNumber);
	bool checkMergedBAMFile();
	void writeBAM();

public:
	TAtlasTest_mergePairs(TParameters & params, TLog* logfile);
	bool run();
};



#endif /* TATLASTESTMERGEPAIRS_H_ */
