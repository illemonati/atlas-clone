/*
 * TVcfTest.h
 *
 *  Created on: Nov 27, 2018
 *      Author: vivian
 */

#ifndef TESTUTILITIES_TVCFTEST_H_
#define TESTUTILITIES_TVCFTEST_H_

#include <vector>
#include <map>
#include "../TestUtilities/TAtlasTest.h"

class TAtlasTest_invariantBed:public TAtlasTest{
private:
	std::string filenameTag;
	std::string vcfFileName;
	std::string bedFileName;
	std::ofstream vcf;

	void setVariables(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TTaskList* TaskList);
	void writeTestVcf();
	void writeLineVcf(std::string chr, std::string pos, std::string gt);
	bool checkBedFile();

public:
	TAtlasTest_invariantBed();
	bool run(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TTaskList* TaskList);

};

#endif /* TESTUTILITIES_TVCFTEST_H_ */
