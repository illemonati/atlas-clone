/*
 * TVcfTest.h
 *
 *  Created on: Nov 27, 2018
 *      Author: vivian
 */

#ifndef TESTS_TVCFTEST_H_
#define TESTS_TVCFTEST_H_

#include "TAtlasTest.h"

#include <vector>
#include <map>

class TAtlasTest_invariantBed:public TAtlasTest{
private:
	std::string filenameTag;
	std::string vcfFileName;
	std::string bedFileName;
	std::ofstream vcf;

	void setVariables(TParameters & params, TLog* Logfile, TTaskList* TaskList);
	void writeTestVcf();
	void writeLineVcf(std::string chr, std::string pos, std::string gt);
	bool checkBedFile();

public:
	TAtlasTest_invariantBed();
	bool run(TParameters & params, TLog* Logfile, TTaskList* TaskList);

};

#endif /* TESTS_TVCFTEST_H_ */
