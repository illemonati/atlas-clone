/*
 * TVcfTest.h
 *
 *  Created on: Nov 27, 2018
 *      Author: vivian
 */

#ifndef TVCFTEST_H_
#define TVCFTEST_H_

#include "atlasTaskSwitcher.h"
#include "TAtlasTest.h"

#include <vector>
#include <map>

class TAtlasTest_invariantBed:public TAtlasTest{
private:
	std::string filenameTag;
	std::string vcfFileName;
	std::string bedFileName;
	std::ofstream vcf;

	void writeTestVcf();
	void writeLineVcf(std::string chr, std::string pos, std::string gt);
	bool checkBedFile();

public:
	TAtlasTest_invariantBed(TParameters & params, TLog* logfile);
	bool run();

};

#endif /* TVCFTEST_H_ */
