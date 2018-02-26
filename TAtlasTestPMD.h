/*
 * TAtlasTestPMD.h
 *
 *  Created on: Feb 22, 2018
 *      Author: vivian
 */

#ifndef TATLASTESTPMD_H_
#define TATLASTESTPMD_H_

#include "TAtlasTest.h"


class TAtlasTest_PMDEmpiric:public TAtlasTest{
private:
	std::string filenameTag;
	std::string bamFileName;
	std::string fastaFileName;
	std::string pmdEmpiricFileName;
	float alpha, beta;
	int minReadLength, maxReadLength;
	std::string firstPMDStringCT, firstPMDStringGA, secondPMDStringCT, secondPMDStringGA, thirdPMDStringCT, thirdPMDStringGA;
	std::string CTpatterns[3];
	std::string GApatterns[3];
	std::ofstream out;

	bool checkPMDEmpiricFile();

public:
	TAtlasTest_PMDEmpiric(TParameters & params, TLog* logfile);
	bool run();

};



#endif /* TATLASTESTPMD_H_ */
