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
	std::string poolRGFileName;
	std::ofstream outPMD;
	std::ofstream outPool;
	float alpha, beta;
	int minReadLength, maxReadLength;
	std::string firstPMDStringCT, firstPMDStringGA, secondPMDStringCT, secondPMDStringGA, thirdPMDStringCT, thirdPMDStringGA;
	std::string CTpatterns[3];
	std::string GApatterns[3];

	void defineVariables(TParameters & params, TLog* Logfile);
	bool checkPMDEmpiricFile();

public:
	TAtlasTest_PMDEmpiric();
	bool run(TParameters & params, TLog* Logfile);

};



#endif /* TATLASTESTPMD_H_ */
