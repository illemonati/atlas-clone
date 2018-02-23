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
public:
	TAtlasTest_PMDEmpiric(TParameters & params, TLog* logfile);

};



#endif /* TATLASTESTPMD_H_ */
