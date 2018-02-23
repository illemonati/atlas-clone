/*
 * TAtlasTestPMD.cpp
 *
 *  Created on: Feb 22, 2018
 *      Author: vivian
 */

#include "TAtlasTestPMD.h"


TAtlasTest_PMDEmpiric::TAtlasTest_PMDEmpiric(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "recalSimulation";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";

}
