/*
 * TAtlasTest.h
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */

#ifndef TATLASTEST_H_
#define TATLASTEST_H_

#include "TTaskList.h"
#include "../commonutilities/IntegrationTests/TTest.h"
#include "../commonutilities/gzstream.h"
#include "../bamtools/api/BamWriter.h"
#include "../bamtools/api/BamReader.h"
#include "../bamtools/api/SamHeader.h"
#include "../bamtools/api/BamAlignment.h"
#include "../TGenotypeMap.h"
#include "../TQualityMap.h"
#include <vector>
#include <map>
//------------------------------------------
//TAtlasTest
//------------------------------------------
//Base class for individual tests.
//Tests can be combined to suites in TATlasTesting

//class TAtlasTest{
//protected:
//	TLog* logfile;
//	TParameters _testParams;
//	std::string _testingPrefix;
//	std::string _name;
//
//	bool runTGenomeFromInputfile(std::string task);
//
//public:
//	TAtlasTest(TParameters & params, TLog* Logfile);
//	virtual ~TAtlasTest(){};
//
//	std::string name(){return _name;};
//
//	virtual bool run(){
//		return true;
//	};
//};

//------------------------------------------
//TAtlasTest_pileup
//------------------------------------------
class TAtlasTest_pileup:public TTest{
private:
	TLog* logfile;
	std::string _testingPrefix;

	int phredError;
	int readLength;
	int chrLength;
	std::vector<int> depths;
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	std::string filenameTag;
	std::string bamFileName;
	std::string fastaName;
	std::string readGroupName;
	double emissionTolerance; //relative error allowed to accommodate rounding issues when reading numbers from file

	void defineVariables(TParameters & params, TLog* Logfile);
	void writeFasta();
	void writeBAM();
	bool checkPileupFile();

public:
	TAtlasTest_pileup();
	~TAtlasTest_pileup(){};
	bool run(TParameters & parameters, TLog* Logfile, TTaskList * TaskList);
};

//------------------------------------------
//TAtlasTest_allelicDepth
//------------------------------------------
class TAtlasTest_allelicDepth:public TTest{
private:
	TLog* logfile;
	std::string _testingPrefix;

	int phredError;
	TQualityMap qualMap;
	std::string filenameTag;
	std::string bamFileName;
	std::string readGroupName;

	void defineVariables(TParameters & params, TLog* Logfile);
	void writeBAM();
	bool checkAllelicDepthTable();

public:
	TAtlasTest_allelicDepth();
	~TAtlasTest_allelicDepth(){};
	bool run(TParameters & params, TLog* Logfile, TTaskList * TaskList);
};

//------------------------------------------
//TAtlasTest_theta
//------------------------------------------
class TAtlasTest_theta:public TTest{
private:
	TLog* logfile;
	std::string _testingPrefix;

	std::string filenameTag;
	std::string bamFileName;
	double simTheta;

	void defineVariables(TParameters & params, TLog* Logfile);
	bool checkThetaFile();

public:
	TAtlasTest_theta();
	~TAtlasTest_theta(){};
	bool run(TParameters & params, TLog* Logfile, TTaskList * TaskList);
};

#endif /* TATLASTEST_H_ */
