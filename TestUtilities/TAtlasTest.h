/*
 * TAtlasTest.h
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */

#ifndef TATLASTEST_H_
#define TATLASTEST_H_

#include "TTaskList.h"
#include "TTest.h"
#include "gzstream.h"
#include "../bamtools/api/BamWriter.h"
#include "../bamtools/api/BamReader.h"
#include "../bamtools/api/SamHeader.h"
#include "../bamtools/api/BamAlignment.h"
#include <vector>
#include <map>


//------------------------------------------
//TAtlasTest
//------------------------------------------
//Base class for individual tests.
//Tests can be combined to suites in TATlasTesting

class TAtlasTest : public coretools::TTest{
protected:
	std::string _testingPrefix;

public:
	TAtlasTest();
	virtual ~TAtlasTest(){};

	std::string name(){return _name;};
};

//------------------------------------------
//TAtlasTest_pileup
//------------------------------------------
class TAtlasTest_pileup:public TAtlasTest{
private:
	int phredError;
	int readLength;
	int chrLength;
	std::vector<int> depths;
	std::string filenameTag;
	std::string bamFileName;
	std::string fastaName;
	std::string readGroupName;
	double emissionTolerance; //relative error allowed to accommodate rounding issues when reading numbers from file

	void setVariables(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TTaskList* TaskList);
	void writeFasta();
	void writeBAM();
	bool checkPileupFile();

public:
	TAtlasTest_pileup();
	~TAtlasTest_pileup(){};
	bool run(coretools::TParameters & parameters, coretools::TLog* Logfile, coretools::TTaskList* TaskList);
};

//------------------------------------------
//TAtlasTest_allelicDepth
//------------------------------------------
class TAtlasTest_allelicDepth:public TAtlasTest{
private:
	int phredError;
	std::string filenameTag;
	std::string bamFileName;
	std::string readGroupName;

	void setVariables(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TTaskList* TaskList);
	void writeBAM();
	bool checkAllelicDepthTable();

public:
	TAtlasTest_allelicDepth();
	~TAtlasTest_allelicDepth(){};
	bool run(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TTaskList * TaskList);
};

//------------------------------------------
//TAtlasTest_theta
//------------------------------------------
class TAtlasTest_theta:public TAtlasTest{
private:
	std::string filenameTag;
	std::string bamFileName;
	double simTheta;

	void defineVariables(coretools::TParameters & params, coretools::TLog* Logfile);
	bool checkThetaFile();

public:
	TAtlasTest_theta();
	~TAtlasTest_theta(){};
	bool run(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TTaskList * TaskList);
};

#endif /* TATLASTEST_H_ */
