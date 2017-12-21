/*
 * TAtlasTest.h
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */

#ifndef TATLASTEST_H_
#define TATLASTEST_H_

#include "atlasTaskSwitcher.h"

#include <vector>
#include <map>

#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamHeader.h"
#include "bamtools/api/BamAlignment.h"

//------------------------------------------
//TAtlasTest
//------------------------------------------
//Base class for individual tests.
//Tests can be combined to suites in TATlasTesting

class TAtlasTest{
protected:
	TLog* logfile;
	TParameters _testParams;
	std::string _testingPrefix;
	std::string _name;

	bool runTGenomeFromInputfile(std::string task);

public:
	TAtlasTest(TParameters & params, TLog* Logfile);
	virtual ~TAtlasTest(){};

	std::string name(){return _name;};

	virtual bool run(){
		return true;
	};
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
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	std::string filenameTag;
	std::string bamFileName;
	std::string readGroupName;
	double emissionTolerance; //relative error allowed to accommodate rounding issues when reading numbers from file

	void writeBAM();
	bool checkPileupFile();

public:
	TAtlasTest_pileup(TParameters & params, TLog* logfile);
	~TAtlasTest_pileup(){};
	bool run();
};

//------------------------------------------
//TAtlasTest_BQSRSimulation
//------------------------------------------
class TAtlasTest_BQSRSimulation:public TAtlasTest{
private:
	std::string filenameTag;
	std::string bamFileName;
	std::string fastaFileName;
	std::string qualityDist;
	int meanQual;
	double sdQual;
	int minQual;
	int maxQual;
	int phi1;
	double phi2;
	double revIntercept;
	double acceptedDelta;


public:
	TAtlasTest_BQSRSimulation(TParameters & params, TLog* logfile);
	~TAtlasTest_BQSRSimulation(){};
	bool run();
	int trueQual(int & phi1, double & phi2, int & fakeQual);
	bool checkBQSRFile();

};

#endif /* TATLASTEST_H_ */
