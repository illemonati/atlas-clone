/*
 * TAtlasTestRecalibration.h
 *
 *  Created on: Jan 22, 2018
 *      Author: linkv
 */

#ifndef TATLASTESTRECALIBRATION_H_
#define TATLASTESTRECALIBRATION_H_

#include "TAtlasTest.h"
#include "../simulation/TSimulatorQualityTransformation.h"

//------------------------------------------
//TAtlasTest_recalSimulation
//------------------------------------------
class TAtlasTest_recalSimulation:public TAtlasTest{
private:
	std::string filenameTag;
	std::string bamFileName;
	std::string fastaFileName;
	int meanQual;
	double sdphredInt;
	int minPhredInt;
	int maxPhredInt;
	std::string qualityDist;
	std::string recalParamString;
	std::string recalParamsFileName;
	std::string poolRGFileName;
	std::ofstream outRecalParams;
	std::ofstream outRecalPool;

public:
	TAtlasTest_recalSimulation(TParameters & params, TLog* logfile);
	~TAtlasTest_recalSimulation(){};
	bool run();
	bool checkRecalFile();
};

//------------------------------------------
//TAtlasTest_BQSRSimulation
//------------------------------------------
class TAtlasTest_BQSRSimulation:public TAtlasTest{
private:
	std::string filenameTag;
	std::string bamFileName;
	std::string fastaFileName;
	int meanQual;
	double sdphredInt;
	int minPhredInt;
	int maxPhredInt;
	int phi1;
	double phi2;
	double revIntercept;
	double acceptedDelta;
	std::string qualityDist;
//	double alpha;
//	double beta;
	int minReadLen;
	int maxReadLen;
//	std::string readLengthDist;
	double positionEffectSlope;
	double positionEffectIntercept;


public:
	TAtlasTest_BQSRSimulation(TParameters & params, TLog* logfile);
	~TAtlasTest_BQSRSimulation(){};
	bool run();
	double trueQual(int & phi1, double & phi2, int & fakeQual);
	bool checkBQSRQualityFile();
	bool checkBQSRPositionFile();
	void calculateSlopeIntercept();
	double trueScaling(int & pos);

};

//------------------------------------------
//TAtlasTest_qualityTransformation
//------------------------------------------
class TAtlasTest_qualityTransformationRecalPlain:public TAtlasTest{
private:
	std::vector<double> trueParams;
	std::vector<std::vector<double>> qualTransTable;

protected:
	std::string filenameTag;
	std::string bamFileName;
	std::string recalParamString;
	int maxReadLength;
	TRandomGenerator* randomGenerator;
	std::string qualDistString;
	std::vector<int> qualDistVec;
	TSimulatorQualityDist* qualityDist;
	TSimulatorQualityTransformationRecal* recalObject;
	bool readTransformationFile();
	bool checkTransformation(std::vector<int> trueQualScores);

public:
	virtual bool run();
	TAtlasTest_qualityTransformationRecalPlain(TParameters & params, TLog* logfile);
	virtual ~TAtlasTest_qualityTransformationRecalPlain(){
		delete randomGenerator;
		delete qualityDist;
		delete recalObject;
	};
};

//------------------------------------------

class TAtlasTest_qualityTransformationRecalBinned:public TAtlasTest_qualityTransformationRecalPlain{
public:
	bool run();
	TAtlasTest_qualityTransformationRecalBinned(TParameters & params, TLog* logfile);
	~TAtlasTest_qualityTransformationRecalBinned(){};
};


#endif /* TATLASTESTRECALIBRATION_H_ */
