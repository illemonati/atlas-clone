/*
 * TAtlasTestRecalibration.h
 *
 *  Created on: Jan 22, 2018
 *      Author: linkv
 */

#ifndef TATLASTESTRECALIBRATION_H_
#define TATLASTESTRECALIBRATION_H_

#include "TAtlasTest.h"
#include "../Simulations/TSimulatorQualityTransformation.h"

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
	std::vector<std::string> tmpVec;
	std::vector<double> trueParams;
	std::string recalParamsFileName;
	std::string poolRGFileName;
	std::ofstream outRecalParams;
	std::ofstream outRecalPool;

	void defineVariables(TParameters & params, TLog* Logfile);

public:
	TAtlasTest_recalSimulation();
	~TAtlasTest_recalSimulation(){};
	bool run(TParameters & params, TLog* logfile);
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

	void setVariables(TParameters & params, TLog* logfile);


public:
	TAtlasTest_BQSRSimulation();
	~TAtlasTest_BQSRSimulation(){};
	bool run(TParameters & params, TLog* logfile);
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
	void setVariables(TParameters & params, TLog* Logfile);

public:
	TAtlasTest_qualityTransformationRecalPlain();
	virtual ~TAtlasTest_qualityTransformationRecalPlain(){
		delete randomGenerator;
		delete qualityDist;
		delete recalObject;
	};
	virtual bool run(TParameters & params, TLog* logfile);
};

//------------------------------------------

class TAtlasTest_qualityTransformationRecalBinned:public TAtlasTest_qualityTransformationRecalPlain{
private:
	void setVariables(TParameters & params, TLog* Logfile);

public:
	TAtlasTest_qualityTransformationRecalBinned();
	~TAtlasTest_qualityTransformationRecalBinned(){};

	bool run(TParameters & params, TLog* logfile);

};


#endif /* TATLASTESTRECALIBRATION_H_ */
