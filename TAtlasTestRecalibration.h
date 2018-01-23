/*
 * TAtlasTestRecalibration.h
 *
 *  Created on: Jan 22, 2018
 *      Author: linkv
 */

#ifndef TATLASTESTRECALIBRATION_H_
#define TATLASTESTRECALIBRATION_H_

#include "TAtlasTest.h"

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


#endif /* TATLASTESTRECALIBRATION_H_ */
