/*
 * TDistanceCalculator.h
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#ifndef TDISTANCEESTIMATOR_H_
#define TDISTANCEESTIMATOR_H_

#include "TParameters.h"
#include "TGLF.h"
#include <math.h>
#include "TGenotypeMap.h"
#include "TQualityMap.h"

//--------------------------------------------
//TPhiToGenoMap
//--------------------------------------------
class TGenoToPhiMap{
public:
	TGenotypeMap genoMap;
	int**   genoToPhiMap;

	TGenoToPhiMap();
	~TGenoToPhiMap(){
		for(int i=0; i<10; ++i)
			delete[] genoToPhiMap[i];
		delete[] genoToPhiMap;
	};

	int& operator()(int & g1, int & g2){
		return genoToPhiMap[g1][g2];
	};
	int& operator()(Genotype & g1, Genotype & g2){
		return genoToPhiMap[g1][g2];
	};
};

//--------------------------------------------
//TGenocombinationToBaseMap
//--------------------------------------------
class TGenocombinationToBaseMap{
public:
	TGenotypeMap genoMap;
	bool*** genotypeCombinationHasBase;

	TGenocombinationToBaseMap();
	~TGenocombinationToBaseMap(){
		for(int i=0; i<10; ++i){
			for(int j=0; j<10; ++j)
				delete[] genotypeCombinationHasBase[i][j];
			delete[] genotypeCombinationHasBase[i];
		}
		delete[] genotypeCombinationHasBase;
	};

	bool& operator()(int & g1, int & g2, int & base){
		return genotypeCombinationHasBase[g1][g2][base];
	};
	bool& operator()(Genotype & g1, Genotype & g2, Base & b){
		return genotypeCombinationHasBase[g1][g2][b];
	};
};

//----------------------------------------------------
//TDistanceClass
//----------------------------------------------------
class TDistance{
public:
	TDistance();
	virtual ~TDistance(){
		delete[] distanceWeight;
	}
	double* distanceWeight; //weight for each phi class towards the distance.
	virtual double calculateDistance(double* phi);
};

class TDistanceProbMismatch:public TDistance{
public:
	TDistanceProbMismatch();
};

class TDistanceEuclidian:public TDistance{
public:
	double calculateDistance(double* phi);
};

class TDistanceUser:public TDistance{
public:
	TDistanceUser(std::vector<double> vec);
};

//--------------------------------------------
//TDistanceEstimate
//--------------------------------------------
class TEMforDistanceEstimation{
private:
	TLog* logfile;
	TGenotypeMap genoMap;
	TGenoToPhiMap genoToPhiMap;
	TGenocombinationToBaseMap genoToBaseMap;
	//TQualityMap phredToLik;
	TGlfConverter glfConverter;

	//settings
	int maxNumEMIterations;
	double epsilonForEM;

	//tmp variables
	double old_LL;
	double* K; //normalizing constant
	double** probGeno;
	double** P_G;
	double** P_G_one_site;
//	double* distanceWeight; //weight for each phi class towards the distance.
	TDistance* distanceObject;

//	void calculateDistance();
	void guessPi(std::vector<uint16_t*> & genoQual1, std::vector<uint16_t*> & genoQual2);
	void guessPhi(std::vector<uint16_t*> & genoQual1, std::vector<uint16_t*> & genoQual2);
	void fill_K(TBaseFrequencies  & thesePi);
	void fill_P_g_given_phi_pi(double* phi, TBaseFrequencies & pi);

public:
	TBaseFrequencies pi;
	double* phi;
	double LL;
	double distance;

	TEMforDistanceEstimation(TLog* Logfile, TParameters & params);
	~TEMforDistanceEstimation(){
		delete[] phi;
		delete[] K;
		for(int g1=0; g1<10; ++g1){
			delete[] probGeno[g1];
			delete[] P_G[g1];
			delete[] P_G_one_site[g1];
		}
		delete[] probGeno;
		delete[] P_G;
		delete[] P_G_one_site;
		delete distanceObject;
//		delete[] distanceWeight;
	};

	bool estimatePhiWithEM(std::vector<uint16_t*> & genoQual1, std::vector<uint16_t*> & genoQual2);
};

//--------------------------------------------
//TDistanceEstimator
//--------------------------------------------
class TDistanceEstimator{
private:
	TLog* logfile;
	int maxNumEMIterations;
	double epsilonForEM;
	std::string outputName;

	//GLF files
	int numGLFs;
	std::vector<std::string> GLFNames;
	TGlfReader* glfs;
	bool readersOpened;

	void openGLF(TParameters & params);
	void closeGLF();


	void estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object);
	bool moveToNextCommonChr(TGlfReader & g1, TGlfReader & g2);
	bool advance(TGlfReader & g1, TGlfReader & g2);
	void readCommonSites(std::vector<uint16_t*> & genoQual1, std::vector<uint16_t*> & genoQual2, TGlfReader & g1, TGlfReader & g2);
	void estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, TGlfReader & g1, TGlfReader & g2, gz::ogzstream & out);

	void estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, long windowLen);
	void estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, std::string filename, TGlfReader & g1, TGlfReader & g2, long windowLen);

	void writeDistanceEstimates(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd, int numsitesWithData, TEMforDistanceEstimation & EM_object);
	void writeDistanceEstimates(gz::ogzstream & out, int numsitesWithData, TEMforDistanceEstimation & EM_object);
	void writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd);
	void writeDistanceEstimatesNoData(gz::ogzstream & out);

public:
	TDistanceEstimator(TLog* Logfile, TParameters & params);
	~TDistanceEstimator(){
		closeGLF();
	};

	void printGLF(TParameters & params);
	void estimateDistances(TParameters & params);

};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimateDist:public TTask{
public:
	TTask_estimateDist(){ _explanation = "Estimating the genetic distance between individuals"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TDistanceEstimator distEst(Logfile, Parameters);
		distEst.estimateDistances(Parameters);
	};
};


#endif /* TDISTANCEESTIMATOR_H_ */
