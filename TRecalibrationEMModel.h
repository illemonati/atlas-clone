/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include "TQualityMap.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include "TBase.h"


//--------------------------------------------------------------------
// TRecalibrationEMQualityPositionMap
// Look-up tables for position and quality. Only indexes will be stored.
//--------------------------------------------------------------------
class TRecalibrationEMQualityPositionMap{
public:
	int maxQual;
	int maxPos;
	double* eta;
	double* etaSquared;
	double* position;
	double* positionSquared;
	bool initialized;

	TRecalibrationEMQualityPositionMap(){
		initialized = false;
		initialize(500, 255); //TODO: think about default!
	};

	~TRecalibrationEMQualityPositionMap(){
		clear();
	};

	void clear(){
		if(initialized){
			delete[] eta;
			delete[] etaSquared;
			delete[] position;
			delete[] positionSquared;
			initialized = false;
		}
	};

	void initialize(int MaxQual, int MaxPos){
		clear();
		maxQual = MaxQual;
		maxPos = MaxPos;
		eta = new double[maxQual+1];
		etaSquared = new double[maxQual+1];
		position = new double[maxPos+1];
		positionSquared = new double[maxPos+1];
		initialized = true;

		//fill qualities. Use TQualityMap for conversion
		TQualityMap qualiMap;
		for(int q=0; q<=maxQual; q++){
			double eps = qualiMap.qualityToError(q);
			if(eps < 0.0000000001) eps = 0.0000000001;
			else if(eps > 0.9999999999) eps = 0.9999999999;

			eta[q] = log(eps / (1.0 - eps));
			etaSquared[q] = eta[q] * eta[q];
		}

		//fill positions
		for(int p = 0; p<=maxPos; p++){
			position[p] = p;
			positionSquared[p] = p * p;
		}
	};
};


//--------------------------------------------------------------------
// TRecalibrationEMReadData
// Per site data storage
//--------------------------------------------------------------------
struct TRecalibrationEMReadData{
	uint8_t quality;
	uint8_t position;
	float D[4];
	uint8_t context;
	uint16_t readGroup;
	bool isSecond;

	void setD(Base base, double PMD_CT, double PMD_GA){
		switch(base){
			case A: D[0] = 0.0; //geno = AA
					D[1] = 1.0; //geno = CC
					D[2] = 1.0 - PMD_GA; //geno = GG
					D[3] = 1.0; //geno = TT
					break;
			case C: D[0] = 1.0; //geno = AA
					D[1] = PMD_CT; //geno = CC
					D[2] = 1.0; //geno = GG
					D[3] = 1.0; //geno = TT
					break;
			case G: D[0] = 1.0; //geno = AA
					D[1] = 1.0; //geno = CC
					D[2] = PMD_GA; //geno = GG
					D[3] = 1.0; //geno = TT
					break;
			case T: D[0] = 1.0; //geno = AA
					D[1] = 1.0 - PMD_CT; //geno = CC
					D[2] = 1.0; //geno = GG
					D[3] = 0.0; //geno = TT
					break;
			case N:
					D[0] = 0.0;
					D[1] = 0.0;
					D[2] = 0.0;
					D[3] = 0.0;
					break;
		}
	};
};

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
// Object to store for which qualities and positions data is available.
//--------------------------------------------------------------------
class TRecalibrationEMDataTable{
public:
	int numReadGroups;
	int maxQual;
	bool*** qualities; //qualities[readGroup][first/second][quality]
	unsigned int** maxPos; //maxPos[readGroup][first/second]

	TRecalibrationEMDataTable(int NumReadGroups, int MaxQual){
		numReadGroups = NumReadGroups;
		maxQual = MaxQual;

		qualities = new bool**[numReadGroups];
		maxPos = new unsigned int*[numReadGroups];
		for(int rg = 0; rg<numReadGroups; ++rg){
			qualities[rg] = new bool*[2];
			qualities[rg][0] = new bool[maxQual];
			qualities[rg][1] = new bool[maxQual];

			maxPos[rg] = new unsigned int[2];
		}

		clear();
	};

	~TRecalibrationEMDataTable(){
		for(int rg = 0; rg<numReadGroups; ++rg){
			delete[] qualities[rg][0];
			delete[] qualities[rg][1];
			delete[] qualities[rg];

			delete[] maxPos[rg];
		}

		delete[] qualities;
		delete[] maxPos;
	};

	void clear(){
		for(int rg = 0; rg<numReadGroups; ++rg){
			for(int i=0; i<maxQual; ++i){
				qualities[rg][0][i] = 0;
				qualities[rg][1][i] = 0;
			}
			maxPos[rg][0] = 0;
			maxPos[rg][1] = 0;
		}
	};

	void add(TRecalibrationEMReadData & data){
		++qualities[data.readGroup][(int) data.isSecond][data.quality];
		if(maxPos[data.readGroup][data.isSecond] < data.position)
				maxPos[data.readGroup][data.isSecond] = data.position;
	};

};

//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel_Base{
protected:
	int numReadGroups;
	uint8_t* numParamsPerReadGroup;
	//int numParameterSets; //is equal to numReadGroups * 2 if all read groups have parameters for first and second mates
	unsigned int totNumParams;
	int* readGroupShifts;
	TRecalibrationEMQualityPositionMap qualPosMap;

	double** betas; //betas of the model
	double** oldBetas; //use during estimation
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	bool initialized;
	bool EMParamsInitialized;

	void _allocateBetaMemory();
	void _freeBetaMemory();
	double _calcEpsilon(double & eta);

public:

	long numSitesAdded;

	TRecalibrationEMModel_Base();
	virtual ~TRecalibrationEMModel_Base(){ _freeBetaMemory(); };

	virtual void initializeWithValues(std::vector<std::string> & vec, int NumReadGroups){ throw "void initializeWithValues(...) not defined for TRecalibrationEMModel_Base!"; };
	virtual void initializeWithHeader(std::vector<std::string> & vec, int NumReadGroups){ throw "void initializeWithHeader(...) not defined for TRecalibrationEMModel_Base!"; };
	virtual void initializeWithDataTable(TRecalibrationEMDataTable & dataTable){ throw "void initializeWithDataTable(...) not defined for TRecalibrationEMModel_Base!"; };
	void initializeEMParams();
	bool setParams(std::vector<std::string> & vec, int & rg);
	void setEMParamsToZero();
	virtual double calcEpsilon(TRecalibrationEMReadData & data){ throw "double calcEpsilon(TRecalibrationEMReadData & data) not defined for TRecalibrationEMModel_Base!"; };
	virtual void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data){ throw "void addToFandJacobian(...) not defined for TRecalibrationEMModel_Base!"; };
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	void rejectProposedParameters();
	double getSteepestGradient();
	virtual void writeHeader(std::ofstream & out){ throw "void writeHeader(std::ofstream & out) not defined for TRecalibrationEMModel_Base!"; };
	void writeParametersToFile(std::ofstream & out, const uint8_t & readGroup);
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();
	virtual double getErrorRate(TBase & base){ throw "double getErrorRate(TBase & base) not defined for TRecalibrationEMModel_Base!"; };
};

class TRecalibrationEMModel:public TRecalibrationEMModel_Base{
protected:
	void _initialize(int NumReadGroups);

	uint8_t numParams;

public:
	TRecalibrationEMModel();

	void initializeWithValues(std::vector<std::string> & vec, int NumReadGroups);
	void initializeWithHeader(std::vector<std::string> & vec, int NumReadGroups);
	void initializeWithDataTable(TRecalibrationEMDataTable & dataTable);

	virtual double calcEpsilon(TRecalibrationEMReadData & data);
	virtual void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data);
	virtual void writeHeader(std::ofstream & out);
	virtual double getErrorRate(TBase & base);
};

class TRecalibrationEMModelNoContext:public TRecalibrationEMModel{
public:
	TRecalibrationEMModelNoContext();

	double calcEpsilon(TRecalibrationEMReadData & data);
	void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data);
	void writeHeader(std::ofstream & out);
	double getErrorRate(TBase & base);
};

class TRecalibrationEMModelPositionSpecific:public TRecalibrationEMModel_Base{
private:
	int maxPos;
	unsigned int* maxPosPerReadGroup;
	unsigned int numParamsWithoutPositions;

	void _initializeFromVectorSize(std::vector<std::string> & vec, int NumReadGroups);

public:
	TRecalibrationEMModelPositionSpecific();

	void initializeWithValues(std::vector<std::string> & vec, int NumReadGroups);
	void initializeWithHeader(std::vector<std::string> & vec, int NumReadGroups);
	void initializeWithDataTable(TRecalibrationEMDataTable & dataTable);
	double calcEpsilon(TRecalibrationEMReadData & data);
	void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, TRecalibrationEMReadData* & data);
	void writeHeader(std::ofstream & out);
	double getErrorRate(int rg, double originalErrorRate, const int & posInRead, const uint8_t & context);
};



#endif /* TRECALIBRATIONEMMODEL_H_ */
