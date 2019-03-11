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
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel_Base{
protected:
	int numParameters;
	int myShift;
	TRecalibrationEMQualityPositionMap qualPosMap;
	long numSitesAdded;

	double* betas; //betas of the model
	double* oldBetas; //use during estimation
	bool initialized;

	void _allocateBetaMemory();
	void _freeBetaMemory();
	double _calcEpsilon(double & eta);

public:

	TRecalibrationEMModel_Base(int Shift);
	virtual ~TRecalibrationEMModel_Base(){ _freeBetaMemory(); };

	int getNextShift(){ return myShift + numParameters; };
	bool setParams(std::vector<std::string> & vec);
	void proposeNewParameters(double & lambda, arma::mat & JxF);
	void rejectProposedParameters();
	virtual double calcEpsilon(TRecalibrationEMReadData & data){ throw "double calcEpsilon(TRecalibrationEMReadData & data) not defined for TRecalibrationEMModel_Base!"; };
	virtual void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMReadData* & data, const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta){ throw "void addToFandJacobian(...) not defined for TRecalibrationEMModel_Base!"; };
	virtual void writeHeader(std::ofstream & out){ throw "void writeHeader(std::ofstream & out) not defined for TRecalibrationEMModel_Base!"; };
	void writeParametersToFile(std::ofstream & out);
	virtual double getErrorRate(TBase & base){ throw "double getErrorRate(TBase & base) not defined for TRecalibrationEMModel_Base!"; };
};

class TRecalibrationEMModel_noRecal:public TRecalibrationEMModel_Base{
public:
	TRecalibrationEMModel_noRecal(int Shift);
	virtual double getErrorRate(TBase & base);
};

class TRecalibrationEMModel_qualFuncPosFunc:public TRecalibrationEMModel_Base{
protected:
	void _initialize(int NumReadGroups);

public:
	TRecalibrationEMModel_qualFuncPosFunc(int Shift);
	TRecalibrationEMModel_qualFuncPosFunc::TRecalibrationEMModel_qualFuncPosFunc(std::vector<std::string> & vec, int Shift);

	virtual double calcEpsilon(TRecalibrationEMReadData & data);
	virtual void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMReadData* & data, const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta);
	virtual void writeHeader(std::ofstream & out);
	virtual double getErrorRate(TBase & base);
};

class TRecalibrationEMModel_qualFuncPosFuncNoContext:public TRecalibrationEMModel_qualFuncPosFunc{
public:
	TRecalibrationEMModel_qualFuncPosFuncNoContext(int Shift);
	TRecalibrationEMModel_qualFuncPosFuncNoContext(std::vector<std::string> & vec, int Shift);

	double calcEpsilon(TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMReadData* & data, const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta);
	void writeHeader(std::ofstream & out);
	double getErrorRate(TBase & base);
};

class TRecalibrationEMModel_qualFuncPosSpecific:public TRecalibrationEMModel_Base{
private:
	int maxPos;
	int numParamsWithoutPositions;

public:
	TRecalibrationEMModel_qualFuncPosSpecific(int Shift, int MaxPos);
	TRecalibrationEMModel_qualFuncPosSpecific(std::vector<std::string> & vec, int Shift);

	double calcEpsilon(TRecalibrationEMReadData & data);
	void addToFandJacobian(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMReadData* & data, const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta);
	void writeHeader(std::ofstream & out);
	double getErrorRate(TBase & base);
};



#endif /* TRECALIBRATIONEMMODEL_H_ */
