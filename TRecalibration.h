/*
 * TRecalibration.h
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#ifndef TRecalIBRATION_H_
#define TRecalIBRATION_H_

#include "bamtools/api/BamReader.h"
#include "TSite.h"
//#include "TLog.h"

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
class TBase;
class TRecalibrationEM{
public:
	int numParams;
	double* params;
	double* newParams; //used during EM
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	double maxF; //largest change during Newton-Ralphson
	long numSitesAdded;
	double logLikelihood;

	TRecalibrationEM(TParameters* arguments, TLog* logfile);
	~TRecalibrationEM(){
		delete[] params;
		delete[] newParams;
	};
	TRecalibrationEM(TLog* logfile);
	void setParams(double* Params);
	double calcEta(TBase* base);
	double calcEta(TBase* base, double* theseParams);
	double calcEpsilon(const double & eta);
	double calcEpsilon(TBase* base);
	double calcEpsilon(TBase* base, double* theseParams);
	void initEMStep();
	void initNetwonRalphsonStep();
	void saveParams();
	void addSiteToJacobianAndF(std::vector<TBase*> & bases, TBaseFrequencies* freqs);
	void runNewtonRalphson();
	void writeHeader(std::ofstream & out);
	void writeParams(std::ofstream & out);
	void resetLikelihood();
	void addSiteToLikelihood(std::vector<TBase*> & bases, TBaseFrequencies* freqs);
};

//---------------------------------------------------------------
//RecalibrationBQSR
//---------------------------------------------------------------
//covariates to take into account:
// - read base (A, G, C, T)
// -

class TRecalibrationBQSR_cell{
public:
	double curEpsilon;
	bool estimationConverged;
	double firstDerivative, secondDerivative;
	TPMD* pmdObject;
	long numObservations;


	TRecalibrationBQSR_cell();
	~TRecalibrationBQSR_cell(){};
	void init(double initialError, TPMD* PmdObject);
	void addBase(TBase* base, Base & RefBase);
	bool estimateEpsilon(double & convergenceThreshold);
};

class TRecalibrationBQSR{
public:
	TLog* logfile;
	BamTools::SamHeader* bamHeader;
	TPMD* pmdObject;
	int numReadGroups;
	int minQ, maxQ, numQ;
	double initialError;
	double convergenceThreshold;
	bool estimationConverged;
	TRecalibrationBQSR_cell** BQSR_cells; //read group x quality

	TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TPMD* PmdObject, TLog* Logfile);
	~TRecalibrationBQSR(){};

	void addSite(TSite & site);
	bool estimateEpsilon();
	void writeToFile(std::string filename);
};


/*
class TRecalibrationBQSR_cellA:public TRecalibrationBQSR_cell{
public:
	TPMD* pmdObject;



	TRecalibrationBQSR_cellA(double Error, TPMD* PmdObject, double ConvergenceThreshold);
	~TRecalibrationBQSR_cellA(){};

	double getD(TBase* base, char & RefBase);
	void addBase(TBase* base, char & RefBase);
	bool estimateEpsilon();
};


class TRecalibrationBQSR_cellC:public TRecalibrationBQSR_cell{
public:
	long N_1, N_2;
	double D;
	TRecalibrationBQSR_cellC(double Error, int pos, TPMD* PmdObject);
	~TRecalibrationBQSR_cellC(){};

	void addBase(TBase* base, char & RefBase);
	bool estimateEpsilon();
};

class TRecalibrationBQSR_cellT:public TRecalibrationBQSR_cellC{
public:
	long N_3;

	TRecalibrationBQSR_cellT(double Error, int pos, TPMD* PmdObject);
	~TRecalibrationBQSR_cellT(){};

	void addBase(TBase* base, char & RefBase);
	bool estimateEpsilon();
};


class TRecalibrationBQSR_cellG:public TRecalibrationBQSR_cellA{
public:
	TRecalibrationBQSR_cellG(double Error, TPMD* PmdObject, double ConvergenceThreshold);
	~TRecalibrationBQSR_cellG(){};

	double getD(TBase* base, char & RefBase);
};

*/


#endif /* TRecalIBRATION_H_ */
