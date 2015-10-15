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

class TQualityIndex{
public:
	int minQ, maxQ, numQ, last;
	int* index;

	TQualityIndex(int MinQ, int MaxQ){
		minQ= MinQ;
		maxQ = MaxQ;
		numQ = maxQ - minQ + 1;
		last = numQ - 1;

		//fill index
		index = new int[maxQ + 1];
		for(int i=0; i < maxQ + 1; ++i){
			if(i < minQ) index[i] = 0;
			else index[i] = i - minQ;
		}
	};

	int& getIndex(int & quality){
		if(quality > maxQ) return last;
		return index[quality];
	};
};

class TBQSR_cellQuality{
public:
	double curEstimate;
	bool estimationConverged;
	double firstDerivative, secondDerivative;
	TPMD* pmdObject;
	long numObservations;
	long numMatches;
	double F;

	TBQSR_cellQuality();
	virtual ~TBQSR_cellQuality(){};
	void empty();
	void init(double initialError, TPMD* PmdObject);
	void addToDerivatives(TBase* base, Base & RefBase, double & epsilon);
	virtual void addBase(TBase* base, Base & RefBase);
	bool estimate(double & convergenceThreshold);
};

class TBQSR_cellPosition:public TBQSR_cellQuality{
public:
	TBQSR_cellQuality** BQSR_cells_readGroup_quality; //read group x quality
	TQualityIndex* qualityIndex;

	TBQSR_cellPosition();
	~TBQSR_cellPosition(){};

	void init(TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, TPMD* PmdObject);
	void addBase(TBase* base, Base & RefBase);
};


class TBQSR_cellContext:public TBQSR_cellPosition{
public:
	TBQSR_cellPosition** BQSR_quality_position; //read group x quality
	bool considerPosition;

	TBQSR_cellContext();
	~TBQSR_cellContext(){};

	void init(TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_quality_position, TQualityIndex* QualityIndex, TPMD* PmdObject);
	void init(TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, TPMD* PmdObject);
	void addBase(TBase* base, Base & RefBase);
};


class TRecalibrationBQSR{
public:
	TLog* logfile;
	BamTools::SamHeader* bamHeader;
	TPMD* pmdObject;
	int numReadGroups;
	TQualityIndex* qualityIndex;
	double initialError;
	double convergenceThreshold;
	bool estimationConverged;
	int maxPos;
	int numContexts;

	//recal tables
	bool qualityConverged;
	TBQSR_cellQuality** BQSR_cells_readGroup_quality; //read group x quality
	bool considerPosition, positionConverged;
	TBQSR_cellPosition** BQSR_cells_readGroup_position; //read group x position
	bool considerContext, contextConverged;
	TBQSR_cellContext** BQSR_cells_quality_context; //quality x context

	TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TPMD* PmdObject, TLog* Logfile);
	~TRecalibrationBQSR(){
		for(int i=0; i<numReadGroups; ++i){
			delete[] BQSR_cells_readGroup_quality[i];
		}
		delete[] BQSR_cells_readGroup_quality;

		if(considerPosition){
			for(int i=0; i<numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_position[i];
			}
			delete[] BQSR_cells_readGroup_position;
		}

		if(considerContext){
			for(int i=0; i<qualityIndex->numQ; ++i){
				delete[] BQSR_cells_quality_context[i];
			}
			delete[] BQSR_cells_quality_context;
		}
	};

	void addSite(TSite & site);
	bool estimateEpsilon();
	void writeToFile(std::string filenameTag);
};



#endif /* TRecalIBRATION_H_ */
