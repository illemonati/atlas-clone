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

//---------------------------------------------------------------
//TQualityIndex
//---------------------------------------------------------------
class TQualityIndex{
public:
	int minQ, maxQ, numQ, last, first;
	int* index;

	TQualityIndex(int MinQ, int MaxQ){
		minQ= MinQ;
		maxQ = MaxQ;
		numQ = maxQ - minQ + 1;
		last = numQ - 1;
		first = 0;

		//fill index
		index = new int[maxQ + 1];
		for(int i=0; i < maxQ + 1; ++i){
			if(i < minQ) index[i] = 0;
			else index[i] = i - minQ;
		}
	};

	~TQualityIndex(){
		delete[] index;
	};

	int& getIndex(int & quality){
		if(quality < 0) throw "Quality is negative!";
		if(quality > maxQ) return last;
		return index[quality];
	};

	int getQuality(int & index){
		if(index < 0) throw "Quality index is negative!";
		if(index > numQ) return maxQ;
		return minQ + index;
	};
};

//---------------------------------------------------------------
//TRecalibration: default = no recalibration
//---------------------------------------------------------------
class TRecalibration{
public:
	TRecalibration(){};
	virtual ~TRecalibration(){};

	int makePhredInt(double & epsilon){
		return round(-10.0 * log10(epsilon));
	};

	double makePhred(double & epsilon){
		if(epsilon < 0.0000000001) return 100.0;
		return -10.0 * log10(epsilon);
	};

	double dePhred(int quality){
		return pow(10.0, (double) quality / -10.0);
	};
	virtual double getErrorRate(TBase* base){
		return dePhred(base->quality);
	};
	void calcEmissionProbabilities(TSite & site, TPMD & pmdObject){
		for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
			(*it)->fillEmissionProbabilitiesCore(pmdObject, getErrorRate(*it));
		}
		site.calcEmissionProbabilities();
	};
};

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
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

class TBQSR_cell{
private:
	double numMatches;
	void addToDerivatives(TBase* base, Base & RefBase);

public:
	double curEstimate;
	bool estimationConverged;
	double firstDerivative, secondDerivative;
	TPMD* pmdObject;
	double numObservations;

	double F;
	double LL;

	TBQSR_cell();
	virtual ~TBQSR_cell(){};
	void empty();
	void init(double initialError, TPMD* PmdObject);
	void set(double error){curEstimate = error;};
	void set(double error, long NumObservations){curEstimate = error; numObservations=NumObservations;};
	double getD(TBase* base, Base & RefBase);
	virtual void addBase(TBase* base, Base & RefBase);
	void runNewtonRalphson(double & convergenceThreshold);
	virtual bool estimate(double & convergenceThreshold);
};

class TBQSR_cellPosition:public TBQSR_cell{
public:
	TBQSR_cell** BQSR_cells_readGroup_quality; //read group x quality
	TQualityIndex* qualityIndex;

	TBQSR_cellPosition();
	~TBQSR_cellPosition(){};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, TPMD* PmdObject);
	void addBase(TBase* base, Base & RefBase);
	void addToDerivatives(TBase* base, Base & RefBase, double & epsilon);
	bool estimate(double & convergenceThreshold);
};


class TBQSR_cellContext:public TBQSR_cellPosition{
public:
	TBQSR_cellPosition** BQSR_quality_position; //read group x quality
	bool considerPosition;

	TBQSR_cellContext();
	~TBQSR_cellContext(){};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_quality_position, TQualityIndex* QualityIndex, TPMD* PmdObject);
	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, TPMD* PmdObject);
	void addBase(TBase* base, Base & RefBase);
};


class TRecalibrationBQSR:public TRecalibration{
private:
	TQualityIndex* qualityIndex;
	BamTools::SamHeader* bamHeader;
	TLog* logfile;
	TPMD* pmdObject;
	TGenotypeMap genoMap;
	int numReadGroups;
	bool estimatetionRequired;
	double convergenceThreshold;
	bool estimationConverged;
	int maxPos;
	int numContexts;

	//recal tables
	bool qualityConverged;
	TBQSR_cell** BQSR_cells_readGroup_quality; //read group x quality
	bool considerPosition, positionConverged;
	TBQSR_cellPosition** BQSR_cells_readGroup_position; //read group x position
	bool considerContext, contextConverged;
	TBQSR_cellContext** BQSR_cells_readGroup_context; //quality x context

	//-------------------------
	/*
	TBQSR_cellPosition* LLSurface;
	int numLLSurfacePoints;
	bool surfaceCalculated;
	*/
	//-------------------------

	int findReadGroupIndex(std::string & name);
	void initializeBQSRReadGroupQualityTable(std::string filename);
	void initializeBQSRReadGroupQualityTable(TParameters & params);
	void initializeBQSRReadGroupPositionTable(std::string filename);
	void initializeBQSRReadGroupPositionTable(TParameters & params);
	void initializeBQSRReadGroupContextTable(std::string filename);
	void initializeBQSRReadGroupContextTable(TParameters & params);

public:
	TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TPMD* PmdObject, TLog* Logfile);
	~TRecalibrationBQSR(){
		for(int i=0; i<numReadGroups; ++i){
			delete[] BQSR_cells_readGroup_quality[i];
		}
		delete[] BQSR_cells_readGroup_quality;
		delete qualityIndex;

		if(considerPosition){
			for(int i=0; i<numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_position[i];
			}
			delete[] BQSR_cells_readGroup_position;
		}

		if(considerContext){
			for(int i=0; i<numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_context[i];
			}
			delete[] BQSR_cells_readGroup_context;
		}
	};

	void addSite(TSite & site);
	bool estimateEpsilon();
	void writeToFile(std::string filenameTag);
	bool allConverged();
	double getErrorRate(TBase* base);
};



#endif /* TRecalIBRATION_H_ */
