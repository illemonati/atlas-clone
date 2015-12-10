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
#include <omp.h>

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
//TQualityTransformTable
//---------------------------------------------------------------
class TQualityTransformTable{
public:
	int maxQ;
	double** table; //old qual / new qual

	TQualityTransformTable(int MaxQ){
		maxQ = MaxQ;
		table = new double*[maxQ];
		for(int i=0; i<maxQ; ++i){
			table[i] = new double[maxQ];
			for(int j=0; j<maxQ; ++j){
				table[i][j] = 0;
			}
		}
	};

	~TQualityTransformTable(){
		for(int i=0; i<maxQ; ++i){
			delete[] table[i];
		}
		delete[] table;
	};

	void add(int oldQ, int newQ){
		if(oldQ < maxQ && newQ < maxQ){
			table[oldQ][newQ] += 1.0;
		}
	};

	double size(){
		double size = 0;
		for(int i=0; i<maxQ; ++i){
			for(int j=0; j<maxQ; ++j){
				size += table[i][j];
			}
		}
		return size;
	};

	void printTable(std::ofstream & out){
		//print header
		out << "oldQ/newQ";
		for(int i=0; i<maxQ; ++i) out << "\t" << i;
		out << "\n";

		//get total
		double sum = size();

		//print rows
		for(int i=0; i<maxQ; ++i){
			out << i;
			for(int j=0; j<maxQ; ++j){
				out << "\t" << table[i][j] / sum;
			}
			out << "\n";
		}
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

	double dePhred(double quality){
		return pow(10.0, quality / -10.0);
	};

	virtual double getErrorRate(TBase* base){
		return dePhred(base->quality);
	};

	virtual int getQuality(TBase* base){
		return base->quality;
	};

	char getQualityAsChar(TBase* base){
		return getQuality(base) + 33;
	};

	void calcEmissionProbabilities(TSite & site){
		//first calculate for each base
		for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
			(*it)->fillEmissionProbabilitiesCore(getErrorRate(*it));
		}
		//then for the site
		site.calcEmissionProbabilities();
	};
	void addSiteToQualityTransformTable(TSite & site, TQualityTransformTable & QT){
		for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
			QT.add((*it)->quality, getQuality(*it));
		}
	};

	virtual bool requiresEstimation(){ return false;};
};

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
class TRecalibrationEMSite{
public:
	double** q; //covariates such as quality, position etc.
	double** D; //D for the emission probabilities: depends on genotype and base!
	double** B; //B = 4/3 D - 1
	int* context;
	int* readGroup;
	int* readGroupShifts;
	double* epsilon;
	double* P_g_given_d_oldBeta;
	int numReads;
	bool initialized;

	TRecalibrationEMSite();
	TRecalibrationEMSite(TSite & site);
	double dePhred(double quality){
		return pow(10.0, quality / -10.0);
	};
	~TRecalibrationEMSite();
	void calcEpsilon(double** params);
	double fill_P_g_given_d_beta_AND_calcLL(double** oldParams, double* freqs);
	double calcLL(double** oldParams, double* freqs);
	double calcQ(double** newParams);
	void addToJacobianAndF(arma::mat & Jacobian, arma::vec & F, double** params);
};

class TRecalibrationEMWindow{
public:
	std::vector<TRecalibrationEMSite*> sites;
	double* freqs; //base frequencies

	TRecalibrationEMWindow(TBaseFrequencies* baseFreqs);
	~TRecalibrationEMWindow(){
		delete[] freqs;
		for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
			delete *site;
		}
		sites.clear();
	};
	void addSite(TSite & site);
	double fill_P_g_given_d_beta_AND_calcLL(double** oldParams);
	double calcLL(double** oldParams);
	double calcQ(double** newParams);
	void addToJacobianAndF(arma::mat & Jacobian, arma::vec & F, double** params);
};

class TRecalibrationEM:public TRecalibration{
public:
	TLog* logfile;
	BamTools::SamHeader* bamHeader;
	int numReadGroups;
	std::string* readGroupNames;
	int numParams;
	int totNumParams;
	double** params;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;

	//variables for EM
	bool estimatetionRequired;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	double** newParams; //used during EM
	double** tmpParams; //used during NR
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	long numSitesAdded;

	TRecalibrationEM(BamTools::SamHeader* BamHeader, TParameters & params, TLog* Logfile);
	~TRecalibrationEM(){
		for(int i=0; i<numReadGroups; ++i){
			delete[] params[i];
			delete[] newParams[i];
			delete[] tmpParams[i];
		}
		delete[] params;
		delete[] newParams;
		delete[] tmpParams;
		delete[] readGroupNames;
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
	};
	bool requiresEstimation(){ return estimatetionRequired;};
	void addNewWindow(TBaseFrequencies* freqs);
	void addSite(TSite & site);
	double getErrorRate(TBase* base, double** theseParams);
	double getErrorRate(TBase* base);
	void runNewtonRaphson(double** theseParams, int & maxNewtonraphsonIteratios, double & maxFThreshold, TLog* logfile, std::string debugFilename);
	void runEM(std::string outputName);
	void writeCurrentEstimates(std::string filename, double & LL);
	void writeHeader(std::ofstream & out);
	void writeParams(std::ofstream & out, double & LL);
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
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

	double firstDerivativeSave, secondDerivativeSave;

	//TPMD* pmdObject;
	double numObservations;
	double numObservationsTmp;

	double F;
	double LL;

	TBQSR_cell();
	virtual ~TBQSR_cell(){};
	void empty();
	void reopenEstimation();
	void init(double initialError);
	void set(double error){curEstimate = error;};
	void set(double error, std::string & NumObservations);
	double getD(TBase* base, Base & RefBase);
	virtual void addBase(TBase* base, Base & RefBase);
	void runNewtonRaphson(double & convergenceThreshold);
	virtual bool estimate(double & convergenceThreshold, long & minObservations);
	std::string getNumObsForPrinting();
};

class TBQSR_cellPosition:public TBQSR_cell{
public:
	TBQSR_cell** BQSR_cells_readGroup_quality; //read group x quality
	TQualityIndex* qualityIndex;

	TBQSR_cellPosition();
	~TBQSR_cellPosition(){};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex);
	void addBase(TBase* base, Base & RefBase);
	void addToDerivatives(TBase* base, Base & RefBase, double & epsilon);
	bool estimate(double & convergenceThreshold, long & minObservations);
};

class TBQSR_cellPositionRev:public TBQSR_cellPosition{
public:
	TBQSR_cellPosition ** BQSR_cells_readGroup_position; //read group x position
	bool considerPosition;

	TBQSR_cellPositionRev();
	~TBQSR_cellPositionRev(){};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup, TQualityIndex* QualityIndex);
	void init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex);
	void addBase(TBase* base, Base & RefBase);
};

class TBQSR_cellContext:public TBQSR_cellPositionRev{
public:
	TBQSR_cellPositionRev** BQSR_cells_readGroup_position_rev; //read group x quality
	bool considerPositionRev;

	TBQSR_cellContext();
	~TBQSR_cellContext(){};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_quality_position, TBQSR_cellPositionRev** gotBQSR_cells_quality_position_rev, TQualityIndex* QualityIndex);
	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_quality_position, TQualityIndex* QualityIndex);
	void init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex);
	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex);
	void addBase(TBase* base, Base & RefBase);
};

//TODO: make class TBQSR_table!

class TRecalibrationBQSR:public TRecalibration{
private:
	TQualityIndex* qualityIndex;
	BamTools::SamHeader* bamHeader;
	TLog* logfile;
	TGenotypeMap genoMap;
	int numReadGroups;
	bool estimatetionRequired;
	double convergenceThreshold;
	bool estimationConverged;
	int maxPos;
	int numContexts;
	long minObservations;

	//recal tables
	bool qualityConverged;
	TBQSR_cell** BQSR_cells_readGroup_quality; //read group x quality
	bool considerPosition, positionConverged;
	TBQSR_cellPosition** BQSR_cells_readGroup_position; //read group x position
	bool considerPositionReverse, positionReverseConverged;
	TBQSR_cellPositionRev** BQSR_cells_readGroup_position_reverse; //read group x position
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
	void initializeBQSRReadGroupQualityTable(TParameters & params);
	void initializeBQSRReadGroupQualityTableFromFile(TParameters & params);
	void initializeBQSRReadGroupPositionTable(TParameters & params);
	void initializeBQSRReadGroupPositionTableFromFile(TParameters & params);

	void initializeBQSRReadGroupPositionReverseTable(TParameters & params);
	void initializeBQSRReadGroupPositionReverseTableFromFile(TParameters & params);

	void initializeBQSRReadGroupContextTable(TParameters & params);
	void initializeBQSRReadGroupContextTableFromFile(TParameters & params);

	bool requiresEstimation(){ return estimatetionRequired;};

public:
	TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TLog* Logfile);
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

		if(considerPositionReverse){
			for(int i=0; i<numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_position_reverse[i];
			}
			delete[] BQSR_cells_readGroup_position_reverse;
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
	void reopenEstimation();
	double getErrorRate(TBase* base);
	int getQuality(TBase* base);
};



#endif /* TRecalIBRATION_H_ */
