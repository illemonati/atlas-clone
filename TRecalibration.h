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
#include "TReadGroups.h"
#include "bamtools/api/SamHeader.h"

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

	int& getIndex(const int & quality){
		if(quality < 0) throw "Quality is negative!";
		if(quality > maxQ) return last;
		return index[quality];
	};

	int getQuality(const int & index){
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
protected:
	bool mergedInd;
	int* readGroupMap;
	int numReadGroups;
	int origNumReadGroups;
	bool readGroupMapInitialized;

	TQualityMap qualityMap;

public:
	TRecalibration();

	virtual ~TRecalibration(){
		if(readGroupMapInitialized) delete[] readGroupMap;
	};

	virtual bool recalibrationChangesQualities(){
		return false;
	};

	void initializeReadGroupMap(BamTools::SamHeader* bamHeader, TParameters & params, TLog* logfile);

	int makePhredInt(double & epsilon){
		return round(-10.0 * log10(epsilon));
	};

	double makePhred(double & epsilon){
		if(epsilon < 0.0000000001) return 100.0;
		return -10.0 * log10(epsilon);
	};

	double makePhred(float & epsilon){
		if(epsilon < 0.0000000001) return 100.0;
		return -10.0 * log10(epsilon);
	};

	double dePhred(double quality){
		return pow(10.0, quality / -10.0);
	};

	char getQualityAsChar(const TBase & base, int & minOutQuality, int & maxOutQuality){
		int qual = getQuality(base) + 33;
		if(qual > maxOutQuality) qual = maxOutQuality;
		if(qual < minOutQuality) qual = minOutQuality;
		return qual;
	};

	void calcEmissionProbabilities(TSite & site);
	virtual double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	double getErrorRateFromBase(const TBase & base);
	virtual int getQuality(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context){
		return quality;
	}
	virtual int getQuality(const TBase & base){
		return base.quality;
	};

	virtual bool requiresEstimation(){ return false;};
	int findReadGroupIndex(std::string & name, BamTools::SamReadGroupDictionary & readGroups);
};

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
class TRecalibrationEMModel{
protected:
	int numReadGroups;
	int totNumParams;
	int* readGroupShifts;

	double** betas; //betas of the model
	double** oldBetas; //use during estimation
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	bool initialized;
	bool EMParamsInitialized;

	double tmp;
	int tmpIndex;

	void initialize(int NumReadGroups);

public:
	int numParams;
	long numSitesAdded;

	TRecalibrationEMModel();
	TRecalibrationEMModel(int NumReadGroups);
	virtual ~TRecalibrationEMModel(){
		delete[] readGroupShifts;
		for(int r=0; r<numReadGroups; ++r){
			delete[] betas[r];
			delete[] oldBetas[r];
		}
		delete[] betas;
		delete[] oldBetas;
	};

	bool setParams(std::vector<std::string> & vec, int & rg);
	void initializeEMParams();
	void setEMParamsToZero();
	virtual double calcEpsilon(const short & readGroup, float* & q, const short & context);
	virtual void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, float** & q, short* & readGroup, short* & context);
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	void rejectProposedParameters();
	double getSteepestGradient();
	virtual void writeParametersToFile(std::ofstream & out, const int & readGroup);
	void printJacobianToStdOut();
	virtual double getErrorRate(int rg, double originalErrorRate, const int & posInRead, const int & context);
};

class TRecalibrationEMModelNoContext:public TRecalibrationEMModel{
public:
	TRecalibrationEMModelNoContext(int NumReadGroups);

	double calcEpsilon(const short & readGroup, float* & q, const short & context);
	void addToFandJacobian(const int & numReads, double* & weights, double* & weightsJacobian, const float & P_g_given_d_oldBeta, float** & q, short* & readGroup, short* & context);
	void writeParametersToFile(std::ofstream & out, const int & readGroup);
	double getErrorRate(int rg, double originalErrorRate, const int & posInRead, const int & context);
};

class TRecalibrationEMSite{
public:
	float** q; //covariates such as quality, position etc.
	float** D; //D for the emission probabilities: depends on genotype and base!
	float** B; //B = 4/3 D - 1
	short* context;
	short* readGroup;
	short* readGroupShifts;
	float* epsilon;
	float* P_g_given_d_oldBeta;
	int numReads;
	bool initialized;

	TRecalibrationEMSite();
	TRecalibrationEMSite(TSite & site, int* readGroupMap);
	double dePhred(double quality){
		double tmp = pow(10.0, quality / -10.0);
		if(tmp < 0.0000000001) return 0.0000000001;
		if(tmp > 0.9999999999) return 0.9999999999;
		return tmp;
	};
	virtual ~TRecalibrationEMSite();
	virtual void calcEpsilon(TRecalibrationEMModel* & model);
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel* & model, double* freqs);
	double calcLL(TRecalibrationEMModel* & model, double* freqs);
	double calcQ(TRecalibrationEMModel* & model);
	virtual void addToJacobianAndF(TRecalibrationEMModel* & model);
};

class TRecalibrationEMWindow{
public:
	std::vector<TRecalibrationEMSite*> sites;
	double* freqs; //base frequencies
	int* readGroupMap;

	TRecalibrationEMWindow(TBaseFrequencies* baseFreqs, int* ReadGroupMap);
	virtual ~TRecalibrationEMWindow(){
		delete[] freqs;
		for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
			delete *site;
		}
		sites.clear();
	};
	virtual void addSite(TSite & site);
	double fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel* & model);
	double calcLL(TRecalibrationEMModel* & model);
	double calcQ(TRecalibrationEMModel* & model);
	void addToJacobianAndF(TRecalibrationEMModel* & model);
	void setEuqalBaseFrequencies();
};

class TRecalibrationEM:public TRecalibration{
public:
	TLog* logfile;
	BamTools::SamHeader* bamHeader;
	std::string* readGroupNames;
	TRecalibrationEMModel* model;
	std::vector<TRecalibrationEMWindow*> windows;
	std::vector<TRecalibrationEMWindow*>::iterator curWindow;

	//variables for EM
	bool estimatetionRequired;
	bool equalBaseFrequencies;
	long numSitesAdded;
	int numEMIterations;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	int maxCoverage; //sites with higher coverage will be ignored

	TRecalibrationEM(BamTools::SamHeader* BamHeader, std::string &name, TParameters & params, TLog* Logfile);
	~TRecalibrationEM(){
		delete[] readGroupNames;
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			delete *curWindow;
		}
		windows.clear();
		delete model;
	};
	bool recalibrationChangesQualities(){ return true; };
	bool requiresEstimation(){ return estimatetionRequired;};
	void addNewWindow(TBaseFrequencies* freqs);
	void addSite(TSite & site);
	void runNewtonRaphson(int & maxNewtonraphsonIteratios, double & maxFThreshold, TLog* logfile, bool & writeTmpTables, std::string debugFilename);
	void runEM(std::string outputName, bool & writeTmpTables);
	void writeCurrentEstimates(std::string filename, double & LL);
	void writeHeader(std::ofstream & out);
	void writeParams(std::ofstream & out, double & LL);
	double calcLL();
	void calcLikelihoodSurface(std::string filename, int numMarginalGridPoints);
	void calcQSurface(std::string filename, int numMarginalGridPoints);
	double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	int getQuality(const TBase & base);
	int getQuality(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
};

//---------------------------------------------------------------
//RecalibrationBQSR
//---------------------------------------------------------------
//covariates to take into account:
// - read base (A, G, C, T)
// -
class TBQSR_cell_base{
public:
 //	TReadGroups newReadGroupObject;
	float curEstimate;
	bool estimationConverged;
	float firstDerivative, secondDerivative;
	float firstDerivativeSave, secondDerivativeSave;
	float numObservations;
	float numObservationsTmp;
	float F;
	float LL;
	int myReadGroup;
	//for storage
	bool store;
	int batchSize;
	int next;

	TBQSR_cell_base();
	virtual ~TBQSR_cell_base(){};
	virtual void empty();
	virtual void init(bool Store, int ReadGroup);
	void reopenEstimation();
	void set(float error){curEstimate = error;};
	void set(float error, std::string & NumObservations);
	float getD(TBase* base, Base & RefBase);
	virtual void addBase(TBase* base, Base & RefBase){throw "TBQSR_cell_base::addBase(TBase* base, Base & RefBase) not defined for base class!";};
	virtual void recalculateDerivativesFromDataInMemory(){throw "TBQSR_cell_base::recalculateDerivativesFromDataInMemory() not defined for base class!";};
	virtual void recalculateLLFromDataInMemory(){throw "TBQSR_cell_base::recalculateLLFromDataInMemory() not defined for base class!";};
	virtual bool estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations){throw "TBQSR_cell_base::estimate(double & convergenceThreshold, double & minEpsilon, long & minObservations) not defined for base class!";};
	void runNewtonRaphson(float & convergenceThreshold);
	std::string getNumObsForPrinting();
	void calcLikelihoodSurfaceAt(int numPositions, double* positions, std::string & tag, std::ofstream & out);
	virtual void calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out){throw "TBQSR_cell_base::calcLikelihoodSurface(int numPositions, std::ofstream & out) not defined for base class!";};;
};


class TBQSR_cell:public TBQSR_cell_base{
public:
	float numMatches;
	//for storage
	std::vector<float*> D_storage;
	std::vector<float*>::reverse_iterator batchIt;
	float* pointerToBatch;

	TBQSR_cell();
	virtual ~TBQSR_cell(){
		clearStorage();
	};
	void empty();
	void clearStorage();
	void init(float initialError, bool Store, int ReadGroup);
	void addBase(TBase* base, Base & RefBase);
	void addToDerivatives(float & D);
	void addToLL(float & D);
	void recalculateDerivativesFromDataInMemory();
	void recalculateLLFromDataInMemory();
	void runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon);
	bool estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations);
	float makePhred(float & epsilon){
		if(epsilon < 0.0000000001) return 100.0;
		return -10.0 * log10(epsilon);
	};

	void calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out);
};

struct BQSRFactorStorage{
	float D;
	float epsilon;
};

class TBQSR_cellPosition:public TBQSR_cell_base{
public:
	TBQSR_cell** BQSR_cells_readGroup_quality; //read group x quality
	TQualityIndex* qualityIndex;
	//for storage
	std::vector<BQSRFactorStorage*> D_storage;
	std::vector<BQSRFactorStorage*>::reverse_iterator batchIt;
	BQSRFactorStorage* pointerToBatch;

	TBQSR_cellPosition();
	~TBQSR_cellPosition(){
		clearStorage();
	};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup);
	void clearStorage();
	virtual float getEpsilon(TBase* base);
	void addBase(TBase* base, Base & RefBase);
	void addToDerivatives(float & D, float & epsilon);
	void addToLL(float & D, float & epsilon);
	void recalculateDerivativesFromDataInMemory();
	void recalculateLLFromDataInMemory();
	void runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon);
	bool estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations);
	void calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out);
};

class TBQSR_cellPositionRev:public TBQSR_cellPosition{
public:
	TBQSR_cellPosition ** BQSR_cells_readGroup_position; //read group x position
	bool considerPosition;

	TBQSR_cellPositionRev();
	~TBQSR_cellPositionRev(){};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup);
	void init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup);
	float getEpsilon(TBase* base);
};

class TBQSR_cellContext:public TBQSR_cellPositionRev{
public:
	TBQSR_cellPositionRev** BQSR_cells_readGroup_position_rev; //read group x quality
	bool considerPositionRev;

	TBQSR_cellContext();
	~TBQSR_cellContext(){};

	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_quality_position, TBQSR_cellPositionRev** gotBQSR_cells_quality_position_rev, TQualityIndex* QualityIndex, bool Store, int ReadGroup);
	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_quality_position, TQualityIndex* QualityIndex, bool Store, int ReadGroup);
	void init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex, bool Store, int ReadGroup);
	void init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup);
	float getEpsilon(TBase* base);
};

//TODO: make class TBQSR_table!

class TRecalibrationBQSR:public TRecalibration{
private:
	TQualityIndex* qualityIndex;
	BamTools::SamHeader* bamHeader;
	BamTools::SamHeader* newHeader;

	TLog* logfile;
	TGenotypeMap genoMap;
	bool estimatetionRequired;
	float convergenceThreshold_F;
	float minEpsilonQuality, minEpsilonFactors;
	bool estimationConverged;
	int maxPos;
	int numContexts;
	long minObservations;
	bool storeDataInMemory;
	bool dataStored;
	bool printLLSurface;
	bool LLSurfacePrinted;
	int numPosLLsurface;

	//recal tables
	bool qualityConverged, estimateQuality;
	TBQSR_cell** BQSR_cells_readGroup_quality; //read group x quality
	bool considerPosition, positionConverged, estimatePosition;
	TBQSR_cellPosition** BQSR_cells_readGroup_position; //read group x position
	bool considerPositionReverse, positionReverseConverged, estimatePositionReverse;
	TBQSR_cellPositionRev** BQSR_cells_readGroup_position_reverse; //read group x position
	bool considerContext, contextConverged, estimateContext;
	TBQSR_cellContext** BQSR_cells_readGroup_context; //quality x context

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

	bool recalibrationChangesQualities(){ return true; };
	bool dataHasBeenStored(){ return dataStored; };
	void addSite(TSite & site);
	void recalculateDerivativesFromDataInMemory();
	bool estimateEpsilon(std::string filenameTag);
	void writeAllToFile(std::string filenameTag);
	void writeCurrentTmpTable(std::string filenameTag);
	void writeQualityToFile(std::string & filenameTag);
	void writePositionToFile(std::string & filenameTag);
	void writePositionReverseToFile(std::string & filenameTag);
	void writeContextToFile(std::string & filenameTag);
	void calculateAndPrintLLSurfaceQuality(std::string & filenameTag);
	void calculateAndPrintLLSurfacePosition(std::string & filenameTag);
	void calculateAndPrintLLSurfaceReversePosition(std::string & filenameTag);
	void calculateAndPrintLLSurfaceContext(std::string & filenameTag);
	bool allConverged();
	void reopenEstimation();
	double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	int getQuality(const TBase & base);
	int getQuality(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
};



#endif /* TRecalIBRATION_H_ */
