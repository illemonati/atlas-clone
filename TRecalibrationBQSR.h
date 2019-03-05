#include "TRecalibration.h"


//---------------------------------------------------------------
//TQualityIndex
//---------------------------------------------------------------
class TQualityIndex{
	//Note: quality as stored in bases ranges from 33 to max!
public:
	int minPhredInt, maxPhredInt, numQ, last, first;
	int* index;

	TQualityIndex(int MinQ, int MaxQ){
		//here minQ and maxQ are not in ascii. by default minQ = 0 and maxQ = 100
		minPhredInt = MinQ;
		maxPhredInt = MaxQ;
		numQ = maxPhredInt - minPhredInt + 1;
		last = numQ - 1;
		first = 0;

		//fill index
		index = new int[numQ + 33];
		for(int i=0; i < maxPhredInt + 1; ++i){
			if(i < minPhredInt) index[i] = 0;
			else index[i] = i - minPhredInt;
		}
	};

	~TQualityIndex(){
		delete[] index;
	};

	int& getIndex(const int & quality){
		if(quality < 33) throw "Quality is negative!";
		if(quality > maxPhredInt + 33) return last;
		return index[quality - 33];
	};

	int getPhredIntFromIndex(const int & index){
		if(index < 0) throw "Quality index is negative!";
		if(index > numQ) return maxPhredInt;
		return minPhredInt + index;
	};
};

//---------------------------------------------------------------
//RecalibrationBQSR
//---------------------------------------------------------------
//covariates to take into account:
// - read base (A, G, C, T)
// -
class TBQSR_cell_base{
public:
	float curEstimate;
	bool estimationConverged, hasData;
	float firstDerivative, secondDerivative;
	float firstDerivativeSave, secondDerivativeSave;
	uint64_t numObservations;
	uint64_t numObservationsTmp;
	double F, oldF;
	double LL;
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
	virtual void addBase(TBase* base, Base & RefBase, TQualityMap & qualiMap){throw "TBQSR_cell_base::addBase(TBase* base, Base & RefBase) not defined for base class!";};
	virtual void recalculateDerivativesFromDataInMemory(){throw "TBQSR_cell_base::recalculateDerivativesFromDataInMemory() not defined for base class!";};
	virtual void recalculateLLFromDataInMemory(){throw "TBQSR_cell_base::recalculateLLFromDataInMemory() not defined for base class!";};
	virtual bool estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations){throw "TBQSR_cell_base::estimate(double & convergenceThreshold, double & minEpsilon, long & minObservations) not defined for base class!";};
	void runNewtonRaphson(float & convergenceThreshold, bool & allowIncreaseInF);
	std::string getNumObsForPrinting();
	void calcLikelihoodSurfaceAt(int numPositions, double* positions, std::string & tag, std::ofstream & out);
	virtual void calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out){throw "TBQSR_cell_base::calcLikelihoodSurface(int numPositions, std::ofstream & out) not defined for base class!";};;
};


class TBQSR_cell:public TBQSR_cell_base{
public:
	uint64_t numMatches;
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
	void addBase(TBase* base, Base & RefBase, TQualityMap & qualiMap);
	void addToDerivatives(float & D);
	void addToLL(float & D);
	void recalculateDerivativesFromDataInMemory();
	void recalculateLLFromDataInMemory();
	void runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon, bool & allowIncreaseInF);
	bool estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations, bool & allowIncreaseInF);
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
	virtual float getEpsilon(TBase* base, TQualityMap & qualiMap);
	void addBase(TBase* base, Base & RefBase, TQualityMap & qualiMap);
	void addToDerivatives(float & D, float & epsilon);
	void addToLL(float & D, float & epsilon);
	void recalculateDerivativesFromDataInMemory();
	void recalculateLLFromDataInMemory();
	void runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon, bool & allowIncreaseInF);
	bool estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations, bool & allowIncreaseInF);
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
	float getEpsilon(TBase* base, TQualityMap & qualiMap);
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
	float getEpsilon(TBase* base, TQualityMap & qualiMap);
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
	int numLoopIncreaseFAllowed;
	int curNewtonRaphsonLoop;
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
	TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TLog* Logfile, TReadGroupMap& ReadGroupMap);
	~TRecalibrationBQSR(){
		for(int i=0; i<readGroupMapObject.numReadGroups; ++i){
			delete[] BQSR_cells_readGroup_quality[i];
		}
		delete[] BQSR_cells_readGroup_quality;
		delete qualityIndex;

		if(considerPosition){
			for(int i=0; i<readGroupMapObject.numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_position[i];
			}
			delete[] BQSR_cells_readGroup_position;
		}

		if(considerPositionReverse){
			for(int i=0; i<readGroupMapObject.numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_position_reverse[i];
			}
			delete[] BQSR_cells_readGroup_position_reverse;
		}

		if(considerContext){
			for(int i=0; i<readGroupMapObject.numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_context[i];
			}
			delete[] BQSR_cells_readGroup_context;
		}
	};

	bool recalibrationChangesQualities(){ return true; };
	bool dataHasBeenStored(){ return dataStored; };
	void addSite(TSite & site, TQualityMap & qualiMap);
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
//	double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	double getErrorRate(TBase & base);
//	int getQuality(TBase & base);
};

