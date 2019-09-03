#include "TRecalibration.h"


//---------------------------------------------------------------
//TQualityIndex
//---------------------------------------------------------------
class TQualityIndex{
	//Note: quality as stored in bases ranges from 33 to max!
	//error is error rate between 0 and 1
	//phred is phred-scaled error as phred = -10 * log10(error)
	//phredInt is (int) phred
	//quality is phredInt + 33
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

	int& getIndexFromQuality(const int & quality){
		if(quality < 33){
			throw "Quality is negative!";
		}
		if(quality > maxPhredInt + 33){
			return last;
		}
		return index[quality - 33];
	};

	int& getIndexFromPhredInt(const int & phredInt){
		if(phredInt < 0){
			throw "Phred int is negative!";
		}
		if(phredInt > maxPhredInt){
			return last;
		}
		return index[phredInt];
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

//TODO: make class TBQSR_table!

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
	virtual void init(int ReadGroup);
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


class TBQSR_cellQuality:public TBQSR_cell_base{
private:
	//for storage
	std::vector<float*> D_storage;
	std::vector<float*>::reverse_iterator batchIt;
	float* pointerToBatch;

public:
	uint64_t numMatches;

	TBQSR_cellQuality();
	virtual ~TBQSR_cellQuality(){
		clearStorage();
	};
	void empty();
	void clearStorage();
	void init(int ReadGroup, float initialError);
	virtual void doStoreData();
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
protected:
	std::vector<BQSRFactorStorage*> D_storage;
	std::vector<BQSRFactorStorage*>::reverse_iterator batchIt;
	BQSRFactorStorage* pointerToBatch;

public:
	TBQSR_cellQuality** BQSR_cells_readGroup_quality; //read group x quality
	TQualityIndex* qualityIndex;

	TBQSR_cellPosition();
	~TBQSR_cellPosition(){
		clearStorage();
	};

	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup);
	virtual void doStoreData();
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

	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup);

	float getEpsilon(TBase* base, TQualityMap & qualiMap);
};

class TBQSR_cellContext:public TBQSR_cellPositionRev{
public:
	TBQSR_cellPositionRev** BQSR_cells_readGroup_position_rev; //read group x quality
	bool considerPositionRev;

	TBQSR_cellContext();
	~TBQSR_cellContext(){};

	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup);
	float getEpsilon(TBase* base, TQualityMap & qualiMap);
};

//-----------------------------------------------------------------
// TRecalibrationBQSRStorage
// Object to store all BQSR tables
//-----------------------------------------------------------------
class TRecalibrationBQSRStorage{
public:
	int numReadGroups;
	int numQuality;
	int maxPos, maxPosReverse;
	int numContexts;

	TGenotypeMap _genoMap;

	//recal tables
	bool considerQuality;
	TBQSR_cellQuality** qualityCells; //read group x quality
	bool considerPosition;
	TBQSR_cellPosition** positionCells; //read group x position
	bool considerPositionReverse;
	TBQSR_cellPositionRev** positionReverseCells; //read group x position
	bool considerContext;
	TBQSR_cellContext** contextCells; //read group x context

	TRecalibrationBQSRStorage();
	~TRecalibrationBQSRStorage();

	void initializeQualityCells(int NumReadGroups, int NumQuality, TQualityMap & qualityMap, TQualityIndex* qualityIndex, const bool & _storeDataInMemory);
	void initializePositionCells(int NumReadGroups, int MaxPos, TQualityIndex* qualityIndex, const bool & _storeDataInMemory);
	void initializePositionReverseCells(int NumReadGroups, int MaxPos, TQualityIndex* qualityIndex, const bool & _storeDataInMemory);
	void initializeContextCells(int NumReadGroups, TQualityIndex* qualityIndex, const bool & _storeDataInMemory);

	int numQualityCells(){ if(considerQuality){ return numReadGroups * numQuality; } else { return 0; } };
	int numPositionCells(){ if(considerPosition){ return numReadGroups * maxPos; } else { return 0; } };
	int numPositionReverseCells(){ if(considerPositionReverse){ return numReadGroups * maxPosReverse; } else { return 0; } };
	int numContextCells(){ if(considerContext){ return numReadGroups * numContexts; } else { return 0; } };

	//HACK! Just need this function in multiple places...
	double getMaxValueFromFile(std::string filename, int col, int numCol);
};

//-----------------------------------------------------------------
// TRecalibrationBQSR
//-----------------------------------------------------------------
class TRecalibrationBQSR:public TRecalibration{
private:
	TLog* _logfile;
	TQualityIndex* qualityIndex;
	TReadGroups* _readGroups;
	TRecalibrationBQSRStorage storage;

	void _initialize(TLog* Logfile, TReadGroups* ReadGroups);

	void _initializeBQSRReadGroupQualityTableFromFile(std::string filename);
	void _initializeBQSRReadGroupPositionTableFromFile(std::string filename);
	void _initializeBQSRReadGroupPositionReverseTableFromFile(std::string filename);
	void _initializeBQSRReadGroupContextTableFromFile(std::string filename);

public:
	TRecalibrationBQSR(TLog* Logfile, TReadGroups* ReadGroups);
	TRecalibrationBQSR(std::string qualityFile, TLog* Logfile, TReadGroups* ReadGroups);
	~TRecalibrationBQSR(){
		if(storage.considerQuality)
			delete qualityIndex;
	};

	void addPositionEffect(std::string filename);
	void addPositionReverseEffect(std::string filename);
	void addContextEffect(std::string filename);

	bool recalibrationChangesQualities(){ return true; };
	double getErrorRate(TBase & base);
};

//-----------------------------------------------------------------
// TRecalibrationBQSREstimator
//-----------------------------------------------------------------
class TRecalibrationBQSREstimator{
private:
	TLog* _logfile;
	TQualityMap _qualityMap;
	TQualityIndex* qualityIndex;
	TReadGroups* _readGroups;
	TReadGroupMap* _readGroupMap;
	TRecalibrationBQSRStorage storage;

	//estimation
	bool _estimatetionRequired;
	float _convergenceThreshold_F;
	float _minEpsilonQuality, _minEpsilonFactors;
	int _numLoopIncreaseFAllowed;
	int _curNewtonRaphsonLoop;
	bool _estimationConverged;

	long _minObservations;
	bool _storeDataInMemory;
	bool _dataStored;

	//LL surfaces
	bool _printLLSurface;
	bool _LLSurfacePrinted;
	int _numPosLLsurface;

	//recal tables
	bool _qualityConverged, _estimateQuality;
	bool _positionConverged, _estimatePosition;
	bool _positionReverseConverged, _estimatePositionReverse;
	bool _contextConverged, _estimateContext;

	void _initializeBQSRReadGroupQualityTable(TParameters & params);
	void _fillBQSRReadGroupQualityTableFromFile(std::string filename);

	void _initializeBQSRReadGroupPositionTable(TParameters & params);
	void _fillBQSRReadGroupPositionTableFromFile(std::string filename);

	void _initializeBQSRReadGroupPositionReverseTable(TParameters & params);
	void _fillBQSRReadGroupPositionReverseTableFromFile(std::string filename);

	void _initializeBQSRReadGroupContextTable(TParameters & params);
	void _fillBQSRReadGroupContextTableFromFile(std::string filename);

	bool _requiresEstimation(){ return _estimatetionRequired; };

public:
	TRecalibrationBQSREstimator(TParameters & params, TLog* Logfile, TReadGroups* ReadGroups, TReadGroupMap* ReadGroupMap);
	~TRecalibrationBQSREstimator(){
		if(storage.considerQuality)
			delete qualityIndex;
	};

	bool dataHasBeenStored(){ return _dataStored; };
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
};


