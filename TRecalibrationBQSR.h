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


class TBQSR_cell:public TBQSR_cell_base{
private:
	//for storage
	std::vector<float*> D_storage;
	std::vector<float*>::reverse_iterator batchIt;
	float* pointerToBatch;

public:
	uint64_t numMatches;

	TBQSR_cell();
	virtual ~TBQSR_cell(){
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
	TBQSR_cell** BQSR_cells_readGroup_quality; //read group x quality
	TQualityIndex* qualityIndex;

	TBQSR_cellPosition();
	~TBQSR_cellPosition(){
		clearStorage();
	};

	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup);
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

	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup);

	float getEpsilon(TBase* base, TQualityMap & qualiMap);
};

class TBQSR_cellContext:public TBQSR_cellPositionRev{
public:
	TBQSR_cellPositionRev** BQSR_cells_readGroup_position_rev; //read group x quality
	bool considerPositionRev;

	TBQSR_cellContext();
	~TBQSR_cellContext(){};

	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup);
	void init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup);
	float getEpsilon(TBase* base, TQualityMap & qualiMap);
};


//-----------------------------------------------------------------

//-----------------------------------------------------------------
class TRecalibrationBQSR:public TRecalibration{
private:
	TQualityIndex* qualityIndex;
	TReadGroups* _readGroups;

	TLog* _logfile;
	TGenotypeMap _genoMap;

	int _maxPos;
	int _numContexts;

	//recal tables
	TBQSR_cell** _BQSR_cells_readGroup_quality; //read group x quality
	bool _considerPosition;
	TBQSR_cellPosition** BQSR_cells_readGroup_position; //read group x position
	bool _considerPositionReverse;
	TBQSR_cellPositionRev** _BQSR_cells_readGroup_position_reverse; //read group x position
	bool _considerContext;
	TBQSR_cellContext** _BQSR_cells_readGroup_context; //quality x context

	void _initializeBQSRReadGroupQualityTableFromFile(std::string filename);
	void _initializeBQSRReadGroupPositionTableFromFile(std::string filename);
	void _initializeBQSRReadGroupPositionReverseTableFromFile(std::string filename);
	void _initializeBQSRReadGroupContextTableFromFile(std::string filename);

public:
	TRecalibrationBQSR(TParameters & params, TLog* Logfile, TReadGroups* ReadGroups);
	~TRecalibrationBQSR(){
		for(int i=0; i<_readGroups->size(); ++i){
			delete[] _BQSR_cells_readGroup_quality[i];
		}
		delete[] _BQSR_cells_readGroup_quality;
		delete qualityIndex;

		if(_considerPosition){
			for(int i=0; i<_readGroups->size(); ++i){
				delete[] BQSR_cells_readGroup_position[i];
			}
			delete[] BQSR_cells_readGroup_position;
		}

		if(_considerPositionReverse){
			for(int i=0; i<_readGroups->size(); ++i){
				delete[] _BQSR_cells_readGroup_position_reverse[i];
			}
			delete[] _BQSR_cells_readGroup_position_reverse;
		}

		if(_considerContext){
			for(int i=0; i<_readGroups->size(); ++i){
				delete[] _BQSR_cells_readGroup_context[i];
			}
			delete[] _BQSR_cells_readGroup_context;
		}
	};

	bool recalibrationChangesQualities(){ return true; };
	double getErrorRate(TBase & base);
};

/*
class TRecalibrationBQSREstimator:public TRecalibration{
private:
	TQualityIndex* qualityIndex;
	TReadGroups* _readGroups;
	TReadGroupMap* _readGroupMap;

	TLog* _logfile;
	TGenotypeMap _genoMap;
	bool _estimatetionRequired;
	float _convergenceThreshold_F;
	float _minEpsilonQuality, _minEpsilonFactors;
	int _numLoopIncreaseFAllowed;
	int _curNewtonRaphsonLoop;
	bool _estimationConverged;
	int _maxPos;
	int _numContexts;
	long _minObservations;
	bool _storeDataInMemory;
	bool _dataStored;
	bool _printLLSurface;
	bool _LLSurfacePrinted;
	int _numPosLLsurface;

	//recal tables
	bool _qualityConverged, _estimateQuality;
	TBQSR_cell** _BQSR_cells_readGroup_quality; //read group x quality
	bool _considerPosition, _positionConverged, _estimatePosition;
	TBQSR_cellPosition** BQSR_cells_readGroup_position; //read group x position
	bool _considerPositionReverse, _positionReverseConverged, _estimatePositionReverse;
	TBQSR_cellPositionRev** _BQSR_cells_readGroup_position_reverse; //read group x position
	bool _considerContext, _contextConverged, _estimateContext;
	TBQSR_cellContext** _BQSR_cells_readGroup_context; //quality x context

	void _initializeBQSRReadGroupQualityTable(TParameters & params);
	void _initializeBQSRReadGroupQualityTableFromFile(TParameters & params);
	void _initializeBQSRReadGroupPositionTable(TParameters & params);
	void _initializeBQSRReadGroupPositionTableFromFile(TParameters & params);

	void _initializeBQSRReadGroupPositionReverseTable(TParameters & params);
	void _initializeBQSRReadGroupPositionReverseTableFromFile(TParameters & params);

	void _initializeBQSRReadGroupContextTable(TParameters & params);
	void _initializeBQSRReadGroupContextTableFromFile(TParameters & params);

	bool _requiresEstimation(){ return _estimatetionRequired; };

public:
	TRecalibrationBQSR(TParameters & params, TLog* Logfile, TReadGroups* ReadGroups, TReadGroupMap* ReadGroupMap);
	~TRecalibrationBQSR(){
		for(int i=0; i<_readGroupMap->numReadGroups; ++i){
			delete[] _BQSR_cells_readGroup_quality[i];
		}
		delete[] _BQSR_cells_readGroup_quality;
		delete qualityIndex;

		if(_considerPosition){
			for(int i=0; i<_readGroupMap->numReadGroups; ++i){
				delete[] BQSR_cells_readGroup_position[i];
			}
			delete[] BQSR_cells_readGroup_position;
		}

		if(_considerPositionReverse){
			for(int i=0; i<_readGroupMap->numReadGroups; ++i){
				delete[] _BQSR_cells_readGroup_position_reverse[i];
			}
			delete[] _BQSR_cells_readGroup_position_reverse;
		}

		if(_considerContext){
			for(int i=0; i<_readGroupMap->numReadGroups; ++i){
				delete[] _BQSR_cells_readGroup_context[i];
			}
			delete[] _BQSR_cells_readGroup_context;
		}
	};

	bool recalibrationChangesQualities(){ return true; };
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
//	double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	double getErrorRate(TBase & base);
//	int getQuality(TBase & base);
};

*/

