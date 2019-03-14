#include "TRecalibrationBQSR.h"
//---------------------------------------------------------------
//TBQSR_cell_base BQSR
//---------------------------------------------------------------
TBQSR_cell_base::TBQSR_cell_base(){
	curEstimate = 0.01;
	estimationConverged = false;
	hasData = false;
	firstDerivative = 0.0;
	firstDerivativeSave = 0.0;
	secondDerivative = 0.0;
	secondDerivativeSave = 0.0;
	numObservations = 0.0;
	numObservationsTmp = 0.0;
	F = 999999999999.0;
	oldF = 999999999999.0;
	LL = 0.0;
	myReadGroup = -1;
	store = false;
	batchSize = 100000;
	next = 0;
}

void TBQSR_cell_base::empty(){
	if(!estimationConverged){
		numObservationsTmp = 0;
		firstDerivativeSave = firstDerivative;
		secondDerivativeSave = secondDerivative;
		firstDerivative = 0.0;
		secondDerivative = 0.0;
		LL = 0.0;
	}
}

void TBQSR_cell_base::init(int ReadGroup){
	store = false;
	myReadGroup = ReadGroup;
}

void TBQSR_cell_base::reopenEstimation(){
	estimationConverged = false;
	empty();
}

void TBQSR_cell_base::set(float error, std::string & NumObservations){
	curEstimate = error;
	if(curEstimate <= 0.0) curEstimate = 0.000000001;
	if(curEstimate >= 1.0) curEstimate = 0.9;
	if(NumObservations == "-") numObservations = 0;
	else numObservations = pow(10.0, stringToDouble(NumObservations));
};

float TBQSR_cell_base::getD(TBase* base, Base & RefBase){
	float D = 0.0;
	switch(base->getBaseAsEnum()){
		case A: if(RefBase == A){
					D = 1.0;
					break;
				}
				if(RefBase == G) D = base->PMD_GA;
				break;
		case C: if(RefBase == C) D = 1.0 - base->PMD_CT;
				break;
		case G: if(RefBase == G) D = 1.0 - base->PMD_GA;
				break;
		case T: if(RefBase == C) D = base->PMD_CT;
		        else if(RefBase == T) D = 1.0;
				break;
		case N: throw "Can not add base with unknown reference to BQSR cell!";
	}
	return D;
}

void TBQSR_cell_base::runNewtonRaphson(float & convergenceThreshold, bool & allowIncreaseInF){
	if(F != F) throw "F is not a number!";
	oldF = F;

	curEstimate = curEstimate - firstDerivative / secondDerivative;
	//decide on convergence
	F = fabs(firstDerivative / (float) numObservations);

	if(F < convergenceThreshold) estimationConverged = true;
	if(oldF < F && !allowIncreaseInF) estimationConverged = true;
}


std::string TBQSR_cell_base::getNumObsForPrinting(){
	if(numObservations == 0) return "-";
	else return toString(log10((double) numObservations));
}

void TBQSR_cell_base::calcLikelihoodSurfaceAt(int numPositions, double* positions, std::string & tag, std::ofstream & out){
	bool estimationConvergedTmp = estimationConverged;
	estimationConverged = false;
	float curEstimateTmp = curEstimate;

	for(int i=0; i<numPositions; ++i){
		curEstimate = positions[i];
		recalculateDerivativesFromDataInMemory();
		recalculateLLFromDataInMemory();
		out << tag << "\t" << positions[i] << "\t" << LL << " \t" << firstDerivative << "\t" << secondDerivative << std::endl;
	}

	curEstimate = curEstimateTmp;
	estimationConverged = estimationConvergedTmp;
}

//---------------------------------------------------------------
//TBQSR_cell BQSR
//---------------------------------------------------------------
TBQSR_cellQuality::TBQSR_cellQuality():TBQSR_cell_base(){
	numMatches = 0;
	pointerToBatch = NULL;
}

void TBQSR_cellQuality::init(int ReadGroup, float initialError){
	TBQSR_cell_base::init(ReadGroup);

	curEstimate = initialError;
	if(curEstimate <= 0.0) curEstimate = 0.000000001;
	if(curEstimate >= 1.0) curEstimate = 0.9;

	//storage
	store = false;
}

void TBQSR_cellQuality::doStoreData(){
	store = true;
	D_storage.push_back(new float[batchSize]);
	batchIt = D_storage.rbegin();
	pointerToBatch = *batchIt;
};

void TBQSR_cellQuality::empty(){
	if(!estimationConverged){
		TBQSR_cell_base::empty();
		if(!store) numMatches = 0;
	} else {
		clearStorage();
	}
}

void TBQSR_cellQuality::clearStorage(){
	if(store){
		for(batchIt = D_storage.rbegin(); batchIt != D_storage.rend(); ++batchIt)
			delete[] *batchIt;
		D_storage.clear();
		next = 0;
	}
}

void TBQSR_cellQuality::addBase(TBase* base, Base & RefBase, TQualityMap & qualiMap){
	if(!estimationConverged){
		if(store){
			if(next == batchSize){
				//add new batch
				D_storage.push_back(new float[batchSize]);
				batchIt = D_storage.rbegin();
				pointerToBatch = *batchIt;
				next = 0;
			}

			//add D to batch
			pointerToBatch[next] = getD(base, RefBase);
			addToDerivatives(pointerToBatch[next]);
			++next;
		} else {
			float D = getD(base, RefBase);
			addToDerivatives(D);
		}
		++numObservationsTmp;
		if(base->getBaseAsEnum() == RefBase) ++numMatches;
	}
}

void TBQSR_cellQuality::addToDerivatives(float & D){
	float oneMinus4D = 1.0 - 4.0 * D;
	if(curEstimate >= 1.0) curEstimate = 0.999999;
	firstDerivative += oneMinus4D / (-4.0*D*curEstimate + 3.0*D + curEstimate);
	float tmpF = oneMinus4D / (D*(3.0-4.0*curEstimate) + curEstimate);
	secondDerivative -= tmpF * tmpF;
}

void TBQSR_cellQuality::addToLL(float & D){
	LL += log((1.0-D)*curEstimate/3.0 + D*(1.0-curEstimate));
}

void TBQSR_cellQuality::recalculateDerivativesFromDataInMemory(){
	if(!estimationConverged){
		//set to zero
		empty();

		//first the last batch, which is not filled to the end
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
		for(int i=0; i<next; ++i){ //next is set when adding sites
			addToDerivatives(pointerToBatch[i]);
		}

		//and now the other batches
		++batchIt;
		for(;batchIt != D_storage.rend(); ++batchIt){
			pointerToBatch = *batchIt;
			for(int i=0; i<batchSize; ++i){
				addToDerivatives(pointerToBatch[i]);
			}
		}
	}
}

void TBQSR_cellQuality::recalculateLLFromDataInMemory(){
	LL = 0.0;

	//first the last batch, which is not filled to the end
	batchIt = D_storage.rbegin();
	pointerToBatch = *batchIt;
	for(int i=0; i<next; ++i){ //next is set when adding sites
		addToLL(pointerToBatch[i]);
	}

	//and now the other batches
	++batchIt;
	for(;batchIt != D_storage.rend(); ++batchIt){
		pointerToBatch = *batchIt;
		for(int i=0; i<batchSize; ++i){
			addToLL(pointerToBatch[i]);
		}
	}
}

void TBQSR_cellQuality::runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon, bool & allowIncreaseInF){
	//need Newton-Raphson to estimate tmpEpsilon
	float oldEstimate = curEstimate;

	runNewtonRaphson(convergenceThreshold, allowIncreaseInF);

	//check boundaries
	if(curEstimate <= 0.0){
		curEstimate = 0.000000001;
		if(oldEstimate == 0.00000001)
			estimationConverged = true; //if estimate is repeatedly below, accept
	} else if(curEstimate >= 1.0){
		curEstimate = 0.999999999;
		if(oldEstimate == 0.999999999)
			estimationConverged = true; //if estimate is repeatedly above, accept
	}

	//do not allow big jump in quality -> max +/- 10!
	if(curEstimate / oldEstimate > 10.0){
		curEstimate = oldEstimate * 10.0;
	} else if(oldEstimate / curEstimate > 10.0){
		curEstimate = oldEstimate / 10.0;
	}

	//check if quality did not change
	if(fabs(makePhred(oldEstimate) - makePhred(curEstimate)) < minEpsilon) estimationConverged = true;
}

bool TBQSR_cellQuality::estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations, bool & allowIncreaseInF){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		if(store){
			numObservations = (D_storage.size() - 1) * batchSize + next;
		} else numObservations = numObservationsTmp;

		//check if there is sufficient data
		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
			hasData = false;
			return estimationConverged;
		} else hasData = true;

		//estimate
		if(numMatches >= numObservations){ //tmpEpsilon = 0
			curEstimate = 1.0 / (double) (numObservations + 1.0);
			estimationConverged = true;
		} else if(numMatches < 1.0){ // tmpEpsilon = 1.0
			std::cout << "numMatches < 1" << std::endl;
			curEstimate = 1.0;
			estimationConverged = true;
		} else {
			//need Newton-Raphson to estimate tmpEpsilon
			runNewtonRaphsonAndCheck(convergenceThreshold, minEpsilon, allowIncreaseInF);
		}
	}

	return estimationConverged;
}

void TBQSR_cellQuality::calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out){
	double* positions = new double[numPositions];
	//now fill between 0.000000001 and 0.999999999
	double delta = (0.999999999 - 0.000000001) / (numPositions - 1.0);
	for(int i=0; i<numPositions; ++i){
		positions[i] = 0.000000001 + delta * (double) i;
	}

	//calc surface!
	calcLikelihoodSurfaceAt(numPositions, positions, tag, out);
}

//---------------------------------------------------------------
//TBQSR_cellPosition BQSR
//---------------------------------------------------------------
TBQSR_cellPosition::TBQSR_cellPosition():TBQSR_cell_base(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
	pointerToBatch = NULL;
}

void TBQSR_cellPosition::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup){
	TBQSR_cell_base::init(ReadGroup);
	BQSR_cells_readGroup_quality = gotBQSR_cells_quality_readGroup;
	qualityIndex = QualityIndex;
}

void TBQSR_cellPosition::doStoreData(){
	D_storage.push_back(new BQSRFactorStorage[batchSize]);
	batchIt = D_storage.rbegin();
	pointerToBatch = *batchIt;
};

void TBQSR_cellPosition::clearStorage(){
	if(store){
		for(batchIt = D_storage.rbegin(); batchIt != D_storage.rend(); ++batchIt)
			delete[] *batchIt;
		D_storage.clear();
		next = 0;
	}
}

void TBQSR_cellPosition::addToDerivatives(float & D, float & epsilon){
	double epsMinus4Deps = epsilon - 4.0 * D * epsilon;
	firstDerivative += epsMinus4Deps / (-4.0*D*epsilon*curEstimate + 3.0*D + epsilon*curEstimate);
	double tmpF = epsMinus4Deps / (D*(3.0-4.0*epsilon*curEstimate) + epsilon*curEstimate);
	secondDerivative -= tmpF * tmpF;
}

void TBQSR_cellPosition::addToLL(float & D, float & epsilon){
	LL += log((1.0-D)*epsilon*curEstimate/3.0 + D*(1.0-epsilon*curEstimate));
}

float TBQSR_cellPosition::getEpsilon(TBase* base, TQualityMap & qualMap){
	return  BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndexFromQuality(qualMap.errorToPhredInt(base->errorRate))].curEstimate;
}

void TBQSR_cellPosition::addBase(TBase* base, Base & RefBase, TQualityMap & qualiMap){
	if(!estimationConverged){
		if(store){
			if(next == batchSize){
				//add new batch
				D_storage.push_back(new BQSRFactorStorage[batchSize]);
				batchIt = D_storage.rbegin();
				pointerToBatch = *batchIt;
				next = 0;
			}

			//add D to batch
			pointerToBatch[next].D = getD(base, RefBase);
			pointerToBatch[next].epsilon = getEpsilon(base, qualiMap);
			addToDerivatives(pointerToBatch[next].D, pointerToBatch[next].epsilon);
			++next;
		} else {
			float D = getD(base, RefBase);
			float eps = getEpsilon(base, qualiMap);
			addToDerivatives(D, eps);
		}
		++numObservationsTmp;
	}
};

void TBQSR_cellPosition::recalculateDerivativesFromDataInMemory(){
	if(!estimationConverged){
		//set to zero
		empty();

		//first the last batch, which is not filled to the end
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
		for(int i=0; i<next; ++i){
			addToDerivatives(pointerToBatch[i].D, pointerToBatch[i].epsilon);
		}

		//and now the other batches
		++batchIt;
		for(;batchIt != D_storage.rend(); ++batchIt){
			pointerToBatch = *batchIt;
			for(int i=0; i<batchSize; ++i){
				addToDerivatives(pointerToBatch[i].D, pointerToBatch[i].epsilon);
			}
		}
	}
};

void TBQSR_cellPosition::recalculateLLFromDataInMemory(){
	LL = 0.0;

	//first the last batch, which is not filled to the end
	batchIt = D_storage.rbegin();
	pointerToBatch = *batchIt;
	for(int i=0; i<next; ++i){
		addToLL(pointerToBatch[i].D, pointerToBatch[i].epsilon);
	}

	//and now the other batches
	++batchIt;
	for(;batchIt != D_storage.rend(); ++batchIt){
		pointerToBatch = *batchIt;
		for(int i=0; i<batchSize; ++i){
			addToLL(pointerToBatch[i].D, pointerToBatch[i].epsilon);
		}
	}
}

void TBQSR_cellPosition::runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon, bool & allowIncreaseInF){
	//need Newton-Raphson to estimate tmpEpsilon
	float oldEstimate = curEstimate;
	runNewtonRaphson(convergenceThreshold, allowIncreaseInF);

	//check boundaries
	if(curEstimate < 0.0){
		curEstimate = 0.01;
		if(oldEstimate == 0.01)
			estimationConverged = true; //if estimate is repeatedly below, accept
	} else if(curEstimate > 10000.0){
		curEstimate = 100.0;
		if(oldEstimate == 100.0)
			estimationConverged = true; //if estimate is repeatedly above, accept
	}

	//check if quality did not change
	if(fabs(oldEstimate - curEstimate) < minEpsilon) estimationConverged = true;
}


bool TBQSR_cellPosition::estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations, bool & allowIncreaseInF){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		if(store){
			numObservations = (D_storage.size() - 1) * batchSize + next;
		} else numObservations = numObservationsTmp;

		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
		} else {
			//need Newton-Raphson to estimate tmpEpsilon
			runNewtonRaphsonAndCheck(convergenceThreshold, minEpsilon, allowIncreaseInF);
		}
	}

	return estimationConverged;
}

void TBQSR_cellPosition::calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out){
	double* positions = new double[numPositions];
	//now fill between 0.01 and 100 -> use log10
	double delta = 4.0 / (numPositions - 1.0);
	for(int i=0; i<numPositions; ++i){
		positions[i] = pow(10.0, -2.0 + delta * (double) i);
	}

	//calc surface!
	calcLikelihoodSurfaceAt(numPositions, positions, tag, out);
}

//---------------------------------------------------------------
TBQSR_cellPositionRev::TBQSR_cellPositionRev():TBQSR_cellPosition(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
	BQSR_cells_readGroup_position = NULL;
	considerPosition = false;
}

void TBQSR_cellPositionRev::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup){
	TBQSR_cellPosition::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position = gotBQSR_cells_position_readGroup;
	considerPosition = true;
}

void TBQSR_cellPositionRev::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup){
	TBQSR_cellPosition::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position = NULL;
}

float TBQSR_cellPositionRev::getEpsilon(TBase* base, TQualityMap & qualMap){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndexFromQuality(qualMap.errorToPhredInt(base->errorRate))].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->distFrom5Prime].curEstimate;
	return  epsilonAlpha;
}

//---------------------------------------------------------------
TBQSR_cellContext::TBQSR_cellContext():TBQSR_cellPositionRev(){
	BQSR_cells_readGroup_position_rev = NULL;
	considerPositionRev = false;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup, gotBQSR_cells_position_readGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_cells_positionRev_readGroup;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup, gotBQSR_cells_position_readGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_cells_positionRev_readGroup;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cellQuality** gotBQSR_cells_quality_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

float TBQSR_cellContext::getEpsilon(TBase* base, TQualityMap & qualMap){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndexFromQuality(qualMap.errorToPhredInt(base->errorRate))].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->distFrom5Prime].curEstimate;
	if(considerPositionRev) epsilonAlpha *= BQSR_cells_readGroup_position_rev[myReadGroup][base->distFrom3Prime].curEstimate;
	return  epsilonAlpha;
};

//---------------------------------------------------------------
//TRecalibrationBQSRStorage
//---------------------------------------------------------------
TRecalibrationBQSRStorage::TRecalibrationBQSRStorage(){
	numReadGroups = 0;
	numQuality = 0;
	maxPos = 0;
	maxPosReverse = 0;
	numContexts = 20;

	considerQuality = false;
	considerPosition = false;
	considerPositionReverse = false;
	considerContext = false;

	qualityCells = NULL;
	positionCells = NULL;
	positionReverseCells = NULL;
	contextCells = NULL;
};

TRecalibrationBQSRStorage::~TRecalibrationBQSRStorage(){
	if(considerQuality){
		for(int i=0; i<numReadGroups; ++i){
			delete[] qualityCells[i];
		}
		delete[] qualityCells;
	}

	if(considerPosition){
		for(int i=0; i<numReadGroups; ++i){
			delete[] positionCells[i];
		}
		delete[] positionCells;
	}

	if(considerPositionReverse){
		for(int i=0; i<numReadGroups; ++i){
			delete[] positionReverseCells[i];
		}
		delete[] positionReverseCells;
	}

	if(considerContext){
		for(int i=0; i<numReadGroups; ++i){
			delete[] contextCells[i];
		}
		delete[] contextCells;
	}
};

void TRecalibrationBQSRStorage::initializeQualityCells(int NumReadGroups, int NumQuality, TQualityMap & qualityMap, TQualityIndex* qualityIndex){
	numReadGroups = NumReadGroups;
	numQuality = NumQuality;
	for(int i=0; i<numReadGroups; ++i){
		qualityCells[i] = new TBQSR_cellQuality[NumQuality];
		for(int q=0; q<NumQuality; ++q)
			qualityCells[i][q].init(i, qualityMap.qualityToError(qualityIndex->getPhredIntFromIndex(q)));
	}
	considerQuality = true;
};

void TRecalibrationBQSRStorage::initializePositionCells(int NumReadGroups, int MaxPos, TQualityIndex* qualityIndex){
	if(NumReadGroups != numReadGroups)
		throw "Unequal number of read groups in TRecalibrationBQSRStorage!";
	maxPos = MaxPos;

	for(int r=0; r<numReadGroups; ++r){
		positionCells[r] = new TBQSR_cellPosition[maxPos];
		for(int p=0; p<maxPos; ++p)
			positionCells[r][p].init(r, qualityIndex, qualityCells);
	}
	considerPosition = true;
};

void TRecalibrationBQSRStorage::initializePositionReverseCells(int NumReadGroups, int MaxPos, TQualityIndex* qualityIndex){
	if(NumReadGroups != numReadGroups)
		throw "Unequal number of read groups in TRecalibrationBQSRStorage!";
	maxPosReverse = MaxPos;

	for(int r=0; r<numReadGroups; ++r){
		positionReverseCells[r] = new TBQSR_cellPositionRev[maxPosReverse];
		for(int p=0; p< maxPosReverse; ++p){
			if(considerPosition) positionReverseCells[r][p].init(r, qualityIndex, qualityCells, positionCells);
			else positionReverseCells[r][p].init(r, qualityIndex, qualityCells);
		}
	}

	considerPositionReverse = true;
};

void TRecalibrationBQSRStorage::initializeContextCells(int NumReadGroups, TQualityIndex* qualityIndex){
	if(NumReadGroups != numReadGroups)
		throw "Unequal number of read groups in TRecalibrationBQSRStorage!";
	contextCells = new TBQSR_cellContext*[numReadGroups];
	for(int r=0; r<numReadGroups; ++r){
		contextCells[r] = new TBQSR_cellContext[numContexts];
		for(int c=0; c<numContexts; ++c){
			if(considerPosition && considerPositionReverse) contextCells[r][c].init(r, qualityIndex, qualityCells, positionCells, positionReverseCells);
			else if(considerPosition && !considerPositionReverse) contextCells[r][c].init(r, qualityIndex, qualityCells, positionCells);
			else if(!considerPosition && considerPositionReverse) contextCells[r][c].init(r, qualityIndex, qualityCells, positionReverseCells);
			else contextCells[r][c].init(r, qualityIndex, qualityCells);
		}
	}
};

double TRecalibrationBQSRStorage::getMaxValueFromFile(std::string filename, int col, int numCol){
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "' for reading!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;

	//skip header
	std::string tmpF;
	std::getline(file, tmpF);
	++lineNum;

	double max = 0.0;
	bool maxInitialized = false;

	//parse file to get max position
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != numCol) throw "Found " + toString(vec.size()) + " instead of " + toString(numCol) + " columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			double m = stringToDouble(vec[col]);
			if(!maxInitialized)
				max = m;
			else if(m > max) max = m;
		}
	}

	return max;
};

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(TLog* Logfile, TReadGroups* ReadGroups){
	_initialize(Logfile, ReadGroups);
};

TRecalibrationBQSR::TRecalibrationBQSR(TParameters & params, TLog* Logfile, TReadGroups* ReadGroups):TRecalibration(){
	_initialize(Logfile, ReadGroups);

	//initialize quality effect
	std::string filename = params.getParameterString("BQSRQuality");
	_initializeBQSRReadGroupQualityTableFromFile(filename);

	//Do we consider the effect of the position in read (cycle)?
	if(params.parameterExists("BQSRPosition")){
		filename = params.getParameterString("BQSRPosition");
		_initializeBQSRReadGroupPositionTableFromFile(filename);
	}

	if(params.parameterExists("BQSRPositionReverse")){
		filename = params.getParameterString("BQSRPositionReverse");
		_initializeBQSRReadGroupPositionReverseTableFromFile(filename);
	}

	//Do we consider the context (dinucleotide)?
	if(params.parameterExists("BQSRContext")){
		filename = params.getParameterString("BQSRContext");
		_initializeBQSRReadGroupContextTableFromFile(filename);
	}
};

void TRecalibrationBQSR::_initialize(TLog* Logfile, TReadGroups* ReadGroups){
	_logfile = Logfile;
	_readGroups = ReadGroups;
	qualityIndex = NULL;
	_type = "BQSR";

};

void TRecalibrationBQSR::_initializeBQSRReadGroupQualityTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x quality table from file '" + filename + "' ...");

	//get maxQ from file
	int maxQ = storage.getMaxValueFromFile(filename, 1, 5);

	//initialize quality index
	qualityIndex = new TQualityIndex(0, maxQ);

	//create corresponding objects
	storage.initializeQualityCells(_readGroups->size(), qualityIndex->numQ, _qualityMap, qualityIndex);

	//open file
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//tmp variables for file parsing
	std::vector<std::string> vec;
	long lineNum = 0;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0 && _readGroups->readGroupExists(vec[0])){
			//set quality and empirical error rate
			int readGroup = _readGroups->find(vec[0]);

			if(readGroup >= 0){ //returns -1 if read group does not exist
				int q = stringToInt(vec[1]);
				double phredEmpiric = stringToDouble(vec[3]);
				storage.qualityCells[readGroup][qualityIndex->getIndexFromQuality(q+33)].set(_qualityMap.phredToError(phredEmpiric), vec[4]);
			} else throw "readGroup " + vec[0] + " does not exist in BAM file header!";
		}
	}

	//done!
	_logfile->done();
	_logfile->conclude("Considering qualities up to " + toString(maxQ));
}

void TRecalibrationBQSR::_initializeBQSRReadGroupPositionTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x position table from file '" + filename + "' ...");

	//get maxPos from file
	int maxPos = storage.getMaxValueFromFile(filename, 1, 5);

	//create corresponding objects
	storage.initializePositionCells(_readGroups->size(), maxPos, qualityIndex);

	//open file
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position table from file '" + filename + "'!";

	//tmp variables for file parsing
	std::vector<std::string> vec;
	long lineNum = 0;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0 && _readGroups->readGroupExists(vec[0])){
			//set quality and empirical error rate
			int readGroup = _readGroups->find(vec[0]);
			int p = stringToInt(vec[1]);
			double alpha = stringToDouble(vec[3]);
			storage.positionCells[readGroup][p-1].set(alpha, vec[4]);
		}
	}

	//done!
	_logfile->done();
	_logfile->conclude("Considering positions up to " + toString(storage.maxPos));
}

void TRecalibrationBQSR::_initializeBQSRReadGroupPositionReverseTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x position reverse table from file '" + filename + "' ...");

	//get maxPos from file
	int maxPos = storage.getMaxValueFromFile(filename, 1, 5);

	//create corresponding objects
	storage.initializePositionReverseCells(_readGroups->size(), maxPos, qualityIndex);

	//open file
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position reverse table from file '" + filename + "'!";

	//tmp variables for file parsing
	std::vector<std::string> vec;
	long lineNum = 0;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";

			if(_readGroups->readGroupExists(vec[0])){
				//set quality and empirical error rate
				int readGroup = _readGroups->find(vec[0]);
				int p = stringToInt(vec[1]);
				double alpha = stringToDouble(vec[3]);
				storage.positionReverseCells[readGroup][p-1].set(alpha, vec[4]);
			}
		}
	}

	//done!
	_logfile->done();
	_logfile->conclude("Considering positions reverse up to " + toString(storage.maxPos));
}

void TRecalibrationBQSR::_initializeBQSRReadGroupContextTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x context table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x context table from file '" + filename + "'!";

	//construct for each read group in bam file
	storage.initializeContextCells(_readGroups->size(), qualityIndex);

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			if(_readGroups->readGroupExists(vec[0])){
				//set quality and empirical error rate
				int readGroup = _readGroups->find(vec[0]);
				int context = storage._genoMap.getContext(vec[1][0], vec[1][1]);
				double alpha = stringToDouble(vec[3]);
				storage.contextCells[readGroup][context].set(alpha, vec[4]);
			}
		}
	}

	//done!
	_logfile->done();
	_logfile->conclude("Considering context");
}

double TRecalibrationBQSR::getErrorRate(TBase & base){
	double q = storage.qualityCells[base.readGroup][qualityIndex->getIndexFromQuality(_qualityMap.errorToQuality(base.errorRate))].curEstimate;

	if(storage.considerPosition){
		if(base.distFrom5Prime > storage.maxPos)
			q *= storage.positionCells[base.readGroup][storage.maxPos].curEstimate;
		else
			q *= storage.positionCells[base.readGroup][base.distFrom5Prime].curEstimate;
	}

	if(storage.considerPositionReverse) q *= storage.positionReverseCells[base.readGroup][base.distFrom3Prime].curEstimate;

	if(storage.considerContext) q *= storage.contextCells[base.readGroup][base.context].curEstimate;
	if(q > 1.0) q = 1.0; //make sure the scaling does not lead to errors > 1.0!
	return q;
};


//---------------------------------------------------------------
//RecalibrationBQSREstimator
//---------------------------------------------------------------
TRecalibrationBQSREstimator::TRecalibrationBQSREstimator(TParameters & params, TLog* Logfile, TReadGroups* ReadGroups, TReadGroupMap* ReadGroupMap){
	_readGroupMap = ReadGroupMap;
	_estimatetionRequired = false;
	_estimationConverged = false;
	_curNewtonRaphsonLoop = 0;
	qualityIndex = NULL;
	_logfile = Logfile;
	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMap;

	_storeDataInMemory = params.parameterExists("storeInMemory");
	if(_storeDataInMemory) _logfile->list("Will store D in memory to iterate Newton-Raphson faster");
	_dataStored = false;

	//check if BQSR table readGroup x Quality is given, or has to be estimated
	_initializeBQSRReadGroupQualityTable(params);

	//Do we also consider the effect of the position in read (cycle)?
	_initializeBQSRReadGroupPositionTable(params);
	_initializeBQSRReadGroupPositionReverseTable(params);

	//Do we also consider the context (dinucleotide)?
	_initializeBQSRReadGroupContextTable(params);

	//read Newton-Raphson arguments from user
	_convergenceThreshold_F = params.getParameterDoubleWithDefault("maxF", 0.0000001);
	_minEpsilonQuality = params.getParameterDoubleWithDefault("minEpsQuality", 0.000001);
	_minEpsilonFactors = params.getParameterDoubleWithDefault("minEpsFactors", 0.0001);
	_numLoopIncreaseFAllowed = params.getParameterIntWithDefault("numLoopIncreaseFAllowed", 3);
	if(_estimatetionRequired){
		_logfile->startIndent("Conditions to stop Newton-Raphson algorithm:");
		_logfile->list("Stopping Newton-Raphson if F < " + toString(_convergenceThreshold_F));
		_logfile->list("Allow F to increase in first " + toString(_numLoopIncreaseFAllowed) + " loops.");
		_logfile->list("Stopping Newton-Raphson if the change in quality is < " + toString(_minEpsilonQuality));
		_logfile->list("Stopping Newton-Raphson if the change in a factor (e.g. position) is < " + toString(_minEpsilonFactors));
		_logfile->endIndent();
	}

	//get minimal number of observations to conduct estimation
	_minObservations = params.getParameterLongWithDefault("minObservations", 32000);

	//do we print LL surfaces?
	_numPosLLsurface = params.getParameterIntWithDefault("LLSurface", 0);
	if(_numPosLLsurface > 0) _printLLSurface = true;
	else _printLLSurface = false;
	_LLSurfacePrinted = false;
}

void TRecalibrationBQSREstimator::_initializeBQSRReadGroupQualityTable(TParameters & params){

	if(params.parameterExists("knownBQSRQuality")){
		_qualityConverged = true;
		_estimateQuality = false;

		//initialize quality index
		std::string filename = params.getParameterString("knownQuality");
		int maxQ = storage.getMaxValueFromFile(filename, 1, 5);
		qualityIndex = new TQualityIndex(0, maxQ);

		//initialize BQSR table
		storage.initializeQualityCells(_readGroupMap->getNumReadGroups(), qualityIndex->numQ, _qualityMap, qualityIndex);
		_fillBQSRReadGroupQualityTableFromFile(filename);
	} else {
		_logfile->list("Will estimate BQSR quality table.");
		_qualityConverged = false;
		_estimateQuality = true;
		_estimatetionRequired = true;

		//initialize quality index
		int minQ = params.getParameterIntWithDefault("minQ", 0);
		int maxQ = params.getParameterIntWithDefault("maxQ", 100);
		_logfile->list("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
		qualityIndex = new TQualityIndex(minQ, maxQ);

		//initialize BQSR table
		storage.initializeQualityCells(_readGroupMap->getNumReadGroups(), qualityIndex->numQ, _qualityMap, qualityIndex);

		//check if we set initial values from file
		if(params.parameterExists("initialBQSRQuality")){
			_fillBQSRReadGroupQualityTableFromFile(params.getParameterString("initialQuality"));
			_qualityConverged = false;
			_estimateQuality = true;
		}
	}
};

void TRecalibrationBQSREstimator::_fillBQSRReadGroupQualityTableFromFile(std::string filename){
	_logfile->listFlush("Reading BQSR readGroup x quality table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			if(_readGroups->readGroupExists(vec[0])){
				//get read group and match it to internal read group ID
				int readGroupIndex = _readGroupMap->getIndex(_readGroups->find(vec[0]));
				int q = stringToInt(vec[1]);
				float phredEmpiric = _qualityMap.phredToError(stringToDouble(vec[3]));
				storage.qualityCells[readGroupIndex][qualityIndex->getIndexFromPhredInt(q)].set(phredEmpiric);
			}
		}
	}

	//done!
	_logfile->done();
};

void TRecalibrationBQSREstimator::_initializeBQSRReadGroupPositionTable(TParameters & params){
	if(params.parameterExists("estimateBQSRPosition") || params.parameterExists("initialBQSRPosition")){
		_logfile->list("Will estimate BQSR position table.");
		_estimatetionRequired = true;
		_estimatePosition = true;
		_positionConverged = false;

		//get max pos
		int maxPos = params.getParameterInt("maxPos");
		if(maxPos < 1) throw "Max position has to be larger than zero!";
		_logfile->list("Considering positions up to " + toString(maxPos));
		storage.initializePositionCells(_readGroupMap->getNumReadGroups(), maxPos, qualityIndex);

		//check if we set initial values from file
		if(params.parameterExists("initialBQSRPosition"))
			_fillBQSRReadGroupPositionTableFromFile(params.getParameterString("initialBQSRPosition"));

	} else {
		_estimatePosition = false;
		_positionConverged = true;
		if(params.parameterExists("knownBQSRPosition")){
			//get maxPos from file
			std::string filename = params.getParameterString("knownBQSRPosition");
			int maxPos = storage.getMaxValueFromFile(filename, 1, 5);

			//read file
			_fillBQSRReadGroupPositionTableFromFile(filename);
		}
	}
};

void TRecalibrationBQSREstimator::_fillBQSRReadGroupPositionTableFromFile(std::string filename){
	_logfile->listFlush("Reading BQSR readGroup x position table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position table from file '" + filename + "'!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			if(_readGroups->readGroupExists(vec[0])){
				//set quality and empirical error rate
				int readGroup = _readGroups->find(vec[0]);
				int p = stringToInt(vec[1]);
				if(p <= storage.maxPos){
					double alpha = stringToDouble(vec[3]);
					storage.positionCells[readGroup][p-1].set(alpha, vec[4]);
				}
			}
		}
	}

	//done!
	_logfile->done();
};

void TRecalibrationBQSREstimator::_initializeBQSRReadGroupPositionReverseTable(TParameters & params){
	if(params.parameterExists("estimateBQSRPositionReverse") || params.parameterExists("initialBQSRPosition")){
		_logfile->list("Will estimate BQSR position reverse table.");
		_estimatetionRequired = true;
		_estimatePositionReverse = true;
		_positionReverseConverged = false;

		//get max pos
		int maxPos = 0;
		if(params.parameterExists("maxPosReverse"))
			maxPos = params.getParameterInt("maxPosReverse");
		else
			maxPos = params.getParameterInt("maxPos");
		if(maxPos < 1) throw "Max position has to be larger than zero!";

		_logfile->list("Considering reverse positions up to " + toString(maxPos));
		storage.initializePositionReverseCells(_readGroupMap->getNumReadGroups(), maxPos, qualityIndex);

		//check if we set initial values from file
		if(params.parameterExists("initialBQSRPosition"))
			_fillBQSRReadGroupPositionTableFromFile(params.getParameterString("initialBQSRPositionReverse"));

	} else {
		_estimatePositionReverse = false;
		_positionReverseConverged = true;
		if(params.parameterExists("knownBQSRPositionReverse")){
			//get maxPos from file
			std::string filename = params.getParameterString("knownBQSRPositionReverse");
			int maxPos = storage.getMaxValueFromFile(filename, 1, 5);

			//read file
			_fillBQSRReadGroupPositionReverseTableFromFile(filename);
		}
	}
};

void TRecalibrationBQSREstimator::_fillBQSRReadGroupPositionReverseTableFromFile(std::string filename){
	_logfile->listFlush("Reading BQSR readGroup x reverse position table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x reverse position table from file '" + filename + "'!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			if(_readGroups->readGroupExists(vec[0])){
				//set quality and empirical error rate
				int readGroup = _readGroups->find(vec[0]);
				int p = stringToInt(vec[1]);
				if(p <= storage.maxPosReverse){
					double alpha = stringToDouble(vec[3]);
					storage.positionReverseCells[readGroup][p-1].set(alpha, vec[4]);
				}
			}
		}
	}

	//done!
	_logfile->done();
};

void TRecalibrationBQSREstimator::_initializeBQSRReadGroupContextTable(TParameters & params){
	if(params.parameterExists("estimateBQSRContext") || params.parameterExists("initialBQSRContext")){
		_logfile->list("Will estimate BQSR context table.");
		_estimatetionRequired = true;
		_estimateContext = true;
		_contextConverged = false;

		storage.initializeContextCells(_readGroupMap->getNumReadGroups(), qualityIndex);

		//check if we set initial values from file
		if(params.parameterExists("initialBQSRContext"))
			_fillBQSRReadGroupPositionTableFromFile(params.getParameterString("initialBQSRContext"));

	} else {
		_estimateContext = false;
		_contextConverged = true;
		if(params.parameterExists("knownBQSRContext"))
			_fillBQSRReadGroupPositionTableFromFile(params.getParameterString("knownBQSRContext"));
	}
};

void TRecalibrationBQSREstimator::_fillBQSRReadGroupContextTableFromFile(std::string filename){
	_logfile->listFlush("Reading BQSR readGroup x context table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x context table from file '" + filename + "'!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//set quality and empirical error rate

			if(_readGroups->readGroupExists(vec[0])){
				//set quality and empirical error rate
				int readGroup = _readGroups->find(vec[0]);
				int context = storage._genoMap.getContext(vec[1][0], vec[1][1]);
				double alpha = stringToDouble(vec[3]);
				storage.contextCells[readGroup][context].set(alpha, vec[4]);
			}
		}
	}

	//done!
	_logfile->done();
	_logfile->conclude("Considering context");
};

void TRecalibrationBQSREstimator::addSite(TSite & site, TQualityMap & qualMap){
	if(site.referenceBase != 'N'){
		Base refBase = site.getRefBaseAsEnum();
		if(!_qualityConverged){
			for(TBase* it : site.bases){
				storage.qualityCells[_readGroupMap->getIndex(it->readGroup)][qualityIndex->getIndexFromQuality(qualMap.errorToPhredInt(it->errorRate))].addBase(it, refBase, qualMap);
			}
		}
		else if(storage.considerPosition && !_positionConverged){
			for(TBase* it : site.bases){
				if(it->distFrom5Prime >= storage.maxPos) throw "Position of base is > maxPos specified!";
				storage.positionCells[_readGroupMap->getIndex(it->readGroup)][it->distFrom5Prime].addBase(it, refBase, qualMap);
			}
		}
		else if(storage.considerPositionReverse && !_positionReverseConverged){
			for(TBase* it : site.bases){
				if(it->distFrom3Prime >= storage.maxPosReverse) throw "Position of base is > maxPos specified!";
				storage.positionReverseCells[_readGroupMap->getIndex(it->readGroup)][it->distFrom3Prime].addBase(it, refBase, qualMap);
			}
		} else if(storage.considerContext && !_contextConverged){
			for(TBase* it : site.bases)
				storage.contextCells[_readGroupMap->getIndex(it->readGroup)][it->context].addBase(it, refBase, qualMap);
		}
	}
};

void TRecalibrationBQSREstimator::recalculateDerivativesFromDataInMemory(){
	if(!_qualityConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int j=0; j<qualityIndex->numQ; ++j){
				storage.qualityCells[r][j].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(storage.considerPosition && !_positionConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int p=0; p<storage.maxPos; ++p){
				storage.positionCells[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(storage.considerPositionReverse && !_positionReverseConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int p=0; p<storage.maxPosReverse; ++p){
				storage.positionReverseCells[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	} else if(storage.considerContext && !_contextConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int c=0; c<storage.numContexts; ++c){
				storage.contextCells[r][c].recalculateDerivativesFromDataInMemory();
			}
		}
	}
};

bool TRecalibrationBQSREstimator::estimateEpsilon(std::string filenameTag){
	++_curNewtonRaphsonLoop;
	bool allowIncreaseInF = false;
	if(_curNewtonRaphsonLoop <= _numLoopIncreaseFAllowed) allowIncreaseInF = true;

	//recalc derivatives if data is in memory. Otherwise, derivatives were calculated when data was added.
	if(_dataStored) recalculateDerivativesFromDataInMemory();

	//estimate tmpEpsilon, if not yet done
	_estimationConverged = true;
	int numCellsNotConverged = 0;
	double maxF = 0.0;

	//readGroup x quality
	//-------------------------------------------------------
	if(!_qualityConverged){
		//do we print LL surface? Only print once!
		if(_printLLSurface && !_LLSurfacePrinted){
			calculateAndPrintLLSurfaceQuality(filenameTag);
			_LLSurfacePrinted = true;
		}
		//now do estimation
		_logfile->listFlush("Estimating tmpEpsilon for readGroup x quality table ...");
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int j=0; j<qualityIndex->numQ; ++j){
				if(!storage.qualityCells[i][j].estimate(_convergenceThreshold_F, _minEpsilonQuality, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(storage.qualityCells[i][j].hasData && storage.qualityCells[i][j].F > maxF) maxF = storage.qualityCells[i][j].F;
			}
		}

		//report
		_logfile->done();
		if(numCellsNotConverged == 0){
			_qualityConverged = true;
			_logfile->list("Estimation converged in all cells!");
		} else {
			_qualityConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (_readGroupMap->getNumReadGroups() * qualityIndex->numQ));
			_logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		_logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!_qualityConverged){
			//empty all cells
			for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
				for(int j=0; j<qualityIndex->numQ; ++j){
					storage.qualityCells[i][j].empty();
				}
			}
			_estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(_storeDataInMemory) _dataStored = true;
		} else {
			//write to file
			writeQualityToFile(filenameTag);

			//empty storage
			if(_storeDataInMemory){
				for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
					for(int j=0; j<qualityIndex->numQ; ++j){
						storage.qualityCells[i][j].clearStorage();
					}
				}
			}
			_dataStored = false;

			//what's next?
			if(!storage.considerPosition && !storage.considerPositionReverse && !storage.considerContext) _estimationConverged = true;
			else _estimationConverged = false;
			_LLSurfacePrinted = false;
		}
		return _estimationConverged;
	}

	//estimate tmpEpsilon for position, if not yet done
	//-------------------------------------------------------
	if(storage.considerPosition && !_positionConverged){
		//do we print LL surface? Only print once!
		if(_printLLSurface && !_LLSurfacePrinted){
			calculateAndPrintLLSurfacePosition(filenameTag);
			_LLSurfacePrinted = true;
		}

		//now do estimation
		_logfile->listFlush("Estimating tmpEpsilon for readGroup x position table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<storage.maxPos; ++p){
				if(!storage.positionCells[i][p].estimate(_convergenceThreshold_F, _minEpsilonFactors, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(storage.positionCells[i][p].hasData && storage.positionCells[i][p].F > maxF) maxF = storage.positionCells[i][p].F;
			}
		}

		//report
		_logfile->done();
		if(numCellsNotConverged == 0){
			_positionConverged = true;
			_logfile->list("Estimation converged in all cells!");
		} else {
			_positionConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (_readGroupMap->getNumReadGroups() * storage.maxPos));
			_logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		_logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!_positionConverged){
			//empty all cells
			for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
				for(int p=0; p<storage.maxPos; ++p){
					storage.positionCells[i][p].empty();
				}
			}
			_estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(_storeDataInMemory) _dataStored = true;
		} else {
			//write to file
			writePositionToFile(filenameTag);

			//empty storage
			if(_storeDataInMemory){
				for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
					for(int p=0; p<storage.maxPos; ++p){
						storage.positionCells[i][p].clearStorage();
					}
				}
			}
			_dataStored = false;

			//what's next?
			if(!storage.considerPositionReverse && !storage.considerContext) _estimationConverged = true;
			else _estimationConverged = false;
			_LLSurfacePrinted = false;
		}
		return _estimationConverged;
	}

	//estimate tmpEpsilon for position reverse, if not yet done
	//-------------------------------------------------------
	if(storage.considerPositionReverse && !_positionReverseConverged){
		//do we print LL surface? Only print once!
		if(_printLLSurface && !_LLSurfacePrinted){
			calculateAndPrintLLSurfaceReversePosition(filenameTag);
			_LLSurfacePrinted = true;
		}

		//now do estimation
		_logfile->listFlush("Estimating tmpEpsilon for readGroup x position reverse table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<storage.maxPosReverse; ++p){
				if(!storage.positionReverseCells[i][p].estimate(_convergenceThreshold_F, _minEpsilonFactors, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(storage.positionReverseCells[i][p].hasData && storage.positionReverseCells[i][p].F > maxF) maxF = storage.positionReverseCells[i][p].F;
			}
		}

		//report
		_logfile->done();
		if(numCellsNotConverged == 0){
			_positionReverseConverged = true;
			_logfile->list("Estimation converged in all cells!");
		} else {
			_positionReverseConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (_readGroupMap->getNumReadGroups() * storage.maxPosReverse));
			_logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		_logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!_positionReverseConverged){
			//empty all cells
			for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
				for(int p=0; p<storage.maxPosReverse; ++p){
					storage.positionReverseCells[i][p].empty();
				}
			}
			_estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(_storeDataInMemory) _dataStored = true;
		} else {
			//write to file
			writePositionReverseToFile(filenameTag);

			//empty storage
			if(_storeDataInMemory){
				for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
					for(int p=0; p<storage.maxPosReverse; ++p){
						storage.positionReverseCells[i][p].clearStorage();
					}
				}
			}
			_dataStored = false;

			//what's next?
			if(!storage.considerContext) _estimationConverged = true;
			else _estimationConverged = false;
			_LLSurfacePrinted = false;
		}
		return _estimationConverged;
	}

	//estimate tmpEpsilon for context
	//-------------------------------------------------------
	if(storage.considerContext && !_contextConverged){
		//do we print LL surface? Only print once!
		if(_printLLSurface && !_LLSurfacePrinted){
			calculateAndPrintLLSurfaceContext(filenameTag);
			_LLSurfacePrinted = true;
		}

		//now do estimation
		_logfile->listFlush("Estimating tmpEpsilon for quality x context table ...");
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int c=0; c<storage.numContexts; ++c){
				if(!storage.contextCells[r][c].estimate(_convergenceThreshold_F, _minEpsilonFactors, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(storage.contextCells[r][c].hasData && storage.contextCells[r][c].F > maxF) maxF = storage.contextCells[r][c].F;
			}
		}

		//report
		_logfile->done();
		if(numCellsNotConverged == 0){
			_contextConverged = true;
			_logfile->list("Estimation converged in all cells!");
		} else {
			_contextConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (_readGroupMap->getNumReadGroups() * storage.numContexts));
			_logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		_logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!_contextConverged){
			//empty all cells
			for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
				for(int c=0; c<storage.numContexts; ++c){
					storage.contextCells[r][c].empty();
				}
			}
			_estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(_storeDataInMemory) _dataStored = true;
		} else {
			//write to file
			writeContextToFile(filenameTag);

			//empty storage
			if(_storeDataInMemory){
				for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
					for(int c=0; c<storage.numContexts; ++c){
						storage.contextCells[i][c].clearStorage();
					}
				}
			}
			_dataStored = false;
			_estimationConverged = true;
		}
		return _estimationConverged;
	}

	//return true on final convergence
	return _estimationConverged;
}

void TRecalibrationBQSREstimator::writeAllToFile(std::string filenameTag){
	//write readGroup x Quality table
	writeQualityToFile(filenameTag);
	//write readGroup x position table
	if(storage.considerPosition){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	if(storage.considerPositionReverse){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	if(storage.considerContext){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSREstimator::writeCurrentTmpTable(std::string filenameTag){
	//write readGroup x Quality table
	if(!_qualityConverged) writeQualityToFile(filenameTag);

	//write readGroup x position table
	else if(storage.considerPosition && !_positionConverged){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	else if(storage.considerPositionReverse && !_positionReverseConverged){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	else if(storage.considerContext && !_contextConverged){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSREstimator::writeQualityToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x quality table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int i=0; i<_readGroups->size(); ++i){
		int rg_index = _readGroupMap->getIndex(i);
		for(int q=0; q<qualityIndex->numQ; ++q){
			out << _readGroups->getName(i) << "\t" << qualityIndex->getPhredIntFromIndex(q)
				<< "\tM\t" << _qualityMap.errorToPhred(storage.qualityCells[rg_index][q].curEstimate)
				<< "\t" << storage.qualityCells[rg_index][q].getNumObsForPrinting()
				<< "\t" << storage.qualityCells[rg_index][q].firstDerivativeSave
				<< "\t" << storage.qualityCells[rg_index][q].secondDerivativeSave
				<< "\t" << storage.qualityCells[rg_index][q].F << "\t"
				<< storage.qualityCells[rg_index][q].estimationConverged << "\n";
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSREstimator::writePositionToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x position table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int i=0; i<_readGroups->size(); ++i){
		int rg_index = _readGroupMap->getIndex(i);
		for(int p=0; p<storage.maxPos; ++p){
			out << _readGroups->getName(i) << "\t" << p+1 << "\tM\t" << storage.positionCells[rg_index][p].curEstimate
				<< "\t" << storage.positionCells[rg_index][p].getNumObsForPrinting()
				<< "\t" << storage.positionCells[rg_index][p].firstDerivativeSave
				<< "\t" << storage.positionCells[rg_index][p].secondDerivativeSave
				<< "\t" << storage.positionCells[rg_index][p].F
				<< "\t" << storage.positionCells[rg_index][p].estimationConverged << "\n";

		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSREstimator::writePositionReverseToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Reverse_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x position reverse table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int i=0; i<_readGroups->size(); ++i){
		int rg_index = _readGroupMap->getIndex(i);
		for(int p=0; p<storage.maxPosReverse; ++p){
			out << _readGroups->getName(i) << "\t" << p+1 << "\tM\t" << storage.positionReverseCells[rg_index][p].curEstimate
				<< "\t" << storage.positionReverseCells[rg_index][p].getNumObsForPrinting()
				<< "\t" << storage.positionReverseCells[rg_index][p].firstDerivativeSave
				<< "\t" << storage.positionReverseCells[rg_index][p].secondDerivativeSave
				<< "\t" << storage.positionReverseCells[rg_index][p].F
				<< "\t" << storage.positionReverseCells[rg_index][p].estimationConverged << "\n";
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSREstimator::writeContextToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x context table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tContext\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int r=0; r<_readGroups->size(); ++r){
		int rg_index = _readGroupMap->getIndex(r);
		for(int c=0; c<storage.numContexts; ++c){
			out << _readGroups->getName(r) << "\t" << storage._genoMap.getContextString(c) << "\tM\t" << storage.contextCells[rg_index][c].curEstimate
				<< "\t" << storage.contextCells[rg_index][c].getNumObsForPrinting()
				<< "\t" << storage.contextCells[rg_index][c].firstDerivativeSave
				<< "\t" << storage.contextCells[rg_index][c].secondDerivativeSave
				<< "\t" << storage.contextCells[rg_index][c].F
				<< "\t" << storage.contextCells[rg_index][c].estimationConverged << "\n";
		}
	}
	out.close();
	_logfile->done();
}


void TRecalibrationBQSREstimator::calculateAndPrintLLSurfaceQuality(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x quality and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tQuality\terrorRate\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int q=0; q<qualityIndex->numQ; ++q){
			storage.qualityCells[_readGroupMap->getIndex(r)][q].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + toString(qualityIndex->getPhredIntFromIndex(q)), out);
		}
	}
	out.close();
		_logfile->done();
}

void TRecalibrationBQSREstimator::calculateAndPrintLLSurfacePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tPosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int p=0; p<storage.maxPos; ++p){
			storage.positionCells[_readGroupMap->getIndex(r)][p].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + toString(p+1), out);
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSREstimator::calculateAndPrintLLSurfaceReversePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_ReversePosition_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x reverse position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tReversePosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int p=0; p<storage.maxPosReverse; ++p){
			storage.positionReverseCells[_readGroupMap->getIndex(r)][p].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + toString(p+1), out);
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSREstimator::calculateAndPrintLLSurfaceContext(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x context and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tContext\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int c=0; c<storage.numContexts; ++c){
			storage.contextCells[_readGroupMap->getIndex(r)][c].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + storage._genoMap.getContextString(c), out);
		}
	}
	out.close();
	_logfile->done();
}

bool TRecalibrationBQSREstimator::allConverged(){
	if(!_qualityConverged) return false;
	if(storage.considerPosition && !_positionConverged) return false;
	if(storage.considerPositionReverse && !_positionReverseConverged) return false;
	if(storage.considerContext && !_contextConverged) return false;
	return true;
}

void TRecalibrationBQSREstimator::reopenEstimation(){
	//resets all cells not to have converged
	if(_estimateQuality){
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int q=0; q<qualityIndex->numQ; ++q){
				storage.qualityCells[i][q].reopenEstimation();
			}
		}
		_qualityConverged = false;
	}

	//also for position
	if(storage.considerPosition && _estimatePosition){
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<storage.maxPos; ++p){
				storage.positionCells[i][p].reopenEstimation();
			}
		}
		_positionConverged = false;
	}

	//reverse position
	if(storage.considerPositionReverse && _estimatePositionReverse){
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<storage.maxPosReverse; ++p){
				storage.positionReverseCells[i][p].reopenEstimation();
			}
		}
		_positionReverseConverged = false;
	}

	//and context
	if(storage.considerContext && _estimateContext){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int c=0; c<storage.numContexts; ++c){
				storage.contextCells[r][c].reopenEstimation();
			}
		}
		_contextConverged = false;
	}
}

