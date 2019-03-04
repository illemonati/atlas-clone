/*
 * TRecalibrationBQSR.cpp
 *
 *  Created on: Mar 4, 2019
 *      Author: phaentu
 */


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

void TBQSR_cell_base::init(bool Store, int ReadGroup){
	store = Store;
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
TBQSR_cell::TBQSR_cell():TBQSR_cell_base(){
	numMatches = 0;
	pointerToBatch = NULL;
}

void TBQSR_cell::init(float initialError, bool Store, int ReadGroup){
	TBQSR_cell_base::init(Store, ReadGroup);

	curEstimate = initialError;
	if(curEstimate <= 0.0) curEstimate = 0.000000001;
	if(curEstimate >= 1.0) curEstimate = 0.9;

	//storage
	if(store){
		D_storage.push_back(new float[batchSize]);
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
	}
}

void TBQSR_cell::empty(){
	if(!estimationConverged){
		TBQSR_cell_base::empty();
		if(!store) numMatches = 0;
	} else {
		clearStorage();
	}
}

void TBQSR_cell::clearStorage(){
	if(store){
		for(batchIt = D_storage.rbegin(); batchIt != D_storage.rend(); ++batchIt)
			delete[] *batchIt;
		D_storage.clear();
		next = 0;
	}
}

void TBQSR_cell::addBase(TBase* base, Base & RefBase){
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

void TBQSR_cell::addToDerivatives(float & D){
	float oneMinus4D = 1.0 - 4.0 * D;
	if(curEstimate >= 1.0) curEstimate = 0.999999;
	firstDerivative += oneMinus4D / (-4.0*D*curEstimate + 3.0*D + curEstimate);
	float tmpF = oneMinus4D / (D*(3.0-4.0*curEstimate) + curEstimate);
	secondDerivative -= tmpF * tmpF;
}

void TBQSR_cell::addToLL(float & D){
	LL += log((1.0-D)*curEstimate/3.0 + D*(1.0-curEstimate));
}

void TBQSR_cell::recalculateDerivativesFromDataInMemory(){
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

void TBQSR_cell::recalculateLLFromDataInMemory(){
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

void TBQSR_cell::runNewtonRaphsonAndCheck(float & convergenceThreshold, float & minEpsilon, bool & allowIncreaseInF){
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

bool TBQSR_cell::estimate(float & convergenceThreshold, float & minEpsilon, long & minObservations, bool & allowIncreaseInF){
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

void TBQSR_cell::calcLikelihoodSurface(int numPositions, std::string tag, std::ofstream & out){
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

void TBQSR_cellPosition::init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cell_base::init(Store, ReadGroup);
	BQSR_cells_readGroup_quality = gotBQSR_cells_quality_readGroup;
	qualityIndex = QualityIndex;

	//storage
	if(store){
		D_storage.push_back(new BQSRFactorStorage[batchSize]);
		batchIt = D_storage.rbegin();
		pointerToBatch = *batchIt;
	}
}

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

float TBQSR_cellPosition::getEpsilon(TBase* base){
	return  BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(base->quality)].curEstimate;
}

void TBQSR_cellPosition::addBase(TBase* base, Base & RefBase){
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
			pointerToBatch[next].epsilon = getEpsilon(base);
			addToDerivatives(pointerToBatch[next].D, pointerToBatch[next].epsilon);
			++next;
		} else {
			float D = getD(base, RefBase);
			float eps = getEpsilon(base);
			addToDerivatives(D, eps);
		}
		++numObservationsTmp;
	}
}

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
}

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

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position = gotBQSR_quality_position;
	considerPosition = true;
}

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position = NULL;
}

float TBQSR_cellPositionRev::getEpsilon(TBase* base){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(base->quality)].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->posInRead].curEstimate;
	return  epsilonAlpha;
}

//---------------------------------------------------------------
TBQSR_cellContext::TBQSR_cellContext():TBQSR_cellPositionRev(){
	BQSR_cells_readGroup_position_rev = NULL;
	considerPositionRev = false;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex, bool Store, int ReadGroup){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex, Store, ReadGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

float TBQSR_cellContext::getEpsilon(TBase* base){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(base->quality)].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->posInRead].curEstimate;
	if(considerPositionRev) epsilonAlpha *= BQSR_cells_readGroup_position_rev[myReadGroup][base->posInReadRev].curEstimate;
	return  epsilonAlpha;
}

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TLog* Logfile, TReadGroupMap& ReadGroupMap):TRecalibration(ReadGroupMap){
	logfile = Logfile;
	bamHeader = BamHeader;
	estimatetionRequired = false;
	estimationConverged = false;
	curNewtonRaphsonLoop = 0;
	numContexts = 20;
	qualityIndex = NULL;
	maxPos = 0;

	storeDataInMemory = params.parameterExists("storeInMemory");
	if(storeDataInMemory) logfile->list("Will store D in memory to iterate Newton-Raphson faster");
	//if(mergeReadGroupsRecalibration) logfile->list("Pooling all read groups for BQSR recalibration");
	dataStored = false;

	//check if BQSR table readGroup x Quality is given, or has to be estimated
	initializeBQSRReadGroupQualityTable(params);

	//Do we also consider the effect of the position in read (cycle)?
	initializeBQSRReadGroupPositionTable(params);
	initializeBQSRReadGroupPositionReverseTable(params);

	//Do we also consider the context (dinucleotide)?
	initializeBQSRReadGroupContextTable(params);

	//read Newton-Raphson arguments from user
	convergenceThreshold_F = params.getParameterDoubleWithDefault("maxF", 0.0000001);
	minEpsilonQuality = params.getParameterDoubleWithDefault("minEpsQuality", 0.000001);
	minEpsilonFactors = params.getParameterDoubleWithDefault("minEpsFactors", 0.0001);
	numLoopIncreaseFAllowed = params.getParameterIntWithDefault("numLoopIncreaseFAllowed", 3);
	if(estimatetionRequired){
		logfile->startIndent("Conditions to stop Newton-Raphson algorithm:");
		logfile->list("Stopping Newton-Raphson if F < " + toString(convergenceThreshold_F));
		logfile->list("Allow F to increase in first " + toString(numLoopIncreaseFAllowed) + " loops.");
		logfile->list("Stopping Newton-Raphson if the change in quality is < " + toString(minEpsilonQuality));
		logfile->list("Stopping Newton-Raphson if the change in a factor (e.g. position) is < " + toString(minEpsilonFactors));
		logfile->endIndent();
	}

	//get minimal number of observations to conduct estimation
	minObservations = params.getParameterLongWithDefault("minObservations", 32000);

	//do we print LL surfaces?
	numPosLLsurface = params.getParameterIntWithDefault("LLSurface", 0);
	if(numPosLLsurface > 0) printLLSurface = true;
	else printLLSurface = false;
	LLSurfacePrinted = false;
}

void TRecalibrationBQSR::initializeBQSRReadGroupQualityTable(TParameters & params){
	if(params.parameterExists("BQSRQuality")) initializeBQSRReadGroupQualityTableFromFile(params);
	else {
		qualityConverged = false;
		estimateQuality = true;
		estimatetionRequired = true;
		int minQ = params.getParameterIntWithDefault("minQ", 0);
		int maxQ = params.getParameterIntWithDefault("maxQ", 100);
		logfile->list("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
		qualityIndex = new TQualityIndex(minQ, maxQ);

		//initialize BQSR table
		BQSR_cells_readGroup_quality = new TBQSR_cell*[readGroupMapObject.getNumReadGroups()];
		for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
			BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
			for(int q=0; q<qualityIndex->numQ; ++q){
				BQSR_cells_readGroup_quality[i][q].init(qualityMap.phredIntToErrorMap[qualityIndex->getPhredIntFromIndex(q)], storeDataInMemory, i);
			}
		}
	}
}
//[qualityIndex->getIndex(base->quality)]
void TRecalibrationBQSR::initializeBQSRReadGroupQualityTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRQuality");
	logfile->listFlush("Constructing BQSR readGroup x quality table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_quality = new TBQSR_cell*[readGroupMapObject.getOrigNumReadGroups()];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	int minQ = 33;
	int maxQ = 133;
	int q;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get min and max quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			q = stringToInt(vec[1]);
			if(q > maxQ)
				maxQ = q;
			if(q < minQ)
				minQ = q;
		}
	}

	//initialize quality index
	qualityIndex = new TQualityIndex(minQ, maxQ);

	//create corresponding objects
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i){
		BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q)
			BQSR_cells_readGroup_quality[i][q].init(qualityMap.qualityToError(qualityIndex->getPhredIntFromIndex(q)), storeDataInMemory, i);
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double phredEmpiric;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				q = stringToInt(vec[1]);
				phredEmpiric = stringToDouble(vec[3]);
				BQSR_cells_readGroup_quality[readGroup][qualityIndex->getIndex(q+33)].set(qualityMap.phredToError(phredEmpiric), vec[4]);
			} else throw "readGroup " + vec[0] + " does not exist in BAM file header!";
		}
	}

	//set that no estimation is not required, unless asked for
	if(params.parameterExists("estimateBQSRQuality")){
		qualityConverged = false;
		estimateQuality = true;
	} else {
		qualityConverged = true;
		estimateQuality = false;
	}

	//done!
	logfile->done();
	logfile->conclude("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
}


void TRecalibrationBQSR::initializeBQSRReadGroupPositionTable(TParameters & params){
	if(params.parameterExists("BQSRPosition")) initializeBQSRReadGroupPositionTableFromFile(params);
	else {
		positionConverged = false;
		considerPosition = params.parameterExists("estimateBQSRPosition");
		if(considerPosition){
			estimatetionRequired = true;
			estimatePosition = true;
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions up to " + toString(maxPos));
			BQSR_cells_readGroup_position = new TBQSR_cellPosition*[readGroupMapObject.getNumReadGroups()];
			for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
				BQSR_cells_readGroup_position[r] = new TBQSR_cellPosition[maxPos];
				for(int p=0; p<maxPos; ++p) BQSR_cells_readGroup_position[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
			}
		} else {
			BQSR_cells_readGroup_position = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupPositionTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPosition");
	logfile->listFlush("Constructing BQSR readGroup x position table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position = new TBQSR_cellPosition*[readGroupMapObject.getOrigNumReadGroups()];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get max position
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			p = stringToInt(vec[1]);
			if(p > maxPos) maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[readGroupMapObject.getOrigNumReadGroups()];
	for(int r=0; r<readGroupMapObject.getOrigNumReadGroups(); ++r){
		BQSR_cells_readGroup_position[r] = new TBQSR_cellPosition[maxPos];
		isListed[r] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				p = stringToInt(vec[1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_position[readGroup][p-1].set(alpha, vec[4]);
				isListed[readGroup][p-1] = true;
			}
		}
	}

	//check if we miss positions
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPosition")){
		positionConverged = false;
		estimatePosition = true;
	} else {
		positionConverged = true;
		estimatePosition = false;
	}
	considerPosition = true;

	//done!
	logfile->done();
	logfile->conclude("Considering positions up to " + toString(maxPos));
}


//the functions are almost identical to the other position -> put in class!
void TRecalibrationBQSR::initializeBQSRReadGroupPositionReverseTable(TParameters & params){
	if(params.parameterExists("BQSRPositionReverse")) initializeBQSRReadGroupPositionReverseTableFromFile(params);
	else {
		positionReverseConverged = false;
		considerPositionReverse = params.parameterExists("estimateBQSRPositionReverse");
		if(considerPositionReverse){
			estimatePositionReverse = true;
			estimatetionRequired = true;
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions reverse up to " + toString(maxPos));
			BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[readGroupMapObject.getNumReadGroups()];
			for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
				BQSR_cells_readGroup_position_reverse[r] = new TBQSR_cellPositionRev[maxPos];
				for(int p=0; p<maxPos; ++p){
					if(considerPosition) BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
					else BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
				}
			}
		} else {
			BQSR_cells_readGroup_position_reverse = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupPositionReverseTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPositionReverse");
	logfile->listFlush("Constructing BQSR readGroup x position reverse table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position reverse table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[readGroupMapObject.getOrigNumReadGroups()];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get max position
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			p = stringToInt(vec[1]);
			if(p > maxPos) maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[readGroupMapObject.getOrigNumReadGroups()];
	for(int r=0; r<readGroupMapObject.getOrigNumReadGroups(); ++r){
		BQSR_cells_readGroup_position_reverse[r] = new TBQSR_cellPositionRev[maxPos];
		isListed[r] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			if(considerPosition) BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
			else BQSR_cells_readGroup_position_reverse[r][p].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				p = stringToInt(vec[1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_position_reverse[readGroup][p-1].set(alpha, vec[4]);
				isListed[readGroup][p-1] = true;
			}
		}
	}

	//check if we miss positions
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPositionReverse")){
		positionReverseConverged = false;
		estimatePositionReverse = true;
	} else {
		positionReverseConverged = true;
		estimatePositionReverse = false;
	}
	considerPositionReverse = true;

	//done!
	logfile->done();
	logfile->conclude("Considering positions reverse up to " + toString(maxPos));
}

void TRecalibrationBQSR::initializeBQSRReadGroupContextTable(TParameters & params){
	if(params.parameterExists("BQSRContext")) initializeBQSRReadGroupContextTableFromFile(params);
	else {
		contextConverged = false;
		considerContext = params.parameterExists("estimateBQSRContext");
		if(considerContext){
			estimateContext = true;
			estimatetionRequired = true;
			logfile->list("Considering context");
			BQSR_cells_readGroup_context = new TBQSR_cellContext*[readGroupMapObject.getNumReadGroups()];
			for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
				BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
				for(int c=0; c<numContexts; ++c){
					if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
					else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
					else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
					else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
				}
			}
		} else {
			BQSR_cells_readGroup_context = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupContextTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRContext");
	logfile->listFlush("Constructing BQSR readGroup x context table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x context table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_context = new TBQSR_cellContext*[readGroupMapObject.getOrigNumReadGroups()];
	for(int r=0; r<readGroupMapObject.getOrigNumReadGroups(); ++r){
		BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
		for(int c=0; c<numContexts; ++c){
			if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
			else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, storeDataInMemory, r);
			else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex, storeDataInMemory, r);
			else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex, storeDataInMemory, r);
		}
	}

	//create object to check of all contexts have been initialized!
	bool** isListed = new bool*[readGroupMapObject.getOrigNumReadGroups()];
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i){
		isListed[i] = new bool[numContexts];
		for(int c=0; c<numContexts; ++c){
			isListed[i][c] = false;
		}
	}

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header
	int context;
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0], bamHeader->ReadGroups);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				context = genoMap.getContext(vec[1][0], vec[1][1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_context[readGroup][context].set(alpha, vec[4]);
				isListed[readGroup][context] = true;
			}
		}
	}

	//check if we miss contexts
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i, ++it){
		for(int c=0; c<numContexts; ++c){
			if(!isListed[i][c]) throw "Context " + genoMap.getContextString(c) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRContext")){
		contextConverged = false;
		estimateContext = true;
	} else {
		contextConverged = true;
		estimateContext = false;
	}
	considerContext = true;

	//done!
	logfile->done();
	logfile->conclude("Considering context");
}

void TRecalibrationBQSR::addSite(TSite & site){
	if(site.referenceBase != 'N'){
		Base refBase = site.getRefBaseAsEnum();
		if(!qualityConverged){
			for(TBase* it : site.bases){
				BQSR_cells_readGroup_quality[readGroupMapObject[it->readGroup]][qualityIndex->getIndex(it->quality)].addBase(it, refBase);
			}
		}
		else if(considerPosition && !positionConverged){
			for(TBase* it : site.bases){
				if(it->posInRead >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position[readGroupMapObject[it->readGroup]][it->posInRead].addBase(it, refBase);
			}
		}
		else if(considerPositionReverse && !positionReverseConverged){
			for(TBase* it : site.bases){
				if(it->posInReadRev >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position_reverse[readGroupMapObject[it->readGroup]][it->posInReadRev].addBase(it, refBase);
			}
		} else if(considerContext && !contextConverged){
			for(TBase* it : site.bases)
				BQSR_cells_readGroup_context[readGroupMapObject[it->readGroup]][it->context].addBase(it, refBase);
		}
	}
}

void TRecalibrationBQSR::recalculateDerivativesFromDataInMemory(){
	if(!qualityConverged){
		for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
			for(int j=0; j<qualityIndex->numQ; ++j){
				BQSR_cells_readGroup_quality[r][j].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(considerPosition && !positionConverged){
		for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(considerPositionReverse && !positionReverseConverged){
		for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position_reverse[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	} else if(considerContext && !contextConverged){
		for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
			for(int c=0; c<numContexts; ++c){
				BQSR_cells_readGroup_context[r][c].recalculateDerivativesFromDataInMemory();
			}
		}
	}
}

bool TRecalibrationBQSR::estimateEpsilon(std::string filenameTag){
	++curNewtonRaphsonLoop;
	bool allowIncreaseInF = false;
	if(curNewtonRaphsonLoop <= numLoopIncreaseFAllowed) allowIncreaseInF = true;

	//recalc derivatives if data is in memory. Otherwise, derivatives were calculated when data was added.
	if(dataStored) recalculateDerivativesFromDataInMemory();

	//estimate tmpEpsilon, if not yet done
	estimationConverged = true;
	int numCellsNotConverged = 0;
	double maxF = 0.0;

	//readGroup x quality
	//-------------------------------------------------------
	if(!qualityConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfaceQuality(filenameTag);
			LLSurfacePrinted = true;
		}
		//now do estimation
		logfile->listFlush("Estimating tmpEpsilon for readGroup x quality table ...");
		for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
			for(int j=0; j<qualityIndex->numQ; ++j){
				if(!BQSR_cells_readGroup_quality[i][j].estimate(convergenceThreshold_F, minEpsilonQuality, minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(BQSR_cells_readGroup_quality[i][j].hasData && BQSR_cells_readGroup_quality[i][j].F > maxF) maxF = BQSR_cells_readGroup_quality[i][j].F;
			}
		}

		//report
		logfile->done();
		if(numCellsNotConverged == 0){
			qualityConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			qualityConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (readGroupMapObject.getNumReadGroups() * qualityIndex->numQ));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!qualityConverged){
			//empty all cells
			for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
				for(int j=0; j<qualityIndex->numQ; ++j){
					BQSR_cells_readGroup_quality[i][j].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writeQualityToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
					for(int j=0; j<qualityIndex->numQ; ++j){
						BQSR_cells_readGroup_quality[i][j].clearStorage();
					}
				}
			}
			dataStored = false;

			//what's next?
			if(!considerPosition && !considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
			LLSurfacePrinted = false;
		}
		return estimationConverged;
	}

	//estimate tmpEpsilon for position, if not yet done
	//-------------------------------------------------------
	if(considerPosition && !positionConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfacePosition(filenameTag);
			LLSurfacePrinted = true;
		}

		//now do estimation
		logfile->listFlush("Estimating tmpEpsilon for readGroup x position table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position[i][p].estimate(convergenceThreshold_F, minEpsilonFactors, minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(BQSR_cells_readGroup_position[i][p].hasData && BQSR_cells_readGroup_position[i][p].F > maxF) maxF = BQSR_cells_readGroup_position[i][p].F;
			}
		}

		//report
		logfile->done();
		if(numCellsNotConverged == 0){
			positionConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			positionConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (readGroupMapObject.getNumReadGroups() * maxPos));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!positionConverged){
			//empty all cells
			for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
				for(int p=0; p<maxPos; ++p){
					BQSR_cells_readGroup_position[i][p].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writePositionToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
					for(int p=0; p<maxPos; ++p){
						BQSR_cells_readGroup_position[i][p].clearStorage();
					}
				}
			}
			dataStored = false;

			//what's next?
			if(!considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
			LLSurfacePrinted = false;
		}
		return estimationConverged;
	}

	//estimate tmpEpsilon for position reverse, if not yet done
	//-------------------------------------------------------
	if(considerPositionReverse && !positionReverseConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfaceReversePosition(filenameTag);
			LLSurfacePrinted = true;
		}

		//now do estimation
		logfile->listFlush("Estimating tmpEpsilon for readGroup x position reverse table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position_reverse[i][p].estimate(convergenceThreshold_F, minEpsilonFactors, minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(BQSR_cells_readGroup_position_reverse[i][p].hasData && BQSR_cells_readGroup_position_reverse[i][p].F > maxF) maxF = BQSR_cells_readGroup_position_reverse[i][p].F;
			}
		}

		//report
		logfile->done();
		if(numCellsNotConverged == 0){
			positionReverseConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			positionReverseConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (readGroupMapObject.getNumReadGroups() * maxPos));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!positionReverseConverged){
			//empty all cells
			for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
				for(int p=0; p<maxPos; ++p){
					BQSR_cells_readGroup_position_reverse[i][p].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writePositionReverseToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
					for(int p=0; p<maxPos; ++p){
						BQSR_cells_readGroup_position_reverse[i][p].clearStorage();
					}
				}
			}
			dataStored = false;

			//what's next?
			if(!considerContext) estimationConverged = true;
			else estimationConverged = false;
			LLSurfacePrinted = false;
		}
		return estimationConverged;
	}

	//estimate tmpEpsilon for context
	//-------------------------------------------------------
	if(considerContext && !contextConverged){
		//do we print LL surface? Only print once!
		if(printLLSurface && !LLSurfacePrinted){
			calculateAndPrintLLSurfaceContext(filenameTag);
			LLSurfacePrinted = true;
		}

		//now do estimation
		logfile->listFlush("Estimating tmpEpsilon for quality x context table ...");
		for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
			for(int c=0; c<numContexts; ++c){
				if(!BQSR_cells_readGroup_context[r][c].estimate(convergenceThreshold_F, minEpsilonFactors, minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(BQSR_cells_readGroup_context[r][c].hasData && BQSR_cells_readGroup_context[r][c].F > maxF) maxF = BQSR_cells_readGroup_context[r][c].F;
			}
		}

		//report
		logfile->done();
		if(numCellsNotConverged == 0){
			contextConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			contextConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (readGroupMapObject.getNumReadGroups() * numContexts));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!contextConverged){
			//empty all cells
			for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
				for(int c=0; c<numContexts; ++c){
					BQSR_cells_readGroup_context[r][c].empty();
				}
			}
			estimationConverged = false;

			//does data need to be added again? Not if stored!
			if(storeDataInMemory) dataStored = true;
		} else {
			//write to file
			writeContextToFile(filenameTag);

			//empty storage
			if(storeDataInMemory){
				for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
					for(int c=0; c<numContexts; ++c){
						BQSR_cells_readGroup_context[i][c].clearStorage();
					}
				}
			}
			dataStored = false;
			estimationConverged = true;
		}
		return estimationConverged;
	}

	//return true on final convergence
	return estimationConverged;
}

void TRecalibrationBQSR::writeAllToFile(std::string filenameTag){
	//write readGroup x Quality table
	writeQualityToFile(filenameTag);
	//write readGroup x position table
	if(considerPosition){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	if(considerPositionReverse){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	if(considerContext){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSR::writeCurrentTmpTable(std::string filenameTag){
	//write readGroup x Quality table
	if(!qualityConverged) writeQualityToFile(filenameTag);

	//write readGroup x position table
	else if(considerPosition && !positionConverged){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	else if(considerPositionReverse && !positionReverseConverged){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	else if(considerContext && !contextConverged){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSR::writeQualityToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x quality table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i, ++it){
		for(int q=0; q<qualityIndex->numQ; ++q){
			out << it->ID << "\t" << qualityIndex->getPhredIntFromIndex(q) << "\tM\t" << qualityMap.errorToPhred(BQSR_cells_readGroup_quality[readGroupMapObject[i]][q].curEstimate) << "\t" << BQSR_cells_readGroup_quality[readGroupMapObject[i]][q].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_quality[readGroupMapObject[i]][q].firstDerivativeSave << "\t" << BQSR_cells_readGroup_quality[readGroupMapObject[i]][q].secondDerivativeSave << "\t" << BQSR_cells_readGroup_quality[readGroupMapObject[i]][q].F << "\t" << BQSR_cells_readGroup_quality[readGroupMapObject[i]][q].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->done();
}

void TRecalibrationBQSR::writePositionToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x position table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position[readGroupMapObject[i]][p].curEstimate << "\t" << BQSR_cells_readGroup_position[readGroupMapObject[i]][p].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_position[readGroupMapObject[i]][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position[readGroupMapObject[i]][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position[readGroupMapObject[i]][p].F << "\t" << BQSR_cells_readGroup_position[readGroupMapObject[i]][p].estimationConverged;
			out << "\n";

		}
	}
	out.close();
	logfile->done();
}

void TRecalibrationBQSR::writePositionReverseToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Reverse_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x position reverse table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<readGroupMapObject.getOrigNumReadGroups(); ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position_reverse[readGroupMapObject[i]][p].curEstimate << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMapObject[i]][p].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMapObject[i]][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMapObject[i]][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMapObject[i]][p].F << "\t" << BQSR_cells_readGroup_position_reverse[readGroupMapObject[i]][p].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->done();
}

void TRecalibrationBQSR::writeContextToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x context table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tContext\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<readGroupMapObject.getOrigNumReadGroups(); ++r, ++it){
		for(int c=0; c<numContexts; ++c){
			out << it->ID << "\t" << genoMap.getContextString(c) << "\tM\t" << BQSR_cells_readGroup_context[readGroupMapObject[r]][c].curEstimate << "\t" << BQSR_cells_readGroup_context[readGroupMapObject[r]][c].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_context[readGroupMapObject[r]][c].firstDerivativeSave << "\t" << BQSR_cells_readGroup_context[readGroupMapObject[r]][c].secondDerivativeSave << "\t" << BQSR_cells_readGroup_context[readGroupMapObject[r]][c].F << "\t" << BQSR_cells_readGroup_context[readGroupMapObject[r]][c].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->done();
}


void TRecalibrationBQSR::calculateAndPrintLLSurfaceQuality(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x quality and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tQuality\terrorRate\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r, ++it){
		for(int q=0; q<qualityIndex->numQ; ++q){
			BQSR_cells_readGroup_quality[r][q].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + toString(qualityIndex->getPhredIntFromIndex(q)), out);
		}
	}
	out.close();
		logfile->done();
}

void TRecalibrationBQSR::calculateAndPrintLLSurfacePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tPosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r, ++it){
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position[r][p].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + toString(p+1), out);
		}
	}
	out.close();
	logfile->done();
}

void TRecalibrationBQSR::calculateAndPrintLLSurfaceReversePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_ReversePosition_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x reverse position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tReversePosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r, ++it){
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position_reverse[r][p].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + toString(p+1), out);
		}
	}
	out.close();
	logfile->done();
}

void TRecalibrationBQSR::calculateAndPrintLLSurfaceContext(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_LLSurface.txt";
	logfile->listFlush("Calculating LL surface for readGroup x context and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tContext\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r, ++it){
		for(int c=0; c<numContexts; ++c){
			BQSR_cells_readGroup_context[r][c].calcLikelihoodSurface(numPosLLsurface, it->ID + "\t" + genoMap.getContextString(c), out);
		}
	}
	out.close();
	logfile->done();
}

bool TRecalibrationBQSR::allConverged(){
	if(!qualityConverged) return false;
	if(considerPosition && !positionConverged) return false;
	if(considerPositionReverse && !positionReverseConverged) return false;
	if(considerContext && !contextConverged) return false;
	return true;
}

void TRecalibrationBQSR::reopenEstimation(){
	//resets all cells not to have converged
	if(estimateQuality){
		for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
			for(int q=0; q<qualityIndex->numQ; ++q){
				BQSR_cells_readGroup_quality[i][q].reopenEstimation();
			}
		}
		qualityConverged = false;
	}

	//also for position
	if(considerPosition && estimatePosition){
		for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position[i][p].reopenEstimation();
			}
		}
		positionConverged = false;
	}

	//reverse position
	if(considerPositionReverse && estimatePositionReverse){
		for(int i=0; i<readGroupMapObject.getNumReadGroups(); ++i){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position_reverse[i][p].reopenEstimation();
			}
		}
		positionReverseConverged = false;
	}

	//and context
	if(considerContext && estimateContext){
		for(int r=0; r<readGroupMapObject.getNumReadGroups(); ++r){
			for(int c=0; c<numContexts; ++c){
				BQSR_cells_readGroup_context[r][c].reopenEstimation();
			}
		}
		contextConverged = false;
	}
}
