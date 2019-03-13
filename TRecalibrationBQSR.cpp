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
TBQSR_cell::TBQSR_cell():TBQSR_cell_base(){
	numMatches = 0;
	pointerToBatch = NULL;
}

void TBQSR_cell::init(int ReadGroup, float initialError){
	TBQSR_cell_base::init(ReadGroup);

	curEstimate = initialError;
	if(curEstimate <= 0.0) curEstimate = 0.000000001;
	if(curEstimate >= 1.0) curEstimate = 0.9;

	//storage
	store = false;
}

void TBQSR_cell::doStoreData(){
	store = true;
	D_storage.push_back(new float[batchSize]);
	batchIt = D_storage.rbegin();
	pointerToBatch = *batchIt;
};

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

void TBQSR_cell::addBase(TBase* base, Base & RefBase, TQualityMap & qualiMap){
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

void TBQSR_cellPosition::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup){
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
	return  BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(qualMap.errorToPhredInt(base->errorRate))].curEstimate;
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

void TBQSR_cellPositionRev::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup){
	TBQSR_cellPosition::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position = gotBQSR_cells_position_readGroup;
	considerPosition = true;
}

void TBQSR_cellPositionRev::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup){
	TBQSR_cellPosition::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position = NULL;
}

float TBQSR_cellPositionRev::getEpsilon(TBase* base, TQualityMap & qualMap){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(qualMap.errorToPhredInt(base->errorRate))].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->distFrom5Prime].curEstimate;
	return  epsilonAlpha;
}

//---------------------------------------------------------------
TBQSR_cellContext::TBQSR_cellContext():TBQSR_cellPositionRev(){
	BQSR_cells_readGroup_position_rev = NULL;
	considerPositionRev = false;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup, gotBQSR_cells_position_readGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_cells_positionRev_readGroup;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPosition** gotBQSR_cells_position_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup, gotBQSR_cells_position_readGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_cells_positionRev_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position_rev = gotBQSR_cells_positionRev_readGroup;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(int ReadGroup, TQualityIndex* QualityIndex, TBQSR_cell** gotBQSR_cells_quality_readGroup){
	TBQSR_cellPositionRev::init(ReadGroup, QualityIndex, gotBQSR_cells_quality_readGroup);
	BQSR_cells_readGroup_position_rev = NULL;
}

float TBQSR_cellContext::getEpsilon(TBase* base, TQualityMap & qualMap){
	float epsilonAlpha = BQSR_cells_readGroup_quality[myReadGroup][qualityIndex->getIndex(qualMap.errorToPhredInt(base->errorRate))].curEstimate;
	if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[myReadGroup][base->distFrom5Prime].curEstimate;
	if(considerPositionRev) epsilonAlpha *= BQSR_cells_readGroup_position_rev[myReadGroup][base->distFrom3Prime].curEstimate;
	return  epsilonAlpha;
}

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(TParameters & params, TLog* Logfile, TReadGroups* ReadGroups):TRecalibration(){
	_logfile = Logfile;
	_readGroups = ReadGroups;
	_numContexts = 20;
	qualityIndex = NULL;
	_maxPos = 0;

	_type = "BQSR";

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
}

void TRecalibrationBQSR::_initializeBQSRReadGroupQualityTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x quality table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//construct for each read group in bam file
	_BQSR_cells_readGroup_quality = new TBQSR_cell*[_readGroups->size()];

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
	for(int i=0; i<_readGroups->size(); ++i){
		_BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q)
			_BQSR_cells_readGroup_quality[i][q].init(i, _qualityMap.qualityToError(qualityIndex->getPhredIntFromIndex(q)));
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0 && _readGroups->readGroupExists(vec[0])){
			//set quality and empirical error rate
			int readGroup = _readGroups->find(vec[0]);

			if(readGroup >= 0){ //returns -1 if read group does not exist
				q = stringToInt(vec[1]);
				double phredEmpiric = stringToDouble(vec[3]);
				_BQSR_cells_readGroup_quality[readGroup][qualityIndex->getIndex(q+33)].set(_qualityMap.phredToError(phredEmpiric), vec[4]);
			} else throw "readGroup " + vec[0] + " does not exist in BAM file header!";
		}
	}

	//done!
	_logfile->done();
	_logfile->conclude("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
}

void TRecalibrationBQSR::_initializeBQSRReadGroupPositionTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x position table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position = new TBQSR_cellPosition*[_readGroups->size()];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	_maxPos = 0;
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
			int p = stringToInt(vec[1]);
			if(p > _maxPos) _maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[_readGroups->size()];
	for(int r=0; r<_readGroups->size(); ++r){
		BQSR_cells_readGroup_position[r] = new TBQSR_cellPosition[_maxPos];
		isListed[r] = new bool[_maxPos];
		for(int p=0; p<_maxPos; ++p){
			BQSR_cells_readGroup_position[r][p].init(r, qualityIndex, _BQSR_cells_readGroup_quality);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
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
			BQSR_cells_readGroup_position[readGroup][p-1].set(alpha, vec[4]);
			isListed[readGroup][p-1] = true;
		}
	}

	//check if we miss positions
	for(int i=0; i<_readGroups->size(); ++i){
		for(int p=0; p<_maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + _readGroups->getName(i) + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//done!
	_considerPosition = true;
	_logfile->done();
	_logfile->conclude("Considering positions up to " + toString(_maxPos));
}

void TRecalibrationBQSR::_initializeBQSRReadGroupPositionReverseTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x position reverse table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position reverse table from file '" + filename + "'!";

	//construct for each read group in bam file
	_BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[_readGroups->size()];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	_maxPos = 0;
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
			int p = stringToInt(vec[1]);
			if(p > _maxPos) _maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[_readGroups->size()];
	for(int r=0; r<_readGroups->size(); ++r){
		_BQSR_cells_readGroup_position_reverse[r] = new TBQSR_cellPositionRev[_maxPos];
		isListed[r] = new bool[_maxPos];
		for(int p=0; p<_maxPos; ++p){
			if(_considerPosition) _BQSR_cells_readGroup_position_reverse[r][p].init(r, qualityIndex, _BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position);
			else _BQSR_cells_readGroup_position_reverse[r][p].init(r, qualityIndex, _BQSR_cells_readGroup_quality);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
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
			_BQSR_cells_readGroup_position_reverse[readGroup][p-1].set(alpha, vec[4]);
			isListed[readGroup][p-1] = true;
		}
	}

	//check if we miss positions
	for(int i=0; i<_readGroups->size(); ++i){
		for(int p=0; p<_maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + _readGroups->getName(i) + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//done!
	_considerPositionReverse = true;
	_logfile->done();
	_logfile->conclude("Considering positions reverse up to " + toString(_maxPos));
}

void TRecalibrationBQSR::_initializeBQSRReadGroupContextTableFromFile(std::string filename){
	_logfile->listFlush("Constructing BQSR readGroup x context table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x context table from file '" + filename + "'!";

	//construct for each read group in bam file
	_BQSR_cells_readGroup_context = new TBQSR_cellContext*[_readGroups->size()];
	for(int r=0; r<_readGroups->size(); ++r){
		_BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[_numContexts];
		for(int c=0; c<_numContexts; ++c){
			if(_considerPosition && _considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(r, qualityIndex, _BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, _BQSR_cells_readGroup_position_reverse);
			else if(_considerPosition && !_considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(r, qualityIndex, _BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position);
			else if(!_considerPosition && _considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(r, qualityIndex, _BQSR_cells_readGroup_quality, _BQSR_cells_readGroup_position_reverse);
			else _BQSR_cells_readGroup_context[r][c].init(r, qualityIndex, _BQSR_cells_readGroup_quality);
		}
	}

	//create object to check of all contexts have been initialized!
	bool** isListed = new bool*[_readGroups->size()];
	for(int i=0; i<_readGroups->size(); ++i){
		isListed[i] = new bool[_numContexts];
		for(int c=0; c<_numContexts; ++c){
			isListed[i][c] = false;
		}
	}

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
			//set quality and empirical error rate

			if(_readGroups->readGroupExists(vec[0])){
				//set quality and empirical error rate
				int readGroup = _readGroups->find(vec[0]);
				int context = _genoMap.getContext(vec[1][0], vec[1][1]);
				double alpha = stringToDouble(vec[3]);
				_BQSR_cells_readGroup_context[readGroup][context].set(alpha, vec[4]);
				isListed[readGroup][context] = true;
			}
		}
	}

	//check if we miss contexts
	for(int i=0; i<_readGroups->size(); ++i){
		for(int c=0; c<_numContexts; ++c){
			if(!isListed[i][c]) throw "Context " + _genoMap.getContextString(c) + " is not listed for read group '" + _readGroups->getName(i) + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//done!
	_considerContext = true;
	_logfile->done();
	_logfile->conclude("Considering context");
}

double TRecalibrationBQSR::getErrorRate(TBase & base){
	double q = _BQSR_cells_readGroup_quality[base.readGroup][qualityIndex->getIndex(_qualityMap.errorToQuality(base.errorRate))].curEstimate;
	if(_considerPosition) q *= BQSR_cells_readGroup_position[base.readGroup][base.distFrom5Prime].curEstimate;
	if(_considerPositionReverse) q *= _BQSR_cells_readGroup_position_reverse[base.readGroup][base.distFrom3Prime].curEstimate;
	if(_considerContext) q *= _BQSR_cells_readGroup_context[base.readGroup][base.context].curEstimate;
	if(q > 1.0) q = 1.0; //make sure the scaling does not lead to errors > 1.0!
	return q;
}

/*
//---------------------------------------------------------------
//RecalibrationBQSREstimator
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(TParameters & params, TLog* Logfile, TReadGroups* ReadGroups, TReadGroupMap* ReadGroupMap):TRecalibration(){
	_logfile = Logfile;
	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMap;
	_estimatetionRequired = false;
	_estimationConverged = false;
	_curNewtonRaphsonLoop = 0;
	_numContexts = 20;
	qualityIndex = NULL;
	_maxPos = 0;

	_type = "BQSR";

	_storeDataInMemory = params.parameterExists("storeInMemory");
	if(_storeDataInMemory) _logfile->list("Will store D in memory to iterate Newton-Raphson faster");
	//if(mergeReadGroupsRecalibration) logfile->list("Pooling all read groups for BQSR recalibration");
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

void TRecalibrationBQSR::_initializeBQSRReadGroupQualityTable(TParameters & params){
	if(params.parameterExists("BQSRQuality")) _initializeBQSRReadGroupQualityTableFromFile(params);
	else {
		_qualityConverged = false;
		_estimateQuality = true;
		_estimatetionRequired = true;
		int minQ = params.getParameterIntWithDefault("minQ", 0);
		int maxQ = params.getParameterIntWithDefault("maxQ", 100);
		_logfile->list("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
		qualityIndex = new TQualityIndex(minQ, maxQ);

		//initialize BQSR table
		_BQSR_cells_readGroup_quality = new TBQSR_cell*[_readGroupMap->getNumReadGroups()];
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			_BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
			for(int q=0; q<qualityIndex->numQ; ++q){
				_BQSR_cells_readGroup_quality[i][q].init(_qualityMap.phredIntToErrorMap[qualityIndex->getPhredIntFromIndex(q)], _storeDataInMemory, i);
			}
		}
	}
}

void TRecalibrationBQSR::_initializeBQSRReadGroupQualityTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRQuality");
	_logfile->listFlush("Constructing BQSR readGroup x quality table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//construct for each read group in bam file
	_BQSR_cells_readGroup_quality = new TBQSR_cell*[_readGroupMap->getOrigNumReadGroups()];

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
	for(int i=0; i<_readGroupMap->getOrigNumReadGroups(); ++i){
		_BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q)
			_BQSR_cells_readGroup_quality[i][q].init(_qualityMap.qualityToError(qualityIndex->getPhredIntFromIndex(q)), _storeDataInMemory, i);
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0 && _readGroups->readGroupExists(vec[0])){
			//set quality and empirical error rate
			int readGroup = _readGroups->find(vec[0]);

			if(readGroup >= 0){ //returns -1 if read group does not exist
				q = stringToInt(vec[1]);
				double phredEmpiric = stringToDouble(vec[3]);
				_BQSR_cells_readGroup_quality[readGroup][qualityIndex->getIndex(q+33)].set(_qualityMap.phredToError(phredEmpiric), vec[4]);
			} else throw "readGroup " + vec[0] + " does not exist in BAM file header!";
		}
	}

	//set that no estimation is not required, unless asked for
	if(params.parameterExists("estimateBQSRQuality")){
		_qualityConverged = false;
		_estimateQuality = true;
	} else {
		_qualityConverged = true;
		_estimateQuality = false;
	}

	//done!
	_logfile->done();
	_logfile->conclude("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
}


void TRecalibrationBQSR::_initializeBQSRReadGroupPositionTable(TParameters & params){
	if(params.parameterExists("BQSRPosition")) _initializeBQSRReadGroupPositionTableFromFile(params);
	else {
		_positionConverged = false;
		_considerPosition = params.parameterExists("estimateBQSRPosition");
		if(_considerPosition){
			_estimatetionRequired = true;
			_estimatePosition = true;
			_maxPos = params.getParameterInt("maxPos");
			if(_maxPos < 1) throw "Max position has to be larger than zero!";
			_logfile->list("Considering positions up to " + toString(_maxPos));
			BQSR_cells_readGroup_position = new TBQSR_cellPosition*[_readGroupMap->getNumReadGroups()];
			for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
				BQSR_cells_readGroup_position[r] = new TBQSR_cellPosition[_maxPos];
				for(int p=0; p<_maxPos; ++p) BQSR_cells_readGroup_position[r][p].init(_BQSR_cells_readGroup_quality, qualityIndex, _storeDataInMemory, r);
			}
		} else {
			BQSR_cells_readGroup_position = NULL;
		}
	}
}

void TRecalibrationBQSR::_initializeBQSRReadGroupPositionTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPosition");
	_logfile->listFlush("Constructing BQSR readGroup x position table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position = new TBQSR_cellPosition*[_readGroupMap->getOrigNumReadGroups()];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	_maxPos = 0;
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
			int p = stringToInt(vec[1]);
			if(p > _maxPos) _maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[_readGroupMap->getOrigNumReadGroups()];
	for(int r=0; r<_readGroupMap->getOrigNumReadGroups(); ++r){
		BQSR_cells_readGroup_position[r] = new TBQSR_cellPosition[_maxPos];
		isListed[r] = new bool[_maxPos];
		for(int p=0; p<_maxPos; ++p){
			BQSR_cells_readGroup_position[r][p].init(_BQSR_cells_readGroup_quality, qualityIndex, _storeDataInMemory, r);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
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
			BQSR_cells_readGroup_position[readGroup][p-1].set(alpha, vec[4]);
			isListed[readGroup][p-1] = true;
		}
	}

	//check if we miss positions
	for(int i=0; i<_readGroups->size(); ++i){
		for(int p=0; p<_maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + _readGroups->getName(i) + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPosition")){
		_positionConverged = false;
		_estimatePosition = true;
	} else {
		_positionConverged = true;
		_estimatePosition = false;
	}
	_considerPosition = true;

	//done!
	_logfile->done();
	_logfile->conclude("Considering positions up to " + toString(_maxPos));
}


//the functions are almost identical to the other position -> put in class!
void TRecalibrationBQSR::_initializeBQSRReadGroupPositionReverseTable(TParameters & params){
	if(params.parameterExists("BQSRPositionReverse")) _initializeBQSRReadGroupPositionReverseTableFromFile(params);
	else {
		_positionReverseConverged = false;
		_considerPositionReverse = params.parameterExists("estimateBQSRPositionReverse");
		if(_considerPositionReverse){
			_estimatePositionReverse = true;
			_estimatetionRequired = true;
			_maxPos = params.getParameterInt("maxPos");
			if(_maxPos < 1) throw "Max position has to be larger than zero!";
			_logfile->list("Considering positions reverse up to " + toString(_maxPos));
			_BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[_readGroupMap->getNumReadGroups()];
			for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
				_BQSR_cells_readGroup_position_reverse[r] = new TBQSR_cellPositionRev[_maxPos];
				for(int p=0; p<_maxPos; ++p){
					if(_considerPosition) _BQSR_cells_readGroup_position_reverse[r][p].init(_BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, _storeDataInMemory, r);
					else _BQSR_cells_readGroup_position_reverse[r][p].init(_BQSR_cells_readGroup_quality, qualityIndex, _storeDataInMemory, r);
				}
			}
		} else {
			_BQSR_cells_readGroup_position_reverse = NULL;
		}
	}
}

void TRecalibrationBQSR::_initializeBQSRReadGroupPositionReverseTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPositionReverse");
	_logfile->listFlush("Constructing BQSR readGroup x position reverse table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position reverse table from file '" + filename + "'!";

	//construct for each read group in bam file
	_BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[_readGroupMap->getOrigNumReadGroups()];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	_maxPos = 0;
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
			int p = stringToInt(vec[1]);
			if(p > _maxPos) _maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[_readGroupMap->getOrigNumReadGroups()];
	for(int r=0; r<_readGroupMap->getOrigNumReadGroups(); ++r){
		_BQSR_cells_readGroup_position_reverse[r] = new TBQSR_cellPositionRev[_maxPos];
		isListed[r] = new bool[_maxPos];
		for(int p=0; p<_maxPos; ++p){
			if(_considerPosition) _BQSR_cells_readGroup_position_reverse[r][p].init(_BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, _storeDataInMemory, r);
			else _BQSR_cells_readGroup_position_reverse[r][p].init(_BQSR_cells_readGroup_quality, qualityIndex, _storeDataInMemory, r);
			isListed[r][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
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
			_BQSR_cells_readGroup_position_reverse[readGroup][p-1].set(alpha, vec[4]);
			isListed[readGroup][p-1] = true;
		}
	}

	//check if we miss positions
	for(int i=0; i<_readGroups->size(); ++i){
		for(int p=0; p<_maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + _readGroups->getName(i) + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPositionReverse")){
		_positionReverseConverged = false;
		_estimatePositionReverse = true;
	} else {
		_positionReverseConverged = true;
		_estimatePositionReverse = false;
	}
	_considerPositionReverse = true;

	//done!
	_logfile->done();
	_logfile->conclude("Considering positions reverse up to " + toString(_maxPos));
}

void TRecalibrationBQSR::_initializeBQSRReadGroupContextTable(TParameters & params){
	if(params.parameterExists("BQSRContext")) _initializeBQSRReadGroupContextTableFromFile(params);
	else {
		_contextConverged = false;
		_considerContext = params.parameterExists("estimateBQSRContext");
		if(_considerContext){
			_estimateContext = true;
			_estimatetionRequired = true;
			_logfile->list("Considering context");
			_BQSR_cells_readGroup_context = new TBQSR_cellContext*[_readGroupMap->getNumReadGroups()];
			for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
				_BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[_numContexts];
				for(int c=0; c<_numContexts; ++c){
					if(_considerPosition && _considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, _BQSR_cells_readGroup_position_reverse, qualityIndex, _storeDataInMemory, r);
					else if(_considerPosition && !_considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, _storeDataInMemory, r);
					else if(!_considerPosition && _considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, _BQSR_cells_readGroup_position_reverse, qualityIndex, _storeDataInMemory, r);
					else _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, qualityIndex, _storeDataInMemory, r);
				}
			}
		} else {
			_BQSR_cells_readGroup_context = NULL;
		}
	}
}

void TRecalibrationBQSR::_initializeBQSRReadGroupContextTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRContext");
	_logfile->listFlush("Constructing BQSR readGroup x context table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x context table from file '" + filename + "'!";

	//construct for each read group in bam file
	_BQSR_cells_readGroup_context = new TBQSR_cellContext*[_readGroupMap->getOrigNumReadGroups()];
	for(int r=0; r<_readGroupMap->getOrigNumReadGroups(); ++r){
		_BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[_numContexts];
		for(int c=0; c<_numContexts; ++c){
			if(_considerPosition && _considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, _BQSR_cells_readGroup_position_reverse, qualityIndex, _storeDataInMemory, r);
			else if(_considerPosition && !_considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, _storeDataInMemory, r);
			else if(!_considerPosition && _considerPositionReverse) _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, _BQSR_cells_readGroup_position_reverse, qualityIndex, _storeDataInMemory, r);
			else _BQSR_cells_readGroup_context[r][c].init(_BQSR_cells_readGroup_quality, qualityIndex, _storeDataInMemory, r);
		}
	}

	//create object to check of all contexts have been initialized!
	bool** isListed = new bool*[_readGroupMap->getOrigNumReadGroups()];
	for(int i=0; i<_readGroupMap->getOrigNumReadGroups(); ++i){
		isListed[i] = new bool[_numContexts];
		for(int c=0; c<_numContexts; ++c){
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

			if(_readGroups->readGroupExists(vec[0])){
				//set quality and empirical error rate
				int readGroup = _readGroups->find(vec[0]);
				context = _genoMap.getContext(vec[1][0], vec[1][1]);
				alpha = stringToDouble(vec[3]);
				_BQSR_cells_readGroup_context[readGroup][context].set(alpha, vec[4]);
				isListed[readGroup][context] = true;
			}
		}
	}

	//check if we miss contexts
	for(int i=0; i<_readGroups->size(); ++i){
		for(int c=0; c<_numContexts; ++c){
			if(!isListed[i][c]) throw "Context " + _genoMap.getContextString(c) + " is not listed for read group '" + _readGroups->getName(i) + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRContext")){
		_contextConverged = false;
		_estimateContext = true;
	} else {
		_contextConverged = true;
		_estimateContext = false;
	}
	_considerContext = true;

	//done!
	_logfile->done();
	_logfile->conclude("Considering context");
}

void TRecalibrationBQSR::addSite(TSite & site, TQualityMap & qualMap){
	if(site.referenceBase != 'N'){
		Base refBase = site.getRefBaseAsEnum();
		if(!_qualityConverged){
			for(TBase* it : site.bases){
				_BQSR_cells_readGroup_quality[_readGroupMap[it->readGroup]][qualityIndex->getIndex(qualMap.errorToPhredInt(it->errorRate))].addBase(it, refBase, qualMap);
			}
		}
		else if(_considerPosition && !_positionConverged){
			for(TBase* it : site.bases){
				if(it->distFrom5Prime >= _maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position[_readGroupMap[it->readGroup]][it->distFrom5Prime].addBase(it, refBase, qualMap);
			}
		}
		else if(_considerPositionReverse && !_positionReverseConverged){
			for(TBase* it : site.bases){
				if(it->distFrom3Prime >= _maxPos) throw "Position of base is > maxPos specified!";
				_BQSR_cells_readGroup_position_reverse[_readGroupMap[it->readGroup]][it->distFrom3Prime].addBase(it, refBase, qualMap);
			}
		} else if(_considerContext && !_contextConverged){
			for(TBase* it : site.bases)
				_BQSR_cells_readGroup_context[_readGroupMap[it->readGroup]][it->context].addBase(it, refBase, qualMap);
		}
	}
}

void TRecalibrationBQSR::recalculateDerivativesFromDataInMemory(){
	if(!_qualityConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int j=0; j<qualityIndex->numQ; ++j){
				_BQSR_cells_readGroup_quality[r][j].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(_considerPosition && !_positionConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int p=0; p<_maxPos; ++p){
				BQSR_cells_readGroup_position[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	}
	else if(_considerPositionReverse && !_positionReverseConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int p=0; p<_maxPos; ++p){
				_BQSR_cells_readGroup_position_reverse[r][p].recalculateDerivativesFromDataInMemory();
			}
		}
	} else if(_considerContext && !_contextConverged){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int c=0; c<_numContexts; ++c){
				_BQSR_cells_readGroup_context[r][c].recalculateDerivativesFromDataInMemory();
			}
		}
	}
}

bool TRecalibrationBQSR::estimateEpsilon(std::string filenameTag){
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
				if(!_BQSR_cells_readGroup_quality[i][j].estimate(_convergenceThreshold_F, _minEpsilonQuality, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(_BQSR_cells_readGroup_quality[i][j].hasData && _BQSR_cells_readGroup_quality[i][j].F > maxF) maxF = _BQSR_cells_readGroup_quality[i][j].F;
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
					_BQSR_cells_readGroup_quality[i][j].empty();
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
						_BQSR_cells_readGroup_quality[i][j].clearStorage();
					}
				}
			}
			_dataStored = false;

			//what's next?
			if(!_considerPosition && !_considerPositionReverse && !_considerContext) _estimationConverged = true;
			else _estimationConverged = false;
			_LLSurfacePrinted = false;
		}
		return _estimationConverged;
	}

	//estimate tmpEpsilon for position, if not yet done
	//-------------------------------------------------------
	if(_considerPosition && !_positionConverged){
		//do we print LL surface? Only print once!
		if(_printLLSurface && !_LLSurfacePrinted){
			calculateAndPrintLLSurfacePosition(filenameTag);
			_LLSurfacePrinted = true;
		}

		//now do estimation
		_logfile->listFlush("Estimating tmpEpsilon for readGroup x position table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<_maxPos; ++p){
				if(!BQSR_cells_readGroup_position[i][p].estimate(_convergenceThreshold_F, _minEpsilonFactors, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(BQSR_cells_readGroup_position[i][p].hasData && BQSR_cells_readGroup_position[i][p].F > maxF) maxF = BQSR_cells_readGroup_position[i][p].F;
			}
		}

		//report
		_logfile->done();
		if(numCellsNotConverged == 0){
			_positionConverged = true;
			_logfile->list("Estimation converged in all cells!");
		} else {
			_positionConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (_readGroupMap->getNumReadGroups() * _maxPos));
			_logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		_logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!_positionConverged){
			//empty all cells
			for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
				for(int p=0; p<_maxPos; ++p){
					BQSR_cells_readGroup_position[i][p].empty();
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
					for(int p=0; p<_maxPos; ++p){
						BQSR_cells_readGroup_position[i][p].clearStorage();
					}
				}
			}
			_dataStored = false;

			//what's next?
			if(!_considerPositionReverse && !_considerContext) _estimationConverged = true;
			else _estimationConverged = false;
			_LLSurfacePrinted = false;
		}
		return _estimationConverged;
	}

	//estimate tmpEpsilon for position reverse, if not yet done
	//-------------------------------------------------------
	if(_considerPositionReverse && !_positionReverseConverged){
		//do we print LL surface? Only print once!
		if(_printLLSurface && !_LLSurfacePrinted){
			calculateAndPrintLLSurfaceReversePosition(filenameTag);
			_LLSurfacePrinted = true;
		}

		//now do estimation
		_logfile->listFlush("Estimating tmpEpsilon for readGroup x position reverse table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<_maxPos; ++p){
				if(!_BQSR_cells_readGroup_position_reverse[i][p].estimate(_convergenceThreshold_F, _minEpsilonFactors, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(_BQSR_cells_readGroup_position_reverse[i][p].hasData && _BQSR_cells_readGroup_position_reverse[i][p].F > maxF) maxF = _BQSR_cells_readGroup_position_reverse[i][p].F;
			}
		}

		//report
		_logfile->done();
		if(numCellsNotConverged == 0){
			_positionReverseConverged = true;
			_logfile->list("Estimation converged in all cells!");
		} else {
			_positionReverseConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (_readGroupMap->getNumReadGroups() * _maxPos));
			_logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		_logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!_positionReverseConverged){
			//empty all cells
			for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
				for(int p=0; p<_maxPos; ++p){
					_BQSR_cells_readGroup_position_reverse[i][p].empty();
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
					for(int p=0; p<_maxPos; ++p){
						_BQSR_cells_readGroup_position_reverse[i][p].clearStorage();
					}
				}
			}
			_dataStored = false;

			//what's next?
			if(!_considerContext) _estimationConverged = true;
			else _estimationConverged = false;
			_LLSurfacePrinted = false;
		}
		return _estimationConverged;
	}

	//estimate tmpEpsilon for context
	//-------------------------------------------------------
	if(_considerContext && !_contextConverged){
		//do we print LL surface? Only print once!
		if(_printLLSurface && !_LLSurfacePrinted){
			calculateAndPrintLLSurfaceContext(filenameTag);
			_LLSurfacePrinted = true;
		}

		//now do estimation
		_logfile->listFlush("Estimating tmpEpsilon for quality x context table ...");
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int c=0; c<_numContexts; ++c){
				if(!_BQSR_cells_readGroup_context[r][c].estimate(_convergenceThreshold_F, _minEpsilonFactors, _minObservations, allowIncreaseInF))
					++numCellsNotConverged;
				if(_BQSR_cells_readGroup_context[r][c].hasData && _BQSR_cells_readGroup_context[r][c].F > maxF) maxF = _BQSR_cells_readGroup_context[r][c].F;
			}
		}

		//report
		_logfile->done();
		if(numCellsNotConverged == 0){
			_contextConverged = true;
			_logfile->list("Estimation converged in all cells!");
		} else {
			_contextConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (_readGroupMap->getNumReadGroups() * _numContexts));
			_logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		_logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!_contextConverged){
			//empty all cells
			for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
				for(int c=0; c<_numContexts; ++c){
					_BQSR_cells_readGroup_context[r][c].empty();
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
					for(int c=0; c<_numContexts; ++c){
						_BQSR_cells_readGroup_context[i][c].clearStorage();
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

void TRecalibrationBQSR::writeAllToFile(std::string filenameTag){
	//write readGroup x Quality table
	writeQualityToFile(filenameTag);
	//write readGroup x position table
	if(_considerPosition){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	if(_considerPositionReverse){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	if(_considerContext){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSR::writeCurrentTmpTable(std::string filenameTag){
	//write readGroup x Quality table
	if(!_qualityConverged) writeQualityToFile(filenameTag);

	//write readGroup x position table
	else if(_considerPosition && !_positionConverged){
		writePositionToFile(filenameTag);
	}

	//write readGroup x position_rev table
	else if(_considerPositionReverse && !_positionReverseConverged){
		writePositionReverseToFile(filenameTag);
	}

	//write readGroup x context table
	else if(_considerContext && !_contextConverged){
		writeContextToFile(filenameTag);
	}
}

void TRecalibrationBQSR::writeQualityToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x quality table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int i=0; i<_readGroups->size(); ++i){
		for(int q=0; q<qualityIndex->numQ; ++q){
			out << _readGroups->getName(i) << "\t" << qualityIndex->getPhredIntFromIndex(q)
				<< "\tM\t" << _qualityMap.errorToPhred(_BQSR_cells_readGroup_quality[_readGroupMap[i]][q].curEstimate)
				<< "\t" << _BQSR_cells_readGroup_quality[_readGroupMap[i]][q].getNumObsForPrinting()
				<< "\t" << _BQSR_cells_readGroup_quality[_readGroupMap[i]][q].firstDerivativeSave
				<< "\t" << _BQSR_cells_readGroup_quality[_readGroupMap[i]][q].secondDerivativeSave
				<< "\t" << _BQSR_cells_readGroup_quality[_readGroupMap[i]][q].F << "\t"
				<< _BQSR_cells_readGroup_quality[_readGroupMap[i]][q].estimationConverged << "\n";
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSR::writePositionToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x position table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int i=0; i<_readGroups->size(); ++i){
		for(int p=0; p<_maxPos; ++p){
			out << _readGroups->getName(i) << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position[_readGroupMap[i]][p].curEstimate
				<< "\t" << BQSR_cells_readGroup_position[_readGroupMap[i]][p].getNumObsForPrinting()
				<< "\t" << BQSR_cells_readGroup_position[_readGroupMap[i]][p].firstDerivativeSave
				<< "\t" << BQSR_cells_readGroup_position[_readGroupMap[i]][p].secondDerivativeSave
				<< "\t" << BQSR_cells_readGroup_position[_readGroupMap[i]][p].F
				<< "\t" << BQSR_cells_readGroup_position[_readGroupMap[i]][p].estimationConverged << "\n";

		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSR::writePositionReverseToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Reverse_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x position reverse table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tPosition\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int i=0; i<_readGroups->size(); ++i){
		for(int p=0; p<_maxPos; ++p){
			out << _readGroups->getName(i) << "\t" << p+1 << "\tM\t" << _BQSR_cells_readGroup_position_reverse[_readGroupMap[i]][p].curEstimate
				<< "\t" << _BQSR_cells_readGroup_position_reverse[_readGroupMap[i]][p].getNumObsForPrinting()
				<< "\t" << _BQSR_cells_readGroup_position_reverse[_readGroupMap[i]][p].firstDerivativeSave
				<< "\t" << _BQSR_cells_readGroup_position_reverse[_readGroupMap[i]][p].secondDerivativeSave
				<< "\t" << _BQSR_cells_readGroup_position_reverse[_readGroupMap[i]][p].F
				<< "\t" << _BQSR_cells_readGroup_position_reverse[_readGroupMap[i]][p].estimationConverged << "\n";
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSR::writeContextToFile(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_Table.txt";
	_logfile->listFlush("Writing BQSR readGroup x context table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tContext\tEventType\tScaling\tLog10Observations";
	out << "\tFirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int c=0; c<_numContexts; ++c){
			out << _readGroups->getName(r) << "\t" << _genoMap.getContextString(c) << "\tM\t" << _BQSR_cells_readGroup_context[_readGroupMap[r]][c].curEstimate
				<< "\t" << _BQSR_cells_readGroup_context[_readGroupMap[r]][c].getNumObsForPrinting()
				<< "\t" << _BQSR_cells_readGroup_context[_readGroupMap[r]][c].firstDerivativeSave
				<< "\t" << _BQSR_cells_readGroup_context[_readGroupMap[r]][c].secondDerivativeSave
				<< "\t" << _BQSR_cells_readGroup_context[_readGroupMap[r]][c].F
				<< "\t" << _BQSR_cells_readGroup_context[_readGroupMap[r]][c].estimationConverged << "\n";
		}
	}
	out.close();
	_logfile->done();
}


void TRecalibrationBQSR::calculateAndPrintLLSurfaceQuality(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x quality and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tQuality\terrorRate\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int q=0; q<qualityIndex->numQ; ++q){
			_BQSR_cells_readGroup_quality[_readGroupMap[r]][q].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + toString(qualityIndex->getPhredIntFromIndex(q)), out);
		}
	}
	out.close();
		_logfile->done();
}

void TRecalibrationBQSR::calculateAndPrintLLSurfacePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tPosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int p=0; p<_maxPos; ++p){
			BQSR_cells_readGroup_position[_readGroupMap[r]][p].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + toString(p+1), out);
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSR::calculateAndPrintLLSurfaceReversePosition(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_ReversePosition_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x reverse position and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tReversePosition\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int p=0; p<_maxPos; ++p){
			_BQSR_cells_readGroup_position_reverse[_readGroupMap[r]][p].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + toString(p+1), out);
		}
	}
	out.close();
	_logfile->done();
}

void TRecalibrationBQSR::calculateAndPrintLLSurfaceContext(std::string & filenameTag){
	std::string filename = filenameTag + "_BQSR_ReadGroup_Context_LLSurface.txt";
	_logfile->listFlush("Calculating LL surface for readGroup x context and writing to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header
	out << "ReadGroup\tContext\talpha\tLL\tFirstDerivative\tSecondDerivative\n";
	for(int r=0; r<_readGroups->size(); ++r){
		for(int c=0; c<_numContexts; ++c){
			_BQSR_cells_readGroup_context[_readGroupMap[r]][c].calcLikelihoodSurface(_numPosLLsurface, _readGroups->getName(r) + "\t" + _genoMap.getContextString(c), out);
		}
	}
	out.close();
	_logfile->done();
}

bool TRecalibrationBQSR::allConverged(){
	if(!_qualityConverged) return false;
	if(_considerPosition && !_positionConverged) return false;
	if(_considerPositionReverse && !_positionReverseConverged) return false;
	if(_considerContext && !_contextConverged) return false;
	return true;
}

void TRecalibrationBQSR::reopenEstimation(){
	//resets all cells not to have converged
	if(_estimateQuality){
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int q=0; q<qualityIndex->numQ; ++q){
				_BQSR_cells_readGroup_quality[i][q].reopenEstimation();
			}
		}
		_qualityConverged = false;
	}

	//also for position
	if(_considerPosition && _estimatePosition){
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<_maxPos; ++p){
				BQSR_cells_readGroup_position[i][p].reopenEstimation();
			}
		}
		_positionConverged = false;
	}

	//reverse position
	if(_considerPositionReverse && _estimatePositionReverse){
		for(int i=0; i<_readGroupMap->getNumReadGroups(); ++i){
			for(int p=0; p<_maxPos; ++p){
				_BQSR_cells_readGroup_position_reverse[i][p].reopenEstimation();
			}
		}
		_positionReverseConverged = false;
	}

	//and context
	if(_considerContext && _estimateContext){
		for(int r=0; r<_readGroupMap->getNumReadGroups(); ++r){
			for(int c=0; c<_numContexts; ++c){
				_BQSR_cells_readGroup_context[r][c].reopenEstimation();
			}
		}
		_contextConverged = false;
	}
}

double TRecalibrationBQSR::getErrorRate(TBase & base){
	double q = _BQSR_cells_readGroup_quality[base.readGroup][qualityIndex->getIndex(_qualityMap.errorToQuality(base.errorRate))].curEstimate;
	if(_considerPosition) q *= BQSR_cells_readGroup_position[base.readGroup][base.distFrom5Prime].curEstimate;
	if(_considerPositionReverse) q *= _BQSR_cells_readGroup_position_reverse[base.readGroup][base.distFrom3Prime].curEstimate;
	if(_considerContext) q *= _BQSR_cells_readGroup_context[base.readGroup][base.context].curEstimate;
	if(q > 1.0) q = 1.0; //make sure the scaling does not lead to errors > 1.0!
	return q;
}

*/
