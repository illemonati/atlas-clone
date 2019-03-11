/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEM.h"

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
TRecalibrationEMDataTable::TRecalibrationEMDataTable(int NumReadGroups, int MaxQual){
	numReadGroups = NumReadGroups;
	maxQual = MaxQual;

	qualities = new bool**[numReadGroups];
	maxPos = new unsigned int*[numReadGroups];
	countsPerReadGroup = new unsigned int*[numReadGroups];
	for(int rg = 0; rg<numReadGroups; ++rg){
		qualities[rg] = new bool*[2];
		qualities[rg][0] = new bool[maxQual];
		qualities[rg][1] = new bool[maxQual];

		countsPerReadGroup[rg] = new unsigned int[2];
		maxPos[rg] = new unsigned int[2];
	}

	clear();
};

TRecalibrationEMDataTable::~TRecalibrationEMDataTable(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		delete[] qualities[rg][0];
		delete[] qualities[rg][1];
		delete[] qualities[rg];

		delete[] countsPerReadGroup[rg];
		delete[] maxPos[rg];
	}

	delete[] qualities;
	delete[] maxPos;
	delete[] countsPerReadGroup;
};

void TRecalibrationEMDataTable::clear(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		for(int i=0; i<maxQual; ++i){
			qualities[rg][0][i] = 0;
			qualities[rg][1][i] = 0;
		}
		maxPos[rg][0] = 0;
		maxPos[rg][1] = 0;
		countsPerReadGroup[rg][0] = 0;
		countsPerReadGroup[rg][1] = 0;
	}
	totalCounts = 0;
};

void TRecalibrationEMDataTable::add(TRecalibrationEMReadData & data){
	++qualities[data.readGroup][(int) data.isSecond][data.quality];
	if(maxPos[data.readGroup][data.isSecond] < data.position)
		maxPos[data.readGroup][data.isSecond] = data.position;
};

void TRecalibrationEMDataTable::assembleCountsPerReadGroup(){
	totalCounts = 0;
	for(int rg = 0; rg<numReadGroups; ++rg){
		countsPerReadGroup[rg][0] = 0;
		countsPerReadGroup[rg][1] = 0;
		for(int i=0; i<maxQual; ++i){
			countsPerReadGroup[rg][0] += qualities[rg][0][i];
			countsPerReadGroup[rg][1] += qualities[rg][1][i];
		}
		totalCounts += countsPerReadGroup[rg][0] + countsPerReadGroup[rg][1];
	}
};

//--------------------------------------------------------------------
// TRecalibrationEMReadGroupIndex
//--------------------------------------------------------------------
TRecalibrationEMReadGroupIndex::TRecalibrationEMReadGroupIndex(){
	numIndexes = 0;
	numReadGroups = 0;
	initialized = false;
	readGroupInUse = NULL;
	readGroupIndex = NULL;
};

TRecalibrationEMReadGroupIndex::~TRecalibrationEMReadGroupIndex(){
	_free();
}

void TRecalibrationEMReadGroupIndex::_init(int NumReadGroups){
	numReadGroups = NumReadGroups;
	numIndexes = 0;
	readGroupInUse = new bool*[numReadGroups];
	readGroupIndex = new int[numReadGroups];

	for(int rg = 0; rg<numReadGroups; ++rg){
		readGroupInUse[rg] = new bool[2];
		readGroupInUse[rg][0] = false;
		readGroupInUse[rg][1] = false;

		readGroupIndex[rg] = new int[2];
		readGroupIndex[rg][0] = -1;
		readGroupIndex[rg][1] = -1;
	}
};

void TRecalibrationEMReadGroupIndex::initialize(int NumReadGroups){
	_init(NumReadGroups);
};

void TRecalibrationEMReadGroupIndex::initialize(TRecalibrationEMDataTable & dataTable){
	_init(dataTable.numReadGroups);

	//no go through all counts to see if enough data is available.
	dataTable.assembleCountsPerReadGroup();

	for(int mate=0; mate<2; ++mate){
		for(int rg = 0; rg<numReadGroups; ++rg){
			if(dataTable.countsPerReadGroup[rg][mate] > 0){
				readGroupInUse[rg][mate] = true;
				readGroupIndex[rg][mate] = numIndexes;
				++numIndexes;
			} else {
				readGroupInUse[rg][mate] = false;
				readGroupIndex[rg][mate] = -1;
			}
		}
	}
};

void TRecalibrationEMReadGroupIndex::_free(){
	if(initialized){
		for(int rg = 0; rg<numReadGroups; ++rg){
			delete[] readGroupInUse[rg];
			delete[] readGroupIndex[rg];
		}

		delete[] readGroupInUse;
		delete[] readGroupIndex;
	}
	initialized = false;
};

void TRecalibrationEMReadGroupIndex::setAllAsUsed(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		readGroupInUse[rg][0] = true;
		readGroupInUse[rg][1] = true;
		readGroupIndex[rg][0] = 2*rg;
		readGroupIndex[rg][1] = 2*rg + 1;
	}
	numIndexes = 2*numReadGroups;
};

void TRecalibrationEMReadGroupIndex::setAllToSingleIndex(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		readGroupInUse[rg][0] = true;
		readGroupInUse[rg][1] = true;
		readGroupIndex[rg][0] = 0;
		readGroupIndex[rg][1] = 0;
	}
	numIndexes = 1;
};

int TRecalibrationEMReadGroupIndex::setAsUsed(int readGroup, bool isSecondMate){
	readGroupInUse[readGroup][isSecondMate] = true;
	readGroupIndex[readGroup][isSecondMate] = numIndexes;
	++numIndexes;
	return readGroupIndex[readGroup][isSecondMate];
};

bool TRecalibrationEMReadGroupIndex::nextNotInUse(std::pair<int, bool> & pair){
	//check if there are read groups not in use. If so, return true and fill pair. Otherwise, return false
	for(int rg = 0; rg<numReadGroups; ++rg){
		for(int j=0; j<2; ++j){
			if(!readGroupInUse[rg][j]){
				pair.first = rg;
				pair.second = j;
				return true;
			}
		}
	}
	return false;
};

//---------------------------------------------------------------
//TRecalibrationEMEstimationParameters
//---------------------------------------------------------------

TRecalibrationEMEstimationParameters::TRecalibrationEMEstimationParameters(TParameters & args, TLog* logfile){
	logfile->startIndent("Settings regarding the EM algorithm:");
	minRequiredObservationsPerReadGroup = args.getParameterIntWithDefault("minOberservations", 100);
	logfile->list("Will only consider read groups with at least " + toString(minRequiredObservationsPerReadGroup) + " observations.");
	numEMIterations = args.getParameterIntWithDefault("iterations", 100);
	logfile->list("Will perform at max " + toString(numEMIterations) + " EM iterations.");
	maxEpsilon = args.getParameterDoubleWithDefault("maxEps", 0.000001);
	logfile->list("Will stop EM when deltaLL < " + toString(maxEpsilon));
	NewtonRaphsonNumIterations = args.getParameterIntWithDefault("NRiterations", 10);
	logfile->list("Will conduct at max " + toString(NewtonRaphsonNumIterations) + " Newton-Raphson iterations");
	NewtonRaphsonMaxF = args.getParameterDoubleWithDefault("maxF", 0.0001);
	logfile->list("Will stop Newton-Raphson when F < " + toString(NewtonRaphsonMaxF));
	equalBaseFrequencies = args.parameterExists("equalBaseFreq");
	if(equalBaseFrequencies) logfile->list("Will assume equal base frequencies {0.25, 0.25, 0.25, 0.25}");
	logfile->endIndent();
};


//---------------------------------------------------------------
//RecalibrationEMSite
//---------------------------------------------------------------
TRecalibrationEMSite::TRecalibrationEMSite(){
	numReads = 0;
	data = NULL;
};

TRecalibrationEMSite::~TRecalibrationEMSite(){
	if(numReads > 0)
		delete[] data;
};

TRecalibrationEMSite::TRecalibrationEMSite(TSite & site, TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap){
	numReads = site.bases.size();
	data = new TRecalibrationEMReadData[numReads];
	int k = 0;
	for(TBase* it : site.bases){
		//read group. Note: also encodes whether it is first or second mate
		data[k].readGroup = ReadGroupMap[it->readGroup]; // + ReadGroupMap.numReadGroups * it->isSecondMate;

		//quality
		data[k].quality = qualiMap.errorToQuality(it->errorRate);

		//position
		data[k].position = it->posInRead;

		//context
		data[k].context = it->context;

		//isSecond
		data[k].isSecond = it->isSecondMate;

		//now also store D: D[ref][read]
		data[k].setD(it->getBaseAsEnum(), it->PMD_CT, it->PMD_GA);

		++k;
	}
};

void TRecalibrationEMSite::addToDataTable(TRecalibrationEMDataTable & dataTable){
	for(unsigned int k=0; k<numReads; ++k)
		dataTable.add(data[k]);
};

void TRecalibrationEMSite::calcEpsilon(TRecalibrationEMModel_Base* & model, float* & epsilon){
	//calc tmpEpsilon using parameter estimates provided
	for(unsigned int k=0; k<numReads; ++k)
		epsilon[k] = model->calcEpsilon(data[k]);
}

double TRecalibrationEMSite::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel_Base* & model, float* & freqs, float* & epsilon){
	calcEpsilon(model, epsilon);

	//over all genotypes
	double P_g_given_d_theta_denominator = 0.0;
	double tmp;
	float B;
	for(int g=0; g<4; ++g){
		tmp = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			tmp *= B * epsilon[k] - data[k].D[g] + 1.0;
		}
		P_g_given_d_oldBeta[g] = tmp * freqs[g];
		P_g_given_d_theta_denominator += P_g_given_d_oldBeta[g];
	}

	if(P_g_given_d_theta_denominator < 1.0E-25){
		//do again but in log
		double max = 0.0;
		for(int g=0; g<4; ++g){
			tmp = 0.0;
			//loop over all reads
			for(unsigned int k=0; k<numReads; ++k){
				B = 4.0 / 3.0 * data[k].D[g] - 1.0;
				tmp += log(B * epsilon[k] - data[k].D[g] + 1.0);
//				tmp += log(B[g][k] * epsilon[k] - D[g][k] + 1.0);
			}
			P_g_given_d_oldBeta[g] = tmp + log(freqs[g]);
			if(g==0) max = P_g_given_d_oldBeta[g];
			else if(P_g_given_d_oldBeta[g] > max) max = P_g_given_d_oldBeta[g];
		}

		//rescale and delog
		P_g_given_d_theta_denominator = 0.0;
		for(int g=0; g<4; ++g){
			P_g_given_d_oldBeta[g] = exp(P_g_given_d_oldBeta[g] - max);
			P_g_given_d_theta_denominator += P_g_given_d_oldBeta[g];
		}
	}

	//calculate P(g|d, theta)
	for(int g=0; g<4; ++g){
		P_g_given_d_oldBeta[g] = P_g_given_d_oldBeta[g] / P_g_given_d_theta_denominator;
	}

	//return LL = P_g_given_d_theta_denominator
	return log(P_g_given_d_theta_denominator);
}

double TRecalibrationEMSite::calcLL(TRecalibrationEMModel_Base* & model, float* & freqs, float* & epsilon){
	calcEpsilon(model, epsilon);

	//over all genotypes
	double LL = 0.0;
	double tmp;
	float B;
	for(int g=0; g<4; ++g){
		tmp = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			tmp *= B * epsilon[k] - data[k].D[g] + 1.0;
//			tmp *= B[g][k] * epsilon[k] - D[g][k] + 1.0;
		}
		LL += tmp * freqs[g];
	}

	//return LL = P_g_given_d_theta_denominator
	return log(LL);
}

double TRecalibrationEMSite::calcQ(TRecalibrationEMModel_Base* & model, float* & epsilon){
	calcEpsilon(model, epsilon);

	//now calculate P(d, g, new params)
	double P_d_given_g_beta;
	double Q = 0.0;
	float B;

	for(int g=0; g<4; ++g){
		P_d_given_g_beta = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			P_d_given_g_beta *= B * epsilon[k] - data[k].D[g] + 1;
//			P_d_given_g_beta *= B[g][k] * epsilon[k] - D[g][k] + 1;
		}

		if(P_d_given_g_beta < 1.0E-50) P_d_given_g_beta = 1.0E-50;
		Q += P_g_given_d_oldBeta[g] * log(P_d_given_g_beta);
	}

	return Q;
}

void TRecalibrationEMSite::addToJacobianAndF(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMModel_Base* & model, float* & epsilon){
	//calculate tmpEpsilon with current parameters
	calcEpsilon(model, epsilon);

	//tmp variables
	double* weights = new double[numReads];
	double* eps1MinusEps = new double[numReads];
	double* oneMinus2Eps = new double[numReads];
	double* weightJacobian = new double[numReads];
	for(unsigned int k=0; k<numReads; ++k){
		eps1MinusEps[k] = epsilon[k] * (1.0 - epsilon[k]);
		oneMinus2Eps[k] = 1.0 - 2.0 * epsilon[k];
	}

	//fill F and Jacobian
	for(int g=0; g<4; ++g){
		//calc weights
		//------------
		for(unsigned int k=0; k<numReads; ++k){
			double B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			weights[k] = B / (1.0 - data[k].D[g] + B * epsilon[k]) * eps1MinusEps[k];
			weightJacobian[k] = P_g_given_d_oldBeta[g] * weights[k] * (oneMinus2Eps[k] - weights[k]);
		}

		model->addToFandJacobian(F, Jacobian, data, numReads, weights, weightJacobian, P_g_given_d_oldBeta[g]);
	}

	//delete tmp variables
	delete[] weights;
	delete[] eps1MinusEps;
	delete[] oneMinus2Eps;
	delete[] weightJacobian;
}

//---------------------------------------------------------------
//TRecalibrationEMWindow
//---------------------------------------------------------------
TRecalibrationEMWindow::TRecalibrationEMWindow(TBaseFrequencies* baseFreqs, TReadGroupMap* ReadGroupMap){
	freqs = new float[4];
	for(int i=0; i<4; ++i) freqs[i] = (*baseFreqs)[i];
	readGroupMapObject = ReadGroupMap;
}

unsigned int TRecalibrationEMWindow::getMaxDepth(){
	unsigned int maxDepth = 0;
	for(TRecalibrationEMSite* site : sites){
		if(maxDepth < site->numReads)
			maxDepth = site->numReads;
	}
	return maxDepth;
};

void TRecalibrationEMWindow::addSite(TSite & site, TQualityMap & qualiMap){
	if(site.hasData)
		sites.push_back(new TRecalibrationEMSite(site, *readGroupMapObject, qualiMap));
}

long TRecalibrationEMWindow::numSites(){
	return sites.size();
}

long TRecalibrationEMWindow::numSitesDepthTwoOrMore(){
	long _numSites = 0;
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		if((*site)->numReads > 1)
		++_numSites;
	}
	return _numSites;
}

void TRecalibrationEMWindow::addToDataTable(TRecalibrationEMDataTable & dataTable){
	for(TRecalibrationEMSite* site : sites)
		site->addToDataTable(dataTable);
}

long TRecalibrationEMWindow::cumulativeDepth(){
	long cumulDepth = 0;
	for(TRecalibrationEMSite* site : sites){
		cumulDepth += site->numReads;
	}
	return cumulDepth;
}

double TRecalibrationEMWindow::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModel_Base* & model, float* & tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites){
		LL += site->fill_P_g_given_d_beta_AND_calcLL(model, freqs, tmpEpsilon);
	}
	return LL;
}

double TRecalibrationEMWindow::calcLL(TRecalibrationEMModel_Base* & model, float* & tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites)
		LL += site->calcLL(model, freqs, tmpEpsilon);
	return LL;
}

double TRecalibrationEMWindow::calcQ(TRecalibrationEMModel_Base* & model, float* & tmpEpsilon){
	double Q = 0.0;
	for(TRecalibrationEMSite* site : sites)
		Q += site->calcQ(model, tmpEpsilon);
	return Q;
}

void TRecalibrationEMWindow::addToJacobianAndF(arma::vec & F, arma::mat & Jacobian, TRecalibrationEMModel_Base* & model, float* & tmpEpsilon){
	for(TRecalibrationEMSite* site : sites)
		site->addToJacobianAndF(F, Jacobian, model, tmpEpsilon);
}

void TRecalibrationEMWindow::setEuqalBaseFrequencies(){
	for(int i=0; i<4; ++i) freqs[i] = 0.25;
}

//---------------------------------------------------------------
//TRecalibrationEM
//---------------------------------------------------------------
TRecalibrationEM::TRecalibrationEM(BamTools::SamHeader* BamHeader, TLog* Logfile, TReadGroupMap& ReadGroupMap):TRecalibration(ReadGroupMap){
	logfile = Logfile;
	tmpEpsilon = NULL;
	tmpEpsilonInitialized = false;

	//read groups
	bamHeader = BamHeader;
	readGroupNames = new std::string[readGroupMapObject.origNumReadGroups];
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<readGroupMapObject.origNumReadGroups; ++r, ++it){
		readGroupNames[r] = it->ID;
	}

	maxDepth = -1;
	models = NULL;
	modelsInitialized = false;

	estimationParametersInitialized = false;
	EMParameters = NULL;
};

void TRecalibrationEM::_addModel(std::string modelTag, std::vector<std::string> & values, bool verbose){
	trimString(modelTag);

	//get shift of next model
	int shift = 0;
	if(models.size() > 0)
		shift = (*models.rbegin())->getNextShift();

	//add model according to tag
	if(modelTag == "full"){
		if(verbose) logfile->list("Will use full model with quality, quality squared, position, position squared and 20 context specific intercepts.");
		models.push_back(new TRecalibrationEMModel_qualFuncPosFunc(shift, values));
	} else if(modelTag == "noContext"){
		if(verbose) logfile->list("Will use simplified model with only quality, quality squared, position, position squared and one intercept.");
		models.push_back(new TRecalibrationEMModel_qualFuncPosFuncNoContext(shift, values));
	} else if(modelTag == "positionSpecific"){
		if(verbose) logfile->list("Will use a model with quality, quality squared, 20 context and each position.");
		models.push_back(new TRecalibrationEMModel_qualFuncPosSpecific(shift, values));
	} else if(modelTag == "noRecal"){
		if(verbose) logfile->list("Will use a model that does not recalibrate.");
		models.push_back(new TRecalibrationEMModel_noRecal(shift));
	} else throw "Unknown recalibration model '" + modelTag + "'!";
};

void TRecalibrationEM::_addModelForEstimation(std::string modelTag, int maxPos){
	trimString(modelTag);

	//get shift of next model
	int shift = 0;
	if(models.size() > 0)
		shift = (*models.rbegin())->getNextShift();

	//add model according to tag
	if(modelTag == "full"){
		models.push_back(new TRecalibrationEMModel_qualFuncPosFunc(shift));
	} else if(modelTag == "noContext"){
		models.push_back(new TRecalibrationEMModel_qualFuncPosFuncNoContext(shift));
	} else if(modelTag == "positionSpecific"){
		models.push_back(new TRecalibrationEMModel_qualFuncPosSpecific(shift, maxPos));
	} else if(modelTag == "noRecal"){
		models.push_back(new TRecalibrationEMModel_noRecal(shift));
	} else throw "Unknown recalibration model '" + modelTag + "'!";
};

void TRecalibrationEM::initialize(std::string string){
	//initialize from string or file
	int pos = string.find_first_of('[');
	if(pos == std::string::npos)
		initializeRecalibrationParametersFromFile(string);
	else
		initializeRecalibrationParametersFromString(string);
};

void TRecalibrationEM::initializeRecalibrationParametersFromString(std::string & string){
	//string has format model[param1, param2, param3, ...]
	logfile->startIndent("Initializing recal with string '" + string + "' for all read groups:");

	//fill internal read group index
	readGroupIndex.initialize(readGroupMapObject.getNumReadGroups());
	readGroupIndex.setAllToSingleIndex();

	//read model tag
	size_t pos = string.find_first_of('[');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing '['!\nEither provide a valid file name or a string of format 'modelTag[Parameter1, Parameter2, ...]'.";
	std::string modelTag = string.substr(0,pos);
	string.erase(0, pos+1);

	//read parameters
	pos = string.find_first_of(']');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing ']'!\nEither provide a valid file name or a string of format 'modelTag[Parameter1, Parameter2, ...]'.";
	std::vector<std::string> tmpVec, vec;
	fillVectorFromString(string.substr(0, pos), tmpVec, ",");
	repeatIndexes(tmpVec, vec);

	//initialize model
	_addModel(modelTag, vec, true);
	modelsInitialized = true;

	logfile->endIndent();
};

void TRecalibrationEM::initializeRecalibrationParametersFromFile(std::string filename){
	//read parameters from file
	logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "' for reading!";

	//read model on first line
	std::string tmp;
	std::getline(file, tmp);
	size_t pos = tmp.find_first_of('=');
	if(pos == std::string::npos)
		throw "Unable to read recal file: model not provided on first line!";
	std::string modelTag = tmp.substr(pos);

	//fill internal read group index and create mode array
	readGroupIndex.initialize(readGroupMapObject.getNumReadGroups());

	//tmp variables for reading
	int lineNum = 0;
	std::vector<std::string> vec;

	//parse file to read details for each read group
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			//check if read group exists
			if(bamHeader->ReadGroups.Contains(vec[0])) {
				//read read group, mate and model
				int rg = readGroupMapObject[findReadGroupIndex(vec[0], bamHeader->ReadGroups)];
				bool isSecondMate = vec[1];
				std::string modelTag = vec[2];

				//clean up vec to only contain parameters (remove read group, mate, model and LL)
				vec.erase(vec.begin(), vec.begin() + 3);

				//create model
				readGroupIndex.setAsUsed(rg, isSecondMate);
				_addModel(vec[2], vec, false);
			} else {
				logfile->warning("Read group '" + vec[0] + "' does not exist in the BAM header! Are you using the correct recal file?");
			}
		}
	}

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	std::pair<int, bool> missingReadGroupInfo;
	vec.clear();
	bool foundOne;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo)){
		if(!foundOne){
			logfile->warning("Missing read groups in file '" + filename + "'!");
			logfile->startIndent("Will assume no recalibration for the following read groups:");
			foundOne = true;
		}

		_printReadGroupNameToLogfile(missingReadGroupInfo);

		//create a model without recalibration
		readGroupIndex.setAsUsed(missingReadGroupInfo.first, missingReadGroupInfo.second);
		_addModel("noRecal", vec, false);
	}
	if(foundOne)
		logfile->endIndent();

	logfile->done();
};

void TRecalibrationEM::_printReadGroupNameToLogfile(std::pair<int, bool> & missingReadGroupInfo){
	if(missingReadGroupInfo.second)
		logfile->list(readGroupNames[missingReadGroupInfo.first] + "(second mate)");
	else
		logfile->list(readGroupNames[missingReadGroupInfo.first] + "(first mate)");
};

void TRecalibrationEM::initializeForParameterEstimation(TParameters & args){
	EMParameters = new TRecalibrationEMEstimationParameters(args);
	modelTagForEstimation = args.getParameterStringWithDefault("model", "full");
	estimationParametersInitialized = true;
};

void TRecalibrationEM::performEstimation(std::string outputName, bool & writeTmpTables){
	if(!estimationParametersInitialized)
		throw "Can not run estimation: EM parameters were not initialized!";

	//count data available for recal
	logfile->listFlush("Counting data available for recal ...");
	TRecalibrationEMDataTable dataTable(readGroupMapObject.numReadGroups, 500);
	addToDataTable(dataTable);
	long _numSites = numSites();
	logfile->done();

	logfile->conclude("Number of sites with data: " + toString(_numSites));
	logfile->conclude("Number of sites with depth > 1: " + toString(numSitesDepthTwoOrMore()));
	logfile->conclude("Number of observations: " + toString(dataTable.totalCounts));
	if(_numSites < 100) throw "Less than 100 sites available for recalibration - aborting estimation!";


	//initialize models based on data tables?
	/*
	THINK: about how to initialize models:
	- currently: a single tag used for all
	- but maybe we should use a file with per read group tags (how to deal with merged read groups?)
	- or choose automatically based on available data?
	*/

	//initialize models
	readGroupIndex.initialize(readGroupMapObject.getNumReadGroups());
	for(int rg = 0; rg < readGroupMapObject.getNumReadGroups(); ++rg){
		for(int mate = 0; mate < 2; ++mate){
			if(dataTable.countsPerReadGroup[rg][0] > 0){
				readGroupIndex.setAsUsed(rg, mate);
				_addModelForEstimation(modelTagForEstimation, dataTable.maxPos[rg][mate]);
			}
		}
	}

	CHECK: report read groups with limited data as warning
	Also report read group swithout data, but not as warning.
	EMParameters->minRequiredObservationsPerReadGroup

	//report models not initialized
	std::pair<int, bool> missingReadGroupInfo;
	bool foundOne;
	while(readGroupIndex.nextNotInUse(missingReadGroupInfo)){
		if(!foundOne){
			logfile->warning("Lack of data for some read groups!");
			logfile->startIndent("The following read groups will have less than " + toString() + " observations and are ignored:");
			foundOne = true;
		}

		_printReadGroupNameToLogfile(missingReadGroupInfo);
	}
	if(foundOne)
		logfile->endIndent();

	//run EM
	runEM();

	//writing final estimates
	filename = outputName + "_recalibrationEM.txt";
	logfile->listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename, LL);
	logfile->done();
};

void TRecalibrationEM::runEM(std::string outputName, bool & writeTmpTables){
	//run EM
	logfile->startNumbering("Running EM algorithm to find MLE recalibration parameters:");

	//initialize tmp variables in windows and model
	prepareWindowsforEM();
	models->initializeEMParams();

	double LL, deltaLL, oldLL = 0.0;
	std::ofstream out;
	std::string filename;

	//running iterations
	for(int iter = 0; iter < numEMIterations; ++iter){
		logfile->number("EM Iteration:"); logfile->addIndent();

		//calculate P(g|d, oldbeta) for all sites and calculate LL
		LL = 0.0;
		logfile->listFlush("Calculating P(g|d, beta') ...");
		for(TRecalibrationEMWindow* curWindow : windows)
			LL += curWindow->fill_P_g_given_d_beta_AND_calcLL(models, tmpEpsilon);
		logfile->done();
		logfile->conclude("Current Log Likelihood = " + toString(LL));

		//DEBUG--------------------------------------------------------
		//calc Q surface for current old params
		//calcQSurface(outputName + "_Qsurface_EMiteration_" + toString(iter) + ".txt", 21);
		//DEBUG--------------------------------------------------------

		//check if we break based on LL
		if(iter > 0){
			deltaLL = LL - oldLL;
			logfile->conclude("Epsilon = " + toString(deltaLL));
			if(iter > 0 && deltaLL < maxEpsilon){
				logfile->conclude("EM has converged (tmpEpsilon < " + toString(maxEpsilon) + ")");
				break;
			} else oldLL = LL;
		} else oldLL = LL;


		//run NewtonRaphson until convergence
		runNewtonRaphson(NewtonRaphsonNumIterations, NewtonRaphsonMaxF, logfile);

		//write current estimates to file
		if(writeTmpTables){
			filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
			logfile->listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename, LL);
			logfile->done();
		}

		//end loop
		logfile->endIndent();
		if(iter == numEMIterations - 1) logfile->warning("EM has not converged after maximum number of iterations!");
	}

	//finalize
	logfile->endNumbering();
};

void TRecalibrationEM::addNewWindow(TBaseFrequencies* freqs){
	windows.push_back(new TRecalibrationEMWindow(freqs, &readGroupMapObject));
	//set iterator
	curWindow = windows.end(); --curWindow;
	if(equalBaseFrequencies) (*curWindow)->setEuqalBaseFrequencies();
}

void TRecalibrationEM::addSite(TSite & site, TQualityMap & qualiMap){
	(*curWindow)->addSite(site, qualityMap);
}

long TRecalibrationEM::numSites(){
	long _numSites = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		_numSites += curWindow->numSites();
	return _numSites;
}

long TRecalibrationEM::numSitesDepthTwoOrMore(){
	long _numSites = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		_numSites += curWindow->numSitesDepthTwoOrMore();
	return _numSites;
}

void TRecalibrationEM::addToDataTable(TRecalibrationEMDataTable & dataTable){
	for(TRecalibrationEMWindow* curWindow : windows)
		curWindow->addToDataTable(dataTable);
	dataTable.assembleCountsPerReadGroup();
};

long TRecalibrationEM::cumulativeDepth(){
	long cumulDepth = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		cumulDepth += curWindow->cumulativeDepth();
	return cumulDepth;
}

void TRecalibrationEM::prepareWindowsforEM(){
	if(tmpEpsilonInitialized) delete[] tmpEpsilon;

	int maxDepth = 0;
	for(TRecalibrationEMWindow* curWindow : windows){
		int tmp = curWindow->getMaxDepth();
		if(tmp > maxDepth)
			maxDepth = tmp;
	}

	//now crate array
	tmpEpsilon = new float[maxDepth];
	tmpEpsilonInitialized = true;
}

void TRecalibrationEM::runNewtonRaphson(int & maxNewtonRaphsonIteratios, double & maxFThreshold, TLog* logfile){
	//variables
	double maxF;
	double lambda; //used in backtracking
	bool acceptMove;
	bool NRconverged = false;

	//calculate Q at current location
	double Q;
	double curQ = 0.0;
	for(TRecalibrationEMWindow* curWindow : windows)
		curQ += curWindow->calcQ(models, tmpEpsilon);

	//run up to maxNewtonRaphsonIteratios iterations, but stop if max(F) < maxFThreshold
	logfile->startIndent("Running Newton-Raphson optimization:");
	for(int i=0; i<maxNewtonRaphsonIteratios; ++i){
		logfile->startIndent("Running iteration " + toString(i+1) + ":");
		logfile->listFlush("Calculating Jacobian and gradient ...");

		//set to zero
		models->setEMParamsToZero();

		//fill Jacobin and F: loop over all sites
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			(*curWindow)->addToJacobianAndF(models, tmpEpsilon);
		}

		//now solve J^-1 x F
		if(models->solveJxF()){
			logfile->done();
/*
			std::cout << "----------------------------------------------" << std::endl;
			std::cout << "JxF " << JxF << std::endl;
			std::cout << "----------------------------------------------" << std::endl;
*/

			//update params for each read group using backtracking
			lambda = 1.0;
			acceptMove = false;
			while(!acceptMove){
				logfile->listFlush("Proposing move with lambda = " + toString(lambda) + " ...");
				models->proposeNewParameters(lambda);

				//calculate Q at new location
				Q = 0.0;
				for(TRecalibrationEMWindow* curWindow : windows)
					Q += curWindow->calcQ(models, tmpEpsilon);

				//check if we accept or backtrack
				if(Q > curQ){
					acceptMove = true; //accept
					logfile->write(" accepting move!");
					logfile->conclude("Q was reduced from " + toString(curQ) + " to " + toString(Q));
					curQ = Q;
				} else {
					lambda = lambda / 2.0; //backtrack;
					logfile->write(" rejecting move!");
					models->rejectProposedParameters();
					if(lambda < 0.000000001){
						acceptMove = true; //accept
						NRconverged = true;
						logfile->conclude("No improvement even with lambda = " + toString(lambda) + ", aborting Newton-Raphson.");
					}
				}
			}
		} else {
			models->printJacobianToStdOut();
			throw "Issue solving JxF in TRecalibrationEM::runNewtonRalphson()! This may be due to a lack of data. Consider adding more sites.";
		}

		//get largest gradient (F) to check if we break
		maxF = models->getSteepestGradient();
		logfile->conclude("max(F) = " + toString(maxF));
		logfile->endIndent();
		if(maxF < maxFThreshold || NRconverged) break;
	}
	logfile->endIndent();
}




void TRecalibrationEM::writeCurrentEstimates(std::string filename, double & LL){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	writeHeader(out);
	writeParams(out, LL);
	out.close();
}

void TRecalibrationEM::writeHeader(std::ofstream & out){
	out << "readGroup\t";
	models->writeHeader(out);
	out << "\tLL" << std::endl;
}

void TRecalibrationEM::writeParams(std::ofstream & out, double & LL){
	for(int r=0; r<readGroupMapObject.getOrigNumReadGroups(); ++r){
		out << readGroupNames[r];
		models->writeParametersToFile(out, readGroupMapObject[r]);
		out << "\t" << LL;
		out << std::endl;
	}
}

double TRecalibrationEM::calcLL(){
	double LL = 0.0;
	for(TRecalibrationEMWindow* curWindow : windows)
		LL += curWindow->calcLL(models, tmpEpsilon);
	return LL;
}

/*
void TRecalibrationEM::calcLikelihoodSurface(std::string filename, int numMarginalGridPoints){
	double LL;

	//open outputfile
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	out << "beta0\tbeta1\tLL" << std::endl;

	//set min, max and step for each parameter
	double min[5];
	min[0] = -5.0;
	min[1] = -5.0;
	min[2] = -1.0;
	min[3] = -1.0;
	min[4] = -1.0;


	double max[5];
	max[0] = 10.0;
	max[1] = 10.0;
	max[2] = 1.0;
	max[3] = 1.0;
	max[4] = 1.0;

	double step[5];
	for(int i=0; i<5; ++i){
		step[i] = (max[i] - min[i]) / (numMarginalGridPoints - 1.0);
	}

	//without last two
	for(int r=0; r<numReadGroups; ++r){
		params[r][3] = 0.0;
		params[r][4] = 0.0;
	}

	//Loop over parameters
	for(int p1=0; p1<numMarginalGridPoints; ++p1){
		//for(int r=0; r<numReadGroups; ++r) params[r][0] = min[0] + p1 * step[0];
		params[0][0] = min[0] + p1 * step[0];
		for(int p2=0; p2<numMarginalGridPoints; ++p2){
			//for(int r=0; r<numReadGroups; ++r) params[r][1] = min[1] + p2 * step[1];
			params[0][1] = min[1] + p2 * step[1];
			for(int p3=0; p3<numMarginalGridPoints; ++p3){
				//for(int r=0; r<numReadGroups; ++r) params[r][2] = min[2] + p3 * step[2];
				params[0][2] = min[2] + p3 * step[2];
				//for(int p4=0; p4<numMarginalGridPoints; ++p4){
					//for(int r=0; r<numReadGroups; ++r) params[r][3] = min[3] + p4 * step[3];
					//for(int p5=0; p5<numMarginalGridPoints; ++p5){
						//for(int r=0; r<numReadGroups; ++r) params[r][4] = min[4] + p5 * step[4];


						//calculate LL
						LL = 0.0;
						for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
							LL += (*curWindow)->calcLL(params);
						}

						//write to file
						for(int i=0; i<5; ++i) out << params[0][i] << "\t";
						out << LL << std::endl;
					//}
				//}
			}
		}
	}

	//close file
	out.close();
}


void TRecalibrationEM::calcQSurface(std::string filename, int numMarginalGridPoints){
	double Q;

	//open outputfile
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	out << "beta0\tbeta1\tQ" << std::endl;

	//set min, max and step for each parameter
	double min[2];
	min[0] = -5.0;
	min[1] = -5.0;

	double max[2];
	max[0] = 10.0;
	max[1] = 10.0;

	double step[2];
	for(int i=0; i<2; ++i){
		step[i] = (max[i] - min[i]) / (numMarginalGridPoints - 1.0);
	}

	//print old params

	//Loop over parameters
	for(int p1=0; p1<numMarginalGridPoints; ++p1){
		for(int r=0; r<numReadGroups; ++r) newParams[r][0] = min[0] + p1 * step[0];
		for(int p2=0; p2<numMarginalGridPoints; ++p2){
			for(int r=0; r<numReadGroups; ++r) newParams[r][1] = min[1] + p2 * step[1];
			//for(int p3=0; p3<numMarginalGridPoints; ++p3){
				//for(int r=0; r<numReadGroups; ++r) params[r][2] = min[2] + p3 * step[2];
				//for(int p4=0; p4<numMarginalGridPoints; ++p4){
					//for(int r=0; r<numReadGroups; ++r) params[r][3] = min[3] + p4 * step[3];
					//for(int p5=0; p5<numMarginalGridPoints; ++p5){
						//for(int r=0; r<numReadGroups; ++r) params[r][4] = min[4] + p5 * step[4];


						//calculate Q
						Q = 0.0;
						//logfile->listFlush("Calculating Q at {" + toString(params[0][0]) + ", " + toString(params[0][1]) + "} ...");
						for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
							Q += (*curWindow)->calcQ(newParams);
						}
						//logfile->done();
						//logfile->conclude("Current Q = " + toString(Q));

						//write to file
						for(int i=0; i<2; ++i) out << newParams[0][i] << "\t";
						out << Q << std::endl;
					//}
				//}
			//}
		}
	}

	//close file
	out.close();
}

*/


double TRecalibrationEM::getErrorRate(TBase & base){
	return models->getErrorRate(base);
}

int TRecalibrationEM::getQuality(TBase & base){
	double q = models->getErrorRate(base);
	//transform to quality
	return qualityMap.errorToQuality(q);
}




