/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEM.h"

//---------------------------------------------------------------
//TRecalibrationEMEstimationParameters
//---------------------------------------------------------------

TRecalibrationEMEstimationParameters::TRecalibrationEMEstimationParameters(TParameters & args, TLog* logfile){
	logfile->startIndent("Settings regarding the EM algorithm:");
	minRequiredObservations = 1000; //constant for reporting
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
		data[k].position = it->distFrom5Prime;

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

double TRecalibrationEMSite::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, float* & freqs, float* & epsilon){
	calcEpsilon(models, epsilon);

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

double TRecalibrationEMSite::calcLL(TRecalibrationEMModels & models, float* & freqs, float* & epsilon){
	calcEpsilon(models, epsilon);

	//over all genotypes
	double LL = 0.0;
	for(int g=0; g<4; ++g){
		double tmp = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			double B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			tmp *= B * epsilon[k] - data[k].D[g] + 1.0;
		}
		LL += tmp * freqs[g];
	}

	//return LL = P_g_given_d_theta_denominator
	return log(LL);
}

double TRecalibrationEMSite::calcQ(TRecalibrationEMModels & models, float* & epsilon){
	calcEpsilon(models, epsilon);

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

void TRecalibrationEMSite::addToJacobianAndF(TRecalibrationEMModels & models, float* & epsilon){
	//calculate tmpEpsilon with current parameters
	calcEpsilon(models, epsilon);

	//tmp variables
	double* eps1MinusEps = new double[numReads];
	double* oneMinus2Eps = new double[numReads];

	for(unsigned int k=0; k<numReads; ++k){
		eps1MinusEps[k] = epsilon[k] * (1.0 - epsilon[k]);
		oneMinus2Eps[k] = 1.0 - 2.0 * epsilon[k];
	}

	//fill F and Jacobian
	for(int g=0; g<4; ++g){
		for(unsigned int k=0; k<numReads; ++k){
			//calc weights
			double B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			double weightF = B / (1.0 - data[k].D[g] + B * epsilon[k]) * eps1MinusEps[k] * P_g_given_d_oldBeta[g];
			double weightJacobian = P_g_given_d_oldBeta[g] * weightF * (oneMinus2Eps[k] - weightF);

			//add to model
			models.addToFandJacobian(data[k], weightF, weightJacobian);
		}
	}

	//delete tmp variables
	delete[] eps1MinusEps;
	delete[] oneMinus2Eps;
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

double TRecalibrationEMWindow::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, float* & tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites){
		LL += site->fill_P_g_given_d_beta_AND_calcLL(models, freqs, tmpEpsilon);
	}
	return LL;
}

double TRecalibrationEMWindow::calcLL(TRecalibrationEMModels & models, float* & tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites)
		LL += site->calcLL(models, freqs, tmpEpsilon);
	return LL;
}

double TRecalibrationEMWindow::calcQ(TRecalibrationEMModels & models, float* & tmpEpsilon){
	double Q = 0.0;
	for(TRecalibrationEMSite* site : sites)
		Q += site->calcQ(models, tmpEpsilon);
	return Q;
}

void TRecalibrationEMWindow::addToJacobianAndF(TRecalibrationEMModels & models, float* & tmpEpsilon){
	for(TRecalibrationEMSite* site : sites)
		site->addToJacobianAndF(models, tmpEpsilon);
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
	estimationParametersInitialized = false;
	EMParameters = NULL;

	models = new TRecalibrationEMModels(readGroupMapObject.getNumReadGroups(), logfile);
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
	models->addSingleModelForAllReadGroups(modelTag, vec, true);

	logfile->endIndent();
};

void TRecalibrationEM::initializeRecalibrationParametersFromFile(std::string filename){
	//read parameters from file
	logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "' for reading!";

	//skip header
	std::string line;
	std::getline(file, line);

	//tmp variables for reading
	int lineNum = 0;
	std::vector<std::string> vec;

	//parse file to read details for each read group
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() != 6)
				throw "Wrong number of entries in file '" + filename + "' on line " + toString(lineNum) + ": expected 6 but found " + toString(vec.size()) + "!";

			//check if read group exists
			if(bamHeader->ReadGroups.Contains(vec[0])) {
				//read read group, mate and model
				int rg = readGroupMapObject[findReadGroupIndex(vec[0], bamHeader->ReadGroups)];
				bool isSecondMate = stringToInt(vec[1]);
				std::string modelTag = vec[2];

				//clean up vec to only contain parameters (remove read group, mate, model and LL)
				vec.erase(vec.begin(), vec.begin() + 3);

				//create model
				models->addModel(rg, isSecondMate, modelTag, vec, false);
			} else {
				logfile->warning("Read group '" + vec[0] + "' does not exist in the BAM header! Are you using the correct recal file?");
			}
		}
	}
	logfile->done();

	//report read groups for which no recal model was given and initialize them as "no_recal" model
	if(models->hasReadGroupsWithoutModel()){
		logfile->warning("Missing read groups in file '" + filename + "'!");
		logfile->startIndent("Will assume no recalibration for the following read groups:");
		models->reportReadGroupsNotUsed(readGroupNames);
		models->addNoRecalModelIfMissing();
		logfile->endIndent();
	}
};

void TRecalibrationEM::initializeForParameterEstimation(TParameters & args){
	EMParameters = new TRecalibrationEMEstimationParameters(args, logfile);
	modelTagForEstimation = args.getParameterStringWithDefault("model", "qualFuncPosFuncContext");
	estimationParametersInitialized = true;
};

void TRecalibrationEM::performEstimation(std::string outputName, bool & writeTmpTables){
	if(!estimationParametersInitialized)
		throw "Can not run estimation: EM parameters were not initialized!";

	//count data available for recal
	logfile->listFlush("Counting data available for recal ...");
	TRecalibrationEMDataTable dataTable(readGroupMapObject.numReadGroups, 500);
	addToDataTable(dataTable);
	int numSitesWithData = numSites();
	logfile->done();

	logfile->conclude("Number of sites with data: " + toString(numSitesWithData));
	logfile->conclude("Number of sites with depth > 1: " + toString(numSitesDepthTwoOrMore()));
	logfile->conclude("Number of observations: " + toString(dataTable.totalCounts));
	if(numSitesWithData < 100) throw "Less than 100 sites available for recalibration - aborting estimation!";


	//initialize models based on data tables?
	/*
	THINK: about how to initialize models:
	- currently: a single tag used for all
	- but maybe we should use a file with per read group tags (how to deal with merged read groups?)
	- or choose automatically based on available data?
	*/

	//initialize models
	//-------------------
	logfile->listFlush("Initializing recalibration models ...");
	int numModelsWithLittleData = 0;
	int numModelsWithoutData = 0;

	for(int rg = 0; rg < readGroupMapObject.getNumReadGroups(); ++rg){
		for(int mate = 0; mate < 2; ++mate){
			if(dataTable.countsPerReadGroup[rg][mate] > 0){
				models->addModel(rg, mate, modelTagForEstimation, dataTable.maxPos[rg][mate]);

				if(dataTable.countsPerReadGroup[rg][mate] < EMParameters->minRequiredObservations)
					++numModelsWithLittleData;
			} else ++numModelsWithoutData;
		}
	}
	logfile->done();
	logfile->conclude("Initialized " + toString(models->numModels()) + " models.");

	//report warning in case of very little data
	if(numModelsWithLittleData > 0){
		logfile->warning("Some read groups have very little data!");
		logfile->startIndent("Consider merging these read groups:");
		for(int rg = 0; rg < readGroupMapObject.getNumReadGroups(); ++rg){
			if(dataTable.countsPerReadGroup[rg][0] > 0 && dataTable.countsPerReadGroup[rg][0] < EMParameters->minRequiredObservations)
				logfile->list(readGroupNames[rg] + " (first mate)");
			if(dataTable.countsPerReadGroup[rg][1] > 0 && dataTable.countsPerReadGroup[rg][1] < EMParameters->minRequiredObservations)
				logfile->list(readGroupNames[rg] + " (second mate)");
		}
		logfile->endIndent();
	}

	//report read groups without data
	if(numModelsWithoutData > 0){
		logfile->startIndent("The following " + toString(numModelsWithoutData) + " read groups do not have data and will not be recalibrated:");
		models->reportReadGroupsNotUsed(readGroupNames);
		logfile->endIndent();
	}

	//run EM
	//-------------------
	runEM(numSitesWithData, outputName, writeTmpTables);

	//writing final estimates
	//-------------------
	std::string filename = outputName + "_recalibrationEM.txt";
	logfile->listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename);
	logfile->done();
};

void TRecalibrationEM::runEM(int numSitesWithData, std::string outputName, bool & writeTmpTables){
	//run EM
	logfile->startNumbering("Running EM algorithm to find MLE recalibration parameters:");

	//initialize tmp variables in windows and model
	prepareWindowsforEM();

	double LL, deltaLL, oldLL = 0.0;
	std::ofstream out;
	std::string filename;

	//running iterations
	for(int iter = 0; iter < EMParameters->numEMIterations; ++iter){
		logfile->number("EM Iteration:"); logfile->addIndent();

		//calculate P(g|d, oldbeta) for all sites and calculate LL
		LL = 0.0;
		logfile->listFlush("Calculating P(g|d, beta') ...");
		for(TRecalibrationEMWindow* curWindow : windows)
			LL += curWindow->fill_P_g_given_d_beta_AND_calcLL(*models, tmpEpsilon);
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
			if(iter > 0 && deltaLL < EMParameters->maxEpsilon){
				logfile->conclude("EM has converged (tmpEpsilon < " + toString(EMParameters->maxEpsilon) + ")");
				break;
			} else oldLL = LL;
		} else oldLL = LL;


		//run NewtonRaphson until convergence
		runNewtonRaphson(numSitesWithData);

		//write current estimates to file
		if(writeTmpTables){
			filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
			logfile->listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename);
			logfile->done();
		}

		//end loop
		logfile->endIndent();
		if(iter == EMParameters->numEMIterations - 1) logfile->warning("EM has not converged after maximum number of iterations!");
	}

	//finalize
	logfile->endNumbering();
};

void TRecalibrationEM::runNewtonRaphson(int numSitesWithData){
	//variables
	double maxF;
	double lambda; //used in backtracking
	bool acceptMove;
	bool NRconverged = false;

	//calculate Q at current location
	double Q;
	double curQ = 0.0;
	for(TRecalibrationEMWindow* curWindow : windows)
		curQ += curWindow->calcQ(*models, tmpEpsilon);

	//run up to maxNewtonRaphsonIteratios iterations, but stop if max(F) < maxFThreshold
	logfile->startIndent("Running Newton-Raphson optimization:");
	for(int i=0; i<EMParameters->NewtonRaphsonNumIterations; ++i){
		logfile->startIndent("Running iteration " + toString(i+1) + ":");
		logfile->listFlush("Calculating Jacobian and gradient ...");

		//set to zero
		models->setEMParamsToZero();

		//fill Jacobin and F: loop over all sites
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			(*curWindow)->addToJacobianAndF(*models, tmpEpsilon);
		}

		//now solve J^-1 x F
		if(models->solveJxF(numSitesWithData)){
			logfile->done();

			//update params for each read group using backtracking
			lambda = 1.0;
			acceptMove = false;
			while(!acceptMove){
				logfile->listFlush("Proposing move with lambda = " + toString(lambda) + " ...");
				models->proposeNewParameters(lambda);

				//calculate Q at new location
				Q = 0.0;
				for(TRecalibrationEMWindow* curWindow : windows)
					Q += curWindow->calcQ(*models, tmpEpsilon);

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
		if(maxF < EMParameters->NewtonRaphsonMaxF || NRconverged) break;
	}
	logfile->endIndent();
};

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
};

void TRecalibrationEM::addNewWindow(TBaseFrequencies* freqs){
	windows.push_back(new TRecalibrationEMWindow(freqs, &readGroupMapObject));
	//set iterator
	curWindow = windows.end(); --curWindow;
	if(EMParameters->equalBaseFrequencies) (*curWindow)->setEuqalBaseFrequencies();
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
};

void TRecalibrationEM::writeCurrentEstimates(std::string filename){
	TOutputFilePlain out(filename);
	models->writeHeader(out);
	models->writeParameters(out, readGroupNames);
};

double TRecalibrationEM::calcLL(){
	double LL = 0.0;
	for(TRecalibrationEMWindow* curWindow : windows)
		LL += curWindow->calcLL(*models, tmpEpsilon);
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





