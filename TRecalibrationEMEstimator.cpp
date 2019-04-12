/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMEstimator.h"

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
		data[k].readGroup = ReadGroupMap[it->readGroup];

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

double TRecalibrationEMSite::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, double* & freqs, double* & epsilon){
	calcEpsilon(models, epsilon);

	//over all genotypes
	double P_g_given_d_theta_denominator = 0.0;
	for(int g=0; g<4; ++g){
		double tmp = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			double B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			tmp *= B * epsilon[k] - data[k].D[g] + 1.0;
		}
		P_g_given_d_oldBeta[g] = tmp * freqs[g];
		P_g_given_d_theta_denominator += P_g_given_d_oldBeta[g];
	}

	double max = 0.0;

	if(P_g_given_d_theta_denominator < 1.0E-20){
		//do again but in log

		for(int g=0; g<4; ++g){
			double tmp = 0.0;
			//loop over all reads
			for(unsigned int k=0; k<numReads; ++k){
				double B = 4.0 / 3.0 * data[k].D[g] - 1.0;
				tmp += log(B * epsilon[k] - data[k].D[g] + 1.0);
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
	return log(P_g_given_d_theta_denominator) + max;
};

double TRecalibrationEMSite::calcLL(TRecalibrationEMModels & models, double* & freqs, double* & epsilon){
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

	if(LL < 1.0E-20){
		//do again but in log
		double max = 0.0;
		double tmp_log[4];
		for(int g=0; g<4; ++g){
			double tmp = 0.0;
			//loop over all reads
			for(unsigned int k=0; k<numReads; ++k){
				double B = 4.0 / 3.0 * data[k].D[g] - 1.0;
				tmp += log(B * epsilon[k] - data[k].D[g] + 1.0);
			}
			tmp_log[g] = tmp + log(freqs[g]);
			if(g==0) max = tmp_log[g];
			else if(tmp_log[g] > max) max = tmp_log[g];
		}

		//rescale and delog
		LL = 0.0;
		for(int g=0; g<4; ++g){
			LL += exp(tmp_log[g] - max);
		}

		return log(LL) + max;
	} else {
		return log(LL);
	}
};

double TRecalibrationEMSite::calcQ(TRecalibrationEMModels & models, double* & epsilon){
	calcEpsilon(models, epsilon);

	//now calculate P(d, g, new params)
	double Q = 0.0;

	for(int g=0; g<4; ++g){
		double P_d_given_g_beta = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			double B = 4.0 / 3.0 * data[k].D[g] - 1.0;
			P_d_given_g_beta *= B * epsilon[k] - data[k].D[g] + 1;
		}

		if(P_d_given_g_beta < 1.0E-50) P_d_given_g_beta = 1.0E-50;
		Q += P_g_given_d_oldBeta[g] * log(P_d_given_g_beta);
	}

	return Q;
};

void TRecalibrationEMSite::addToQ(TRecalibrationEMModels & models){
	for(unsigned int k=0; k<numReads; ++k){
		models.addToQ(data[k], P_g_given_d_oldBeta);
	}
};

void TRecalibrationEMSite::addToJacobianAndF(TRecalibrationEMModels & models, double* & epsilon){
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
	freqs = new double[4];
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

double TRecalibrationEMWindow::fill_P_g_given_d_beta_AND_calcLL(TRecalibrationEMModels & models, double* & tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites){
		LL += site->fill_P_g_given_d_beta_AND_calcLL(models, freqs, tmpEpsilon);
	}
	return LL;
}

double TRecalibrationEMWindow::calcLL(TRecalibrationEMModels & models, double* & tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites)
		LL += site->calcLL(models, freqs, tmpEpsilon);
	return LL;
};

double TRecalibrationEMWindow::calcQ(TRecalibrationEMModels & models, double* & tmpEpsilon){
	double Q = 0.0;
	for(TRecalibrationEMSite* site : sites)
		Q += site->calcQ(models, tmpEpsilon);
	return Q;
};

void TRecalibrationEMWindow::addToQ(TRecalibrationEMModels & models){
	for(TRecalibrationEMSite* site : sites)
		site->addToQ(models);
};

void TRecalibrationEMWindow::addToJacobianAndF(TRecalibrationEMModels & models, double* & tmpEpsilon){
	for(TRecalibrationEMSite* site : sites)
		site->addToJacobianAndF(models, tmpEpsilon);
};

void TRecalibrationEMWindow::setEuqalBaseFrequencies(){
	for(int i=0; i<4; ++i) freqs[i] = 0.25;
};

//---------------------------------------------------------------
//TRecalibrationEMEstimator
//---------------------------------------------------------------
TRecalibrationEMEstimator::TRecalibrationEMEstimator(TParameters & args, TReadGroups* ReadGroups, TLog* Logfile, TReadGroupMap* ReadGroupMap){
	logfile = Logfile;
	tmpEpsilon = NULL;
	tmpEpsilonInitialized = false;
	maxDepth = -1;

	//read groups
	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMap;

	//models
	models = new TRecalibrationEMModels(_readGroupMap->numReadGroups, logfile);

	//estimation parameters
	logfile->startIndent("Settings regarding the EM algorithm:");
	std::string recalFile = args.getParameterString("initialRecalValues", false);
	modelTagForEstimation = args.getParameterStringWithDefault("model", "qualFuncPosFuncContext");
	if(recalFile.empty()){
		logfile->list("Will fit the model '" + modelTagForEstimation + "' for all read groups.");
	} else {
		logfile->list("Will read models and initial values from file '" + recalFile + "'.");
		logfile->list("Will fit the model '" + modelTagForEstimation + "' for read groups not in file '" + recalFile + "'.");
		initializeFromString(recalFile);
	}

	minRequiredObservations = 10000; //constant for reporting
	numEMIterations = args.getParameterIntWithDefault("iterations", 200);
	logfile->list("Will perform at max " + toString(numEMIterations) + " EM iterations.");
	maxEpsilon = args.getParameterDoubleWithDefault("maxEps", 0.000001);
	logfile->list("Will stop EM when deltaLL < " + toString(maxEpsilon));
	NewtonRaphsonNumIterations = args.getParameterIntWithDefault("NRiterations", 20);
	logfile->list("Will conduct at max " + toString(NewtonRaphsonNumIterations) + " Newton-Raphson iterations");
	NewtonRaphsonMaxF = args.getParameterDoubleWithDefault("maxF", 0.0001);
	logfile->list("Will stop Newton-Raphson when F < " + toString(NewtonRaphsonMaxF));
	equalBaseFrequencies = args.parameterExists("equalBaseFreq");
	if(equalBaseFrequencies) logfile->list("Will assume equal base frequencies {0.25, 0.25, 0.25, 0.25}");
	logfile->endIndent();
};

void TRecalibrationEMEstimator::initializeFromString(const std::string string){
	models->createModels(string, *_readGroups, *_readGroupMap);
}

void TRecalibrationEMEstimator::performEstimation(std::string outputName, bool & writeTmpTables){
	//count data available for recal
	logfile->listFlush("Counting data available for recal ...");
	TRecalibrationEMDataTable dataTable(_readGroups->size(), 500);
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
	std::vector<int> qualities;

	for(int rg = 0; rg < _readGroupMap->numReadGroups; ++rg){
		for(int mate = 0; mate < 2; ++mate){
			if(dataTable.countsPerReadGroup[rg][mate] > 0){
				dataTable.fillVectorWithUsedQualities(rg, mate, qualities);
				models->addModelIfItDoesNotExist(rg, mate, modelTagForEstimation, qualities, dataTable.maxPos[rg][mate]);
				if(dataTable.countsPerReadGroup[rg][mate] < minRequiredObservations)
					++numModelsWithLittleData;
			} else {
				if(models->modelExists(rg, mate)){
					models->removeModel(rg, mate);
				}
				++numModelsWithoutData;
			}
		}
	}
	logfile->done();
	logfile->conclude("Initialized " + toString(models->numModels()) + " models.");

	//report warning in case of very little data
	if(numModelsWithLittleData > 0){
		logfile->warning("Some read groups have very little data!");
		logfile->startIndent("Consider merging these read groups:");

		for(int rg = 0; rg < _readGroups->size(); ++rg){
			int index = _readGroupMap->getIndex(rg);
			if(dataTable.countsPerReadGroup[index][0] > 0 && dataTable.countsPerReadGroup[index][0] < minRequiredObservations)
				logfile->list(_readGroups->getName(rg)  + " (first mate): only " + toString(dataTable.countsPerReadGroup[index][0]) + " observations.");
			if(dataTable.countsPerReadGroup[index][1] > 0 && dataTable.countsPerReadGroup[index][1] < minRequiredObservations)
				logfile->list(_readGroups->getName(rg) + " (second mate): only " + toString(dataTable.countsPerReadGroup[index][0]) + " observations.");
		}
		logfile->endIndent();
	}

	//report read groups without data
	if(numModelsWithoutData > 0){
		logfile->startIndent("The following " + toString(numModelsWithoutData) + " read groups do not have data and will not be recalibrated:");
		models->reportReadGroupsNotUsed(*_readGroups, *_readGroupMap);
		logfile->endIndent();
	}

	//run EM
	//-------------------
	_runEM(numSitesWithData, outputName, writeTmpTables);

	//writing final estimates
	//-------------------
	std::string filename = outputName + "_recalibrationEM.txt";
	logfile->listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename);
	logfile->done();
};

void TRecalibrationEMEstimator::_runEM(int numSitesWithData, std::string outputName, bool & writeTmpTables){
	//run EM
	logfile->startNumbering("Running EM algorithm to find MLE recalibration parameters:");

	//initialize tmp variables in windows and model
	_initializeTmpEpsilon();

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
			LL += curWindow->fill_P_g_given_d_beta_AND_calcLL(*models, tmpEpsilon);
		logfile->done();
		logfile->conclude("Current Log Likelihood = " + toString(LL));

		//DEBUG--------------------------------------------------------
		//calc and print LL
		//logfile->conclude("DEBUG LL = " + toString(calcLL()));
		//DEBUG--------------------------------------------------------

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
		_runNewtonRaphson(numSitesWithData);

		//write current estimates to file
		if(writeTmpTables){
			filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
			logfile->listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename);
			logfile->done();
		}

		//end loop
		logfile->endIndent();
		if(iter == numEMIterations - 1) logfile->warning("EM has not converged after maximum number of iterations!");
	}

	//finalize
	logfile->endNumbering();
};

void TRecalibrationEMEstimator::_runNewtonRaphson(int numSitesWithData){
	//calculate Q at current location
	models->setQToZero();
	for(TRecalibrationEMWindow* curWindow : windows)
		curWindow->addToQ(*models);
	double curQ = models->curQ();

	//run up to maxNewtonRaphsonIteratios iterations, but stop if max(F) < maxFThreshold
	logfile->startIndent("Running Newton-Raphson optimization:");
	for(int i=0; i<NewtonRaphsonNumIterations; ++i){
		logfile->startIndent("Running Newton-Raphson iteration " + toString(i+1) + ":");

		//fill Jacobin and F: loop over all sites
		logfile->listFlush("Calculating Jacobian and gradient ...");
		models->setEMParamsToZero();
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow)
			(*curWindow)->addToJacobianAndF(*models, tmpEpsilon);

		//now solve J^-1 x F
		models->solveJxF();
		logfile->done();

		//update params for each read group using backtracking
		double lambda = 1.0;
		int log2_lambda = 0;
		int numUpdatedModels = 0;
		int numUpdatedModels_old;

		while(numUpdatedModels < models->numModels() && lambda > 1.0E-20){
			//propose move
			logfile->listFlush("Proposing move with log2(lambda) = " + toString(log2_lambda) + " ... ");
			models->proposeNewParameters(lambda);

			//calculate Q at new location
			models->setQToZero();
			for(TRecalibrationEMWindow* curWindow : windows)
				curWindow->addToQ(*models);

			//check if we accept or backtrack
			numUpdatedModels_old = numUpdatedModels;
			numUpdatedModels = models->acceptProposedParametersBasedOnQ();
			double Q = models->curQ();

			logfile->write(toString(numUpdatedModels) + "/" + toString(models->numModels()) + " models converged.");

			if(numUpdatedModels > numUpdatedModels_old){
				logfile->conclude("Q was increased from " + toString(curQ) + " to " + toString(Q));
			}
			curQ = Q;

			//backtrack
			lambda = lambda / 2.0; //backtrack;
			--log2_lambda;
		}

		if(numUpdatedModels < models->numModels()){
			logfile->conclude("Some models did not improve even with log2(lambda) = " + toString(log2_lambda) + ", aborting Newton-Raphson.");
		}

		//get largest gradient (F) to check if we break NR optimization
		double maxF = models->getSteepestGradient();
		logfile->conclude("max(F) = " + toString(maxF));
		logfile->endIndent();
		if(maxF < NewtonRaphsonMaxF || numUpdatedModels == 0) break;
	}
	logfile->endIndent();
};

void TRecalibrationEMEstimator::_initializeTmpEpsilon(){
	if(!tmpEpsilonInitialized){
		int maxDepth = 0;
		for(TRecalibrationEMWindow* curWindow : windows){
			int tmp = curWindow->getMaxDepth();
			if(tmp > maxDepth)
				maxDepth = tmp;
		}

		//now crate array
		tmpEpsilon = new double[maxDepth];
		tmpEpsilonInitialized = true;
	}
};

void TRecalibrationEMEstimator::addNewWindow(TBaseFrequencies* freqs){
	windows.push_back(new TRecalibrationEMWindow(freqs, _readGroupMap));
	//set iterator
	curWindow = windows.end(); --curWindow;
	if(equalBaseFrequencies) (*curWindow)->setEuqalBaseFrequencies();
}

void TRecalibrationEMEstimator::addSite(TSite & site, TQualityMap & qualiMap){
	(*curWindow)->addSite(site, qualiMap);
}

long TRecalibrationEMEstimator::numSites(){
	long _numSites = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		_numSites += curWindow->numSites();
	return _numSites;
}

long TRecalibrationEMEstimator::numSitesDepthTwoOrMore(){
	long _numSites = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		_numSites += curWindow->numSitesDepthTwoOrMore();
	return _numSites;
}

void TRecalibrationEMEstimator::addToDataTable(TRecalibrationEMDataTable & dataTable){
	for(TRecalibrationEMWindow* curWindow : windows)
		curWindow->addToDataTable(dataTable);
	dataTable.assembleCountsPerReadGroup();
};

long TRecalibrationEMEstimator::cumulativeDepth(){
	long cumulDepth = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		cumulDepth += curWindow->cumulativeDepth();
	return cumulDepth;
};

void TRecalibrationEMEstimator::writeCurrentEstimates(std::string filename){
	TOutputFilePlain out(filename);
	models->writeHeader(out);
	models->writeParameters(out, *_readGroups, *_readGroupMap);
};

double TRecalibrationEMEstimator::calcLL(){
	_initializeTmpEpsilon();
	double LL = 0.0;
	for(TRecalibrationEMWindow* curWindow : windows)
		LL += curWindow->calcLL(*models, tmpEpsilon);
	return LL;
};

/*
void TRecalibrationEMEstimator::calcLikelihoodSurface(std::string filename, int numMarginalGridPoints){
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




