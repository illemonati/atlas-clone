/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "../GenotypeLikelihoods/TRecalibrationEMEstimator.h"

namespace GenotypeLikelihoods{

/*
//---------------------------------------------------------------
//RecalibrationEMSite
//---------------------------------------------------------------
TRecalibrationEMSite::TRecalibrationEMSite(){
	trueBase = N;
};

TRecalibrationEMSite::TRecalibrationEMSite(TSite & site, BAM::TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap){
	trueBase = N;
	_save(site, ReadGroupMap, qualiMap);
};

TRecalibrationEMSite::TRecalibrationEMSite(TSite & site, BAM::TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap, const Base TrueBase){
	trueBase = TrueBase;
	_save(site, ReadGroupMap, qualiMap);
};

void TRecalibrationEMSite::_save(TSite & site, BAM::TReadGroupMap & ReadGroupMap, TQualityMap & qualiMap){
	data.reserve(site.depth());

	for(auto& b : site){
		data.push_back(*b);
	}
};

void TRecalibrationEMSite::addToDataTable(TRecalibrationEMDataTables & dataTable){
	for(auto& b : data){
		dataTable.add(b);
	}
};

double TRecalibrationEMSite::fill_P_g_given_d_beta_AND_calcLL(TSequencingErrorModels & models, TPostMortemDamage _pmd, double* freqs, double* epsilon){
	calcEpsilon(models, epsilon);

	_pmd.

	for(uint8_t b = 0; b<4; ++b){
		P_bbar_given_d_oldTheta

	}



	//------ OLD ------------- ///
	//over all genotypes
	double P_g_given_d_theta_denominator = 0.0;
	for(int g=0; g<4; ++g){
		double tmp = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			tmp *= calcB(data[k].D[g]) * epsilon[k] - data[k].D[g] + 1.0;
		}
		P_bbar_given_d_oldTheta[g] = tmp * freqs[g];
		P_g_given_d_theta_denominator += P_bbar_given_d_oldTheta[g];
	}

	double max = 0.0;

	if(P_g_given_d_theta_denominator < 1.0E-20){
		//do again but in log

		for(int g=0; g<4; ++g){
			double tmp = 0.0;
			//loop over all reads
			for(unsigned int k=0; k<numReads; ++k){
				tmp += log( calcB(data[k].D[g]) * epsilon[k] - data[k].D[g] + 1.0 );
			}
			P_bbar_given_d_oldTheta[g] = tmp + log(freqs[g]);
			if(g==0) max = P_bbar_given_d_oldTheta[g];
			else if(P_bbar_given_d_oldTheta[g] > max) max = P_bbar_given_d_oldTheta[g];
		}

		//rescale and delog
		P_g_given_d_theta_denominator = 0.0;
		for(int g=0; g<4; ++g){
			P_bbar_given_d_oldTheta[g] = exp(P_bbar_given_d_oldTheta[g] - max);
			P_g_given_d_theta_denominator += P_bbar_given_d_oldTheta[g];
		}
	}

	//calculate P(g|d, theta)
	for(int g=0; g<4; ++g){
		P_bbar_given_d_oldTheta[g] = P_bbar_given_d_oldTheta[g] / P_g_given_d_theta_denominator;
	}

	//return LL = P_g_given_d_theta_denominator
	return log(P_g_given_d_theta_denominator) + max;
};

double TRecalibrationEMSite::calcLL(TSequencingErrorModels & models, double* freqs, double* epsilon){
	calcEpsilon(models, epsilon);

	//over all genotypes
	double LL = 0.0;
	for(int g=0; g<4; ++g){
		double tmp = 1.0;
		//loop over all reads
		for(unsigned int k=0; k<numReads; ++k){
			tmp *= calcB(data[k].D[g]) * epsilon[k] - data[k].D[g] + 1.0;
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
				tmp += log( calcB(data[k].D[g]) * epsilon[k] - data[k].D[g] + 1.0 );
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

void TRecalibrationEMSite::addToQ(TSequencingErrorModels & models){
	if(trueBase == N){
		for(unsigned int k=0; k<numReads; ++k){
			models.addToQ(data[k], P_bbar_given_d_oldTheta);
		}
	} else {
		for(unsigned int k=0; k<numReads; ++k){
			models.addToQ(data[k], trueBase);
		}
	}
};

void TRecalibrationEMSite::_addToJacobianAndF(TSequencingErrorModels & models, double* epsilon){
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
			double B = calcB(data[k].D[g]);
			double weightF = B / (1.0 - data[k].D[g] + B * epsilon[k]) * eps1MinusEps[k];
			double weightJacobian = weightF * (oneMinus2Eps[k] - weightF);
			weightF *=  P_bbar_given_d_oldTheta[g];
			weightJacobian *=  P_bbar_given_d_oldTheta[g];

			//add to model
			models.addToFandJacobian(data[k], weightF, weightJacobian);
		}
	}

	//delete tmp variables
	delete[] eps1MinusEps;
	delete[] oneMinus2Eps;
};

void TRecalibrationEMSite::_addToJacobianAndFKnownGenotype(TSequencingErrorModels & models, double* epsilon){
	//fill F and Jacobian
	for(unsigned int k=0; k<numReads; ++k){
		//tmp variables
		double eps1MinusEps = epsilon[k] * (1.0 - epsilon[k]);
		double oneMinus2Eps = 1.0 - 2.0 * epsilon[k];

		//calc weights
		double B = calcB(data[k].D[trueBase]);
		double weightF = B / (1.0 - data[k].D[trueBase] + B * epsilon[k]) * eps1MinusEps;
		double weightJacobian = weightF * (oneMinus2Eps - weightF);

		//add to model
		models.addToFandJacobian(data[k], weightF, weightJacobian);
	}
};

void TRecalibrationEMSite::addToJacobianAndF(TSequencingErrorModels & models, double* epsilon){
	//calculate tmpEpsilon with current parameters
	calcEpsilon(models, epsilon);

	if(trueBase == N){
		_addToJacobianAndF(models, epsilon);
	} else {
		_addToJacobianAndFKnownGenotype(models, epsilon);
	}
};
*/


//---------------------------------------------------------------
//TRecalibrationEMWindow
//---------------------------------------------------------------
TRecalibrationEMWindow::TRecalibrationEMWindow(const TBaseData & baseFreqs, BAM::TReadGroupMap* ReadGroupMap){
	freqs[A] = baseFreqs.at(A);
	freqs[C] = baseFreqs.at(C);
	freqs[G] = baseFreqs.at(G);
	freqs[T] = baseFreqs.at(T);

	readGroupMapObject = ReadGroupMap;
}

uint32_t TRecalibrationEMWindow::getMaxDepth(){
	uint32_t maxDepth = 0;
	for(auto& s : sites){
		if(maxDepth < s.depth())
			maxDepth = s.depth();
	}
	return maxDepth;
};

void TRecalibrationEMWindow::addSite(TSite & site){
	if(!site.empty())
		sites.emplace_back(site);
};

size_t TRecalibrationEMWindow::size(){
	return sites.size();
};

size_t TRecalibrationEMWindow::numSitesDepthTwoOrMore(){
	size_t _numSites = 0;
	for(auto& s : sites){
		if(s.depth() > 1)
		++_numSites;
	}
	return _numSites;
};

void TRecalibrationEMWindow::addToDataTable(TRecalibrationEMDataTables & dataTables){
	for(auto& s : sites)
		dataTables.add(s);
};

uint64_t TRecalibrationEMWindow::cumulativeDepth(){
	uint64_t cumulDepth = 0;
	for(auto& s : sites){
		cumulDepth += s.depth();
	}
	return cumulDepth;
};

/*
double TRecalibrationEMWindow::fill_P_g_given_d_beta_AND_calcLL(TSequencingErrorModels & models, double* tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites){
		LL += site->fill_P_g_given_d_beta_AND_calcLL(models, freqs, tmpEpsilon);
	}
	return LL;
};

double TRecalibrationEMWindow::calcLL(TSequencingErrorModels & models, double* tmpEpsilon){
	double LL = 0.0;
	for(TRecalibrationEMSite* site : sites)
		LL += site->calcLL(models, freqs, tmpEpsilon);
	return LL;
};

void TRecalibrationEMWindow::addToQ(TSequencingErrorModels & models){
	for(TRecalibrationEMSite* site : sites)
		site->addToQ(models);
};

void TRecalibrationEMWindow::addToJacobianAndF(TSequencingErrorModels & models, double* tmpEpsilon){
	for(TRecalibrationEMSite* site : sites)
		site->addToJacobianAndF(models, tmpEpsilon);
};

*/

void TRecalibrationEMWindow::setEuqalBaseFrequencies(){
	for(int i=0; i<4; ++i) freqs[i] = 0.25;
};

//---------------------------------------------------------------
//TRecalibrationEMEstimator
//---------------------------------------------------------------
TRecalibrationEMEstimator::TRecalibrationEMEstimator(TParameters & args, BAM::TReadGroups* ReadGroups, TLog* Logfile, BAM::TReadGroupMap* ReadGroupMap){
	_logfile = Logfile;
	maxDepth = -1;

	//read groups
	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMap;

	//recal models
	_logfile->startIndent("Settings regarding the EM algorithm:");
	std::string recalFile = args.getParameterString("recal", false);
	std::string modelTagForEstimation = args.getParameterStringWithDefault("model", "quality=polynomial(2);position=polynomial(2);context=specific");

	//initialize from file?
	if(recalFile.empty()){
		_logfile->list("Will fit the model '" + modelTagForEstimation + "' for all read groups.");
	} else {
		_logfile->list("Will read models and initial values from file '" + recalFile + "'.");
		_logfile->list("Will fit the model '" + modelTagForEstimation + "' for read groups not in file '" + recalFile + "'.");
	}

	//parse model string
	std::string error;
	if(!covariateDefitionForEstimation.parse(modelTagForEstimation, error)){
		throw error + "!";
	}

	//estimation parameters
	minRequiredObservations = 10000; //constant for reporting
	numEMIterations = args.getParameterIntWithDefault("iterations", 200);
	_logfile->list("Will perform at max " + toString(numEMIterations) + " EM iterations.");
	maxEpsilon = args.getParameterDoubleWithDefault("maxEps", 0.000001);
	_logfile->list("Will stop EM when deltaLL < " + toString(maxEpsilon));
	NewtonRaphsonNumIterations = args.getParameterIntWithDefault("NRiterations", 20);
	_logfile->list("Will conduct at max " + toString(NewtonRaphsonNumIterations) + " Newton-Raphson iterations");
	NewtonRaphsonMaxF = args.getParameterDoubleWithDefault("maxF", 0.0001);
	_logfile->list("Will stop Newton-Raphson when F < " + toString(NewtonRaphsonMaxF));
	equalBaseFrequencies = true;
	if(args.parameterExists("estimateBaseFrequencies")){
		equalBaseFrequencies = false;
		_logfile->list("Will estimate the base frequencies. (parameter ''estimateBaseFrequencies)");
	} else if(equalBaseFrequencies)
		_logfile->list("Will assume equal base frequencies {0.25, 0.25, 0.25, 0.25}. (use 'estimateBaseFrequencies' to estimate them)");
	_logfile->endIndent();
};


void TRecalibrationEMEstimator::_initializeModels(){
	//count data available for recal
	_logfile->listFlush("Counting data available for recal ...");
	dataTables.init(_readGroups->size(), 500, 1000, 500);
	addToDataTable(dataTables);
	int numSitesWithData = numSites();
	_logfile->done();

	_logfile->conclude("Number of sites with data: " + toString(numSitesWithData));
	_logfile->conclude("Number of sites with depth > 1: " + toString(numSitesDepthTwoOrMore()));
	_logfile->conclude("Number of observations: " + toString(dataTables._totalCounts));
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
	_logfile->startIndent("Initializing recalibration models:");
	_sequencingErrorModels.prepareModelsForEstimation(_readGroups, _readGroupMap, _logfile);

	//first initialize from file, if given
	if(!recalFile.empty()){
		_sequencingErrorModels.addModelsFromFile(recalFile, &dataTables);
	}

	//now add default model for all others and check if existing (from file) match data range
	_logfile->listFlush("Initializing default models ...");
	int numModelsWithLittleData = 0;
	int numModelsWithoutData = 0;
	for(uint16_t rg = 0; rg < _readGroupMap->getNumReadGroups(); ++rg){
		for(int mate = 0; mate < 2; ++mate){
			TRecalibrationEMDataTable* table = dataTables.table(rg, mate);
			if(table->size() > 0){
				_sequencingErrorModels.addModel(rg, mate, covariateDefitionForEstimation, table);
				if(table->size() < minRequiredObservations)
					++numModelsWithLittleData;
			} else {
				if(_sequencingErrorModels.modelExists(rg, mate)){
					_sequencingErrorModels.removeModel(rg, mate);
				}
				++numModelsWithoutData;
			}
		}
	}
	_logfile->done();
	_logfile->conclude("Initialized " + toString(_sequencingErrorModels.numModels()) + " models.");

	//report warning in case of very little data
	if(numModelsWithLittleData > 0){
		_logfile->warning("Some read groups have very little data!");
		_logfile->startIndent("Consider merging these read groups:");

		for(size_t rg = 0; rg < _readGroups->size(); ++rg){
			int index = _readGroupMap->getIndex(rg);
			TRecalibrationEMDataTable* table = dataTables.table(index, false);
			if(table->size() > 0 && table->size() < minRequiredObservations)
				_logfile->list(_readGroups->getName(rg)  + " (first mate): only " + toString(table->size()) + " observations.");

			table = dataTables.table(index, true);
			if(table->size() > 0 && table->size() < minRequiredObservations)
				_logfile->list(_readGroups->getName(rg)  + " (second mate): only " + toString(table->size()) + " observations.");
		}
		_logfile->endIndent();
	}

	//report read groups without data
	if(numModelsWithoutData > 0){
		_logfile->startIndent("The following " + toString(numModelsWithoutData) + " read groups do not have data and will not be recalibrated:");
		_sequencingErrorModels.reportReadGroupsNotUsed();
		_logfile->endIndent();
	}
};

void TRecalibrationEMEstimator::initializeFromFile(const std::string string){
	_sequencingErrorModels.createModels(string, _readGroups, _readGroupMap, _logfile);
};

void TRecalibrationEMEstimator::performEstimation(std::string outputName, bool & writeTmpTables){
	//initialize models
	_initializeModels();

	//run EM
	_runEM(outputName, writeTmpTables);

	//writing final estimates
	std::string filename = outputName + "_recalibrationEM.txt";
	_logfile->listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename);
	_logfile->done();
};

void TRecalibrationEMEstimator::performEstimationKnownGenotypes(std::string outputName, bool & writeTmpTables){
	//initialize models
	_initializeModels();

	//run Newton-Raphson optimization
	_runNewtonRaphson();

	//writing final estimates
	std::string filename = outputName + "_recalibrationKnownGenotypes.txt";
	_logfile->listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename);
	_logfile->done();
};

void TRecalibrationEMEstimator::_calculate_EMWeights_epsilon(std::vector<TBaseData> & EMWeights){
	//make sure EM-weight storage is of appropriate size
	EMWeights.resize(dataTables.size());

	//loop over all bases and calculate EM-weights
	size_t index = 0;
	for(auto& w : windows){
		TBaseData& baseFreq = w.expectedBaseFrequencies();
		for(auto& s : w){
			for(auto& b : s){
				TBaseData noPMD;
				_sequencingErrorModels.calculateBaseLikelihoods(b, noPMD);
				_pmd.calculateBaseLikelihoods(b, noPMD, EMWeights[index]);
				EMWeights[index] *= baseFreq;
				EMWeights[index].normalize();

				//increment index
				++index;
			}
		}
	}
};

double TRecalibrationEMEstimator::_calculate_Q_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d){
	_sequencingErrorModels.setQToZero();

	//loop over all bases and add them to Q
	size_t index = 0;
	for(auto& w : windows){
		for(auto& s : w){
			for(auto& b : s){
				_sequencingErrorModels.addToQ(b, EM_weights_bbar_given_d[index]);
				++index;
			}
		}
	}

	//return total Q
	return _sequencingErrorModels.curQ();
};

void TRecalibrationEMEstimator::_calculate_J_F_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d){
	_logfile->listFlush("Calculating Jacobian and gradient ...");
	_sequencingErrorModels.setNewtonRaphsonParamsToZero();

	size_t index = 0;
	for(auto& w : windows){
		for(auto& s : w){
			for(auto& b : s){
				_sequencingErrorModels.addToFandJacobian(b, EM_weights_bbar_given_d[index]);
				++index;
			}
		}
	}

	//solve J^-1 x F
	_sequencingErrorModels.solveJxF();
	_logfile->done();
};



void TRecalibrationEMEstimator::_updateEM_theta_epsilon(){
	_logfile->startIndent("Updating sequencing error models (theta_epsilon):");

	// 1) calculate EM weights
	//-------------------------
	std::vector<TBaseData> EM_weights_bbar_given_d;
	_calculate_EMWeights_epsilon(EM_weights_bbar_given_d);

	// 2) update rho
	size_t index = 0;
	for(auto& w : windows){
		for(auto& s : w){
			for(auto& b : s){
				_sequencingErrorModels.addToFandJacobian(b, EM_weights_bbar_given_d[index]);
				++index;
			}
		}
	}

	// 3) Calculate Q_beta at current location
	double curQ = _calculate_Q_beta(EM_weights_bbar_given_d);

	// 4) Use Newton-Raphson to optimize
	_logfile->startNumbering("Optimizing Q_beta using a Newton-Raphson algorithm:");

	for(int i=0; i<NewtonRaphsonNumIterations; ++i){
		_logfile->startIndent("Running Newton-Raphson iteration " + toString(i+1) + ":");

		// a) fill Jacobin and F
		_calculate_J_F_beta(EM_weights_bbar_given_d);

		// b) update parameters using backtracking
		double lambda = 1.0;
		int log2_lambda = 0;
		int numUpdatedModels = 0;
		int numUpdatedModels_old;

		while(numUpdatedModels < _sequencingErrorModels.numModels() && lambda > 1.0E-20){
			//propose move
			_logfile->listFlush("Proposing move with log2(lambda) = " + toString(log2_lambda) + " ... ");
			_sequencingErrorModels.proposeNewParameters(lambda);

			//calculate Q at new location
			double Q = _calculate_Q_beta(EM_weights_bbar_given_d);

			//check if we accept or backtrack
			numUpdatedModels_old = numUpdatedModels;
			numUpdatedModels = _sequencingErrorModels.acceptProposedParametersBasedOnQ();
			_logfile->write(toString(numUpdatedModels) + "/" + toString(_sequencingErrorModels.numModels()) + " models converged.");

			if(numUpdatedModels > numUpdatedModels_old){
				_logfile->conclude("Q_beta was increased from " + toString(curQ) + " to " + toString(Q));
			}
			curQ = Q;

			//backtrack
			lambda = lambda / 2.0; //backtrack;
			--log2_lambda;
		}

		if(numUpdatedModels < _sequencingErrorModels.numModels()){
			_logfile->conclude("Some models did not improve even with log2(lambda) = " + toString(log2_lambda) + ", aborting Newton-Raphson.");
		}

		// c) adjust parameters post estimation
		_sequencingErrorModels.adjustParametersPostEstimation();

		// d) get largest gradient (F) to check if we break NR optimization
		double maxF = _sequencingErrorModels.getSteepestGradient();
		_logfile->conclude("max(F) = " + toString(maxF));
		_logfile->endIndent();
		if(maxF < NewtonRaphsonMaxF || numUpdatedModels == 0) break;
	}
	_logfile->endIndent();
};

void TRecalibrationEMEstimator::_calculateLL_fullModel(){

	size_t index = 0;
	for(auto& w : windows){
		TGenotypeProbabilities& genoFreq = w.expectedBaseFrequencies();
		for(auto& s : w){
			TGenotypeLikelihoods genotypeLikelihoods;


			for(auto& b : s){
				TBaseData noPMD, withPMD;
				_sequencingErrorModels.calculateBaseLikelihoods(b, noPMD);
				_pmd.calculateBaseLikelihoods(b, noPMD, withPMD);



				//increment index
				++index;
			}
		}
	}

};

void TRecalibrationEMEstimator::_runEM(std::string outputName, bool & writeTmpTables){
	//run EM
	_logfile->startNumbering("Running EM algorithm:");

	double LL, deltaLL, oldLL = 0.0;
	std::ofstream out;
	std::string filename;

	//running iterations
	for(int iter = 0; iter < numEMIterations; ++iter){
		_logfile->number("EM Iteration:"); _logfile->addIndent();

		//update theta_epsilon (sequencing errors)
		_updateEM_theta_epsilon();

		//calculate LL



		//calculate EM weights P(g|d, theta_eps') for all sites and calculate LL
		LL = 0.0;
		_logfile->listFlush("Calculating P(g|d, theta_eps') ...");
		for(TRecalibrationEMWindow* curWindow : windows)
			LL += curWindow->fill_P_g_given_d_beta_AND_calcLL(_sequencingErrorModels, tmpEpsilon);
		_logfile->done();
		_logfile->conclude("Current Log Likelihood = " + toString(LL));

		//check if we break based on LL
		if(iter > 0){
			deltaLL = LL - oldLL;
			_logfile->conclude("Epsilon = " + toString(deltaLL));
			if(iter > 0 && deltaLL < maxEpsilon){
				_logfile->conclude("EM has converged (tmpEpsilon < " + toString(maxEpsilon) + ")");
				break;
			} else oldLL = LL;
		} else oldLL = LL;

		//run NewtonRaphson until convergence
		_runNewtonRaphson();

		//write current estimates to file
		if(writeTmpTables){
			filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
			_logfile->listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename);
			_logfile->done();
		}

		//end loop
		_logfile->endIndent();
		if(iter == numEMIterations - 1) _logfile->warning("EM has not converged after maximum number of iterations!");
	}

	//finalize
	_logfile->endNumbering();
};

void TRecalibrationEMEstimator::_runNewtonRaphson(){
	//calculate Q at current location
	_sequencingErrorModels.setQToZero();
	for(TRecalibrationEMWindow* curWindow : windows)
		curWindow->addToQ(_sequencingErrorModels);
	double curQ = modTSiteStorageTSiteStorageTSiteStorageels.curQ();

	//run up to maxNewtonRaphsonIteratios iterations, but stop if max(F) < maxFThreshold
	_logfile->startIndent("Running Newton-Raphson optimization:");
	for(int i=0; i<NewtonRaphsonNumIterations; ++i){
		_logfile->startIndent("Running Newton-Raphson iteration " + toString(i+1) + ":");

		//fill Jacobin and F: loop over all sites
		_logfile->listFlush("Calculating Jacobian and gradient ...");
		_sequencingErrorModels.setNewtonRaphsonParamsToZero();
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow)
			(*curWindow)->addToJacobianAndF(_sequencingErrorModels, tmpEpsilon);

		//now solve J^-1 x F
		_sequencingErrorModels.solveJxF();
		_logfile->done();

		//update params for each read group using backtracking
		double lambda = 1.0;
		int log2_lambda = 0;
		int numUpdatedModels = 0;
		int numUpdatedModels_old;

		while(numUpdatedModels < _sequencingErrorModels.numModels() && lambda > 1.0E-20){
			//propose move
			_logfile->listFlush("Proposing move with log2(lambda) = " + toString(log2_lambda) + " ... ");
			_sequencingErrorModels.proposeNewParameters(lambda);

			//calculate Q at new location
			_sequencingErrorModels.setQToZero();
			for(TRecalibrationEMWindow* curWindow : windows)
				curWindow->addToQ(_sequencingErrorModels);

			//check if we accept or backtrack
			numUpdatedModels_old = numUpdatedModels;
			numUpdatedModels = _sequencingErrorModels.acceptProposedParametersBasedOnQ();
			double Q = _sequencingErrorModels.curQ();

			_logfile->write(toString(numUpdatedModels) + "/" + toString(_sequencingErrorModels.numModels()) + " models converged.");

			if(numUpdatedModels > numUpdatedModels_old){
				_logfile->conclude("Q was increased from " + toString(curQ) + " to " + toString(Q));
			}
			curQ = Q;

			//backtrack
			lambda = lambda / 2.0; //backtrack;
			--log2_lambda;
		}

		if(numUpdatedModels < _sequencingErrorModels.numModels()){
			_logfile->conclude("Some models did not improve even with log2(lambda) = " + toString(log2_lambda) + ", aborting Newton-Raphson.");
		}

		//adjust parameters post estimation
		_sequencingErrorModels.adjustParametersPostEstimation();

		//get largest gradient (F) to check if we break NR optimization
		double maxF = _sequencingErrorModels.getSteepestGradient();
		_logfile->conclude("max(F) = " + toString(maxF));
		_logfile->endIndent();
		if(maxF < NewtonRaphsonMaxF || numUpdatedModels == 0) break;
	}
	_logfile->endIndent();
};

void TRecalibrationEMEstimator::addNewWindow(const TBaseData & freqs){
	windows.emplace_back(freqs, _readGroupMap);
	//set iterator
	curWindow = windows.end(); --curWindow;
	if(equalBaseFrequencies){
		(*curWindow)->setEuqalBaseFrequencies();
	}
};

void TRecalibrationEMEstimator::addSite(TSite & site){
	(*curWindow)->addSite(site);
};

void TRecalibrationEMEstimator::addSite(TSite & site, TQualityMap & qualiMap, const Base TrueBase){
	(*curWindow)->addSite(site, qualiMap, TrueBase);
};

long TRecalibrationEMEstimator::numSites(){
	long _numSites = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		_numSites += curWindow->size();
	return _numSites;
};

long TRecalibrationEMEstimator::numSitesDepthTwoOrMore(){
	long _numSites = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		_numSites += curWindow->numSitesDepthTwoOrMore();
	return _numSites;
};

void TRecalibrationEMEstimator::addToDataTable(TRecalibrationEMDataTables & dataTable){
	for(TRecalibrationEMWindow* curWindow : windows)
		curWindow->addToDataTable(dataTable);
};

long TRecalibrationEMEstimator::cumulativeDepth(){
	long cumulDepth = 0;
	for(TRecalibrationEMWindow* curWindow : windows)
		cumulDepth += curWindow->cumulativeDepth();
	return cumulDepth;
};

void TRecalibrationEMEstimator::writeCurrentEstimates(const std::string filename){
	_sequencingErrorModels.writeRecalFile(filename);
};

double TRecalibrationEMEstimator::calcLL(){
	_initializTmpVariablesForEstimation();
	double LL = 0.0;
	for(TRecalibrationEMWindow* curWindow : windows)
		LL += curWindow->calcLL(_sequencingErrorModels, tmpEpsilon);
	return LL;
};

}; //end namespace
