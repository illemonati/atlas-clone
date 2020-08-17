/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMEstimator.h"

namespace GenotypeLikelihoods{


//---------------------------------------------------------------
//TRecalibrationEMEstimator
//---------------------------------------------------------------
TRecalibrationEMEstimator::TRecalibrationEMEstimator(TParameters & args, BAM::TReadGroups* ReadGroups, TLog* Logfile, BAM::TReadGroupMap* ReadGroupMap){
	_logfile = Logfile;

	//read groups
	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMap;

	//genotype distribution: currently only allow for haploid
	_genoDist = std::make_unique<TGenotypeDistribution_haploid>();

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
	if(!_modelDefitionForEstimation.parseCovariates(modelTagForEstimation, error)){
		throw error + "!";
	}

	//estimation parameters
	_minRequiredObservations = 10000; //constant for reporting
	_numEMIterations = args.getParameterIntWithDefault("iterations", 200);
	_logfile->list("Will perform at max " + toString(_numEMIterations) + " EM iterations.");
	_minDeltaLL = args.getParameterDoubleWithDefault("minDeltaLL", 0.000001);
	_logfile->list("Will stop EM when deltaLL < " + toString(_minDeltaLL));
	_NewtonRaphsonNumIterations = args.getParameterIntWithDefault("NRiterations", 20);
	_logfile->list("Will conduct at max " + toString(_NewtonRaphsonNumIterations) + " Newton-Raphson iterations");
	_NewtonRaphsonMaxF = args.getParameterDoubleWithDefault("maxF", 0.0001);
	_logfile->list("Will stop Newton-Raphson when F < " + toString(_NewtonRaphsonMaxF));

	//base frequency model
	equalBaseFrequencies = true;
	if(args.parameterExists("estimateBaseFrequencies")){
		equalBaseFrequencies = false;
		_logfile->list("Will estimate the base frequencies. (parameter ''estimateBaseFrequencies)");

		//TODO: implement estimation of genotype distribution!
		throw std::runtime_error("Estimation of genotype distribution not yet implemented!");


	} else if(equalBaseFrequencies){
		_logfile->list("Will assume equal base frequencies {0.25, 0.25, 0.25, 0.25}. (use 'estimateBaseFrequencies' to estimate them)");
		_genoDist->reset();
	}
	_logfile->endIndent();
};

void TRecalibrationEMEstimator::_initializeModels(){
	//count data available for recal
	_logfile->listFlush("Counting data available for recal ...");
	_dataTables.init(_readGroups->size(), 500, 1000, 500);
	addToDataTable(_dataTables);
	int numSitesWithData = numSites();
	_logfile->done();

	_logfile->conclude("Number of sites with data: " + toString(numSitesWithData));
	_logfile->conclude("Number of sites with depth > 1: " + toString(numSitesDepthTwoOrMore()));
	_logfile->conclude("Number of observations: " + toString(_dataTables.size()));
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
	if(!_recalFile.empty()){
		_sequencingErrorModels.addModelsFromFile(_recalFile, &_dataTables);
	}

	//now add default model for all others and check if existing (from file) match data range
	_logfile->listFlush("Initializing default models ...");
	int numModelsWithLittleData = 0;
	int numModelsWithoutData = 0;
	for(uint16_t rg = 0; rg < _readGroupMap->getNumReadGroups(); ++rg){
		for(int mate = 0; mate < 2; ++mate){
			TRecalibrationEMDataTable* table = _dataTables.table(rg, mate);
			if(table->size() > 0){
				_modelDefitionForEstimation.readGroupId = rg;
				_modelDefitionForEstimation.isSecondMate = mate;
				_sequencingErrorModels.addModel(_modelDefitionForEstimation, table);
				if(table->size() < _minRequiredObservations)
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
			TRecalibrationEMDataTable* table = _dataTables.table(index, false);
			if(table->size() > 0 && table->size() < _minRequiredObservations)
				_logfile->list(_readGroups->getName(rg)  + " (first mate): only " + toString(table->size()) + " observations.");

			table = _dataTables.table(index, true);
			if(table->size() > 0 && table->size() < _minRequiredObservations)
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

//----------------------------
//Functions to add data
//----------------------------
void TRecalibrationEMEstimator::addSite(const TSite & site){
	if(!site.empty()){
		_sites.emplace_back(site);
	}
};

size_t TRecalibrationEMEstimator::numSites(){
	return _sites.size();
};

size_t TRecalibrationEMEstimator::numSitesDepthTwoOrMore(){
	size_t _numSites = 0;
	for(auto& s : _sites){
		if(s.depth() > 1)
		++_numSites;
	}
	return _numSites;
};

void TRecalibrationEMEstimator::addToDataTable(TRecalibrationEMDataTables & dataTable){
	for(auto& s : _sites){
		_dataTables.add(s);
	}
};

uint64_t TRecalibrationEMEstimator::cumulativeDepth(){
	uint64_t cumulDepth = 0;
	for(auto& s : _sites){
		cumulDepth += s.depth();
	}
	return cumulDepth;
};

//----------------------------
//Functions for estimation
//----------------------------
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

void TRecalibrationEMEstimator::_fillRelevantBaseFrequencies(TBaseData & baseFreq, const Genotype genotype){
	if(genotype == NN){
		baseFreq = _genoDist->baseFrequencies();
	} else {
		_genoDist->fillBaseFrequences(baseFreq, genotype);
	}
};

void TRecalibrationEMEstimator::_calculate_EMWeights_epsilon(std::vector<TBaseData> & EMWeights){
	//make sure EM-weight storage is of appropriate size
	EMWeights.resize(_dataTables.size());
	//loop over all bases and calculate EM-weights
	size_t index = 0;
	for(auto& s : _sites){
		//get relevant base frequencies: from known genotype or distribution if genotype is unknown
		TBaseData baseFreq;
		_fillRelevantBaseFrequencies(baseFreq, s.genotype());

		//calculate weights per base
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
};

double TRecalibrationEMEstimator::_calculate_Q_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d){
	_sequencingErrorModels.setQToZero();

	//loop over all bases and add them to Q
	size_t index = 0;
	for(auto& s : _sites){
		for(auto& b : s){
			_sequencingErrorModels.addToQ(b, EM_weights_bbar_given_d[index]);
			++index;
		}
	}

	//return total Q
	return _sequencingErrorModels.curQ();
};

void TRecalibrationEMEstimator::_calculate_J_F_beta(const std::vector<TBaseData> & EM_weights_bbar_given_d){
	_logfile->listFlush("Calculating Jacobian and gradient ...");
	_sequencingErrorModels.setNewtonRaphsonParamsToZero();

	size_t index = 0;
	for(auto& s : _sites){
		for(auto& b : s){
			_sequencingErrorModels.addToFandJacobian(b, EM_weights_bbar_given_d[index]);
			++index;
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
	_logfile->listFlushDots("Calculating EM weights");
	std::vector<TBaseData> EM_weights_bbar_given_d;
	_calculate_EMWeights_epsilon(EM_weights_bbar_given_d);
	_logfile->done();

	// 2) update rho
	//-------------------------
	_logfile->listFlushDots("Updating rho");
	size_t index = 0;
	for(auto& s : _sites){
		for(auto& b : s){
			_sequencingErrorModels.addToFandJacobian(b, EM_weights_bbar_given_d[index]);
			++index;
		}
	}
	_logfile->done();

	// 3) Calculate Q_beta at current location
	//-------------------------
	_logfile->listFlushDots("Calculating Q_beta at current location");
	double curQ = _calculate_Q_beta(EM_weights_bbar_given_d);
	_logfile->done();
	_logfile->conclude("Q_beta = ", curQ);

	// 4) Use Newton-Raphson to optimize
	//-------------------------
	_logfile->startNumbering("Optimizing Q_beta using a Newton-Raphson algorithm:");

	for(int i=0; i<_NewtonRaphsonNumIterations; ++i){
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
		if(maxF < _NewtonRaphsonMaxF || numUpdatedModels == 0) break;

		_logfile->endIndent();
	}
	_logfile->endNumbering();
	_logfile->endIndent();
};

double TRecalibrationEMEstimator::_calculateLL_fullModel(){
	double LL = 0.0;
	const TGenotypeProbabilities& genoFreq = _genoDist->genotypeFrequencies();
	for(auto& s : _sites){
		//calculate genotype likelihoods
		TGenotypeLikelihoods genotypeLikelihoods;
		calculateGenotypeLikelihoods(s, genotypeLikelihoods);

		//weight by genotype prior
		if(s.genotype() == NN){
			LL += log(genotypeLikelihoods.weightedSum(genoFreq));
		} else {
			LL += log(genotypeLikelihoods[s.genotype()]);
		}
	}

	return LL;
};

void TRecalibrationEMEstimator::_runEM(std::string outputName, bool & writeTmpTables){
	//run EM
	_logfile->startNumbering("Running EM algorithm:");

	std::ofstream out;
	std::string filename;

	//calculate initial LL
	double oldLL = _calculateLL_fullModel();
	_logfile->conclude("Initial log Likelihood = " + toString(oldLL));

	//running iterations
	for(int iter = 0; iter < _numEMIterations; ++iter){
		_logfile->number("EM Iteration:"); _logfile->addIndent();

		//update theta_epsilon (sequencing errors)
		_updateEM_theta_epsilon();

		//TODO: update other parameters

		//calculate LL
		double LL = _calculateLL_fullModel();
		_logfile->conclude("Current Log Likelihood = " + toString(LL));

		//check if we break based on LL
		double deltaLL = LL - oldLL;
		_logfile->conclude("delta LL = " + toString(deltaLL));
		if(iter > 0 && deltaLL < _minDeltaLL){
			_logfile->conclude("EM has converged (delta LL < " + toString(_minDeltaLL) + ")");
			break;
		}
		oldLL = LL;

		//write current estimates to file
		if(writeTmpTables){
			filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
			_logfile->listFlush("Writing current estimates to file '" + filename + "' ...");
			writeCurrentEstimates(filename);
			_logfile->done();
		}

		//end loop
		_logfile->endIndent();
		if(iter == _numEMIterations - 1) _logfile->warning("EM has not converged after maximum number of iterations!");
	}

	//finalize
	_logfile->endNumbering();
};

//----------------------------
// Other functions
//----------------------------
void TRecalibrationEMEstimator::writeCurrentEstimates(const std::string filename){
	_sequencingErrorModels.writeRecalFile(filename);
};

double TRecalibrationEMEstimator::calcLL(){
	throw std::runtime_error("double TRecalibrationEMEstimator::calcLL() not yet implemented!");
};

void TRecalibrationEMEstimator::calcLikelihoodSurface(std::string filename, int numMarginalGridPoints){
	throw std::runtime_error("void TRecalibrationEMEstimator::calcLikelihoodSurface(std::string filename, int numMarginalGridPoints) not yet implemented!");
};

void TRecalibrationEMEstimator::calcQSurface(std::string filename, int numMarginalGridPoints){
	throw std::runtime_error("TRecalibrationEMEstimator::calcQSurface(std::string filename, int numMarginalGridPoints) not yet implemented!");
};


}; //end namespace
