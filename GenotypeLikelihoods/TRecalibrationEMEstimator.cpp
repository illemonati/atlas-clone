/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMEstimator.h"

namespace GenotypeLikelihoods {

namespace RecalEstimator{

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
//--------------------------------------------------------------------------

TSequencingErrorModelVectorForEstimation::TSequencingErrorModelVectorForEstimation(TSequencingErrorModels & SequencingErrorModels,
																				   const RecalEstimatorTools::TRecalDataTables & DataTables,
																				   const BAM::TReadGroups & ReadGroups,
																				   const BAM::TReadGroupMap & ReadGroupMap,
																				   uint32_t MinRequiredObservations,
																				   TLog* Logfile):
																					   _modelIndex(ReadGroups){
	//Copy models that are 1) in use after pooling and 2) have data.
	//Note: data table is already pooled!

	//prepare storage for reporting
	RecalEstimatorTools::TModelStati modelStati;

	//copy models
	for(auto& r : ReadGroupMap.readGroupsInUse()){
		//add to model stati
		modelStati.add(r);

		//loop over mates
		for(uint8_t mate = 0; mate < 2; ++mate){

			const RecalEstimatorTools::TRecalDataTable& table = DataTables[r][mate];

			//check if there is sufficient data for _first mate
			if(table.size() > 0){
				//check if model is estimatable
				if(SequencingErrorModels[r][mate].estimatable()){

					//copy model and update index
					std::shared_ptr<TSequencingErrorModelRecal>& model = SequencingErrorModels[r].getSharedPointerToRecalModel(mate);
					_models.push_back(model);
					modelStati[r][RecalEstimatorTools::copied].set(mate);
					_modelIndex.set(r, mate, model, ReadGroupMap); //handles pooling

					//check if there is limited data
					if(table.size() < MinRequiredObservations){
						modelStati[r][RecalEstimatorTools::littleData].set(mate);
					}

				} else {
					modelStati[r][RecalEstimatorTools::dataButNoRecal].set(mate);
				}

			} else {
				if(SequencingErrorModels[r][mate].estimatable()){
					modelStati[r][RecalEstimatorTools::noData].set(mate);
				}
			}
		}
	}

	//report models that will be estimated
	modelStati.report(RecalEstimatorTools::copied, "Read groups for which models will be estimated:", ReadGroups, Logfile);
	modelStati.report(RecalEstimatorTools::noData, "Read groups excluded because they have no data:", ReadGroups, Logfile);
	modelStati.report(RecalEstimatorTools::dataButNoRecal, "Read groups with data but no recal model:", ReadGroups, Logfile);
	if(modelStati.num(RecalEstimatorTools::copied) == 0){
		throw "No recal models need estimation!";
	}
	modelStati.report(RecalEstimatorTools::littleData, "Read groups with very little data (consider pooling):", ReadGroups, Logfile);
};

void TSequencingErrorModelVectorForEstimation::fillBaseLikelihoods(const BAM::TSequencedBase & base,  TBaseLikelihoods & baseLikelihoods) const{
	_modelIndex(base).fillBaseLikelihoods(base, baseLikelihoods);
};

//functions to estimate rho
//-------------------------------------------------------------------
//functions to estimate rho
void TSequencingErrorModelVectorForEstimation::prepareRhoEstimationFromEMWeights(){
	for(auto& model : _models){
		model->prepareRhoEstimationFromEMWeights();
	}
};

void TSequencingErrorModelVectorForEstimation::addBaseForRhoEstimation(BAM::TSequencedBase & base, const TBaseLikelihoods & EMWeights){
	_modelIndex(base).addBaseForRhoEstimation(base, EMWeights);
};

void TSequencingErrorModelVectorForEstimation::estimateRho(){
	for(auto& model : _models){
		model->estimateRho();
	}
};

//functions to estimate beta
//-------------------------------------------------------------------

void TSequencingErrorModelVectorForEstimation::addToFandJacobian(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d){
	_modelIndex(base).addToFandJacobian(base, EM_weights_bbar_given_d);
};

void TSequencingErrorModelVectorForEstimation::addToQ(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d){
	_modelIndex(base).addToQ(base, EM_weights_bbar_given_d);
};

void TSequencingErrorModelVectorForEstimation::setNewtonRaphsonParamsToZero(){
	for(auto& model : _models){
		model->setNewtonRaphsonParamsToZero();
	}
};

void TSequencingErrorModelVectorForEstimation::setQToZero(){
	for(auto& model : _models){
		model->setQToZero();
	}
};

double TSequencingErrorModelVectorForEstimation::curQ(){
	double Q = 0.0;
	for(auto& model : _models){
		Q += model->curQ();
	}
	return Q;
};

bool TSequencingErrorModelVectorForEstimation::solveJxF(){
	bool couldSolve = true;
	for(auto& model : _models){
		if(!model->solveJxF()){
			model->printJacobianToStdOut();
			throw "Issue solving JxF! This may be due to a lack of data. Consider adding more sites.";
			couldSolve = false;
		}
	}
	return couldSolve;
};

void TSequencingErrorModelVectorForEstimation::proposeNewParameters(double lambda){
	for(auto& model : _models){
		model->proposeNewParameters(lambda);
	}
};

unsigned int TSequencingErrorModelVectorForEstimation::acceptProposedParametersBasedOnQ(){
	unsigned int numAccepted = 0;
	for(auto& model : _models){
		numAccepted += (unsigned int) model->acceptProposedParametersBasedOnQ();
	}
	return numAccepted;
};

void TSequencingErrorModelVectorForEstimation::adjustParametersPostEstimation(){
	for(auto& model : _models){
		model->adjustParametersPostEstimation();
	}
};

double TSequencingErrorModelVectorForEstimation::getSteepestGradient(){
	double maxF = 0.0;
	for(auto& model : _models){
		double tmp = model->getSteepestGradient();
		if(fabs(tmp) > maxF) maxF = fabs(tmp);
	}
	return maxF;
};

void TSequencingErrorModelVectorForEstimation::print(){
	for(auto& model : _models){
		std::cout << model->getCovariateDefinition() << "\t" << model->getRhoDefinition() << std::endl;
	}
};

//---------------------------------------------------------------
//TRecalibrationEMEstimator
//---------------------------------------------------------------
TRecalibrationEMEstimator::TRecalibrationEMEstimator(TParameters & args, TLog* Logfile, const BAM::TReadGroups* ReadGroups, const BAM::TReadGroupMap* ReadGroupMap){
	_logfile = Logfile;

	//read groups
	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMap;

	//genotype distribution: currently only allow for haploid
	_genoDist = std::make_unique<TGenotypeDistribution_haploid>();

	//estimation parameters
	_logfile->startIndent("Settings regarding the EM algorithm:");

	_minRequiredObservations = 10000; //constant for reporting
	_numEMIterations = args.getParameterWithDefault<int>("iterations", 200);
	_logfile->list("Will perform at max ", _numEMIterations, " EM iterations.");
	_minDeltaLL = args.getParameterWithDefault<double>("minDeltaLL", 0.000001);
	_logfile->list("Will stop EM when deltaLL < ", _minDeltaLL, ".");
	_NewtonRaphsonNumIterations = args.getParameterWithDefault<int>("NRiterations", 20);
	_logfile->list("Will conduct at max ", _NewtonRaphsonNumIterations, " Newton-Raphson iterations");
	_NewtonRaphsonMaxF = args.getParameterWithDefault<double>("maxF", 0.0001);
	_logfile->list("Will stop Newton-Raphson when F < ", _NewtonRaphsonMaxF, ".");

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

	//write tmp tables?
	if(args.parameterExists("writeTmpTables")){
		_writeTmpTables = true;
		_logfile->list("Will write intermediate estimates of EM and Newton-Raphson to file. (parameter 'writeTmpTables')");
	} else {
		_writeTmpTables = false;
		_logfile->list("Will not write intermediate estimates. (use 'writeTmpTables' to get this debug information)");
	}

	_logfile->endIndent();
};

//-----------------------------------------------
// Function to initialize models used for estimation
//-----------------------------------------------
size_t TRecalibrationEMEstimator::_numSitesDepthTwoOrMore(){
	size_t _numSites = 0;
	for(auto& s : _sites){
		if(s.depth() > 1)
		++_numSites;
	}
	return _numSites;
};

void TRecalibrationEMEstimator::_initializeModels(TSequencingErrorModels & SequencingErrorModels){
	//count data available for recal
	_logfile->listFlush("Counting data available for recal ...");
	//Note: data tables pool read groups!
	_dataTables.initialize(_readGroups, _readGroupMap);
	_dataTables.add(_sites);
	_logfile->done();

	_logfile->conclude("Number of sites with data: " + toString(_sites.size()));
	size_t numSitesDepthTwoOrMore = _numSitesDepthTwoOrMore();
	_logfile->conclude("Number of sites with depth > 1: " + toString(numSitesDepthTwoOrMore));
	_logfile->conclude("Number of bases: " + toString(_dataTables.size()));
	if(numSitesDepthTwoOrMore < 100) throw "Less than 100 sites with depth >= 2 available - aborting estimation!";

	//identify models with data that can be estimated
	_logfile->startIndent("Identifying models to estimate:");
	_modelsToEstimate = std::make_unique<TSequencingErrorModelVectorForEstimation>(SequencingErrorModels, _dataTables, *_readGroups, *_readGroupMap, _minRequiredObservations, _logfile);
	_logfile->endIndent();
};

//----------------------------
//Functions to add data
//----------------------------
void TRecalibrationEMEstimator::addSite(const TSite & site){
	if(!site.empty()){
		_sites.emplace_back(site);
	}
};

//----------------------------
//Functions for estimation
//----------------------------
void TRecalibrationEMEstimator::performEstimation(std::string outputName, TSequencingErrorModels & SequencingErrorModels, const TPostMortemDamage & PmdModels){
	//initialize models
	_initializeModels(SequencingErrorModels);

	//run EM
	_runEM(outputName, PmdModels);

	//writing final estimates
	std::string filename = outputName + "_recalibrationEM.txt";
	_logfile->listFlush("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename);
	_logfile->done();
};

void TRecalibrationEMEstimator::_fillRelevantBaseFrequencies(TBaseProbabilities & baseFreq, const genometools::Genotype & genotype){
	if(genotype == genometools::NN){
		baseFreq = _genoDist->baseFrequencies();
	} else {
		_genoDist->fillBaseFrequences(baseFreq, genotype);
	}
};

void TRecalibrationEMEstimator::_calculate_EMWeights_epsilon(std::vector<TBaseLikelihoods> & EMWeights, const TPostMortemDamage & PmdModels){
	//make sure EM-weight storage is of appropriate size
	EMWeights.resize(_dataTables.size());
	//loop over all bases and calculate EM-weights
	size_t index = 0;
	for(auto& s : _sites){
		//get relevant base frequencies P(b): from known genotype or distribution if genotype is unknown
		TBaseProbabilities baseFreq;
		_fillRelevantBaseFrequencies(baseFreq, s.genotype());

		//calculate weights per base
		for(auto& b : s){
			//calculate P(bbar) = \sum_b P(bbar|b)P(b)
			TBaseLikelihoods PMD;
			PmdModels.fillBaseLikelihoods(b, baseFreq, PMD);

			//calculate P(d|bbar)
			_modelsToEstimate->fillBaseLikelihoods(b, EMWeights[index]);

			//calculate P(d|bbar) \propto P(d|bbar)P(bbar)
			EMWeights[index] *= PMD;
			EMWeights[index].normalize();

			//increment index
			++index;
		}
	}
};

double TRecalibrationEMEstimator::_calculate_Q_beta(const std::vector<TBaseLikelihoods> & EM_weights_bbar_given_d){
	_modelsToEstimate->setQToZero();

	//loop over all bases and add them to Q
	size_t index = 0;
	for(auto& s : _sites){
		for(auto& b : s){
			_modelsToEstimate->addToQ(b, EM_weights_bbar_given_d[index]);
			++index;
		}
	}

	//return total Q
	return _modelsToEstimate->curQ();
};

void TRecalibrationEMEstimator::_calculate_J_F_beta(const std::vector<TBaseLikelihoods> & EM_weights_bbar_given_d){
	_logfile->listFlush("Calculating Jacobian and gradient ...");
	_modelsToEstimate->setNewtonRaphsonParamsToZero();

	size_t index = 0;
	for(auto& s : _sites){
		for(auto& b : s){
			_modelsToEstimate->addToFandJacobian(b, EM_weights_bbar_given_d[index]);
			++index;
		}
	}

	//solve J^-1 x F
	_modelsToEstimate->solveJxF();
	_logfile->done();
};

void TRecalibrationEMEstimator::_updateEM_theta_epsilon(const TPostMortemDamage & PmdModels){
	_logfile->startIndent("Updating sequencing error models (theta_epsilon):");

	// 1) calculate EM weights
	//-------------------------
	_logfile->listFlushDots("Calculating EM weights");
	std::vector<TBaseLikelihoods> EM_weights_bbar_given_d;
	_calculate_EMWeights_epsilon(EM_weights_bbar_given_d, PmdModels);
	_logfile->done();

	// 2) update rho
	//-------------------------
	/*
	_logfile->listFlushDots("Updating rho");
	_sequencingErrorModels.prepareRhoEstimationFromEMWeights();
	size_t index = 0;
	for(auto& s : _sites){
		for(auto& b : s){
			_sequencingErrorModels.addBaseForRhoEstimation(b, EM_weights_bbar_given_d[index]);
			++index;
		}
	}
	_sequencingErrorModels.estimateRho();
	_logfile->done();
	*/

	// 3) Calculate Q_beta at current location
	//-------------------------
	_logfile->listFlushDots("Calculating Q_beta at current location");
	double curQ = _calculate_Q_beta(EM_weights_bbar_given_d);
	_logfile->done();
	_logfile->conclude("Q_beta = ", curQ);

	// 4) Use Newton-Raphson to optimize
	//-------------------------
	_logfile->startIndent("Optimizing Q_beta using a Newton-Raphson algorithm:");

	for(int i=0; i<_NewtonRaphsonNumIterations; ++i){
		_logfile->startIndent("Running Newton-Raphson iteration " + toString(i+1) + ":");

		// a) fill Jacobin and F
		_calculate_J_F_beta(EM_weights_bbar_given_d);

		// b) update parameters using backtracking
		double lambda = 1.0;
		int log2_lambda = 0;
		size_t numUpdatedModels = 0;
		size_t numUpdatedModels_old;

		while(numUpdatedModels < _modelsToEstimate->size() && lambda > 1.0E-20){
			//propose move
			_logfile->listFlush("Proposing move with log2(lambda) = " + toString(log2_lambda) + " ... ");
			_modelsToEstimate->proposeNewParameters(lambda);

			std::cout << std::endl;
			_modelsToEstimate->print();
			std::cout << std::endl;

			//calculate Q at new location
			double Q = _calculate_Q_beta(EM_weights_bbar_given_d);

			std::cout << "curQ = " << curQ << ", newQ = " << Q << std::endl;

			//check if we accept or backtrack
			numUpdatedModels_old = numUpdatedModels;
			numUpdatedModels = _modelsToEstimate->acceptProposedParametersBasedOnQ();
			_logfile->write(toString(numUpdatedModels) + "/" + toString(_modelsToEstimate->size()) + " models converged.");

			if(numUpdatedModels > numUpdatedModels_old){
				_logfile->conclude("Q_beta was increased from " + toString(curQ) + " to " + toString(Q));
			}

			//backtrack
			lambda = lambda / 2.0; //backtrack;
			--log2_lambda;
		}

		if(numUpdatedModels < _modelsToEstimate->size()){
			_logfile->conclude("Some models did not improve even with log2(lambda) = " + toString(log2_lambda) + ", aborting Newton-Raphson.");
		}

		// c) adjust parameters post estimation
		_modelsToEstimate->adjustParametersPostEstimation();

		// d) get largest gradient (F) to check if we break NR optimization
		double maxF = _modelsToEstimate->getSteepestGradient();
		_logfile->conclude("max(F) = " + toString(maxF));
		_logfile->endIndent();
		if(maxF < _NewtonRaphsonMaxF || numUpdatedModels == 0) break;

		_logfile->endIndent();
	}
	_logfile->endIndent();
	_logfile->endIndent();
};

double TRecalibrationEMEstimator::_calculateLL_fullModel(const TPostMortemDamage & PmdModels){
	double LL = 0.0;
	const TGenotypeProbabilities& genoFreq = _genoDist->genotypeFrequencies();
	static TGenotypeLikelihoods genotypeLikelihoods;
	for(auto& s : _sites){
		//calculate genotype likelihoods
		_genotypeLikelihoodCalculator.fillGenotypeLikelihoods(s, genotypeLikelihoods, PmdModels, *_modelsToEstimate);

		//weight by genotype prior
		if(s.genotype() == genometools::NN){
			LL += log(genotypeLikelihoods.weightedSum(genoFreq));
		} else {
			LL += log(genotypeLikelihoods[s.genotype()].get());
		}
	}

	return LL;
};

void TRecalibrationEMEstimator::_runEM(std::string outputName, const TPostMortemDamage & PmdModels){
	//run EM
	_logfile->startNumbering("Running EM algorithm:");

	std::ofstream out;
	std::string filename;

	//calculate initial LL
	double oldLL = _calculateLL_fullModel(PmdModels);
	_logfile->conclude("Initial log Likelihood = " + toString(oldLL));

	//running iterations
	for(int iter = 0; iter < _numEMIterations; ++iter){
		_logfile->number("EM Iteration:"); _logfile->addIndent();

		//update theta_epsilon (sequencing errors)
		_updateEM_theta_epsilon(PmdModels);

		//calculate LL
		double LL = _calculateLL_fullModel(PmdModels);
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
		if(_writeTmpTables){
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
void TRecalibrationEMEstimator::writeCurrentEstimates(const std::string){
	//loop over read groups
	for(uint16_t r = 0; r < _readGroups->size(); ++r){
		//loop over mates
		for(uint8_t mate = 0; mate < 2; ++mate){
			//write model

		}
	}
};

double TRecalibrationEMEstimator::calcLL(){
	throw std::runtime_error("double TRecalibrationEMEstimator::calcLL() not yet implemented!");
};

void TRecalibrationEMEstimator::calcLikelihoodSurface(std::string, int){
	throw std::runtime_error("void TRecalibrationEMEstimator::calcLikelihoodSurface(std::string filename, int numMarginalGridPoints) not yet implemented!");
};

void TRecalibrationEMEstimator::calcQSurface(std::string, int){
	throw std::runtime_error("TRecalibrationEMEstimator::calcQSurface(std::string filename, int numMarginalGridPoints) not yet implemented!");
};

}; //end namespace RecalEstimator

}; //end namespace GenotypeLikelihoods
