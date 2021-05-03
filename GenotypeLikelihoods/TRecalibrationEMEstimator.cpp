/*
 * TRecalibrationEM.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMEstimator.h"

namespace GenotypeLikelihoods {

namespace RecalEstimatorTools {

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
TDataTable::TDataTable(){
	maxQual = 0;
	maxFragmentLength = 0;
	maxMQ = 0;
	maxPos = 0;
	counts = 0;
	initialized = false;
	qualities = nullptr;
	fragmentLengths = nullptr;
	MQ = nullptr;
};

TDataTable::TDataTable(const int MaxQual, const int MaxFragmentLength, const int maxMQ){
	initialize(MaxQual, MaxFragmentLength, maxMQ);
};

TDataTable::~TDataTable(){
	delete[] qualities;
	delete[] fragmentLengths;
	delete[] MQ;
};

void TDataTable::initialize(const int MaxQual, const int MaxFragmentLength, const int MaxMQ){
	if(initialized){
		delete[] qualities;
		delete[] fragmentLengths;
		delete[] MQ;
	}
	maxQual = MaxQual;
	qualities = new unsigned int[maxQual];

	maxFragmentLength = MaxFragmentLength;
	fragmentLengths = new unsigned int[maxFragmentLength];

	maxMQ = MaxMQ;
	MQ = new unsigned int[maxMQ];

	initialized = true;
	clear();
};

void TDataTable::add(const BAM::TBase & base){
	++counts;
	++qualities[base.originalQuality_phredInt];
	++fragmentLengths[base.fragmentLength];
	++MQ[base.mappingQuality];
	if(maxPos < base.distFrom5Prime)
		maxPos = base.distFrom5Prime;
};

void TDataTable::clear(){
	maxPos = 0;
	for(int q=0; q<maxQual; ++q)
		qualities[q] = 0;
	counts = 0;
};

uint64_t TDataTable::size() const{
	return counts;
};

void TDataTable::fillVectorWithUsedQualities(std::vector<uint16_t> & Q) const{
	Q.clear();
	for(int i=0; i<maxQual; ++i){
		if(qualities[i] > 0){
			Q.push_back(i);
		}
	};
};

void TDataTable::fillVectorWithUsedFragmentLengths(std::vector<uint16_t> & lengths) const{
	lengths.clear();
	for(int i=0; i<maxFragmentLength; ++i){
		if(fragmentLengths[i] > 0){
			lengths.push_back(i);
		}
	};
};

void TDataTable::fillVectorWithUsedMQ(std::vector<uint16_t> & MQ) const{
	MQ.clear();
	for(int i=0; i<maxMQ; ++i){
		if(MQ[i] > 0){
			MQ.push_back(i);
		}
	};
};

//--------------------------------------------------------------------
TDataTables::TDataTables(){
	_initialized = false;
	_numReadGroups = 0;
	_maxQual = 0;
	_tables = nullptr;
	_totalCounts = 0;
};

TDataTables::TDataTables(const  int NumReadGroups, const int MaxQual, const int MaxFragmentLength, const int MaxMQ){
	_initialized = false;
	init(NumReadGroups, MaxQual, MaxFragmentLength, MaxMQ);
};

TDataTables::~TDataTables(){
	clear();
};

void TDataTables::init(const int NumReadGroups, const int MaxQual, const int MaxFragmentLength, const int MaxMQ){
	//empty if already initialized
	if(_initialized){
		clear();
	}

	//crate storage
	_numReadGroups = NumReadGroups;
	_maxQual = MaxQual;

	_tables = new TDataTable*[_numReadGroups];
	for(int rg = 0; rg<_numReadGroups; ++rg){
		_tables[rg] = new TDataTable[2];
		_tables[rg][0].initialize(_maxQual, MaxFragmentLength, MaxMQ);
		_tables[rg][1].initialize(_maxQual, MaxFragmentLength, MaxMQ) ;
	}

	_initialized = true;

	//set to zero
	reset();
};

void TDataTables::clear(){
	for(int rg = 0; rg<_numReadGroups; ++rg){
		delete[] _tables[rg];
	}

	delete[] _tables;
};

void TDataTables::reset(){
	for(int rg = 0; rg<_numReadGroups; ++rg){
		_tables[rg][0].clear();
		_tables[rg][1].clear();
	}
	_totalCounts = 0;
};

void TDataTables::add(const BAM::TBase & base){
	++_totalCounts;
	_tables[base.readGroupID][(int) base.isSecondMate()].add(base);
};

void TDataTables::add(const TSite & site){
	_totalCounts += site.depth();
	for(std::vector<BAM::TBase>::const_iterator b = site.cbegin(); b != site.cend(); ++b){
		_tables[b->readGroupID][(int) b->isSecondMate()].add(*b);
	}
};

void TDataTables::fillVectorWithUsedQualities(const int readGroupId, const bool isSecondMate, std::vector<uint16_t> & Q){
	_tables[readGroupId][(int) isSecondMate].fillVectorWithUsedQualities(Q);
};

uint64_t TDataTables::size() const{
	return _totalCounts;
};

const TDataTable& TDataTables::table(const int readGroupId, const bool isSecondMate) const{
	return &_tables[readGroupId][(int) isSecondMate];
};

const TDataTable& TDataTables::table(const int readGroupId, const int Mate) const{
	return &_tables[readGroupId][Mate];
};




//--------------------------------------------------------------------------------------
// TModelIndex
//--------------------------------------------------------------------------------------
void TModelIndex::set(const uint16_t & ReadGroupId, const bool & IsSecondMate, const uint16_t & ModelIndex, const BAM::TReadGroupMap & ReadGroupMap){
	//save index for read group AND all read groups pooled with it!

	for(auto& r : ReadGroupMap.readGroupsPooledWith(ReadGroupId)){
		_index[r][IsSecondMate].set(ModelIndex);
	}
};

}; //end namespaceRecal EstimatorTools

//--------------------------------------------------------------------------
// TSequencingErrorModelVectorForEstimation
//--------------------------------------------------------------------------

TSequencingErrorModelVectorForEstimation::TSequencingErrorModelVectorForEstimation(TSequencingErrorModels & SequencingErrorModels,
																				   const RecalEstimatorTools::TDataTables & DataTables,
																				   const BAM::TReadGroups & ReadGroups,
																				   const BAM::TReadGroupMap & ReadGroupMap,
																				   const uint32_t & MinRequiredObservations,
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

			const RecalEstimatorTools::TDataTable& table = DataTables.table(r, mate);

			//check if there is sufficient data for _first mate
			if(table.size() > 0){
				//check if model is estimatable
				if(SequencingErrorModels[r][mate].estimatable()){

					//copy model and update index
					std::shared_ptr<TSequencingErrorModelRecal>& model = SequencingErrorModels[r].getPointerToRecalModel(mate);
					_models.push_back(model);
					modelStati[r][RecalEstimatorTools::copied].set(mate);
					_modelIndex.set(r, mate, model, ReadGroupMap);

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


//functions to estimate rho
//-------------------------------------------------------------------
//functions to estimate rho
void TSequencingErrorModelVectorForEstimation::prepareRhoEstimationFromEMWeights(){
	for(auto& model : _models){
		model->prepareRhoEstimationFromEMWeights();
	}
};

void TSequencingErrorModelVectorForEstimation::addBaseForRhoEstimation(BAM::TBase & base, const TBaseData & EMWeights){
	_modelIndex(base)->addBaseForRhoEstimation(base, EMWeights);
};

void TSequencingErrorModelVectorForEstimation::estimateRho(){
	for(auto& model : _models){
		model->estimateRho();
	}
};

//functions to estimate beta
//-------------------------------------------------------------------

void TSequencingErrorModelVectorForEstimation::addToFandJacobian(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d){
	_modelIndex(base)->addToFandJacobian(base, EM_weights_bbar_given_d);
};

void TSequencingErrorModelVectorForEstimation::addToQ(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d){
	_modelIndex(base)->addToQ(base, EM_weights_bbar_given_d);
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
	std::string modelTagForEstimation = args.getParameterStringWithDefault("model", "quality:polynomial(2);position:polynomial(2);context:specific");


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
	compileDataTable(_dataTables);
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
			TDataTable* table = _dataTables.table(rg, mate);
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
			TDataTable* table = _dataTables.table(index, false);
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

void TRecalibrationEMEstimator::compileDataTable(TDataTables & dataTable){
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
		//get relevant base frequencies P(b): from known genotype or distribution if genotype is unknown
		TBaseData baseFreq;
		_fillRelevantBaseFrequencies(baseFreq, s.genotype());

		//calculate weights per base
		for(auto& b : s){
			//calculate P(bbar) = \sum_b P(bbar|b)P(b)
			TBaseData PMD;
			_pmd.calculateBaseLikelihoods(b, baseFreq, PMD);

			//calculate P(d|bbar)
			_sequencingErrorModels.calculateBaseLikelihoods(b, EMWeights[index]);

			//calculate P(d|bbar) \propto P(d|bbar)P(bbar)
			EMWeights[index] *= PMD;
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
		int numUpdatedModels = 0;
		int numUpdatedModels_old;

		while(numUpdatedModels < _sequencingErrorModels.numModels() && lambda > 1.0E-20){
			//propose move
			_logfile->listFlush("Proposing move with log2(lambda) = " + toString(log2_lambda) + " ... ");
			_sequencingErrorModels.proposeNewParameters(lambda);

			std::cout << std::endl;
			_sequencingErrorModels.print();
			std::cout << std::endl;

			//calculate Q at new location
			double Q = _calculate_Q_beta(EM_weights_bbar_given_d);

			std::cout << "curQ = " << curQ << ", newQ = " << Q << std::endl;

			//check if we accept or backtrack
			numUpdatedModels_old = numUpdatedModels;
			numUpdatedModels = _sequencingErrorModels.acceptProposedParametersBasedOnQ();
			_logfile->write(toString(numUpdatedModels) + "/" + toString(_sequencingErrorModels.numModels()) + " models converged.");

			if(numUpdatedModels > numUpdatedModels_old){
				_logfile->conclude("Q_beta was increased from " + toString(curQ) + " to " + toString(Q));
			}

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
	_logfile->endIndent();
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

		//TODO: update other parameters -> refactor to be modular

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
