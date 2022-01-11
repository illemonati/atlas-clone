/*
 * TThetaEstimator.h
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATOR_H_
#define TTHETAESTIMATOR_H_

#include "TThetaEstimatorData.h"
#include "TRandomGenerator.h"
#include <stdio.h>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include "TFile.h"
#include "TGenotypeData.h"
#include "TLog.h"
#include "TParameters.h"
#include "TWindow.h"

namespace GenotypeLikelihoods{

using coretools::TRandomGenerator;

//---------------------------------------------------------------
//Theta
//---------------------------------------------------------------
struct Theta{
	double theta, expTheta, logTheta, thetaConfidence, LL;
	TBaseProbabilities baseFreq;

	Theta(){
		theta = 0.0;
		thetaConfidence = 0.0;
		expTheta = 0.0;
		logTheta = -9e100;
		LL = -9e100;
	};

	void setTheta(const double val){
		theta = val;
		expTheta = exp(-theta);
		logTheta = log(theta);
		LL = -9e100;
	};

	void setExpTheta(const double val){
		expTheta = val;
		theta = -log(val);
		logTheta = log(theta);
		LL = -9e100;
	};

	void setLogTheta(const double val){
		logTheta = val;
		theta = exp(val);
		expTheta = exp(-theta);
		LL = -9e100;
	};

	void setLogTheta(const double & val, const double & newLL){
		logTheta = val;
		theta = exp(val);
		expTheta = exp(-theta);
		LL = newLL;
	};

	std::string getBaseFrequencyString(){
		return coretools::str::toString("Pi(A) = ", baseFreq[genometools::A], ", Pi(C) = ", baseFreq[genometools::C], ", Pi(G) = ", baseFreq[genometools::G], ", Pi(T) = ", baseFreq[genometools::T]);
	};
};

//---------------------------------------------------------------
//TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator_base{
protected:
	TLog* logfile;
	TRandomGenerator* randomGenerator;

	//data
	TThetaEstimatorData* data;
	bool dataInitialized;
	bool useTmpFile;
	std::string tmpFileName;

	//initial theta
	double initialTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;

	//estimation
	int minSitesWithData;
	GenotypeLikelihoods::TGenotypeProbabilities pGenotype; //P(g|pi, theta)
	Theta theta;
	bool extraVerbose;

	void initDataStorage();
	void readParametersRegardingInitialSearch(TParameters & params);
	void fillPGenotype(GenotypeLikelihoods::TGenotypeProbabilities & pGeno, const double & expTheta, const TBaseProbabilities & baseFrequencies);
	void fillPGenotype(GenotypeLikelihoods::TGenotypeProbabilities & pGeno, const Theta & thisTheta);

	void findGoodStartingTheta(TThetaEstimatorData* thisData, Theta & thisTheta, std::string tag);

public:
	TThetaEstimator_base(TLog* Logfile, TRandomGenerator* randomGenerator);
	TThetaEstimator_base(TParameters & params, TLog* Logfile, TRandomGenerator* randomGenerator);
	TThetaEstimator_base(const TThetaEstimator_base & other);

	virtual ~TThetaEstimator_base(){
		if(dataInitialized)
			delete data;
	};

	TThetaEstimatorData* pointerToDataContainer(){ return data; };

	void fillPGenotype(GenotypeLikelihoods::TGenotypeProbabilities & pGeno){ fillPGenotype(pGeno, theta); };

};

//---------------------------------------------------------------
//TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator:public TThetaEstimator_base{
private:
	//EM parameters
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	bool estimationSuccessful;

	//tmp vectors
	GenotypeLikelihoods::TGenotypeData P_G; // see paper

	double _calcFisherInfo(const TGenotypeProbabilities & _pGenotype, const TGenotypeData deriv_pGenotype);
	bool _NRAllParams();
	void _NROnlyTheta();
	void _runEMForTheta();
	void _estimateConfidenceInterval();

public:
	TThetaEstimator(TParameters & params, TLog* Logfile, TRandomGenerator* randomGenerator);
	TThetaEstimator(TLog* Logfile, TRandomGenerator* randomGenerator);
	TThetaEstimator(const TThetaEstimator & other);

	virtual ~TThetaEstimator(){};

	void clear();
	void add(const TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	void add(const TWindow & window, const TGenotypeLikelihoodCalculator & glCalculator);
	long sizeWithData(){ return data->sizeWithData();};
	bool estimateTheta();
	void setTheta(const double Theta);
	void setBaseFreq(const GenotypeLikelihoods::TBaseProbabilities & BaseFreq);
	void addToHeader(std::vector<std::string> & header, std::string prefix="");
	void writeEstimateFrequenciesAndTheta(coretools::TOutputFile & out);
	void writeResultsToFile(coretools::TOutputFile & out);
	void calcLikelihoodSurface(coretools::TOutputFile & out, uint32_t steps);
	void bootstrapTheta();
};

//---------------------------------------------------------------
//TThetaEstimatorRatio
//---------------------------------------------------------------
class TThetaEstimatorRatio:public TThetaEstimator_base{
private:
	//second data
	TThetaEstimatorData* data2;
	bool data2Initialized;
	Theta theta2;

	//MCMC parameters
	double phiPriorMean;
	double phiPriorVar;
	double phiPriorOneOverTwoVar;
	double sdProposalKernelTheta1;
	double sdProposalKernelTheta2;
	double sdProposalKernelBaseFreq1;
	double sdProposalKernelBaseFreq2;
	int numIterations;
	int burnin;
	int thinning;
	int numAcceptedTheta1;
	int numAcceptedTheta2;
	int numAcceptedBaseFreq1;
	int numAcceptedBaseFreq2;

	void clearCounters();
	void concludeAcceptanceRate(const int & numAccepted, const int & length, std::string name);
	void concludeAcceptanceRateUpdateProposal(const int & numAccepted, const int & length, double & sd, std::string name);
	void concludeAcceptanceRates(const int & length);
	void concludeAcceptanceRatesUpdateProposal(const int & length);
	bool updateTheta(TThetaEstimatorData* thisData, Theta & thisTheta, double otherLogThetaMean, const double & thisSdProposalKernel);
	bool updateBaseFrequencies(TThetaEstimatorData* thisData, Theta & thisTheta, const double & thisSdProposalKernel);
	void oneMCMCIteration();

public:
	TThetaEstimatorRatio(TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TThetaEstimatorRatio(){
		if(data2Initialized)
			delete data2;
	};

	TThetaEstimatorData* pointerToDataContainer2(){ return data2; };

	void estimateRatio(std::string ouputName);
};

//---------------------------------------------------------------
//TThetaOutputFile
//---------------------------------------------------------------
class TThetaOutputFile{
protected:
	coretools::TOutputFile _out;
	std::vector<TThetaEstimator*> _thetaEstimators;
	std::vector<std::string> _prefixes;

	void _writeHeader(){
		std::vector<std::string> header = {"Chr", "start", "end"};

		//add headers of all estimators
		for(size_t i = 0; i < _thetaEstimators.size(); ++i){
			_thetaEstimators[i]->addToHeader(header, _prefixes[i]);
		}
		_out.writeHeader(header);
	};

	void _writeEstimates(){
		for(TThetaEstimator* est: _thetaEstimators){
			est->writeResultsToFile(_out);
		}
		_out << std::endl;
	};

public:
	TThetaOutputFile(){};

	TThetaOutputFile(TThetaEstimator* Estimator, const std::string Prefix){
		addEstimator(Estimator, Prefix);
	};

	TThetaOutputFile(TThetaEstimator* Estimator, const std::string Filename, TLog* logfile){
		addEstimator(Estimator, "");
		open(Filename, logfile);
	};

	~TThetaOutputFile(){
		close();
	};

	coretools::TOutputFile& file(){ return _out; };

	void addEstimator(TThetaEstimator* Estimator, const std::string Prefix){
		if(_out.isOpen()){
			throw "Can not add estimators to an open TThetaOutputFile!";
		}
		_thetaEstimators.push_back(Estimator);
		_prefixes.push_back(Prefix);
	};

	void open(const std::string Filename, TLog* logfile){
		logfile->list("Will write theta estimates to file '" + Filename + "'.");
		_out.open(Filename);
		_writeHeader();
		_out.setPrecision(9);
	};

	void open(TThetaEstimator* Estimator, const std::string Filename, TLog* logfile){
		addEstimator(Estimator, "");
		open(Filename, logfile);
	};

	void close(){
		_out.close();
	};

	void write(const TWindow_base & window){
		_out << window;
		_writeEstimates();
	};

	void write(const std::string & chr, const std::string & start, const std::string & end){
		_out << chr << start << end;
		_writeEstimates();
	};
};

}; //end namespace

#endif /* TTHETAESTIMATOR_H_ */
