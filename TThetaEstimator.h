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

using namespace GenotypeLikelihoods;

//---------------------------------------------------------------
//Theta
//---------------------------------------------------------------
struct Theta{
	double theta, expTheta, logTheta, thetaConfidence, LL;
	double* baseFreq;

	Theta(){
		theta = 0.0;
		thetaConfidence = 0.0;
		expTheta = 0.0;
		logTheta = -9e100;
		LL = -9e100;
		baseFreq = new double[4];
	};

	~Theta(){
		delete[] baseFreq;
	};

	void setTheta(double val){
		theta = val;
		expTheta = exp(-theta);
		logTheta = log(theta);
		LL = -9e100;
	};

	void setExpTheta(double val){
		expTheta = val;
		theta = -log(val);
		logTheta = log(theta);
		LL = -9e100;
	};

	void setLogTheta(double val){
		logTheta = val;
		theta = exp(val);
		expTheta = exp(-theta);
		LL = -9e100;
	};

	void setLogTheta(double & val, double & newLL){
		logTheta = val;
		theta = exp(val);
		expTheta = exp(-theta);
		LL = newLL;
	};

	std::string getBaseFrequencyString(){
		return "Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]);
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
	TGenotypeMap genoMap;

	//initial theta
	double initialTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;

	//estimation
	int minSitesWithData;
	GenotypeLikelihoods::TGenotypeData pGenotype; //P(g|pi, theta)
	Theta theta;
	bool extraVerbose;

	void initDataStorage();
	void readParametersRegardingInitialSearch(TParameters & params);
	void fillPGenotype(TGenotypeData & pGeno, const double & expTheta, const double* baseFrequencies);
	void fillPGenotype(TGenotypeData & pGeno, const Theta & thisTheta);

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

	void fillPGenotype(GenotypeLikelihoods::TGenotypeData & pGeno){ fillPGenotype(pGeno, theta); };

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
	TGenotypeData P_G; // see paper

	double _calcFisherInfo(const TGenotypeData & _pGenotype, const TGenotypeData deriv_pGenotype);
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
	void add(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	long sizeWithData(){ return data->sizeWithData();};
	bool estimateTheta();
	void setTheta(double Theta);
	void setBaseFreq(TBaseFrequencies & BaseFreq);
	void addToHeader(std::vector<std::string> & header, std::string prefix="");
	void writeestimateFrequenciesAndTheta(TOutputFile* out);
	void writeResultsToFile(TOutputFile* out);
	void calcLikelihoodSurface(gz::ogzstream & out, int & steps);
	void bootstrapTheta(TOutputFile* out);
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
	double* tmpBaseFreq;

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

	void initAdditionalTmpStorage();
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
		delete[] tmpBaseFreq;
	};

	TThetaEstimatorData* pointerToDataContainer2(){ return data2; };

	void estimateRatio(std::string ouputName);
};

//---------------------------------------------------------------
//TThetaOutputFile
//---------------------------------------------------------------
class TThetaOutputFile{
protected:
	TOutputFileZipped out;
	std::vector<TThetaEstimator*> thetaEstimators;
	std::vector<std::string> prefixes;

	void writeHeader(){
		std::vector<std::string> header = {"Chr", "start", "end"};

		//add headers of all estimators
		for(size_t i = 0; i < thetaEstimators.size(); ++i){
			thetaEstimators[i]->addToHeader(header, prefixes[i]);
		}
		out.writeHeader(header);
	};

	void writeEstimates(){
		for(TThetaEstimator* est: thetaEstimators){
			est->writeResultsToFile(&out);
		}
		out << std::endl;
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


	void addEstimator(TThetaEstimator* Estimator, const std::string Prefix){
		if(out.isOpen()){
			throw "Can not add estimators to an open TThetaOutputFile!";
		}
		thetaEstimators.push_back(Estimator);
		prefixes.push_back(Prefix);
	};

	void open(const std::string Filename, TLog* logfile){
		logfile->list("Will write theta estimates to file '" + Filename + "'.");
		out.open(Filename);
		writeHeader();
	};

	void open(TThetaEstimator* Estimator, const std::string Filename, TLog* logfile){
		addEstimator(Estimator, "");
		open(Filename, logfile);
	};

	void close(){
		out.close();
	};

	void write(const std::string & chr, const long & start, const long & end){
		out << chr << start << end;
		writeEstimates();
	};

	void write(const std::string & chr, const std::string & start, const std::string & end){
		out << chr << start << end;
		writeEstimates();
	};
};


#endif /* TTHETAESTIMATOR_H_ */
