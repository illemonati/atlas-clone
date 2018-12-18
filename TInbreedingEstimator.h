/*
 * TInbreedingEstimator.h
 *
 *  Created on: Dec 6, 2018
 *      Author: linkv
 */

#ifndef TINBREEDINGESTIMATOR_H_
#define TINBREEDINGESTIMATOR_H_

#include "TParameters.h"
#include "TLog.h"
#include "TQualityMap.h"
#include "TPopulationLikelihoods.h"
#include <limits>

//---------------------------
// allele frequencies p
//---------------------------
class TAlleleFreq{
private:
	float* proposalWidths;
	double* sumIterations;
	double* sumOfSquaresIterations;
	std::vector<double> alleleFreq;

public:
	long numLoci;

	TAlleleFreq();
	TAlleleFreq(std::vector<double> & P, float initialProposalWidth);
	double& operator[](int index){
		//Note: no check on range!
		return alleleFreq[index];
	};
	void adjustProposalWidthAfterBurnin(int* numAcceptedP, int numUpdates);
	double proposeNew(long & locusNum, TRandomGenerator & randomGenerator);
	void update(long & index, double & value);
};

//---------------------------
// alphaOrBeta
//---------------------------

class TAlphaOrBeta{
private:
	double _alphaOrBeta;
	double _logAlphaOrBeta;
public:
	std::string variableName;
	double proposalWidth;

	TAlphaOrBeta();
	TAlphaOrBeta(std::string VariableName, double & ProposalWidth);
	void update(double & newLogValue, double & newNaturalScaleValue);
	void adjustProposalWidthAfterBurnin(int numAccepted, int numUpdates);

	double getLogValue();
	double getNaturalScaleValue();
};

//---------------------------
// TInbreedingEstimator
//---------------------------

class TInbreedingEstimator{
private:
	TRandomGenerator randomGenerator;
	TQualityMap qualMap;

	//log
	TLog* logfile;
	std::string outname;

	//algorithm params
	int numIterations;
	double pi;
	double widthProposalKernelP;
	int numAcceptedF;
	int* numAcceptedP;
	int numAcceptedAlpha;
	int numAcceptedBeta;
	int numBurnins;
	int burninLength;
	int thinning;

	//data
	TPopulationLikelihoods likelihoods;
	unsigned long numLoci;

	//params
	double F;
	//std::vector<double> p;
	TAlleleFreq p;
	TAlphaOrBeta alpha;
	TAlphaOrBeta beta;

//	void initializeAlphaBeta();
	void initializeAlphaAndBeta();
	void initParams(TRandomGenerator & randomGenerator, TParameters & parameters);
	void printTrajectory(gz::ogzstream & tracefile);
	bool updateF();
	bool updateP(uint8_t* data, long & locusNum, int curSampleSize, TAlphaOrBeta & alpha, TAlphaOrBeta & beta);
	bool updateAlphaOrBeta(TAlphaOrBeta & alphaOrBetaToUpdate, TAlphaOrBeta & alphaOrBetaOther);
	double probGenoGivenFAndP(int & genotype, double & F, double & p);
	double logProbPGivenAlphaBeta();
	double logLikelihoodAllInds(uint8_t* data, int curSampleSize, double & thisP, double & thisF, TAlphaOrBeta & alpha, TAlphaOrBeta & beta);
	void wholeLogLikelihood();
	void oneMCMCIteration(int iterationNum);
	void printAcceptanceRates(int numIterations);
	void resetAcceptanceRates();
	void adjustProposalWidths();
	void writeParameterEstimatesOfIteration(gz::ogzstream & out);

public:
	TInbreedingEstimator(TParameters & Parameters, TLog* Logfile);
	~TInbreedingEstimator(){
		delete numAcceptedP;
	}
	void runEstimation();
	void writeLikelihoodForDebuggingAlpha(TParameters & params);

};



#endif /* TINBREEDINGESTIMATOR_H_ */
