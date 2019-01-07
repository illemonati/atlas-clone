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
class TInbreedingF{
private:
	double _F;
	float _probMovingToModelNoF;
	double _sdProposal;
	bool _inModelWithF;
	double _lambda;
	int _posteriorProbModelWithF;

public:
	TInbreedingF();
	TInbreedingF(double F, float & ProbMovingToModelNoF, double & SdProposal, bool InModelWithF, double lambda);
	void adjustProposalWidthAfterBurnin(int numAcceptedF, int numUpdates, TLog* logfile);
	double proposeNew(TRandomGenerator & randomGenerator);
	void update(double value, bool inModelWithF);
	float probMovingToModelNoF();
	double F();
	bool inModelWithF();
	double logPDFExp();
	double PDFExp();
	double lambda();
	int posteriorProbModelWithF();
};

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
	double operator[](long index){
		//Note: no check on range!
		return alleleFreq[index];
	};
	void setSumsToZero();
	void setToValue(double fixedValue);
	void adjustProposalWidthAfterBurnin(int* numAcceptedP, int numUpdates);
	double proposeNew(long & locusNum, TRandomGenerator & randomGenerator);
	void update(long & index, double value);
	double getPosteriorMean(unsigned long & index, int numUpdates);
	double getPosteriorVariance(unsigned long & index, int numUpdates);
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
	double sdF;
	int numAcceptedF;
	int* numAcceptedP;
	int numAcceptedAlpha;
	int numAcceptedBeta;
	int numBurnins;
	int burninLength;
	int thinning;

	//data
	TPopulationLikelihoods likelihoods;
	unsigned int numLoci;

	//params
	TInbreedingF F;
	//std::vector<double> p;
	TAlleleFreq p;
	TAlphaOrBeta alpha;
	TAlphaOrBeta beta;

//	void initializeAlphaBeta();
	void initializeAlphaAndBeta();
	void initParams(TRandomGenerator & randomGenerator, TParameters & parameters);
	bool updateF();
	bool updateP(uint8_t* data, long & locusNum, int curSampleSize, TAlphaOrBeta & alpha, TAlphaOrBeta & beta);
	bool updateAlphaOrBeta(TAlphaOrBeta & alphaOrBetaToUpdate, TAlphaOrBeta & alphaOrBetaOther);
	double logProbPGivenAlphaBeta();
	double logLikelihoodAllInds(uint8_t* data, int curSampleSize, double thisP, double thisF, TAlphaOrBeta & alpha, TAlphaOrBeta & beta);
	void wholeLogLikelihood();
	void oneMCMCIteration(int iterationNum);
	void printAcceptanceRates(int numIterations);
	void resetAcceptanceRates();
	void adjustProposalWidths();
	void writeParameterEstimatesOfIteration(std::ofstream & out);

public:
	TInbreedingEstimator(TParameters & Parameters, TLog* Logfile);
	~TInbreedingEstimator(){
		delete numAcceptedP;
	}
	void runEstimation();
	void writeLikelihoodForDebuggingAlpha(TParameters & params);
	void writeLikelihoodForDebuggingBeta(TParameters & params);
	void writeLikelihoodForDebuggingAlleleFreq(TParameters & params);
	void writeLikelihoodForDebuggingF(TParameters & params);
};



#endif /* TINBREEDINGESTIMATOR_H_ */
