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
	double proposeNew(TRandomGenerator* randomGenerator);
	void updateAndAccept(double value, bool inModelWithF);
	void updateAndReject(bool inModelWithF);
	void resetPosterior();
	float probMovingToModelNoF();
	double F();
	bool inModelWithF();
	double logPDFExp(double & thisF);
	double logPDFExp();
	double PDFExp(double & thisF);
	double PDFExp();
	double lambda();
	int posteriorProbModelWithF();
	double proposalWidth();
};

//---------------------------
// allele frequencies p
//---------------------------
class TAlleleFreq{
private:
	double* sumIterations;
	double* sumOfSquaresIterations;
	std::vector<double> alleleFreq;
	double minAlleleFreq;
	float* proposalWidths;

public:
	long numLoci;


	TAlleleFreq();
	TAlleleFreq(std::vector<double> & P, float initialProposalWidth, int numSamples);
	double operator[](long index){
		//Note: no check on range!
		return alleleFreq[index];
	};
	void setSumsToZero();
	void setToValue(double fixedValue);
	void adjustProposalWidthAfterBurnin(int* numAcceptedP, int numUpdates);
	double proposeNew(long & locusNum, TRandomGenerator* randomGenerator);
	void update(long & index, double value);
	double getPosteriorMean(unsigned long & index, int numUpdates);
	double getPosteriorVariance(unsigned long & index, int numUpdates);
	double getProposalWidth(const unsigned long & index);
	void resetPosterior(const unsigned long & index);
};

//---------------------------
// alphaOrBeta
//---------------------------

class TGamma{
private:
	double _gamma;
	double _logGamma;
public:
	double proposalWidth;

	TGamma();
	TGamma(double & ProposalWidth);
	void update(const double & newLogValue, const double & newNaturalScaleValue);
	void adjustProposalWidthAfterBurnin(int numAccepted, int numUpdates);

	double getLogValue();
	double getNaturalScaleValue();
	double getProposalWidth();
};

//---------------------------
// TInbreedingEstimator
//---------------------------

class TInbreedingEstimator{
private:
	TRandomGenerator* randomGenerator;
	TQualityMap qualMap;

	//log
	TLog* logfile;
	std::string outname;

	//algorithm params
	int numIterations;
	double sdF;
	int numAcceptedF;
	int* numAcceptedP;
	int numAcceptedGamma;
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
	bool trueAlleleFreqProvided;
	std::vector<double> trueAlleleFreq;
	TGamma Gamma;

//	void initializeAlphaBeta();
	void initializeGamma();
	void initParams(TRandomGenerator* randomGenerator, TParameters & parameters);
	bool updateF();
	bool updateP(uint8_t* data, long & locusNum, int curSampleSize, TGamma & Gamma);
	bool updateGamma();
	double logProbPGivenGamma();
	double logLikelihoodAllInds(uint8_t* data, int curSampleSize, double thisP, double thisF);
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
		delete randomGenerator;
	}
	void runEstimation(TParameters & params);
	void writeLikelihoodForDebuggingGamma(TParameters & params);
	void writeLikelihoodForDebuggingAlleleFreq(TParameters & params);
	void writeLikelihoodForDebuggingF(TParameters & params);
};



#endif /* TINBREEDINGESTIMATOR_H_ */
