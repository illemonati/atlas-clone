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
	double logPDFExp(const double & thisF);
	double logPDFExp();
	double PDFExp(const double & thisF);
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
	long numLociModelP;
	bool* modelP;
	float probMovingToModel0;
	double lambda;

	TAlleleFreq();
	TAlleleFreq(std::vector<double> & P, double & initialProposalWidthFactor, const int numSamples, float & ProbMovingToModel0, double & Lambda);
	double operator[](long index){
		//Note: no check on range!
		return alleleFreq[index];
	};
	void initializeModels(const double & cutoff);
	void setSumsToZero();
	void setToValue(double fixedValue);
	void adjustProposalWidthAfterBurnin(int* numAcceptedP, int numUpdates);
	double proposeNew(long & locusNum, TRandomGenerator* randomGenerator);
	void update(long & index, double value);
	double getPosteriorMean(unsigned long & index, int numUpdates);
	double getPosteriorVariance(unsigned long & index, int numUpdates);
	double getProposalWidth(const unsigned long & index);
	long getNumLociInModelP();
	long getNumLociInModel0();
	void resetPosterior(const unsigned long & index);
};

//---------------------------
// gamma
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
// Pi
//---------------------------

class TPi{
private:
	double _pi;
	double _logPi;
public:
	double proposalWidth;

	TPi();
	TPi(double & ProposalWidth, double & initialValue);
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
	int numAcceptedPi;
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
	TPi pi;

//	void initializeAlphaBeta();
	void initializeGamma();
	void initParams(TRandomGenerator* randomGenerator, TParameters & parameters);
	bool updateF();
	bool updateP(uint8_t* data, long & locusNum, int curSampleSize, TGamma & Gamma);
	bool updateGamma();
	bool updatePi();
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
