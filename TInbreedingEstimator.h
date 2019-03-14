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
	void adjustProposalWidthAfterBurnin(int numAcceptedFModelF, int numIterInModelF);
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
	std::vector<double> sumIterations;
	std::vector<double> sumOfSquaresIterations;
	std::vector<double> alleleFreq;
	std::vector<float> proposalWidths;

public:
	std::vector<bool> modelP;
	std::vector<int> posteriorProbModelP;

	long numLoci;
	long numLociModelP;
	float probMovingToModel0;
	double lambda;
	double minAlleleFreq;

	TAlleleFreq();
	TAlleleFreq(std::vector<double> & P, double & initialProposalWidthFactor, const int numSamples, float & ProbMovingToModel0, double & Lambda);

	double operator[](long index){
		//Note: no check on range!
		return alleleFreq[index];
	};
//	void initializeModels();
	void setSumsForPosteriorToZero();
	void setToValue(double fixedValue);
	void adjustProposalWidthAfterBurnin(std::vector<int> & numAcceptedP, std::vector<int> & numUpdates);
	double proposeNew(long & locusNum, TRandomGenerator* randomGenerator);
	void update(long & index, const double & value, const bool modelP);
	double getPosteriorMean(unsigned long & index, int numUpdates);
	double getPosteriorVariance(unsigned long & index, int numUpdates);
	double getProposalWidth(const unsigned long & index);
	long getNumLociInModelP();
	long getNumLociInModel0();
	double logPDFExp(const double & thisP);
	double logPDFExp(const long & thisLocus);
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
	double _minPi;
public:
	double proposalWidth;

	TPi();
	TPi(double & ProposalWidth, double & initialValue);
	void update(const double & newNaturalScaleValue);
	void adjustProposalWidthAfterBurnin(int numAccepted, int numUpdates);

	double getPi();
	double getProposalWidth();
	double proposeNew(TRandomGenerator* randomGenerator);
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
	int numAcceptedFModelF;
	int numIterInModelF;
	std::vector<int> numAcceptedP;
	std::vector<int> numAcceptedPModelP;
	std::vector<int> numIterInModelP;
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
		delete randomGenerator;
	}
	void runEstimation(TParameters & params);
	void writeLikelihoodForDebuggingGamma(TParameters & params);
	void writeLikelihoodForDebuggingAlleleFreq(TParameters & params);
	void writeLikelihoodForDebuggingF(TParameters & params);
};



#endif /* TINBREEDINGESTIMATOR_H_ */
