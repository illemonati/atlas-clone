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
// F
//---------------------------
class TInbreedingF{
private:
	double _F;
	float _probMovingToModelNoF;
	double _sdProposal;
	bool _inModelWithF;
	double _lambda, _logLambda, _expMinusLambda;
	int _posteriorProbModelWithF;
	double _sumIterations, _sumOfSquaresIterations;

public:
	TInbreedingF();
	TInbreedingF(float & ProbMovingToModelNoF, double & SdProposal, bool InModelWithF, double lambda);
	void adjustProposalWidthAfterBurnin(int numAcceptedFModelF, int numIterInModelF);
	double proposeNew(TRandomGenerator* randomGenerator);
	void updatePosteriors(const double & value, const bool & inModelWithF);
	void updateAndAccept(const double & value, const bool & inModelWithF);
	void updateAndReject(bool inModelWithF);
	void resetPosterior();
	double getPosteriorMean(int numUpdates);
	double getPosteriorVariance(int numUpdates);
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
	long _numLociModelP;
	long _numLociWithAcceptanceZero;

public:
	std::vector<bool> modelP;
	std::vector<int> posteriorProbModelP;
	std::vector<double> initialAlleleFreq;
	std::vector<float> proposalWidths;

	long numLoci;
	double probMovingToModel0, probMovingToModelP;
	double logProbMovingToModel0, logProbMovingToModelP;
	double lambda, logLambda, expMinusLambda;
	double minAlleleFreq;

	TAlleleFreq();
	TAlleleFreq(std::vector<double> & P, double & initialProposalWidthFactor, const int numSamples, double & ProbMovingToModel0, double & ProbMovingToModelP, double & Lambda, bool trueAlleleFreqProvided);

	double operator[](long index){
		//Note: no check on range!
		return alleleFreq[index];
	};
//	void initializeModels();
	void setSumsForPosteriorToZero();
	void setToValue(double fixedValue);
	void adjustProposalWidthAfterBurnin(std::vector<int> & numAcceptedP, std::vector<int> & numUpdates);
	double proposeNew(const long & locusNum, TRandomGenerator* randomGenerator);
	void update(const long & index, const double & value, const bool ModelP);
	double getPosteriorMean(unsigned long index, int numUpdates);
	double getPosteriorVariance(unsigned long index, int numUpdates);
	double getProposalWidth(const unsigned long & index);
	long getNumLociInModelP();
	long getNumLociInModel0();
 	long numLociWithAcceptanceZero();
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
	double _logPi;
	double _logOneMinusPi;
public:
	double proposalWidth;

	TPi();
	TPi(double & ProposalWidth, double & initialValue);
	void update(const double & newNaturalScaleValue);
	void adjustProposalWidthAfterBurnin(int numAccepted, int numUpdates);

	double getPi();
	double getLogPi(){ return _logPi; };
	double getLogOneMinusPi(){ return _logOneMinusPi; };
	double getProposalWidth();
	double proposeNew(TRandomGenerator* randomGenerator);
};
//---------------------------
// TInbreedingEstimator
//---------------------------

class TInbreedingEstimator{
private:
	TRandomGenerator* randomGenerator;
//	TQualityMap qualMap;

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

	long numLociToTrace;

	bool shouldUpdateF;
	bool shouldUpdateP;
	bool shouldUpdateGamma;
	bool shouldUpdatePi;

	bool adjustFAfterBurnin;
	bool adjustPAfterBurnin;
	bool adjustGammaAfterBurnin;
	bool adjustPiAfterBurnin;

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
	void initializeGamma(TParameters & parameters);
	void initF(TParameters & parameters);
	void initAlleleFreq(TParameters & parameters);
	void initParams(TRandomGenerator* randomGenerator, TParameters & parameters);
	void resetToInitialValuesDebugging();
	void checkHastingsRatios();
	bool updateF();
	bool updateP(const TSampleLikelihoods* data, const long locusNum, const int curSampleSize, const TGamma Gamma);
	bool updateGamma();
	bool updatePi();
	double logProbPGivenGamma();
	double logLikelihoodAllInds(const TSampleLikelihoods* data, const int curSampleSize, const double thisP, const double thisF, TGlfConverter & glfConverter);
	void wholeLogLikelihood();
	double getLogLikelihoodCurrentParams(const TGlfConverter & glfConverter);
	double compareLikelihoods(int numUpdates, const TGlfConverter & glfConverter);
	void oneMCMCIteration();
	void printAcceptanceRates(int numIterations);
	void resetAcceptanceRates();
	void adjustProposalWidths();
	void writeParameterEstimatesOfIteration(std::ofstream & out);
	void writePosteriors(int i);
	void runBurnins(std::ofstream & out, TParameters & params);
	void runMCMC(std::ofstream & out, TParameters & params);

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
