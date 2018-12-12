#ifndef TRandomGenerator_H_
#define TRandomGenerator_H_

#include <cmath>
//#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <limits>

#include "mathFunctions.h"

//TODO: move math functions to math_functions.h

class TRandomGenerator{
private:
	long _Idum;
	int binomPValueTableSize;

	double ran3();
	double _lowerIncompleteGamma(double alpha, double z);
	double _upperIncompleteGamma(double alpha, double z);
	double _betacf(const double a, const double b, const double x);
	double _betaiapprox(double a, double b, double x);
	double _binomPValue(const int & k, const int & l);
	long get_randomSeedFromCurrentTime(long & addToSeed);
	void init();

public:
	long usedSeed;
	double* factorialTable;
	bool factorialTableInitialized;
	double* factorialTableLn;
	bool factorialTableLnInitialized;
	bool binomPValueTableInitialized;
	double** binomPValueTable;

	TRandomGenerator(long addToSeed){
		setSeed(addToSeed);
		init();
	};
	TRandomGenerator(long addToSeed, bool seedIsFixed){
		setSeed(addToSeed, seedIsFixed);
		init();
	};
	TRandomGenerator(){
		setSeed(0);
		init();
	};
	~TRandomGenerator(){
		if(factorialTableInitialized)
			delete[] factorialTable;
		if(factorialTableLnInitialized)
			delete[] factorialTableLn;
		if(binomPValueTableInitialized){
			for(int i=0; i<binomPValueTableSize; ++i)
				delete[] binomPValueTable[i];
			delete[] binomPValueTable;
		}

	};
	void setSeed(long addToSeed, bool seedIsFixed=false);

	//uniform
	double getRand(){ return ran3(); };
	double getRand(double min, double max);
	int getRand(int min, int maxPlusOne);
	int pickOne(int numElements);
	int pickOne(int numElements, float* probsCumulative);
	int pickOne(int numElements, double* probsCumulative);
	long getRand(long min, long maxPlusOne);

	void shuffle(std::vector<int> & vec);


	//normal
	double getNormalRandom (double dMean, double dStdDev);
	double getLogNormalRandom (double dMean, double dStdDev);
	double normalCumulativeDistributionFunction(double x, double mean, double sigma);
	double normalComplementaryErrorFunction(double x);
	double normalComplementaryErrorCheb(double x);

	//binomial
	double getBiomialRand(double pp, long n);
	double factorial(int n);
	double factorialLn(int n);
	double binomCoeffLn(int n, int k);
	int binomCoeff(int n, int k);
	double binomDensity(int n, int k, double p);
	double binomPValue(int k, int l);

	//gamma
	double gammaln(double z);

	double getGammaRand(double a, double b);
	double getGammaRand(double a);
	double getGammaRand(int ia);

	double gammaLogDensityFunction(double x, double alpha, double beta);
	double gammaCumulativeDistributionFunction(double x, double alpha, double beta);
	double lowerIncompleteGamma(double alpha, double z);
	double upperIncompleteGamma(double alpha, double z);

	//chisquared
	double getChisqRand(double a);

	//beta
	double getBetaDensity(double x, double alpha, double beta);
	double getBetaRandom (double alpha, double beta, double a, double b);
	double getBetaRandom (double alpha, double beta);
	double incompleteBeta(const double a, const double b, const double x);
	double inverseIncompleteBeta(double p, double a, double b);

	//Poisson
	double getPoissonRandom(const double & lambda);

	//exponential
	double getExponentialRandom(const double & lambda);
	double exponentialCumulativeFunction(const double & x, const double & lambda);

	//generalized Pareto
	double getGeneralizedParetoRand(const double & locationMu, const double & scaleSigma, const double & shapeXi);
	double GeneralizedParetoDensityFunction(const double & x, const double & locationMu, const double & scaleSigma, const double & shapeXi);
	double GeneralizedParetoCumulativeFunction(const double & x, const double & locationMu, const double & scaleSigma, const double & shapeXi);

	//geometric
	int getGeometricRandom(double & p);
};

#endif /* TRandomGenerator_H_ */
