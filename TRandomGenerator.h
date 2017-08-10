#ifndef TRandomGenerator_H_
#define TRandomGenerator_H_

#include <cmath>
//#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>

class TRandomGenerator{
private:
	long _Idum;
	double ran3();
	double _lowerIncompleteGamma(double alpha, double z);
	double _upperIncompleteGamma(double alpha, double z);
	long get_randomSeedFromCurrentTime(long & addToSeed);

public:
	long usedSeed;
	double* factorialTable;
	bool factorialTableInitialized;
	double* factorialTableLn;
	bool factorialTableLnInitialized;

	TRandomGenerator(long addToSeed){
		init(addToSeed);
	};
	TRandomGenerator(long addToSeed, bool seedIsFixed){
		if(!seedIsFixed) init(addToSeed);
		else {
			if(addToSeed<0) addToSeed=-addToSeed;
	        usedSeed=addToSeed;
	        _Idum=-addToSeed;
		}
	};
	TRandomGenerator(){
			init(0);
	};
	~TRandomGenerator(){};
	void init(long addToSeed);

	//uniform
	double getRand(){ return ran3(); };
	double getRand(double min, double max);
	int getRand(int min, int maxPlusOne);
	int pickOne(int numElements);
	int pickOne(int numElements, float* probsCumulative);
	int pickOne(int numElements, double* probsCumulative);
	long getRand(long min, long maxPlusOne);

	//normal
	double getNormalRandom (double dMean, double dStdDev);
	double normalCumulativeDistributionFunction(double x, double mean, double sigma);
	double normalComplementaryErrorFunction(double x);
	double normalComplementaryErrorCheb(double x);

	//binomial
	double getBiomialRand(double pp, int n);
	double factorial(int n);
	double factorialLn(int n);
	double binomCoeffLn(int n, int k);
	int binomCoeff(int n, int k);
	double binomDensity(int n, int k, double p);

	//gamma
	double gammaln(double z);

	double getGammaRand(double a, double b);
	double getGammaRand(double a);
	double getGammaRand(int ia);

	double gammaCumulativeDistributionFunction(double x, double alpha, double beta);
	double lowerIncompleteGamma(double alpha, double z);
	double upperIncompleteGamma(double alpha, double z);

	//beta
	double getBetaRandom (double alpha, double beta, double a, double b);
	double getBetaRandom (double alpha, double beta);

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
