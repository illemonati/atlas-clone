/*
 * TRandomGenerator.cpp
 *
 *  Created on: Sep 24, 2009
 *      Author: wegmannd
 */

#include "TRandomGenerator.h"

//---------------------------------------------------------------------------------------
void TRandomGenerator::setSeed(long addToSeed, bool seedIsFixed){
	if(seedIsFixed){
		if(addToSeed<0) addToSeed=-addToSeed;
        usedSeed = addToSeed;
        _Idum = -addToSeed;
	} else {
		_Idum = get_randomSeedFromCurrentTime(addToSeed);
		if(_Idum==0) _Idum=get_randomSeedFromCurrentTime(_Idum);
		if(_Idum < 0) _Idum=-_Idum;
		if(_Idum > 161803398) _Idum = _Idum % 161803397;
		usedSeed = _Idum;
		_Idum = -_Idum;
	}
};

long TRandomGenerator::get_randomSeedFromCurrentTime(long & addToSeed){
	struct timeval time;
	gettimeofday(&time, 0);
	long microseconds = time.tv_sec * 1000 + time.tv_usec / 1000;

	//microseconds = microseconds - 38*365*24*3600; //substract 38 years...
	microseconds = (double) (addToSeed + microseconds);
   return microseconds ;
}

void TRandomGenerator::init(){
	factorialTable = NULL;
	factorialTableInitialized = false;
	factorialTableLn = NULL;
	factorialTableLnInitialized = false;
}

#define MBIG 1000000000L
#define MSEED 161803398L
#define MZ 0
#define FAC (1.0/MBIG)

double TRandomGenerator::ran3(){
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if ((_Idum < 0) || (iff == 0) ) {
                iff=1;
                mj=MSEED-(_Idum < 0 ? -_Idum : _Idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;++i) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;++k)
                        for (i=1;i<=55;++i) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                _Idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

//--------------------------------------------------------
//Binomial Distribution
//--------------------------------------------------------
double TRandomGenerator::getBiomialRand(double pp, long n){
	int j;
	static int nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;

	if (n < 25){
		bnl=0;
		for (j=1;j<=n;j++)
			if (ran3() < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran3();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammaln(en+1.0);
			nold=n;
		}
		if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);

		do {
			do {
				angle=3.141592654*ran3();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));

			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammaln(em+1.0)
					-gammaln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran3() > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}

double TRandomGenerator::factorial(int n){
	if(n < 0 || n > 170){
		std::ostringstream tos;
		tos << n;
		throw "TRandomGenerator::factorial: n = " + tos.str() + " out of range!";
	}

	if(!factorialTableInitialized){
		factorialTable = new double[171];
		factorialTable[0] = 1.0;
		for(int i=1; i<171; i++)
			factorialTable[i] = factorialTable[i-1]*i;
		factorialTableInitialized = true;
	}

	return factorialTable[n];
}

double TRandomGenerator::factorialLn(int n){
	static const int TABLESIZE = 2000;
	if(n < 0){
		std::ostringstream tos;
		tos << n;
		throw "TRandomGenerator::factorial: n = " + tos.str() + "out of range!";
	}

	if(n < TABLESIZE){
		if(!factorialTableLnInitialized){
			factorialTableLn = new double[TABLESIZE];
			factorialTableLn[0] = 0.0;
			for(int i=1; i<TABLESIZE; i++)
				factorialTableLn[i] = gammaln(i+1);
			factorialTableLnInitialized = true;
		}
		return factorialTableLn[n];
	}
	return gammaln(n+1);
}


double TRandomGenerator::binomCoeffLn(int n, int k){
	return factorialLn(n) - factorialLn(k) - factorialLn(n-k);
}

int TRandomGenerator::binomCoeff(int n, int k){
	return factorial(n) - factorial(k) - factorial(n-k);
}

double TRandomGenerator::binomDensity(int n, int k, double p){
	if(p == 0.0){
		if(k == 0) return 1.0;
		else return 0.0;
	}
	if(p == 1.0){
		if(k == n) return 1.0;
		else return 0.0;
	}
	return exp(binomCoeffLn(n,k) + k*log(p) + (n-k)*log(1.0-p));
}

double TRandomGenerator::binomPValue(int k, int l){
	double cumul = 0.0;
	double logHalf = -0.6931472 ;  // = log(0.5);
	int n = l+k;
	if(k < l){
		for(unsigned int i = 0; i <= k; ++i)
			cumul += exp(binomCoeffLn(n, i) + logHalf*n);
	} else {
		for(unsigned int i = k; i <= n; ++i)
			cumul += exp(binomCoeffLn(n, i) + logHalf*n);
	}
	return cumul;
}

//--------------------------------------------------------
//Uniform Distribution
//--------------------------------------------------------
double TRandomGenerator::getRand(double min, double max){
	//return a random number between min and max
	 return (getRand() * (max - min) + min);
}
int TRandomGenerator::getRand(int min, int maxPlusOne){
	//return an random integer between min and maxPlusOne-1
	float r=1.0;
	while(r==1.0) r=getRand(); //we have a number in [0,1[
	return min+floor(r*(maxPlusOne-min));
}

int TRandomGenerator::pickOne(int numElements){
	//if(numElements < 1) throw "TRandomGenerator::pickOne: can not choose an element among less than 1 elements!";
	if(numElements == 1) return 0;
	float r = 1.0;
	while(r == 1.0) r=getRand(); //we have a number in [0,1[
	return floor(r*(numElements));
}

int TRandomGenerator::pickOne(int numElements, float* probsCumulative){
	if(numElements == 1) return 0;
	float r = 1.0;
	//while(r == 1.0) r=getRand(); //we have a number in [0,1[
	r = getRand();
	int i = 0;
	while(r > probsCumulative[i])
		++i;
	return i;
}

int TRandomGenerator::pickOne(int numElements, double* probsCumulative){
	if(numElements == 1) return 0;
	double r = 1.0;
	//while(r == 1.0) r=getRand(); //we have a number in [0,1[
	r = getRand();
	int i = 0;
	while(r > probsCumulative[i])
		++i;
	return i;
}

long TRandomGenerator::getRand(long min, long maxPlusOne){
	//return an random integer between min and maxPlusOne-1
	double r=1.0;
	while(r==1.0) r=getRand(); //we have a number in [0,1[
	return min+floor(r*(maxPlusOne-min));
}

//--------------------------------------------------------
//Normal Random and associated functions
//--------------------------------------------------------
/* Returns a Normal random variate based on a unit variate,
   using a random generator as a source of uniform deviates.
   Adapted from the algorithm described in the book Numerical Recipes by
   Press et al.
*/
double TRandomGenerator::getNormalRandom (double dMean, double dStdDev){
  double w, x1, x2;
   do {
      x1 = 2. * ran3() - 1.;
      x2 = 2. * ran3() - 1.;
      w = x1*x1 + x2*x2;
   } while (w >= 1. || w < 1E-30);

   w = sqrt((-2.*log(w))/w);
   x1 *= w;
   return (x1 * dStdDev + dMean);
}

double TRandomGenerator::normalCumulativeDistributionFunction(double x, double mean, double sigma){
  if(x==0.) return 0.;
  return 0.5*normalComplementaryErrorFunction(-0.707106781186547524*(x-mean)/sigma);
}
double TRandomGenerator::normalComplementaryErrorFunction(double x){
   //see book "numerical recipes"
   if(x>=0.) return normalComplementaryErrorCheb(x);
   else return 2.0- normalComplementaryErrorCheb(-x);
}
double TRandomGenerator::normalComplementaryErrorCheb(double x){
   //see book "numerical recipes"
   double coef[28]={-1.3026537197817094, 6.4196979235649026e-1,
   1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
   3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
   -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
   6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
   9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
   -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};
   int j;
   double t, ty, tmp, d=0., dd=0.;
   t=2./(2.+x);
   ty=4.*t-2;
   for(j=27;j>0;--j){
	  tmp=d;
	  d=ty*d-dd+coef[j];
	  dd=tmp;
   }
   return t*exp(-x*x+0.5*(coef[0]+ty*d)-dd);
}

//--------------------------------------------------------
//Gamma Distribution
//--------------------------------------------------------
/* Log gamma function
* \log{\Gamma(z)}
* AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
*/

double TRandomGenerator::gammaln(double z){
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526 / (z+5);
	x += 12.50734324009056 / (z+4);
	x -= 176.6150291498386 / (z+3);
	x += 771.3234287757674 / (z+2);
	x -= 1259.139216722289 / (z+1);
	x += 676.5203681218835 / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

/*Random deviates from standard gamma distribution with density
         a-1
        x    exp[ -x ]
f(x) = ----------------
         Gamma[a]

where a is the shape parameter.  The algorithm for integer a comes
from numerical recipes, 2nd edition, pp 292-293.  The algorithm for
a<1 uses code from p 213 of Statistical Computing, by Kennedy and
Gentle, 1980 edition.  This algorithm was originally published in:

Ahrens, J.H. and U. Dieter (1974), "Computer methods for sampling from
Gamma, Beta, Poisson, and Binomial Distributions".  COMPUTING
12:223-246.

The mean and variance of these values are both supposed to equal a.
My tests indicate that they do.

This algorithm has problems when a is small.  In single precision, the
problem  arises when a<0.1, roughly.  That is why I have declared
everything as double below.  Trouble is, I still don't know how small
a can be without causing trouble.  Mean and variance are ok at least
down to a=0.01.  f(x) doesn't seem to have a series expansion around
x=0.
****************************************************************/
double TRandomGenerator::getGammaRand(double a, double b) {
   if(b <= 0) {
	   throw "Negative value received in getGammaRand()!";
   }
   return getGammaRand(a)/b;
}


double TRandomGenerator::getGammaRand(double a) {

  int ia;
  double u, b, p, x, y=0.0, recip_a;

  if(a <= 0) {
	  throw "Negative value received in getGammaRand()!";
  }

  ia = (int) floor(a);  /* integer part */
  a -= ia;        /* fractional part */
  if(ia > 0) {
    y = getGammaRand(ia);  /* gamma deviate w/ integer argument ia */
    if(a==0.0) return(y);
  }

  /* get gamma deviate with fractional argument "a" */
  b = (M_E + a)/M_E;
  recip_a = 1.0/a;
  for(;;) {
    u = ran3();
    p = b*u;
    if(p > 1) {
      x = -log((b-p)/a);
      if( ran3() > pow(x, (double) a-1.0)) continue;
      break;
    }
    else {
      x = pow(p, recip_a);
      if( ran3() > exp(-x)) continue;
      break;
    }
  }
  return(x+y);
}

double TRandomGenerator::getGammaRand(int ia){
	/****************************************************************
	gamma deviate for integer shape argument.  Code modified from pp
	292-293 of Numerical Recipes in C, 2nd edition.
	****************************************************************/
  int j;
  double am,e,s,v1,v2,x,y;

  if (ia < 1) throw "Argument below 1 in getGammaRand()!";
  if (ia < 6){
    x=1.0;
    for (j=0; j<ia; j++)
      x *= ran3();
    x = -log(x);
  } else {
    do {
      do {
	    do{                         /* next 4 lines are equivalent */
	    	v1=2.0*ran3()-1.0;       /* to y = tan(Pi * uni()).     */
	    	v2=2.0*ran3()-1.0;
	    } while (v1*v1+v2*v2 > 1.0);
	    y=v2/v1;
	    am=ia-1;
	    s=sqrt(2.0*am+1.0);
	    x=s*y+am;
      } while (x <= 0.0);
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    } while (ran3() > e);
  }
  return(x);
}

//gamma log density function
double TRandomGenerator::gammaLogDensityFunction(double x, double alpha, double beta){
	return alpha * log(beta) - gammaln(alpha) + (alpha-1.0)*log(x) - beta * x;
}

//Functions to calculate cumulative of Gamma
//Adapted from kfunc.c of samtools
double TRandomGenerator::gammaCumulativeDistributionFunction(double x, double alpha, double beta){
	return lowerIncompleteGamma(alpha, beta * x); // do not need to divide by exp(gammaln(alpha)), it is already regularized;
}

/* The following computes regularized incomplete gamma functions.
 * Taken from kfunc.c from samtools. Original comment:
* Formulas are taken from Wiki, with additional input from Numerical
* Recipes in C (for modified Lentz's algorithm) and AS245
* (http://lib.stat.cmu.edu/apstat/245).
*/
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

double TRandomGenerator::lowerIncompleteGamma(double alpha, double z){
	if(z <= 1.0 || z < alpha) return _lowerIncompleteGamma(alpha, z);
	else return 1.0 - _upperIncompleteGamma(alpha, z);
}
double TRandomGenerator::upperIncompleteGamma(double alpha, double z){
	if(z <= 1.0 || z < alpha) return 1.0 - _lowerIncompleteGamma(alpha, z);
	else return _upperIncompleteGamma(alpha, z);
}

// regularized lower incomplete gamma function, by series expansion
double TRandomGenerator::_lowerIncompleteGamma(double alpha, double z){
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (alpha + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return exp(alpha * log(z) - z - gammaln(alpha + 1.0) + log(sum));
}

// regularized upper incomplete gamma function, by continued fraction
double TRandomGenerator::_upperIncompleteGamma(double alpha, double z){
	int j;
	double C, D, f;
	f = 1. + z - alpha; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	// See Numerical Recipes in C, 2nd edition, section 5.2
	for (j = 1; j < 100; ++j) {
		double a = j * (alpha - j), b = (j<<1) + 1 + z - alpha, d;
		D = b + a * D;
		if (D < KF_TINY) D = KF_TINY;
		C = b + a / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(alpha * log(z) - z - gammaln(alpha) - log(f));
}


//--------------------------------------------------------
//Beta Distribution
//--------------------------------------------------------
/*
   returns a variate that is Beta distributed on the interval [a,b]
   with shape parameters alpha and beta.

   The Beta function has two shaping parameters, alpha and beta.
   Setting these parameters to 1.5 and 1.5 yields a normal-like
   distribution, but without tails. If alpha and beta are equal to
   1 it is a uniform distribution.

   If alpha and beta are less than 1, use a rejection algorithm;
   Otherwise use the fact that if x is distributed Gamma(alpha) and y
   Gamma(beta) then x/(x+y) is Beta(alpha, beta).

   The rejection algorithm first a Beta variate is found over the
   interval [0, 1] with not the most efficient algorithm.  This is then
   scaled at the end to desired range.

   It may be tempting to re-use the second number drawn as the first
   random number of the next iteration, and simply draw one more.
   *** Don't do it.  You will produce an incorrect distribution.  You
   must draw two new numbers for the rejection sampling to be correct.

   References:
   - Ripley, Stochastic Simulations, John Wiley and Sons, 1987, p 90.
   - J.H.Maindonald, Statistical Computation, John Wiley and Sons,
     1984, p 370.
*/
double TRandomGenerator::getBetaRandom (double alpha, double beta, double a, double b){
  if (b <= a) throw "Bad shape or range for a beta variate!";
  return (a + getBetaRandom(alpha, beta) * (b-a));   /* Scale to interval [a, b] */
}


double TRandomGenerator::getBetaRandom (double alpha, double beta) {
  double u1, u2, w;

  if (alpha <= 0 || beta <= 0) throw "Bad shape or range for a beta variate!";

  if ((alpha < 1) && (beta < 1))
    /* use rejection */
    do {
      u1 = ran3(); /* Draw two numbers */
      u2 = ran3();

      u1 = pow(u1, (double) 1.0/alpha); /* alpha and beta are > 0 */
      u2 = pow(u2, (double) 1.0/beta);

      w = u1 + u2;

    } while (w > 1.0);

  else {
    /* use relation to Gamma */
    u1 = getGammaRand(alpha);
    u2 = getGammaRand(beta);
    w  = u1 + u2;
  }

  return u1/w;

} /* BetaRandom */

//--------------------------------------------------------
//Poisson Distribution
//--------------------------------------------------------
//Poisson Ranodm from Wikipedia. Works only for small Lambda and may be a bit inefficient.
double TRandomGenerator::getPoissonRandom(const double & lambda){
	double L=exp(-lambda);
	int k=0;
	double p=1;
	do{
		++k;
		p*=getRand();
	} while (p>L);
	return k-1;
}

//--------------------------------------------------------
//Exponential Distribution
//--------------------------------------------------------
//Exponential Ranodm from Wikipedia
double TRandomGenerator::getExponentialRandom(const double & lambda){
	return -log(getRand()) / lambda;
}

double TRandomGenerator::exponentialCumulativeFunction(const double & x, const double & lambda){
	return 1.0 - exp(-lambda * x);
}

//--------------------------------------------------------
//Exponential Distribution
//--------------------------------------------------------
double TRandomGenerator::getGeneralizedParetoRand(const double & locationMu, const double & scaleSigma, const double & shapeXi){
	if(shapeXi == 0.0){
		return locationMu - scaleSigma * log(getRand());
	} else {
		return locationMu + scaleSigma * (pow(getRand(), -shapeXi) - 1.0) / shapeXi;
	}
}

double TRandomGenerator::GeneralizedParetoDensityFunction(const double & x, const double & locationMu, const double & scaleSigma, const double & shapeXi){
	if(locationMu == 0.0){
	   return -log(scaleSigma) - (x - locationMu) / scaleSigma;
	} else {
	   return -log(scaleSigma) + (-1.0-1.0/shapeXi) * log(1.0 + shapeXi * (x - locationMu )/scaleSigma);
	}
}


double TRandomGenerator::GeneralizedParetoCumulativeFunction(const double & x, const double & locationMu, const double & scaleSigma, const double & shapeXi){
	if(x<locationMu) throw "Problem calculating cumulative function for generalized pareto: x < mu (location)!";
	if(shapeXi < 0.0 && x > locationMu - scaleSigma / shapeXi) throw "Problem calculating cumulative function for generalized pareto: x > mu - sigma/xi!";
	if(shapeXi == 0.0){
		return 1.0 - exp((locationMu - x)/scaleSigma);
	} else {
		return 1.0 - pow(1.0 + shapeXi*(x - locationMu)/scaleSigma, -1.0/shapeXi);
	}
}

//--------------------------------------------------------
//Geometric Distribution
//--------------------------------------------------------
int TRandomGenerator::getGeometricRandom(double & p){
	return log(getRand())/log(1.0-p);
}
