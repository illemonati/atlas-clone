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
	binomPValueTableInitialized = false;
	binomPValueTable = NULL;
	binomPValueTableSize = 100;
};

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
	return exp(choose(n,k) + k*log(p) + (n-k)*log(1.0-p));
}

double TRandomGenerator::_binomPValue(const int & k, const int & l){
	double cumul = 0.0;
	int n = k + l;
	for(int i = 0; i <= std::min(k,l); i++){
		cumul += exp(chooseLog(n, i) + -0.6931472*n); // -0.6931472 is log(0.5)
	}
	return cumul;
};

double TRandomGenerator::binomPValue(const int k, const int l){
	int n = l + k;

	if (n < binomPValueTableSize){
		if(!binomPValueTableInitialized){
			binomPValueTable = new double*[binomPValueTableSize];

			//loop over all entries
			for(int i=0; i<binomPValueTableSize; i++) {
				//get row size
				int rowLength = floor(i/2) + 1;

				//allocate memory
				binomPValueTable[i] = new double[rowLength];

				//fill row
				for(int j=0; j<rowLength; j++)
					binomPValueTable[i][j] = _binomPValue(j, i-j);
			}
			binomPValueTableInitialized = true;
		}
		return binomPValueTable[n][std::min(k,l)];
	}
	else return _binomPValue(k, l);
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

char TRandomGenerator::getRandomAlphaNumerichCharacter(){
	//choose a random letter or number
	//numbers are 48 - 57 (10)
	//uppercase letters are 65 - 90 (26)
	//lowercase letters are 97 - 122 (26)
	//in total: 62
	int index = getRand(48, 110);
	if(index > 57) index += 8;
	if(index > 90) index += 7;

	return (char) index;
};

std::string TRandomGenerator::getRandomAlphaNumericString(const int length){
	std::string s;
	for(int i=0; i<length; ++i){
		s+= getRandomAlphaNumerichCharacter();
	}
	return s;
};

void TRandomGenerator::shuffle(std::vector<int> & vec){
	int j;

	for(int i=vec.size()-1; i>0; --i){
		j = pickOne(i+1);
		std::swap(vec[i], vec[j]);
	}
};

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

double TRandomGenerator::getLogNormalRandom (double dMean, double dStdDev){
    double r = getNormalRandom (dMean, dStdDev);
    return exp(r);
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
//Chi-Square Distribution
//--------------------------------------------------------

/* Numerical Recipes 3rd edition, pp. 371: Chisquare(v) = Gamma(v/2, 1/2) = 2Gamma(v/2, 1)
*/

double TRandomGenerator::getChisqRand(double a){
	   if(a <= 0) {
		   throw "Negative value received in getChisqRand()!";
	   }
	   return getGammaRand(a/2)/0.5;
}

//--------------------------------------------------------
//Beta Distribution
//--------------------------------------------------------

double TRandomGenerator::getBetaDensity(double x, double alpha, double beta){
	// this is not an approximation, but the exact beta density function (as defined in https://en.wikipedia.org/wiki/Beta_distribution)
	if (x < 0 || x > 1)
		throw "Input value x should be 0 <= x <= 1!";
	if (alpha < 0 || beta < 0)
		throw "Parameters alpha and beta should be > 0!";
	double tmp1 = tgamma(alpha + beta) / (tgamma(alpha) * tgamma(beta));
	double tmp2 = pow(x, alpha - 1.) * pow (1. - x, beta - 1.);
	return tmp1 * tmp2;
}

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

/* Incomplete Beta functions
 * Reference: Press et al., Numerical Recipes, 3rd edition, pp. 270-273
 */


double TRandomGenerator::incompleteBeta(const double a, const double b, const double x) {
	// returns incomplete beta function I_x(a, b) for positive a and b, and x between 0 and 1
	double bt;
	static const int SWITCH = 3000; // when to switch to quadrature method
	if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
	if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) return x;
	if (a > SWITCH && b > SWITCH) return _betaiapprox(a,b,x);
	bt=exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) return bt*_betacf(a,b,x)/a;
	else return 1.0-bt*_betacf(b,a,1.0-x)/b;
}

double TRandomGenerator::_betacf(const double a, const double b, const double x) {
	// evaluates continued fraction for incomplete beta function by modified Lentz's method
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	const double EPS = std::numeric_limits<double>::epsilon();
	const double FPMIN = std::numeric_limits<double>::min()/EPS;
	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<10000;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	return h;
}

//
double TRandomGenerator::_betaiapprox(double a, double b, double x) {
	// incomplete beta by quadrature; returns I_x(a, b)
	int j;
	double xu,t,sum,ans;
	double a1 = a-1.0, b1 = b-1.0, mu = a/(a+b);
	double lnmu=log(mu),lnmuc=log(1.-mu);
	// y and w are coefficients for Gauss-Legendre quadrature
	static const double y[18] = {0.0021695375159141994,
	0.011413521097787704,0.027972308950302116,0.051727015600492421,
	0.082502225484340941, 0.12007019910960293,0.16415283300752470,
	0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
	0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
	0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
	0.87126389619061517, 0.95698180152629142};
	static const double w[18] = {0.0055657196642445571,
	0.012915947284065419,0.020181515297735382,0.027298621498568734,
	0.034213810770299537,0.040875750923643261,0.047235083490265582,
	0.053244713977759692,0.058860144245324798,0.064039797355015485,
	0.068745323835736408,0.072941885005653087,0.076598410645870640,
	0.079687828912071670,0.082187266704339706,0.084078218979661945,
	0.085346685739338721,0.085983275670394821};

	t = sqrt(a*b/(pow(a+b, 2)*(a+b+1.0)));
	if (x > a/(a+b)) {
		if (x >= 1.0) return 1.0;
		xu = std::min(1.,std::max(mu + 10.*t, x + 5.0*t));
	} else {
		if (x <= 0.0) return 0.0;
		xu = std::max(0.,std::min(mu - 10.*t, x - 5.0*t));
	}
	sum = 0;
	for (j=0;j<18;j++) {
		t = x + (xu-x)*y[j];
		sum += w[j]*exp(a1*(log(t)-lnmu)+b1*(log(1-t)-lnmuc));
	}
	ans = sum*(xu-x)*exp(a1*lnmu-gammaln(a)+b1*lnmuc-gammaln(b)+gammaln(a+b));
	return ans>0.0? 1.0-ans : -ans;
}

double TRandomGenerator::inverseIncompleteBeta(double p, double a, double b){
	// inverse of incomplete beta function
	// Returns x such that I_x(a,b) = p for argument p between 0 and 1. Uses Halley's method
	const double EPS = 1.e-8;
	double pp,t,u,err,x,al,h,w,afac,a1=a-1.,b1=b-1.;
	int j;
	if (p <= 0.) return 0.;
	else if (p >= 1.) return 1.;
	else if (a >= 1. && b >= 1.) {
		pp = (p < 0.5)? p : 1. - p;
		t = sqrt(-2.*log(pp));
		x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
		if (p < 0.5) x = -x;
		al = (pow(x, 2)-3.)/6.;
		h = 2./(1./(2.*a-1.)+1./(2.*b-1.));
		w = (x*sqrt(al+h)/h)-(1./(2.*b-1)-1./(2.*a-1.))*(al+5./6.-2./(3.*h));
		x = a/(a+b*exp(2.*w));
	} else {
		double lna = log(a/(a+b)), lnb = log(b/(a+b));
		t = exp(a*lna)/a;
		u = exp(b*lnb)/b;
		w = t + u;
		if (p < t/w) x = pow(a*w*p,1./a);
		else x = 1. - pow(b*w*(1.-p),1./b);
	}
	afac = -gammaln(a)-gammaln(b)+gammaln(a+b);
	for (j=0;j<10;j++) {
		if (x == 0. || x == 1.) return x;
		err = incompleteBeta(a,b,x) - p;
		t = exp(a1*log(x)+b1*log(1.-x) + afac);
		u = err/t;
		x -= (t = u/(1.-0.5*std::min(1.,u*(a1/x - b1/(1.-x)))));
		if (x <= 0.) x = 0.5*(x + t);
		if (x >= 1.) x = 0.5*(x + t + 1.);
		if (fabs(t) < EPS*x && j > 0) break;
	}
	return x;
}


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

double TRandomGenerator::getExponentialRandomTruncated(const double & lambda, const double & lowerBound, const double & upperBound){
	//copied from stack overflow
	return -log(exp(-lambda * lowerBound) - (exp(-lambda * lowerBound) - exp(-lambda * upperBound)) * getRand()) / lambda;
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
