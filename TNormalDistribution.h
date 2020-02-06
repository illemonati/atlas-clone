/*
 * TNormalDistribution.h
 *
 *  Created on: Feb 6, 2020
 *      Author: wegmannd
 */

#ifndef TNORMALDISTRIBUTION_H_
#define TNORMALDISTRIBUTION_H_

#include <math.h>

//Most methods from "Numerical Recipes"

class TNormalDistribution{
private:

	double _coefficientsChebyshevApproxErrorFunction[28] = {-1.3026537197817094, 6.4196979235649026e-1,
			   1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
			   3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
			   -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
			   6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
			   9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
			   -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};


	double _errorFunctionChebyshev(double x){
		int j;
		double t, ty, tmp, d = 0., dd = 0.;
		t = 2. / (2. + x);
		ty = 4. * t - 2;
		for(j=27; j>0; --j){
		  tmp = d;
		  d = ty * d - dd + _coefficientsChebyshevApproxErrorFunction[j];
		  dd = tmp;
		}
		return t * exp(-x * x + 0.5 * (_coefficientsChebyshevApproxErrorFunction[0] + ty * d) -dd);
	};

public:

	TNormalDistribution(){
		//set coefficients for error function approximation


	};

	double density(double x){
		return 0.398942280401432678 * exp(-0.5 * x * x);
	};

	double density(double x, const double mean, const double sigma){
		double tmp = (x - mean) / sigma;
		return 0.398942280401432678 / sigma * exp(-0.5 * tmp * tmp);
	};

	double cumulativeDensity(const double x){
	  return 0.5 * complementaryErrorFunction(-0.707106781186547524 * x);
	};

	double cumulativeDensity(const double x, const double mean, const double sigma){
	  return 0.5 * complementaryErrorFunction(-0.707106781186547524 * (x - mean) / sigma);
	};

	double errorFunction(const double x){
		if(x >= 0.) return 1.0 - _errorFunctionChebyshev(x);
		else return _errorFunctionChebyshev(-x) - 1.0;
	};

	double complementaryErrorFunction(const double x){
		if(x >= 0.) return _errorFunctionChebyshev(x);
		else return 2.0 - _errorFunctionChebyshev(-x);
	};

};




#endif /* TNORMALDISTRIBUTION_H_ */
