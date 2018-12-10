/*
 * math_functions.h
 *
 *  Created on: Dec 10, 2018
 *      Author: phaentu
 */

#ifndef MATH_FUNCTIONS_H_
#define MATH_FUNCTIONS_H_

#include <math.h>

namespace math{

//------------------------------------------------
// Log gamma function
// AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
//------------------------------------------------
double gammaln(const double z){
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
};

//----------------------------------------------------------------
// Factorials, see Numerical recipes (third edition), chapter 6.1
//----------------------------------------------------------------
double factorial(int n){
	//initialize lookup table
	static double factorialTable[171];
	static bool needsInitialization = true;
	if(needsInitialization){
		factorialTable[0] = 1.0;
		for(int i=1; i<171; i++)
			factorialTable[i] = factorialTable[i-1]*i;
		needsInitialization = false;
	}

	//check range
	if(n < 0 || n > 170){
		std::ostringstream tos;
		tos << n;
		throw "factorial: n = " + tos.str() + " out of range!";
	}

	//return
	return factorialTable[n];
};

double factorialLog(int n){
	//initialize lookup table
	static int tableSize = 2000;
	static double factorialLogTable[tableSize];
	static bool needsInitialization = true;
	if(needsInitialization){
		factorialLogTable[0] = 0.0;
		for(int i=1; i<tableSize; i++)
			factorialLogTable[i] = math::gammaln(i+1);
		needsInitialization = false;
	}

	//check range
	if(n < 0){
		std::ostringstream tos;
		tos << n;
		throw "factorial: n = " + tos.str() + "out of range!";
	}

	//return
	if(n < tableSize) return factorialLogTable[n];
	return math::gammaln(n + 1.0);
};

//------------------------------------------------
// choose function (binomial coefficient)
//------------------------------------------------
double chooseLog(const int n, const int k){
	return factorialLog(n) - factorialLog(k) - factorialLog(n-k);
};

int choose(const int n, const int k){
	return factorial(n) - factorial(k) - factorial(n-k);
};


//------------------------------------------------
//end of name space
}
//------------------------------------------------


#endif /* MATH_FUNCTIONS_H_ */
