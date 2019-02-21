/*
 * math_functions.h
 *
 *  Created on: Dec 10, 2018
 *      Author: phaentu
 */

#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include <math.h>
#include "stringFunctions.h"

//------------------------------------------------
// Log gamma function
// AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
//------------------------------------------------
double gammaln(const double z);

//----------------------------------------------------------------
// Factorials, see Numerical recipes (third edition), chapter 6.1
//----------------------------------------------------------------
double factorial(int n);
double factorialLog(int n);

//------------------------------------------------
// choose function (binomial coefficient)
//------------------------------------------------
double chooseLog(const int n, const int k);
int choose(const int n, const int k);



#endif /* MATHFUNCTIONS_H_ */
