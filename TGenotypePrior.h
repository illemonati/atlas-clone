/*
 * TGenotypePrior.h
 *
 *  Created on: Nov 20, 2018
 *      Author: phaentu
 */

#ifndef TGENOTYPEPRIOR_H_
#define TGENOTYPEPRIOR_H_

#include "TThetaEstimator.h"
#include "TWindow.h"

//---------------------------------------------------------------------------------
// TGenotypePrior
// An class used to serve genotype priors to Bayesian inferences, such as caller
//---------------------------------------------------------------------------------
class TGenotypePrior{
protected:
	double genotypePrior[10];

public:
	TGenotypePrior();
	virtual ~TGenotypePrior();

	virtual void update(TWindow* window, const std::string chrName, TLog* logfile);
	double* getPointerToPrior();
};


class TGenotypePriorUniform:public TGenotypePrior{
public:
	TGenotypePriorUniform();
};

class TGenotypePriorFixedTheta:public TGenotypePrior{
private:
	TThetaEstimator* thetaEstimator;
	bool equalBaseFreq;

public:
	TGenotypePriorFixedTheta(double theta, bool EqualBaseFreq, TLog* logfile);

	~TGenotypePriorFixedTheta();

	void update(TWindow* window, const std::string chrName, TLog* logfile);
};

class TGenotypePriorTheta:public TGenotypePrior{
private:
	TThetaEstimator* thetaEstimator;
	TThetaOutputFile out;
	TLog* logfile;
	double defaultTheta;
	bool hasDefaultTheta;

	void init(TParameters & parameters, std::string & thetaOutputName, TLog* Logfile);

public:
	TGenotypePriorTheta(TParameters & parameters, std::string thetaOutputName, TLog* logfile);
	TGenotypePriorTheta(TParameters & parameters, std::string thetaOutputName, double DefaultTheta, TLog* logfile);
	~TGenotypePriorTheta();

	void update(TWindow* window, const std::string chrName, TLog* logfile);
};




#endif /* TGENOTYPEPRIOR_H_ */
