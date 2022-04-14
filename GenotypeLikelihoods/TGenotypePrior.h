/*
 * TGenotypePrior.h
 *
 *  Created on: May 15, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEPRIOR_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEPRIOR_H_

#include "TThetaEstimator.h"
#include "TWindow.h"

namespace GenotypeLikelihoods{

//---------------------------------------------------------------------------------
// TGenotypePrior
// An class used to serve genotype priors to Bayesian inferences, such as caller
//---------------------------------------------------------------------------------
class TGenotypePrior{
protected:
	TGenotypeProbabilities genotypePrior;

public:
	TGenotypePrior(){};

	virtual ~TGenotypePrior(){};

	virtual void update(const TWindow &, coretools::TLog*, const TGenotypeLikelihoodCalculator &){};
	TGenotypeProbabilities* getPointerToPrior(){ return &genotypePrior; };
	coretools::Probability operator[](const genometools::Genotype & genotype){ return genotypePrior[genotype]; };
};


class TGenotypePriorUniform:public TGenotypePrior{
public:
	TGenotypePriorUniform(){
		reset(genotypePrior);
	};
};

class TGenotypePriorFixedTheta:public TGenotypePrior{
private:
	TThetaEstimator* thetaEstimator;
	bool equalBaseFreq;

public:
	TGenotypePriorFixedTheta(double theta, bool EqualBaseFreq, coretools::TLog* logfile, coretools::TRandomGenerator* randomGenerator){
		thetaEstimator = new TThetaEstimator(logfile, randomGenerator);
		thetaEstimator->setTheta(theta);
		equalBaseFreq = EqualBaseFreq;
		if(equalBaseFreq){
			GenotypeLikelihoods::TBaseProbabilities freq;
			freq.fill(0.25);
			thetaEstimator->setBaseFreq(freq);
		}
		thetaEstimator->fillPGenotype(genotypePrior);
	};

	~TGenotypePriorFixedTheta(){
		delete thetaEstimator;
	};

	void update(const TWindow & window, coretools::TLog* logfile, const TGenotypeLikelihoodCalculator &){
		using genometools::Base;
		if(!equalBaseFreq){
			logfile->listFlush("Estimating base frequencies for prior ...");
			GenotypeLikelihoods::TBaseProbabilities freq = window.estimateBaseFrequencies();
			thetaEstimator->setBaseFreq(freq);
			logfile->done();
			logfile->conclude("Estimated base frequencies: " + coretools::str::toString(freq[Base::A]) + ", " +
					  coretools::str::toString(freq[Base::C]) + ", " + coretools::str::toString(freq[Base::G]) + ", " +
					  coretools::str::toString(freq[Base::T]));
			thetaEstimator->fillPGenotype(genotypePrior);
		}
	};
};

class TGenotypePriorTheta:public TGenotypePrior{
private:
	TThetaEstimator* thetaEstimator;

	TThetaOutputFile out;
	coretools::TLog* logfile;
	double defaultTheta;
	bool hasDefaultTheta;

	void init(coretools::TParameters & parameters, std::string & thetaOutputName, coretools::TLog* Logfile, coretools::TRandomGenerator* randomGenerator){
		logfile = Logfile;
		thetaEstimator = new TThetaEstimator(parameters, logfile, randomGenerator);
		out.open(thetaEstimator, thetaOutputName, logfile);
	};

public:
	TGenotypePriorTheta(coretools::TParameters & parameters, std::string thetaOutputName, coretools::TLog* logfile, coretools::TRandomGenerator* randomGenerator){
		hasDefaultTheta = false;
		defaultTheta = -1.0;

		init(parameters, thetaOutputName, logfile, randomGenerator);
	};

	TGenotypePriorTheta(coretools::TParameters & parameters, std::string thetaOutputName, double DefaultTheta, coretools::TLog* logfile, coretools::TRandomGenerator* randomGenerator){
		hasDefaultTheta = true;
		defaultTheta = DefaultTheta;
		if(defaultTheta < 0.0) throw "Theta must be >= 0.0!";
		init(parameters, thetaOutputName, logfile, randomGenerator);	};

	~TGenotypePriorTheta(){
		out.close();
		delete thetaEstimator;
	};

	void update(const TWindow & window, coretools::TLog* logfile, const TGenotypeLikelihoodCalculator & glCalculator){
		logfile->startIndent("Estimating theta for prior:");
		//clear theta estimator
		(*thetaEstimator).clear();

		//adding sites to estimator
		thetaEstimator->add(window, glCalculator);

		//estimate Theta
		if(!thetaEstimator->estimateTheta()){
			if(hasDefaultTheta){
				logfile->conclude("Will use a default theta of " + coretools::str::toString(defaultTheta) + ".");
				thetaEstimator->setTheta(defaultTheta);
				GenotypeLikelihoods::TBaseProbabilities freq;
				freq.fill(0.25);
				thetaEstimator->setBaseFreq(freq);
			} else
				throw "Please increase window size or provide a default theta!";
		}

		//update prior
		thetaEstimator->fillPGenotype(genotypePrior);
	};
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEPRIOR_H_ */
