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
	TGenotypeData genotypePrior;

public:
	TGenotypePrior(){};

	virtual ~TGenotypePrior(){};

	virtual void update(const TWindow & window, TLog* logfile){};
	TGenotypeData* getPointerToPrior(){ return &genotypePrior; };
	double operator[](const Genotype genotype){ return genotypePrior[genotype]; };
	double operator[](const uint8_t genotype){ return genotypePrior[genotype]; };
};


class TGenotypePriorUniform:public TGenotypePrior{
public:
	TGenotypePriorUniform(){
		for(int g=0; g<10; ++g)
			genotypePrior[g] = 1.0 / 10.0;
	};
};

class TGenotypePriorFixedTheta:public TGenotypePrior{
private:
	TThetaEstimator* thetaEstimator;
	bool equalBaseFreq;

public:
	TGenotypePriorFixedTheta(double theta, bool EqualBaseFreq, TLog* logfile, TRandomGenerator* randomGenerator){
		thetaEstimator = new TThetaEstimator(logfile, randomGenerator);
		thetaEstimator->setTheta(theta);
		equalBaseFreq = EqualBaseFreq;
		if(equalBaseFreq){
			GenotypeLikelihoods::TBaseData freq(0.25);
			thetaEstimator->setBaseFreq(freq);
		}
		thetaEstimator->fillPGenotype(genotypePrior);
	};

	~TGenotypePriorFixedTheta(){
		delete thetaEstimator;
	};

	void update(const TWindow & window, TLog* logfile){
		if(!equalBaseFreq){
			logfile->listFlush("Estimating base frequencies for prior ...");
			GenotypeLikelihoods::TBaseData freq;
			window.estimateBaseFrequencies(freq);
			thetaEstimator->setBaseFreq(freq);
			logfile->done();
			logfile->conclude("Estimated base frequencies: " + toString(freq[A])+ ", " + toString(freq[C]) + ", " + toString(freq[G]) + ", " + toString(freq[T]));
			thetaEstimator->fillPGenotype(genotypePrior);
		}
	};
};

class TGenotypePriorTheta:public TGenotypePrior{
private:
	TThetaEstimator* thetaEstimator;

	TThetaOutputFile out;
	TLog* logfile;
	double defaultTheta;
	bool hasDefaultTheta;

	void init(TParameters & parameters, std::string & thetaOutputName, TLog* Logfile, TRandomGenerator* randomGenerator){
		logfile = Logfile;
		thetaEstimator = new TThetaEstimator(parameters, logfile, randomGenerator);
		out.open(thetaEstimator, thetaOutputName, logfile);
	};

public:
	TGenotypePriorTheta(TParameters & parameters, std::string thetaOutputName, TLog* logfile, TRandomGenerator* randomGenerator){
		hasDefaultTheta = false;
		defaultTheta = -1.0;

		init(parameters, thetaOutputName, logfile, randomGenerator);
	};

	TGenotypePriorTheta(TParameters & parameters, std::string thetaOutputName, double DefaultTheta, TLog* logfile, TRandomGenerator* randomGenerator){
		hasDefaultTheta = true;
		defaultTheta = DefaultTheta;
		if(defaultTheta < 0.0) throw "Theta must be >= 0.0!";
		init(parameters, thetaOutputName, logfile, randomGenerator);	};

	~TGenotypePriorTheta(){
		out.close();
		delete thetaEstimator;
	};

	void update(const TWindow & window, TLog* logfile, const TGenotypeLikelihoodCalculator & glCalculator){
		logfile->startIndent("Estimating theta for prior:");
		//clear theta estimator
		(*thetaEstimator).clear();

		//adding sites to estimator
		thetaEstimator->add(window, glCalculator);

		//estimate Theta
		if(!thetaEstimator->estimateTheta()){
			if(hasDefaultTheta){
				logfile->conclude("Will use a default theta of " + toString(defaultTheta) + ".");
				thetaEstimator->setTheta(defaultTheta);
				GenotypeLikelihoods::TBaseData freq(0.25);
				thetaEstimator->setBaseFreq(freq);
			} else
				throw "Please increase window size or provide a default theta!";
		}

		//update prior
		thetaEstimator->fillPGenotype(genotypePrior);

		//write results to file
		window.writeCoordinates(out, chromosomes);
		out.write()
		out.write(window->_chrName, window->startPos, window->endPos);
		logfile->endIndent();
	};
};


}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEPRIOR_H_ */
