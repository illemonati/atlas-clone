/*
 * TGenotypePrior.h
 *
 *  Created on: Nov 20, 2018
 *      Author: phaentu
 */

#ifndef TGENOTYPEPRIOR_H_
#define TGENOTYPEPRIOR_H_

#include "../TThetaEstimator.h"
#include "../TWindow.h"

//---------------------------------------------------------------------------------
// TGenotypePrior
// An class used to serve genotype priors to Bayesian inferences, such as caller
//---------------------------------------------------------------------------------
class TGenotypePrior{
protected:
	double genotypePrior[10];

public:
	TGenotypePrior(){
		for(int g=0; g<10; ++g)
			genotypePrior[g] = 1.0;
	};

	virtual ~TGenotypePrior(){};

	virtual void update(TWindow* window, const std::string chrName, TLog* logfile){};
	double* getPointerToPrior(){ return genotypePrior; };
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
			TBaseFrequencies freq;
			freq.setEqualBaseFreq();
			thetaEstimator->setBaseFreq(freq);
		}
		thetaEstimator->fillPGenotype(genotypePrior);
	};

	~TGenotypePriorFixedTheta(){
		delete thetaEstimator;
	};

	void update(TWindow* window, const std::string chrName, TLog* logfile){
		if(!equalBaseFreq){
			logfile->listFlush("Estimating base frequencies for prior ...");
			window->estimateBaseFrequencies();
			thetaEstimator->setBaseFreq(window->baseFreq);
			logfile->done();
			logfile->conclude("Estimated base frequencies: " + toString(window->baseFreq.freq[0])+ ", " + toString(window->baseFreq.freq[1]) + ", " + toString(window->baseFreq.freq[2]) + ", " + toString(window->baseFreq.freq[3]));
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



	void update(TWindow* window, const std::string chrName, TLog* logfile){
		logfile->startIndent("Estimating theta and base frequencies:");
		//clear theta estimator
		(*thetaEstimator).clear();

		//adding sites to estimator
		window->addSitesToThetaEstimator(thetaEstimator->pointerToDataContainer());

		//estimate Theta
		if(!thetaEstimator->estimateTheta()){
			if(hasDefaultTheta){
				logfile->conclude("Will use a default theta of " + toString(defaultTheta) + ".");
				thetaEstimator->setTheta(defaultTheta);
				TBaseFrequencies freq;
				freq.setEqualBaseFreq();
				thetaEstimator->setBaseFreq(freq);
			} else
				throw "Please increase window size or provide a default theta!";
		}

		//write results to file
		out.write(window->chrName, window->start, window->end);
		logfile->endIndent();
	};
};




#endif /* TGENOTYPEPRIOR_H_ */
