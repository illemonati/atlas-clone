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

namespace GenotypeLikelihoods{


//--------------------------------------------------------------------
// TBaseLikelihoods
//--------------------------------------------------------------------
class TBaseLikelihoods{
private:
	double likelihoods[5];

public:
	TBaseLikelihoods();
	void operator=(const TBaseLikelihoods & other);
	double& operator[](const Base base){ return likelihoods[base];};
	double at(const Base base) const { return likelihoods[base]; };
	void reset();
};

//--------------------------------------------------------------------
// TGenotypeData
// base class for likelihoods, prior and posterior
//--------------------------------------------------------------------
class TGenotypeData{
protected:
	double data[10];

	void _copyFrom(const TGenotypeData & other);
	void _set(const double val);
	void _normalize();

public:
	TGenotypeData(){};
	virtual ~TGenotypeData(){};

	void operator=(const TGenotypeData & other);
	double operator[](const Genotype genotype){ return data[genotype]; };
	double operator[](const uint8_t genotype){ return data[genotype]; };
	double at(const Genotype genotype) const { return data[genotype]; };
	double at(const uint8_t genotype) const { return data[genotype]; };

	virtual void reset();
	void write(TOutputFileZipped & out) const;
};


//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
class TGenotypeLikelihoods:public TGenotypeData{
public:
	TGenotypeLikelihoods();

	void operator=(const TGenotypeLikelihoods & other);

	void fill(const std::vector<TBaseLikelihoods> & bases);
	void fill(const std::vector<TBaseLikelihoods> & bases, const size_t size);
};

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

//--------------------------------------------------------------------
// TGenotypePosteriorProbabilities
//--------------------------------------------------------------------

class TGenotypePosteriorProbabilities:public TGenotypeData{
public:
	TGenotypePosteriorProbabilities();
	void reset();
	void fill(TGenotypeLikelihoods & likelihoods, TGenotypePrior & prior);
};

}; //end namespace


#endif /* TGENOTYPEPRIOR_H_ */
