/*
 * TGenotypePrior.h
 *
 *  Created on: May 15, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TGENOTYPEPRIOR_H_
#define GENOTYPELIKELIHOODS_TGENOTYPEPRIOR_H_

#include "coretools/Main/TLog.h"

#include "TThetaEstimator.h"
#include "TWindow.h"

namespace GenotypeLikelihoods {

//---------------------------------------------------------------------------------
// TGenotypePrior
// An class used to serve genotype priors to Bayesian inferences, such as caller
//---------------------------------------------------------------------------------
class TGenotypePrior {
protected:
	TGenotypeProbabilities genotypePrior{};
public:
	virtual ~TGenotypePrior() = default;

	virtual void update(const TWindow &window, const TErrorModels &) = 0;
	TGenotypeProbabilities *getPointerToPrior() { return &genotypePrior; };
	coretools::Probability operator[](const genometools::Genotype &genotype) { return genotypePrior[genotype]; };
};

class TGenotypePriorUniform : public TGenotypePrior {
	virtual void update(const TWindow &, const TErrorModels &) override {};
};

class TGenotypePriorFixedTheta : public TGenotypePrior {
private:
	TThetaEstimator *thetaEstimator;
	bool equalBaseFreq;
public:
	TGenotypePriorFixedTheta(double theta, bool EqualBaseFreq) {
		thetaEstimator = new TThetaEstimator;
		thetaEstimator->setTheta(theta);
		equalBaseFreq = EqualBaseFreq;
		if (equalBaseFreq) {
			GenotypeLikelihoods::TBaseProbabilities freq{};
			thetaEstimator->setBaseFreq(freq);
		}
		genotypePrior = thetaEstimator->pGenotype();
	};

	~TGenotypePriorFixedTheta() { delete thetaEstimator; };

	void update(const TWindow &window, const TErrorModels &) override {
		using genometools::Base;
		using coretools::instances::logfile;
		if (!equalBaseFreq) {
			logfile().listFlush("Estimating base frequencies for prior ...");
			GenotypeLikelihoods::TBaseProbabilities freq = window.estimateBaseFrequencies();
			thetaEstimator->setBaseFreq(freq);
			logfile().done();
			logfile().conclude("Estimated base frequencies: " + coretools::str::toString(freq[Base::A]) + ", " +
					  coretools::str::toString(freq[Base::C]) + ", " + coretools::str::toString(freq[Base::G]) +
					  ", " + coretools::str::toString(freq[Base::T]));
			genotypePrior = thetaEstimator->pGenotype();
		}
	};
};

class TGenotypePriorTheta : public TGenotypePrior {
private:
	TThetaEstimator *thetaEstimator;

	TThetaOutputFile out;
	double defaultTheta;
	bool hasDefaultTheta;

	void init(std::string thetaOutputName) {
		thetaEstimator = new TThetaEstimator;
		out.open(thetaEstimator, thetaOutputName);
	};

public:
	TGenotypePriorTheta(const std::string &thetaOutputName) {
		hasDefaultTheta = false;
		defaultTheta    = -1.0;

		init(thetaOutputName);
	};

	TGenotypePriorTheta(std::string thetaOutputName, double DefaultTheta) {
		hasDefaultTheta = true;
		defaultTheta    = DefaultTheta;
		if (defaultTheta < 0.0) UERROR("Theta must be >= 0.0!");
		init(thetaOutputName);
	};

	~TGenotypePriorTheta() {
		out.close();
		delete thetaEstimator;
	};

	void update(const TWindow &window, const TErrorModels &glCalculator) override {
		using coretools::instances::logfile;
		logfile().startIndent("Estimating theta for prior:");
		// clear theta estimator
		(*thetaEstimator).clear();

		// adding sites to estimator
		thetaEstimator->add(window, glCalculator);

		// estimate Theta
		if (!thetaEstimator->estimateTheta()) {
			if (hasDefaultTheta) {
				logfile().conclude("Will use a default theta of " + coretools::str::toString(defaultTheta) + ".");
				thetaEstimator->setTheta(defaultTheta);
				GenotypeLikelihoods::TBaseProbabilities freq{};
				thetaEstimator->setBaseFreq(freq);
			} else
				UERROR("Please increase window size or provide a default theta!");
		}

		// update prior
		genotypePrior = thetaEstimator->pGenotype();
	};
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TGENOTYPEPRIOR_H_ */
