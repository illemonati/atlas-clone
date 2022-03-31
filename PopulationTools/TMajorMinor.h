/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include "TGenotypeData.h"
#include "TGenotypeFrequencies.h"
#include "TGlfMultiReader.h"
#include "TStrongArray.h"
#include <TPopulationLikelihoodLocus.h>
#include <math.h>
#include <vector>

namespace PopulationTools {

using TAlleleicCombinationData =
    coretools::StrongTypes::TStrongArray<coretools::Log10Probability, genometools::AllelicCombination>;

//-----------------------------------------------
// TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase {
protected:
	GLF::TMultiGLFDataOneAllelicCombination genotypeLikelihoods;
	genometools::TGenotypeFrequencies genotypeFrequencies;
	TAlleleicCombinationData L10L_perCombination;
	std::vector<genometools::AllelicCombination> usedAllelicCombinations;

	void useAllAlleleicCombinations();
	void useAllelicCombinationsThatContain(const genometools::Base &base);
	void calculateL10LPerCombination();
	void chooseBestAllelicCombinationAmongThoseWithEqualScores();
	virtual void findMLAllelicCombination(const GLF::TMultiGLFData &data) = 0;
	void _estimateMajorMinor(const GLF::TMultiGLFData &data);
public:
	genometools::Base minor                                = genometools::Base::N;
	genometools::Base major                                = genometools::Base::N;
	genometools::AllelicCombination bestAllelicCombination = genometools::AllelicCombination::NN;
	coretools::Log10Probability L10L                       = 0.0;
	genometools::PhredIntProbability variantQuality{0}; // in phred format for VCF

	TMajorMinorEstimatorBase() = default;
	virtual ~TMajorMinorEstimatorBase() = default;

	void estimateMajorMinor(const GLF::TMultiGLFData &data);
	void estimateMajorMinor(const GLF::TMultiGLFData &data, const genometools::Base &base);
};

class TMajorMinorEstimatorSkotte : public TMajorMinorEstimatorBase {
private:
	double epsilonF;
	genometools::TGenotypeFrequencies priorGenotypeFrequencies;

	void findMLAllelicCombination(const GLF::TMultiGLFData &data);
public:
	TMajorMinorEstimatorSkotte(double EpsilonF);
};

class TMajorMinorEstimatorMLE : public TMajorMinorEstimatorBase {
private:
	double epsilonF;
	std::array<genometools::TGenotypeFrequencies, 6> tmpGenotypeFrequencies;

	coretools::Log10Probability estimateGenotypeFrequencies(const GLF::TMultiGLFData &data,
								const genometools::AllelicCombination &alleleicCombination);
	void findMLAllelicCombination(const GLF::TMultiGLFData &data);
public:
	TMajorMinorEstimatorMLE(double EpsilonF) : epsilonF(EpsilonF) {};
};

void estimateMajorMinor();

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_majorMinor : public coretools::TTask {
public:
	TTask_majorMinor() {
		_citations.insert("Skotte et al. (2012) Genetic Epidemiology");
		_explanation = "Estimating major and minor alles";
	};

	void run() {
		estimateMajorMinor();
	};
};

}; // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
