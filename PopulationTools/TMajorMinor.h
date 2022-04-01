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

namespace PopulationTools {

using TAlleleicCombinationData =
    coretools::StrongTypes::TStrongArray<coretools::Log10Probability, genometools::AllelicCombination>;


//-----------------------------------------------
// TMajorMinorEstimator
//-----------------------------------------------

class TMajorMinorEstimatorBase {
protected:
	using TAlleleicCombinations = std::array<genometools::AllelicCombination, genometools::index(genometools::AllelicCombination::max)>;
	GLF::TMultiGLFDataOneAllelicCombination _genotypeLikelihoods;
	genometools::TGenotypeFrequencies _genotypeFrequencies;
	TAlleleicCombinationData _L10L_perCombination;

	//void useAllAlleleicCombinations();
	//void useAllelicCombinationsThatContain(genometools::Base base);
	void _chooseBestAllelicCombinationAmongThoseWithEqualScores();
	virtual void _findMLAllelicCombination(const GLF::TMultiGLFData &data, const TAlleleicCombinations& used) = 0;
	void _estimateMajorMinor(const GLF::TMultiGLFData &data);
public:
	genometools::Base minor                                = genometools::Base::N;
	genometools::Base major                                = genometools::Base::N;
	genometools::AllelicCombination bestAllelicCombination = genometools::AllelicCombination::NN;
	coretools::Log10Probability L10L                       = 0.0;
	genometools::PhredIntProbability variantQuality{0}; // in phred format for VCF

	virtual ~TMajorMinorEstimatorBase() = default;

	void estimateMajorMinor(const GLF::TMultiGLFData &data);
	void estimateMajorMinor(const GLF::TMultiGLFData &data, genometools::Base base);
};

class TMajorMinorEstimatorSkotte : public TMajorMinorEstimatorBase {
private:
	double epsilonF;
	genometools::TGenotypeFrequencies priorGenotypeFrequencies;

	void _findMLAllelicCombination(const GLF::TMultiGLFData &data, const TAlleleicCombinations& used) override;
public:
	TMajorMinorEstimatorSkotte(double EpsilonF);
};

class TMajorMinorEstimatorMLE : public TMajorMinorEstimatorBase {
private:
	double epsilonF;
	std::array<genometools::TGenotypeFrequencies, 6> tmpGenotypeFrequencies;

	coretools::Log10Probability estimateGenotypeFrequencies(const GLF::TMultiGLFData &data,
								genometools::AllelicCombination alleleicCombination);
	void _findMLAllelicCombination(const GLF::TMultiGLFData &data, const TAlleleicCombinations& used) override;
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
