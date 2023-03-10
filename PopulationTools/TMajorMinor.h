/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include <array>
#include <set>
#include <string>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/TGenotypeFrequencies.h"
#include "TGlfMultiReader.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"

namespace PopulationTools {

using TAlleleicCombinationData =
    coretools::TStrongArray<coretools::Log10Probability, genometools::AllelicCombination>;


//-----------------------------------------------
// TMajorMinorEstimator
//-----------------------------------------------

class TMajorMinorEstimatorBase {
protected:
	GLF::TMultiGLFDataOneAllelicCombination _genotypeLikelihoods;
	genometools::TGenotypeFrequencies _genotypeFrequencies;
	TAlleleicCombinationData _L10L_perCombination;
	genometools::AllelicCombination _bestAllelicCombination = genometools::AllelicCombination::NN;
	genometools::Base _minor                                = genometools::Base::N;
	genometools::Base _major                                = genometools::Base::N;
	genometools::PhredIntProbability _variantQuality{0}; // in phred format for VCF

	virtual void _findMLAllelicCombination(const GLF::TMultiGLFData &data, genometools::Base base) = 0;
public:
	virtual ~TMajorMinorEstimatorBase() = default;

	void estimateMajorMinor(const GLF::TMultiGLFData &data, genometools::Base base = genometools::Base::N);
	constexpr genometools::Base minor() const noexcept {return _minor;}
	constexpr genometools::Base major() const noexcept {return _major;}
	constexpr genometools::PhredIntProbability variantQuality() const noexcept {return _variantQuality;} ;
	const genometools::TGenotypeFrequencies& genotypeFrequencies() const noexcept {return _genotypeFrequencies;}
};

class TMajorMinorEstimatorSkotte : public TMajorMinorEstimatorBase {
private:
	double _epsilonF;
	genometools::TGenotypeFrequencies _priorGenotypeFrequencies;

	void _findMLAllelicCombination(const GLF::TMultiGLFData &data, genometools::Base base) override;
public:
	TMajorMinorEstimatorSkotte(double EpsilonF);
};

class TMajorMinorEstimatorMLE : public TMajorMinorEstimatorBase {
private:
	double _epsilonF;
	std::array<genometools::TGenotypeFrequencies, 6> _tmpGenotypeFrequencies;

	coretools::Log10Probability _estimateGenotypeFrequencies(const GLF::TMultiGLFData &data,
								genometools::AllelicCombination alleleicCombination);
	void _findMLAllelicCombination(const GLF::TMultiGLFData &data, genometools::Base base) override;
public:
	TMajorMinorEstimatorMLE(double EpsilonF) : _epsilonF(EpsilonF) {};
};

struct TMajorMinor {
	void run();
};

} // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
