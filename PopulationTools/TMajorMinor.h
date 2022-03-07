/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include <math.h>
#include <TPopulationLikelihoodLocus.h>
#include "TGlfMultiReader.h"
#include "TRandomGenerator.h"
#include "TGenotypeFrequencies.h"
#include "TGenotypeData.h"
#include <vector>

namespace PopulationTools{

//-------------------------------------
// TAlleleicCombinationData
//-------------------------------------
class TAlleleicCombinationData
	: public GenotypeLikelihoods::TData_base<coretools::Log10Probability, genometools::AllelicCombination,
					     genometools::index(genometools::AllelicCombination::NN)> {
private:
	using TData_base<coretools::Log10Probability, genometools::AllelicCombination,
			 genometools::index(genometools::AllelicCombination::NN)>::_data;

public:
	TAlleleicCombinationData()
		: TData_base<coretools::Log10Probability, genometools::AllelicCombination,
			 genometools::index(genometools::AllelicCombination::NN)>(0.0){};
	~TAlleleicCombinationData() = default;
};

//-----------------------------------------------
//TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase{
protected:
	coretools::TRandomGenerator* randomGenerator;

	GLF::TMultiGLFDataOneAllelicCombination genotypeLikelihoods;
	genometools::TGenotypeFrequencies genotypeFrequencies;
	TAlleleicCombinationData L10L_perCombination;

	std::vector<genometools::AllelicCombination> usedAllelicCombinations;

	void useAllAlleleicCombinations();
	void useAllelicCombinationsThatContain(const genometools::Base & base);
	void calculateL10LPerCombination();
	void chooseBestAllelicCombinationAmongThoseWithEqualScores();
	virtual void findMLAllelicCombination(const GLF::TMultiGLFData & data);

	void _estimateMajorMinor(const GLF::TMultiGLFData & data);

public:
	genometools::Base minor, major;
	genometools::AllelicCombination bestAllelicCombination;
	coretools::Log10Probability L10L;
	genometools::PhredIntProbability variantQuality; // in phred format for VCF

	TMajorMinorEstimatorBase(coretools::TRandomGenerator* RandomGenerator);
	virtual ~TMajorMinorEstimatorBase() = default;

	void estimateMajorMinor(const GLF::TMultiGLFData & data);
	void estimateMajorMinor(const GLF::TMultiGLFData & data, const genometools::Base & base);
};

class TMajorMinorEstimatorSkotte:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	genometools::TGenotypeFrequencies priorGenotypeFrequencies;

	void findMLAllelicCombination(const GLF::TMultiGLFData & data);

public:
	TMajorMinorEstimatorSkotte(coretools::TRandomGenerator* RandomGenerator, double EpsilonF);
	~TMajorMinorEstimatorSkotte() = default;
};

class TMajorMinorEstimatorMLE:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	std::array<genometools::TGenotypeFrequencies, 6> tmpGenotypeFrequencies;

	coretools::Log10Probability estimateGenotypeFrequencies(const GLF::TMultiGLFData & data, const genometools::AllelicCombination & alleleicCombination);
	void findMLAllelicCombination(const GLF::TMultiGLFData & data);

public:
	TMajorMinorEstimatorMLE(coretools::TRandomGenerator* RandomGenerator, double EpsilonF);
	~TMajorMinorEstimatorMLE() = default;
};

//-----------------------------------------------
//TMajorMinor
//-----------------------------------------------
class TMajorMinor{
private:
	coretools::TLog* logfile;
	coretools::TRandomGenerator* randomGenerator;
	bool hasReference;

	//settings
	uint32_t minSamplesWithData;
	genometools::PhredIntProbability minVariantQuality;

public:
	TMajorMinor(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	~TMajorMinor(){};

	void estimateMajorMinor(coretools::TParameters & params);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_majorMinor:public coretools::TTask{
public:
	TTask_majorMinor(){
		_citations.insert("Skotte et al. (2012) Genetic Epidemiology");
		_explanation = "Estimating major and minor alles"; };

	void run(){
		using namespace coretools::instances;
		TMajorMinor majorMinor(parameters(), &logfile(), &randomGenerator());
		majorMinor.estimateMajorMinor(parameters());
	};
};

}; //end namespace

#endif /* TMAJORMINOR_H_ */
