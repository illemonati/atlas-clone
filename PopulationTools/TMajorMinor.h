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
#include "TGenotypeMap.h"
#include "TQualityMap.h"
#include "TRandomGenerator.h"
#include "TGenotypeFrequencies.h"
#include "TGenotypeData.h"
#include <vector>

namespace PopulationTools{

using coretools::TLog;
using coretools::TParameters;
using coretools::TRandomGenerator;
using genometools::Base;
using genometools::Genotype;
using genometools::AllelicCombination;
using genometools::HighPrecisionPhredIntProbability;
using GLF::TMultiGLFData;


//-------------------------------------
// TAlleleicCombinationData
//-------------------------------------
class TAlleleicCombinationData : public GenotypeLikelihoods::TData_base<coretools::Log10Probability, AllelicCombination, genometools::AllelicCombinationEnum, genometools::alleleicCombinationNN>{
private:
	using TData_base<coretools::Log10Probability, AllelicCombination, genometools::AllelicCombinationEnum, genometools::alleleicCombinationNN>::_data;

public:
	TAlleleicCombinationData() : TData_base<coretools::Log10Probability, AllelicCombination, genometools::AllelicCombinationEnum, genometools::alleleicCombinationNN>(0.0) {};
	~TAlleleicCombinationData() = default;
};

//-----------------------------------------------
//TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase{
protected:
	TRandomGenerator* randomGenerator;

	GLF::TMultiGLFDataOneAllelicCombination genotypeLikelihoods;
	TGenotypeFrequencies genotypeFrequencies;
	TAlleleicCombinationData L10L_perCombination;

	std::vector<AllelicCombination> usedAllelicCombinations;

	void useAllAlleleicCombinations();
	void useAllelicCombinationsThatContain(const Base & base);
	void calculateL10LPerCombination();
	void chooseBestAllelicCombinationAmongThoseWithEqualScores();
	virtual void findMLAllelicCombination(const TMultiGLFData & data);

	void _estimateMajorMinor(const TMultiGLFData & data);

public:
	Base minor, major;
	genometools::AllelicCombination bestAllelicCombination;
	coretools::Log10Probability L10L;
	genometools::PhredIntProbability variantQuality; // in phred format for VCF

	TMajorMinorEstimatorBase(TRandomGenerator* RandomGenerator);
	virtual ~TMajorMinorEstimatorBase() = default;

	void estimateMajorMinor(const TMultiGLFData & data);
	void estimateMajorMinor(const TMultiGLFData & data, const Base & base);
};

class TMajorMinorEstimatorSkotte:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	TGenotypeFrequencies priorGenotypeFrequencies;

	void findMLAllelicCombination(const TMultiGLFData & data);

public:
	TMajorMinorEstimatorSkotte(TRandomGenerator* RandomGenerator, const double & EpsilonF);
	~TMajorMinorEstimatorSkotte() = default;
};

class TMajorMinorEstimatorMLE:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	std::array<TGenotypeFrequencies, 6> tmpGenotypeFrequencies;

	coretools::Log10Probability estimateGenotypeFrequencies(const TMultiGLFData & data, const AllelicCombination & alleleicCombination);
	void findMLAllelicCombination(const TMultiGLFData & data);

public:
	TMajorMinorEstimatorMLE(TRandomGenerator* RandomGenerator, const double & EpsilonF);
	~TMajorMinorEstimatorMLE() = default;
};

//-----------------------------------------------
//TMajorMinor
//-----------------------------------------------
class TMajorMinor{
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	bool hasReference;
	gz::ogzstream vcf;
	bool vcfOpened;

	//settings
	uint32_t minSamplesWithData;
	genometools::PhredIntProbability minVariantQuality;

public:
	TMajorMinor(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TMajorMinor(){};

	void estimateMajorMinor(TParameters & params);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_majorMinor:public coretools::TTask{
public:
	TTask_majorMinor(){
		_citations.insert("Skotte et al. (2012) Genetic Epidemiology");
		_explanation = "Estimating major and minor alles"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TMajorMinor majorMinor(Logfile, Parameters, _randomGenerator);
		majorMinor.estimateMajorMinor(Parameters);
	};
};

}; //end namespace

#endif /* TMAJORMINOR_H_ */
