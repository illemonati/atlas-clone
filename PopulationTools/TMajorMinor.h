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

//-----------------------------------------------
//TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase{
protected:
	TGenotypeMap genoMap; //TODO: pass?
	TRandomGenerator* randomGenerator;

	TPopulationLikehoodLocus genotypeLikelihoods;
	TGenotypeFrequencies genotypeFrequencies;
	double* L10L_perCombination;
	std::vector<uint8_t> usedAllelicCombinations;

	void useAllAlleleicCombinations();
	void useAllelicCombinationsThatContain(const Base & base);
	void calculateL10LPerCombination();
	void chooseBestAllelicCombinationAmongThoseWithEqualScores();
	virtual void findMLAllelicCombination(TMultiGLFData & data, TGlfConverter & glfConverter);

	void _estimateMajorMinor(TMultiGLFData & data, TGlfConverter & glfConverter);

public:
	Base minor, major;
	int bestAllelicCombination;
	double L10L;
	int variantQuality; // in phred format for VCF

	TMajorMinorEstimatorBase(TRandomGenerator* RandomGenerator);
	virtual ~TMajorMinorEstimatorBase();

	void estimateMajorMinor(TMultiGLFData & data, TGlfConverter & glfConverter);
	void estimateMajorMinor(TMultiGLFData & data, TGlfConverter & glfConverter, const Base & base);
};

class TMajorMinorEstimatorSkotte:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	TGenotypeFrequencies priorGenotypeFrequencies;

	void findMLAllelicCombination(TMultiGLFData & data, TGlfConverter & glfConverter);

public:
	TMajorMinorEstimatorSkotte(TRandomGenerator* RandomGenerator, double EpsilonF);
	virtual ~TMajorMinorEstimatorSkotte();
};

class TMajorMinorEstimatorMLE:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	TGenotypeFrequencies* tmpGenotypeFrequencies;

	double estimateGenotypeFrequencies(TMultiGLFData & data, const int alleleicCombination, TGlfConverter & glfConverter);
	void findMLAllelicCombination(TMultiGLFData & data, TGlfConverter & glfConverter);

public:
	TMajorMinorEstimatorMLE(TRandomGenerator* RandomGenerator, double EpsilonF);
	~TMajorMinorEstimatorMLE();
};

//-----------------------------------------------
//TMajorMinor
//-----------------------------------------------
class TMajorMinor{
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	BamTools::Fasta reference;
	bool hasReference;
	gz::ogzstream vcf;
	bool vcfOpened;

public:
	TMajorMinor(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TMajorMinor(){};

	void estimateMajorMinor(TParameters & params);
};



#endif /* TMAJORMINOR_H_ */
