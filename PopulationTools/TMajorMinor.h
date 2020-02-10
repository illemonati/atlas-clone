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
#include "TGLF.h"
#include "../TGenotypeMap.h"
#include "../TQualityMap.h"
#include "../TRandomGenerator.h"
#include "TGenotypeFrequencies.h"

//-----------------------------------------------
//TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase{
protected:
	TGenotypeMap genoMap;
	TRandomGenerator* randomGenerator;

	TPopulationLikehoodLocus genotypeLikelihoods;
	TGenotypeFrequencies genotypeFrequencies;
	double* L10L_perCombination;
	std::vector<uint8_t> usedAllelicCombinations;

	void useAllAlleleicCombinations();
	void useAllelicCombinationsThatContain(const Base & base);
	void calculateL10LPerCombination();
	void chooseBestAllelicCombinationAmongThoseWithEqualScores();
	virtual void findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter);

	void _estimateMajorMinor(TGlfMultiReader & glfReader, TGlfConverter & glfConverter);

public:
	Base minor, major;
	int bestAllelicCombination;
	double L10L;
	int variantQuality; // in phred format for VCF

	TMajorMinorEstimatorBase(TRandomGenerator* RandomGenerator);
	virtual ~TMajorMinorEstimatorBase();

	void estimateMajorMinor(TGlfMultiReader & glfReader, TGlfConverter & glfConverter);
	void estimateMajorMinor(TGlfMultiReader & glfReader, TGlfConverter & glfConverter, const Base & base);
};

class TMajorMinorEstimatorSkotte:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	TGenotypeFrequencies priorGenotypeFrequencies;

	void findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter);

public:
	TMajorMinorEstimatorSkotte(TRandomGenerator* RandomGenerator, double EpsilonF);
	virtual ~TMajorMinorEstimatorSkotte();
};

class TMajorMinorEstimatorMLE:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	TGenotypeFrequencies* tmpGenotypeFrequencies;

	double estimateGenotypeFrequencies(TGlfMultiReader & glfReader, const int alleleicCombination, TGlfConverter & glfConverter);
	void findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter);

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
	TGenotypeMap genoMap;
	gz::ogzstream vcf;
	bool vcfOpened;
	TGlfConverter glfConverter;

	void openVCF(std::string filenameTag, TGlfMultiReader & glfReader, bool usePhredLikelihoods);
	void closeVCF();

public:
	TMajorMinor(TParameters & params, TLog* Logfile);
	~TMajorMinor();

	void estimateMajorMinor(TParameters & params);
};



#endif /* TMAJORMINOR_H_ */
