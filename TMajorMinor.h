/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include <math.h>
#include "TGLF.h"
#include "TGenotypeMap.h"
#include "TQualityMap.h"


//-----------------------------------------------
//TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase{
protected:
	TQualityMap qualiMap;
	TGenotypeMap genoMap;
	TRandomGenerator* randomGenerator;

	int numSamples;
	double* genotypeLikelihoods;
	double* L10L_perCombination;

	double calculateLog10Likelihood(double* genotypeFrequencies);
	void fillLikelihoods(uint8_t** phred, Genotype* genotypes);
	void guessGenotypeFrequencies(double* genotypeFrequencies);
	void normalizeGenotypeFrequencies();
	void calculateL10LPerCombination();
	void chooseMLAllelicCombinationAmongThoseWithEqualLikelihood();
	virtual void findMLAllelicCombination(uint8_t** phred);

public:
	bool* isHaploid;
	Base minor, major;
	int allelicCombination;
	double L10L;
	int variantQuality;
	double* genotypeFrequencies;


	TMajorMinorEstimatorBase(int NumSamples, TRandomGenerator* RandomGenerator);
	virtual ~TMajorMinorEstimatorBase();

	void estimateMajorMinor(uint8_t** phred);
};

class TMajorMinorEstimatorSkotte:public TMajorMinorEstimatorBase{
private:
	double* priorGenotypeFrequencies;

	void findMLAllelicCombination(uint8_t** phred);

public:
	TMajorMinorEstimatorSkotte(int NumSamples, TRandomGenerator* RandomGenerator);
	virtual ~TMajorMinorEstimatorSkotte();
};

class TMajorMinorEstimatorMLE:public TMajorMinorEstimatorBase{
private:
	double maxF;
	double** tmpGenotypeFrequencies;

	void estimateGenotypeFrequenciesEM(double* genotypeFrequencies);
	double estimateGenotypeFrequencies(uint8_t** phred, const int alleleicCombination);
	void findMLAllelicCombination(uint8_t** phred);

public:
	TMajorMinorEstimatorMLE(int NumSamples, TRandomGenerator* RandomGenerator, double MaxF);
	~TMajorMinorEstimatorMLE();

};

//-----------------------------------------------
//TMajorMinor
//-----------------------------------------------
class TMajorMinor{
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	TGenotypeMap genoMap;
	gz::ogzstream vcf;
	bool vcfOpened;

	void openVCF(std::string filenameTag, TGlfMultiReader & glfReader, bool usePhredLikelihoods);
	void closeVCF();

public:
	TMajorMinor(TParameters & params, TLog* Logfile);

	void estimateMajorMinor(TParameters & params);

};



#endif /* TMAJORMINOR_H_ */
