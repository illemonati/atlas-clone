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

//-----------------------------------------------
//TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase{
protected:
	TQualityMap qualiMap;
	TGenotypeMap genoMap;

	int numSamples;
	double* genotypeLikelihoods;
	double* mleGenotypeFrequencies;

	double calculateLogLikelihood(double* genotypeFrequencies);
	void fillLikelihoods(uint8_t** phred, Genotype* genotypes);
	void guessGenotypeFrequencies(double* genotypeFrequencies);
	virtual int findMLAllelicCombination(uint8_t** phred);

public:
	TMajorMinorEstimatorBase(int NumSamples);
	virtual ~TMajorMinorEstimatorBase();

	std::pair<Base,Base>  estimateMajorMinor(uint8_t** phred);
};

class TMajorMinorEstimatorSkotte:public TMajorMinorEstimatorBase{
private:
	double* priorGenotypeFrequencies;

	int findMLAllelicCombination(uint8_t** phred);

public:
	TMajorMinorEstimatorSkotte(int NumSamples);
	virtual ~TMajorMinorEstimatorSkotte();
};

class TMajorMinorEstimatorMLE:public TMajorMinorEstimatorBase{
private:
	double maxF;
	double** genotypeFrequencies;

	void estimateGenotypeFrequenciesEM(double* genotypeFrequencies);
	double estimateGenotypeFrequencies(uint8_t** phred, const int alleleicCombination);
	int findMLAllelicCombination(uint8_t** phred);

public:
	TMajorMinorEstimatorMLE(int NumSamples, double MaxF);
	~TMajorMinorEstimatorMLE();

};

//-----------------------------------------------
//TMajorMinor
//-----------------------------------------------
class TMajorMinor{
private:
	TLog* logfile;
	TGenotypeMap genoMap;
	gz::ogzstream vcf;
	bool vcfOpened;

	void openVCF(std::string filenameTag, TGlfMultiReader & glfReader);
	void closeVCF();

public:
	TMajorMinor(TLog* Logfile);

	void estimateMajorMinor(TParameters & params);

};



#endif /* TMAJORMINOR_H_ */
