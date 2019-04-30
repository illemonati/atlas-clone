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
#include "TRandomGenerator.h"
#include "TPopulationLikelihoodStorage.h"
#include "TGenotypeFrequencies.h"

//-----------------------------------------------
//TMajorMinorEstimator
//-----------------------------------------------
class TMajorMinorEstimatorBase{
protected:
	TQualityMap qualiMap;
	TGenotypeMap genoMap;
	TRandomGenerator* randomGenerator;

	TPopulationLikehoodStorage genotypeLikelihoods;
	TGenotypeFrequencies genotypeFrequencies;
	double* L10L_perCombination;

	void calculateL10LPerCombination();
	void chooseMLAllelicCombinationAmongThoseWithEqualLikelihood();
	virtual void findMLAllelicCombination(TGlfMultiReader & glfReader);

public:
	Base minor, major;
	int allelicCombination;
	double L10L;
	int variantQuality;

	TMajorMinorEstimatorBase(TRandomGenerator* RandomGenerator);
	virtual ~TMajorMinorEstimatorBase();

	void estimateMajorMinor(TGlfMultiReader & glfReader);
};

class TMajorMinorEstimatorSkotte:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	TGenotypeFrequencies priorGenotypeFrequencies;

	void findMLAllelicCombination(TGlfMultiReader & glfReader);

public:
	TMajorMinorEstimatorSkotte(TRandomGenerator* RandomGenerator, double EpsilonF);
	virtual ~TMajorMinorEstimatorSkotte();
};

class TMajorMinorEstimatorMLE:public TMajorMinorEstimatorBase{
private:
	double epsilonF;
	TGenotypeFrequencies* tmpGenotypeFrequencies;

	double estimateGenotypeFrequencies(TGlfMultiReader & glfReader, const int alleleicCombination);
	void findMLAllelicCombination(TGlfMultiReader & glfReader);

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
	TGenotypeMap genoMap;
	gz::ogzstream vcf;
	bool vcfOpened;

	void openVCF(std::string filenameTag, TGlfMultiReader & glfReader, bool usePhredLikelihoods);
	void closeVCF();

public:
	TMajorMinor(TParameters & params, TLog* Logfile);
	~TMajorMinor();

	void estimateMajorMinor(TParameters & params);

};



#endif /* TMAJORMINOR_H_ */
