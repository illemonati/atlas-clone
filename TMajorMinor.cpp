/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

//---------------------------------------------------
// TMajorMinorEstimatorBase
//---------------------------------------------------

TMajorMinorEstimatorBase::TMajorMinorEstimatorBase(int NumSamples){
	numSamples = NumSamples;
	genotypeLikelihoods = new double[3*numSamples];
	mleGenotypeFrequencies = new double[3];
	L10L_atMLE = 0.0;
};

TMajorMinorEstimatorBase::~TMajorMinorEstimatorBase(){
	delete[] genotypeLikelihoods;
	delete[] mleGenotypeFrequencies;
};

double TMajorMinorEstimatorBase::calculateLog10Likelihood(double* genotypeFrequencies){
	double LL = 0.0;
	for(int i = 0; i < 3 * numSamples; i += 3){
		LL += log(genotypeLikelihoods[i] * genotypeFrequencies[0]
			    + genotypeLikelihoods[i+1] * genotypeFrequencies[1]
		        + genotypeLikelihoods[i+2] * genotypeFrequencies[2]);
	}
	return LL;
};


void TMajorMinorEstimatorBase::fillLikelihoods(uint8_t** phred, Genotype* genotypes){
	for(int i=0; i<numSamples; ++i){
		int ind = 3*i;
		genotypeLikelihoods[ind]     = qualiMap[ phred[i][genotypes[0]] ];
		genotypeLikelihoods[ind + 1] = qualiMap[ phred[i][genotypes[1]] ];
		genotypeLikelihoods[ind + 2] = qualiMap[ phred[i][genotypes[2]] ];
	}
};

void TMajorMinorEstimatorBase::guessGenotypeFrequencies(double* genotypeFrequencies){
	//calculate by using MLE genotype for each individual
	genotypeFrequencies[0] = 0.0; genotypeFrequencies[1] = 0.0; genotypeFrequencies[2] = 0.0;
	for(int i = 0 ; i < 3 * numSamples; i += 3){
		if(genotypeLikelihoods[i + 1] > genotypeLikelihoods[i]){
			if(genotypeLikelihoods[i + 2] > genotypeLikelihoods[i + 1]) genotypeFrequencies[2] += 1.0;
			else genotypeFrequencies[1] += 1.0;
		} else {
			if(genotypeLikelihoods[i + 2] > genotypeLikelihoods[i]) genotypeFrequencies[2] += 1.0;
			else genotypeFrequencies[0] += 1.0;
		}
	}

	double sum = 0.0;
	for(int g = 0; g<3; ++g){
		genotypeFrequencies[g] /= (double) numSamples;
		if(genotypeFrequencies[g] <= 0.0) genotypeFrequencies[g] = 0.01;
		if(genotypeFrequencies[g] >= 1.0) genotypeFrequencies[g] = 0.99;
		sum += genotypeFrequencies[g];
	}
	for(int g = 0; g<3; ++g){
		genotypeFrequencies[g] /= sum;
	}
};

int TMajorMinorEstimatorBase::findMLAllelicCombination(uint8_t** phred){
	throw "Function TMajorMinorEstimatorBase::findMLAllelicCombination(uint8_t** phred) not implemented for base class!";
};

std::pair<Base,Base>  TMajorMinorEstimatorBase::estimateMajorMinor(uint8_t** phred){
	int mleCombination = findMLAllelicCombination(phred);

	//which one is major?
	if(mleGenotypeFrequencies[0] > mleGenotypeFrequencies[2])
		return std::pair<Base,Base>(genoMap.alleleicCombinationToBase[mleCombination][0], genoMap.alleleicCombinationToBase[mleCombination][1]);
	else
		return std::pair<Base,Base>(genoMap.alleleicCombinationToBase[mleCombination][1], genoMap.alleleicCombinationToBase[mleCombination][0]);
};

//---------------------------------------------------
// TMajorMinorEstimatorSkotte
//---------------------------------------------------
TMajorMinorEstimatorSkotte::TMajorMinorEstimatorSkotte(int NumSamples):TMajorMinorEstimatorBase(NumSamples){
	genotypeLikelihoods = new double[3*numSamples];

	priorGenotypeFrequencies = new double[3];
	priorGenotypeFrequencies[0] = 0.25;
	priorGenotypeFrequencies[1] = 0.50;
	priorGenotypeFrequencies[2] = 0.25;
};

TMajorMinorEstimatorSkotte::~TMajorMinorEstimatorSkotte(){
	delete[] priorGenotypeFrequencies;
};


int TMajorMinorEstimatorSkotte::findMLAllelicCombination(uint8_t** phred){
	int mleCombination = 0;
	fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[0]);
	L10L_atMLE = calculateLog10Likelihood(priorGenotypeFrequencies);
	for(int i=1; i<6; ++i){
		fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[i]);
		double LL = calculateLog10Likelihood(priorGenotypeFrequencies);
		if(LL > L10L_atMLE){
			L10L_atMLE = LL;
			mleCombination = i;
		}
	}

	//now guess genotype frequencies at MLE
	fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[mleCombination]);
	guessGenotypeFrequencies(mleGenotypeFrequencies);

	return mleCombination;
}

//---------------------------------------------------
// TMajorMinorEstimatorMLE
//---------------------------------------------------
TMajorMinorEstimatorMLE::TMajorMinorEstimatorMLE(int NumSamples, double MaxF):TMajorMinorEstimatorBase(NumSamples){
	maxF = MaxF;

	genotypeFrequencies = new double*[6];
	for(int i=0; i<6; ++i) genotypeFrequencies[i] = new double[3];
}

TMajorMinorEstimatorMLE::~TMajorMinorEstimatorMLE(){
	for(int i=0; i<6; ++i) delete[] genotypeFrequencies[i];
	delete[] genotypeFrequencies;
};

void TMajorMinorEstimatorMLE::estimateGenotypeFrequenciesEM(double* genotypeFrequencies){
	//prepare variables
	double weightsNull[3];
	double genoFreq_old[3];

	//estimate initial frequencies from MLEs
	guessGenotypeFrequencies(genotypeFrequencies);

	//run EM for max 1000 steps
	for(int s=0; s<1000; ++s){
		//set genofreq and calc P(g|f)
		for(int g=0; g<3; ++g){
			genoFreq_old[g] = genotypeFrequencies[g];
			genotypeFrequencies[g] = 0.0;
		}

		//estimate new genotype frequencies
		for(int i = 0; i < 3 * numSamples; i += 3){
			double sum = 0.0;
			for(int g = 0; g < 3; ++g){
				weightsNull[g] = genotypeLikelihoods[i + g] * genoFreq_old[g];
				sum += weightsNull[g];
			}
			genotypeFrequencies[0] += weightsNull[0] / sum;
			genotypeFrequencies[2] += weightsNull[2] / sum;
		}

		genotypeFrequencies[0] /= (double) numSamples;
		genotypeFrequencies[2] /= (double) numSamples;
		genotypeFrequencies[1] = 1.0 - genotypeFrequencies[0] - genotypeFrequencies[2];

		//check if we stop
		if(fabs(genotypeFrequencies[0] - genoFreq_old[0]) < maxF && fabs(genotypeFrequencies[2] - genoFreq_old[2]) < maxF) break;
	}
};

double TMajorMinorEstimatorMLE::estimateGenotypeFrequencies(uint8_t** phred, const int alleleicCombination){
	fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[alleleicCombination]);
	estimateGenotypeFrequenciesEM(genotypeFrequencies[alleleicCombination]);
	return calculateLog10Likelihood(genotypeFrequencies[alleleicCombination]);
}

int TMajorMinorEstimatorMLE::findMLAllelicCombination(uint8_t** phred){
	int mleCombination = 0;
	L10L_atMLE = estimateGenotypeFrequencies(phred, 0);
	for(int i=1; i<6; ++i){
		double LL = estimateGenotypeFrequencies(phred, i);
		if(LL > L10L_atMLE){
			L10L_atMLE = LL;
			mleCombination = i;
		}
	}
	mleGenotypeFrequencies = genotypeFrequencies[mleCombination];
	return mleCombination;
}

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
TMajorMinor::TMajorMinor(TLog* Logfile, TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	randomGenerator = RandomGenerator;
	vcfOpened = false;
};

void TMajorMinor::openVCF(std::string filenameTag, TGlfMultiReader & glfReader){
	if(vcfOpened) closeVCF();

	//open vcf file
	filenameTag += "vcf.gz";
	vcf.open(filenameTag.c_str());
	vcfOpened = true;

	//write info
	//TODO: create VCF class to harmonize code across different uses. Also include code in Tiger and other
	vcf << "##fileformat=VCFv4.2\n";
	vcf << "##source=ATLAS_GLF_Caller\n";
	glfReader.writeVCFHeader(vcf);
};

void TMajorMinor::closeVCF(){
	vcf.close();
	vcfOpened = false;
};


void TMajorMinor::estimateMajorMinor(TParameters & params){
	//open GLF files
	TGlfMultiReader glfReader(params, logfile);
	glfReader.setAllActive();

	//estimation method
	std::string method = params.getParameterStringWithDefault("method", "MLE");
	TMajorMinorEstimatorBase* MMEstimator;
	if(method == "Skotte"){
		logfile->list("Will estimate major / minor alleles using the method of Skotte et al. (2012).");
		MMEstimator = new TMajorMinorEstimatorSkotte(glfReader.numSamples());
	} else if(method == "MLE"){
		double maxF = params.getParameterDoubleWithDefault("maxF", 0.00001);
		logfile->list("Will estimate major / minor alleles using the MLE method with maxF " + toString(maxF) + ".");
		MMEstimator = new TMajorMinorEstimatorMLE(glfReader.numSamples(), maxF);
	} else throw "Unknown MajorMinor method '" + method + "'!";

	//TODO: think about filters

	//variables
	double* genotypeLikelihoods = new double[3*glfReader.numSamples()];
	double** genotypeFrequencies = new double*[6];
	for(int i=0; i<6; ++i) genotypeFrequencies[i] = new double[3];

	while(glfReader.readNext()){
		//TODO: think about filters

		//do major minor
		std::pair<Base,Base> majorMinor = MMEstimator->estimateMajorMinor(glfReader.data);

		//write to VCF
		glfReader.writeSiteToVCF(vcf, MMEstimator->L10L_atMLE, genoMap.genotypeMap[majorMinor.first][majorMinor.first], genoMap.genotypeMap[majorMinor.first][majorMinor.second], genoMap.genotypeMap[majorMinor.second][majorMinor.second], randomGenerator);



	}


	//clean storage
	delete[] genotypeLikelihoods;
	for(int i=0; i<6; ++i) delete[] genotypeFrequencies[i];
	delete[] genotypeFrequencies;
};



//---------------------------------------------------
// Skotte
//---------------------------------------------------

