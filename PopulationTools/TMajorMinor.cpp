/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

namespace PopulationTools{

//---------------------------------------------------
// TMajorMinorEstimatorBase
//---------------------------------------------------
TMajorMinorEstimatorBase::TMajorMinorEstimatorBase(TRandomGenerator* RandomGenerator){
	L10L_perCombination = new double[6];
	bestAllelicCombination = 0;
	minor = A;
	major = A;
	L10L = 0.0;
	randomGenerator = RandomGenerator;
	variantQuality = 0;
};

TMajorMinorEstimatorBase::~TMajorMinorEstimatorBase(){
	delete[] L10L_perCombination;
};

void TMajorMinorEstimatorBase::useAllAlleleicCombinations(){
	usedAllelicCombinations = {0, 1, 2, 3, 4, 5};
};

void TMajorMinorEstimatorBase::useAllelicCombinationsThatContain(const Base & base){
	usedAllelicCombinations.clear();

	//get the three alleleic combinations that contain base base
	usedAllelicCombinations.push_back(genoMap.allelicCombinationsMatchingBase[base][0]);
	usedAllelicCombinations.push_back(genoMap.allelicCombinationsMatchingBase[base][1]);
	usedAllelicCombinations.push_back(genoMap.allelicCombinationsMatchingBase[base][2]);
};

void TMajorMinorEstimatorBase::chooseBestAllelicCombinationAmongThoseWithEqualScores(){
	//identify max L10L
	L10L = L10L_perCombination[usedAllelicCombinations[0]];
	for(uint16_t ac : usedAllelicCombinations){
		if(L10L_perCombination[ac] > L10L){
			L10L = L10L_perCombination[ac];
		}
	}

	//select MLE combination
	std::vector<int> best_combinations;
	for(uint16_t i : usedAllelicCombinations){
		if(L10L_perCombination[i] == L10L)
			best_combinations.push_back(i);
	}
	bestAllelicCombination = best_combinations[randomGenerator->sample(best_combinations.size())];
};

void TMajorMinorEstimatorBase::findMLAllelicCombination(TMultiGLFData & data, TGlfConverter & glfConverter){
	throw "Function TMajorMinorEstimatorBase::findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter) not implemented for base class!";
};

void  TMajorMinorEstimatorBase::_estimateMajorMinor(TMultiGLFData & data, TGlfConverter & glfConverter){
	findMLAllelicCombination(data, glfConverter);

	//which one is major?
	if(genotypeFrequencies.alleleFrequency < 0.5){
		major = genoMap.alleleicCombinationToBase[bestAllelicCombination][0];
		minor = genoMap.alleleicCombinationToBase[bestAllelicCombination][1];
	} else {
		major = genoMap.alleleicCombinationToBase[bestAllelicCombination][1];
		minor = genoMap.alleleicCombinationToBase[bestAllelicCombination][0];

		//also flip frequencies
		genotypeFrequencies.flip();
	}

	//calculate variant quality
	int refHomIndex = genoMap.genotypeMap[major][major];
	double LL_fixed_glfPhred = 0.0;
	for(uint32_t i=0; i<data.size; ++i){
		if(data.samples[i].hasData){
			if(data.samples[i].isHaploid)
				LL_fixed_glfPhred += data.samples[i].genotypeLikelihoodsGLF[major];
			else
				LL_fixed_glfPhred += data.samples[i].genotypeLikelihoodsGLF[refHomIndex];
		}
	}

	variantQuality = glfConverter.toPhred(LL_fixed_glfPhred - glfConverter.log10ToGlfFormat(L10L));
};

void  TMajorMinorEstimatorBase::estimateMajorMinor(TMultiGLFData & data, TGlfConverter & glfConverter){
	//use all allelic combinations
	useAllAlleleicCombinations();

	//now estimate
	_estimateMajorMinor(data, glfConverter);
};

void  TMajorMinorEstimatorBase::estimateMajorMinor(TMultiGLFData & data, TGlfConverter & glfConverter, const Base & base){
	//use all allelic combinations
	useAllelicCombinationsThatContain(base);

	//now estimate
	_estimateMajorMinor(data, glfConverter);
};

//---------------------------------------------------
// TMajorMinorEstimatorSkotte
//---------------------------------------------------
TMajorMinorEstimatorSkotte::TMajorMinorEstimatorSkotte(TRandomGenerator* RandomGenerator, double EpsilonF):TMajorMinorEstimatorBase(RandomGenerator){
	 epsilonF = EpsilonF;

	//diploid
	priorGenotypeFrequencies.diploidFrequencies[0] = 0.25;
	priorGenotypeFrequencies.diploidFrequencies[1] = 0.50;
	priorGenotypeFrequencies.diploidFrequencies[2] = 0.25;

	//haploid
	priorGenotypeFrequencies.haploidFrequencies[0] = 0.50;
	priorGenotypeFrequencies.haploidFrequencies[1] = 0.50;
};

TMajorMinorEstimatorSkotte::~TMajorMinorEstimatorSkotte(){};

void TMajorMinorEstimatorSkotte::findMLAllelicCombination(TMultiGLFData & data, TGlfConverter & glfConverter){
	//calculate L10L for each allelic combination used
	for(uint16_t ac : usedAllelicCombinations){
		data.fill(genotypeLikelihoods, ac);
		L10L_perCombination[ac] = priorGenotypeFrequencies.calculateLog10Likelihood(genotypeLikelihoods, glfConverter);
		L10L = L10L_perCombination[ac];
	}

	//pick combination with highest likelihood
	chooseBestAllelicCombinationAmongThoseWithEqualScores();

	//now estimate genotype frequencies at MLE allelic combination
	data.fill(genotypeLikelihoods, bestAllelicCombination);
	genotypeFrequencies.estimate(genotypeLikelihoods, glfConverter, epsilonF);

	//calculate likelihood again with better genotype frequencies
	L10L_perCombination[bestAllelicCombination] = genotypeFrequencies.calculateLog10Likelihood(genotypeLikelihoods, glfConverter);
	L10L = L10L_perCombination[bestAllelicCombination];
};

//---------------------------------------------------
// TMajorMinorEstimatorMLE
//---------------------------------------------------
TMajorMinorEstimatorMLE::TMajorMinorEstimatorMLE(TRandomGenerator* RandomGenerator, double EpsilonF):TMajorMinorEstimatorBase(RandomGenerator){
	epsilonF = EpsilonF;
	tmpGenotypeFrequencies = new TGenotypeFrequencies[6];
};

TMajorMinorEstimatorMLE::~TMajorMinorEstimatorMLE(){
	delete[] tmpGenotypeFrequencies;
};

double TMajorMinorEstimatorMLE::estimateGenotypeFrequencies(TMultiGLFData & data, int thisAlleleicCombination, TGlfConverter & glfConverter){
	data.fill(genotypeLikelihoods, thisAlleleicCombination);
	tmpGenotypeFrequencies[thisAlleleicCombination].estimate(genotypeLikelihoods, glfConverter, epsilonF);
	return tmpGenotypeFrequencies[thisAlleleicCombination].calculateLog10Likelihood(genotypeLikelihoods, glfConverter);
};

void TMajorMinorEstimatorMLE::findMLAllelicCombination(TMultiGLFData & data, TGlfConverter & glfConverter){
	//calculate L10L for each allelic combination
	for(uint16_t ac : usedAllelicCombinations){
		L10L_perCombination[ac] = estimateGenotypeFrequencies(data, ac, glfConverter);
	}

	//pick combination
	chooseBestAllelicCombinationAmongThoseWithEqualScores();
	genotypeFrequencies.set(tmpGenotypeFrequencies[bestAllelicCombination]);
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
TMajorMinor::TMajorMinor(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	vcfOpened = false;
	hasReference = false;
	randomGenerator = RandomGenerator;
};

void TMajorMinor::estimateMajorMinor(TParameters & params){
	//open GLF files
	TGlfMultiReader glfReader(params, logfile);
	glfReader.setAllActive();
	//TODO: add printAll option
	glfReader.onlyJumpToPositionsWithData();

	//add reference, if provided
	if(params.parameterExists("fasta")){
		logfile->list("Will only identify the most likely alternative allele (argument: fasta)");
		std::string fastaFile = params.getParameterString("fasta");
		logfile->list("Reading reference sequence from '" + fastaFile + "'");
		glfReader.addReference(fastaFile);
		hasReference = true;
	}

	//estimation method
	std::string method = params.getParameterStringWithDefault("method", "MLE");
	TMajorMinorEstimatorBase* MMEstimator;
	double maxF = params.getParameterDoubleWithDefault("maxF", 0.0000001);
	if(method == "Skotte"){
		logfile->list("Will estimate major / minor alleles using the Skotte method with maxF " + toString(maxF) + ".");
		MMEstimator = new TMajorMinorEstimatorSkotte(randomGenerator, maxF);
	} else if(method == "MLE"){
		logfile->list("Will estimate major / minor alleles using the MLE method with maxF " + toString(maxF) + ".");
		MMEstimator = new TMajorMinorEstimatorMLE(randomGenerator, maxF);
	} else throw "Unknown MajorMinor method '" + method + "'!";
	bool usePhredLikelihoods = params.parameterExists("phredLik");

	//think about filters
	int minSamplesWithData = params.getParameterIntWithDefault("minSamplesWithData", 1);
	if(minSamplesWithData <= 0)
		throw "minSamplesWithData must be >= 1!";
	if(minSamplesWithData > 0)
		logfile->list("Will only print sites for which at least " + toString(minSamplesWithData) + " samples have data.");
	int minVariantQuality = params.getParameterIntWithDefault("minVariantQual", 0);
	if(minVariantQuality > 0)
		logfile->list("Will only print sites with variant quality >= " + toString(minVariantQuality) + " samples have data.");

	//limit input
	long limitSites = params.getParameterDoubleWithDefault("limitSites", 0);
	if(limitSites > 0)
		logfile->list("Will stop at input position " + toString(limitSites) + ". (parameter 'limitSites')");
	if(limitSites < 0)
		throw "maxPos cannot be negative!";

	//filename tag
	std::string outname = params.getParameterStringWithDefault("out", "ATLAS_majorMinor");
	logfile->list("Will write output files with tag '" + outname + "'. (parameter 'out')");

	//open vcf file
	std::vector<std::string> sampleNames;
	glfReader.fillSampleNamesOfActiveFiles(sampleNames);
	TGlfMultiReaderVcf vcf(outname + ".vcf.gz", "ATLAS_GLF_Caller", sampleNames, randomGenerator);
	if(usePhredLikelihoods){
		vcf.usePhredScaledLikelihoods();
	}

	//vars
	logfile->startIndent("Parsing through glf files:");
	struct timeval start;
    gettimeofday(&start, NULL);
    long counter = 0;

	while(glfReader.readNext()){
		++counter;

		//filter on missingness
		if(glfReader.numActiveSamplesWithData() >= minSamplesWithData){
			//do major minor
			if(hasReference){
				Base ref = glfReader.refBase();
				MMEstimator->estimateMajorMinor(glfReader.data, glfReader.converter, ref);

				//write to VCF
				if(MMEstimator->variantQuality >= minVariantQuality){
					if(MMEstimator->major == ref){
						vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality, glfReader.data, MMEstimator->major, MMEstimator->minor);
					} else {
						vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality, glfReader.data, MMEstimator->minor, MMEstimator->major);
					}
				}
			} else {
				MMEstimator->estimateMajorMinor(glfReader.data, glfReader.converter);

				//write to VCF
				if(MMEstimator->variantQuality >= minVariantQuality){
					vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality, glfReader.data, MMEstimator->major, MMEstimator->minor);
				}
			}
		} //end filter on missingness

		//report progress
		if(counter % 1000000 == 0){
			static struct timeval end;
			gettimeofday(&end, NULL);
			float runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " positions in " + toString(runtime) + " min.");
		}

		//break?
		if(limitSites > 0 && counter == limitSites)
			break;
	}

	logfile->list("Reached end of glf files!");
	logfile->removeIndent();

	//clean storage
	delete MMEstimator;
};

}; //end namespace
