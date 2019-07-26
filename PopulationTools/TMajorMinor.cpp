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

void TMajorMinorEstimatorBase::chooseBestAllelicCombinationAmongThoseWithEqualScores(){
	//select MLE combination
	std::vector<int> best_combinations;
	for(int i=0; i<6; ++i){
		if(L10L_perCombination[i] == L10L)
			best_combinations.push_back(i);
	}
	bestAllelicCombination = best_combinations[randomGenerator->pickOne(best_combinations.size())];
};

void TMajorMinorEstimatorBase::findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter){
	throw "Function TMajorMinorEstimatorBase::findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter) not implemented for base class!";
};

void  TMajorMinorEstimatorBase::estimateMajorMinor(TGlfMultiReader & glfReader, TGlfConverter & glfConverter){
	findMLAllelicCombination(glfReader, glfConverter);

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
	for(int i=0; i<glfReader.numActiveSamples(); ++i){
		if(glfReader.hasData[i]){
			if(glfReader.isHaploid[i])
				LL_fixed_glfPhred += glfReader.data[i][major];
			else
				LL_fixed_glfPhred += glfReader.data[i][refHomIndex];
		}
	}

	variantQuality = glfConverter.toPhred(LL_fixed_glfPhred - glfConverter.log10ToGlfFormat(L10L));
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

void TMajorMinorEstimatorSkotte::findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter){
	//calculate L10L for first allelic combination
	glfReader.fill(genotypeLikelihoods, 0);
	L10L_perCombination[0] = priorGenotypeFrequencies.calculateLog10Likelihood(genotypeLikelihoods, glfConverter);
	L10L = L10L_perCombination[0];

	//calculate L10L for all other allelic combination
	for(int i=1; i<6; ++i){
		glfReader.fill(genotypeLikelihoods, i);
		L10L_perCombination[i] = priorGenotypeFrequencies.calculateLog10Likelihood(genotypeLikelihoods, glfConverter);
		if(L10L_perCombination[i] > L10L){
			L10L = L10L_perCombination[i];
		}
	}

	//pick combination with highest likelihood
	chooseBestAllelicCombinationAmongThoseWithEqualScores();

	//now guess genotype frequencies at MLE
	glfReader.fill(genotypeLikelihoods, bestAllelicCombination);

	genotypeFrequencies.estimate(genotypeLikelihoods, glfConverter, epsilonF);
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

double TMajorMinorEstimatorMLE::estimateGenotypeFrequencies(TGlfMultiReader & glfReader, int thisAlleleicCombination, TGlfConverter & glfConverter){
	glfReader.fill(genotypeLikelihoods, thisAlleleicCombination);
	tmpGenotypeFrequencies[thisAlleleicCombination].estimate(genotypeLikelihoods, glfConverter, epsilonF);
	return tmpGenotypeFrequencies[thisAlleleicCombination].calculateLog10Likelihood(genotypeLikelihoods, glfConverter);
};

void TMajorMinorEstimatorMLE::findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter){
	//calculate L10L for each allelic combination
	L10L_perCombination[0] = estimateGenotypeFrequencies(glfReader, 0, glfConverter);
	L10L = L10L_perCombination[0];

	for(int i=1; i<6; ++i){
		L10L_perCombination[i] = estimateGenotypeFrequencies(glfReader, i, glfConverter);
		if(L10L_perCombination[i] > L10L){
			L10L = L10L_perCombination[i];
		}
	}

	//pick combination
	chooseBestAllelicCombinationAmongThoseWithEqualScores();
	genotypeFrequencies.set(tmpGenotypeFrequencies[bestAllelicCombination]);
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
TMajorMinor::TMajorMinor(TParameters & params, TLog* Logfile){
	file = NULL;
	logfile = Logfile;
	vcfOpened = false;

	//initialize random generator
	//TODO: do the random generator initialization in the task switcher?
	logfile->listFlush("Initializing random generator ...");
	if(params.parameterExists("fixedSeed")){
		randomGenerator = new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator = new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator = new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
};

TMajorMinor::~TMajorMinor(){
	delete randomGenerator;
};

void TMajorMinor::openVCF(std::string filenameTag, TGlfMultiReader & glfReader, bool usePhredLikelihoods){
    if(vcfOpened) closeVCF();

	//open vcf file
    file->open(filenameTag);

	vcfOpened = true;

	//write info
	//TODO: create VCF class to harmonize code across different uses. Also include code in Tiger and other
    glfReader.writeVCFHeader(file, usePhredLikelihoods);
}

void TMajorMinor::closeVCF(){
    file->close();
	vcfOpened = false;
}

void TMajorMinor::estimateMajorMinor(TParameters & params){
	//open GLF files
	TGlfMultiReader glfReader(params, logfile);
	glfReader.setAllActive();

	//estimation method
	std::string method = params.getParameterStringWithDefault("method", "MLE");
	TMajorMinorEstimatorBase* MMEstimator;
	double maxF = params.getParameterDoubleWithDefault("maxF", 0.0000001);
	if(method == "Skotte"){
		logfile->list("Will estimate major / minor alleles using the method of Skotte et al. (2012) with maxF " + toString(maxF) + ".");
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
		logfile->list("Will stop at input position " + toString(limitSites) + ".");
	if(limitSites < 0)
		throw "maxPos cannot be negative!";

	//filename tag
	std::string outname = params.getParameterStringWithDefault("out", "ATLAS_majorMinor");
	logfile->list("Will write output files with tag '" + outname + "'.");

    //create compression file
    std::string compression = params.getParameterStringWithDefault("compression", "zlib");
    int level_compression = params.getParameterIntWithDefault("levelCompression",9);
    if(level_compression < 1 || level_compression > 9){
        throw "Level compression must be a value between [1,9].";
    }

    if(compression=="lz4"){
#ifdef LZ4
        file = new LZ4Adapter();
        logfile->list("Will use LZ4 as compression library.");
        outname += ".glf.lz4";
#else
        throw "ATLAS has been compiled without LZ4. Could not use this type of compression.";
#endif
    }else if(compression=="zstd"){
#ifdef ZSTD
        file= new ZSTDAdapter(level_compression);
        logfile->list("Will use ZStandard as compression library.");
        outname += ".glf.zst";
#else
        throw "ATLAS has been compiled without ZStandard. Could not use this type of compression.";
#endif
    }else if(compression=="none"){
        file = new UncompressedAdapter();
        logfile->list("Will write data uncompressed.");
        outname += ".glf";
    }else{
        file = new GZStreamAdpater();
        logfile->list("Will use zlib as compression library.");
        outname += ".glf.gz";
    }

	//open vcf file    
    openVCF(outname, glfReader, usePhredLikelihoods);

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
			MMEstimator->estimateMajorMinor(glfReader, glfConverter);

			//filter on variant quality
			if(MMEstimator->variantQuality >= minVariantQuality){

				//write to VCF
				//glfReader.writeSiteToVCF(vcf, MMEstimator->variantQuality, genoMap.genotypeMap[MMEstimator->major][MMEstimator->major], genoMap.genotypeMap[MMEstimator->major][MMEstimator->minor], genoMap.genotypeMap[MMEstimator->minor][MMEstimator->minor], randomGenerator, usePhredLikelihoods);
                glfReader.writeSiteToVCF(file, MMEstimator->variantQuality, MMEstimator->major, MMEstimator->minor, randomGenerator, usePhredLikelihoods);
            }
		} //end filter on missngness

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
    closeVCF();

    delete file;
	delete MMEstimator;
};

