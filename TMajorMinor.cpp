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

TMajorMinorEstimatorBase::TMajorMinorEstimatorBase(int NumSamples, TRandomGenerator* RandomGenerator){
	numSamples = NumSamples;
	genotypeLikelihoods = new double[3*numSamples];
	genotypeFrequencies = new double[3];
	L10L_perCombination = new double[6];
	allelicCombination = 0;
	minor = A;
	major = A;
	L10L = 0.0;
	randomGenerator = RandomGenerator;
	variantQuality = 0;
};

TMajorMinorEstimatorBase::~TMajorMinorEstimatorBase(){
	delete[] genotypeLikelihoods;
	delete[] genotypeFrequencies;
};

double TMajorMinorEstimatorBase::calculateLog10Likelihood(double* genotypeFrequencies){
	double LL = 0.0;
	for(int i = 0; i < 3 * numSamples; i += 3){
		LL += log10(genotypeLikelihoods[i] * genotypeFrequencies[0]
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

void TMajorMinorEstimatorBase::chooseMLAllelicCombinationAmongThoseWithEqualLikelihood(){
	//select MLE combination
	std::vector<int> MLE_combinations;
	for(int i=0; i<6; ++i){
		if(L10L_perCombination[i] == L10L)
			MLE_combinations.push_back(i);
	}
	allelicCombination = MLE_combinations[randomGenerator->pickOne(MLE_combinations.size())];
}

void TMajorMinorEstimatorBase::findMLAllelicCombination(uint8_t** phred){
	throw "Function TMajorMinorEstimatorBase::findMLAllelicCombination(uint8_t** phred) not implemented for base class!";
};

void  TMajorMinorEstimatorBase::estimateMajorMinor(uint8_t** phred){
	findMLAllelicCombination(phred);

	//which one is major?
	if(genotypeFrequencies[0] > genotypeFrequencies[2]){
		major = genoMap.alleleicCombinationToBase[allelicCombination][0];
		minor = genoMap.alleleicCombinationToBase[allelicCombination][1];
	} else {
		major = genoMap.alleleicCombinationToBase[allelicCombination][1];
		minor = genoMap.alleleicCombinationToBase[allelicCombination][0];

		//also flip frequencies
		double tmp = genotypeFrequencies[0];
		genotypeFrequencies[0] = genotypeFrequencies[2];
		genotypeFrequencies[2] = tmp;
	}

	//calculate variant quality
	int refHomIndex = genoMap.genotypeMap[major][major];
	double LL_fixed_phred = 0.0;
	for(int i=0; i<numSamples; ++i){
		LL_fixed_phred += phred[i][refHomIndex];
	}

	variantQuality = round(LL_fixed_phred + 10.0 * L10L);
};

//---------------------------------------------------
// TMajorMinorEstimatorSkotte
//---------------------------------------------------
TMajorMinorEstimatorSkotte::TMajorMinorEstimatorSkotte(int NumSamples, TRandomGenerator* RandomGenerator):TMajorMinorEstimatorBase(NumSamples, RandomGenerator){
	genotypeLikelihoods = new double[3*numSamples];

	priorGenotypeFrequencies = new double[3];
	priorGenotypeFrequencies[0] = 0.25;
	priorGenotypeFrequencies[1] = 0.50;
	priorGenotypeFrequencies[2] = 0.25;
};

TMajorMinorEstimatorSkotte::~TMajorMinorEstimatorSkotte(){
	delete[] priorGenotypeFrequencies;
};


void TMajorMinorEstimatorSkotte::findMLAllelicCombination(uint8_t** phred){
	//calculate L10L for each allelic combination
	fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[0]);
	L10L_perCombination[0] = calculateLog10Likelihood(priorGenotypeFrequencies);
	L10L = L10L_perCombination[0];
	for(int i=1; i<6; ++i){
		fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[i]);
		L10L_perCombination[i] = calculateLog10Likelihood(priorGenotypeFrequencies);
		if(L10L_perCombination[i] > L10L){
			L10L = L10L_perCombination[i];
		}
	}

	//pick combination
	chooseMLAllelicCombinationAmongThoseWithEqualLikelihood();

	//now guess genotype frequencies at MLE
	fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[allelicCombination]);
	guessGenotypeFrequencies(genotypeFrequencies);
}

//---------------------------------------------------
// TMajorMinorEstimatorMLE
//---------------------------------------------------
TMajorMinorEstimatorMLE::TMajorMinorEstimatorMLE(int NumSamples, TRandomGenerator* RandomGenerator, double MaxF):TMajorMinorEstimatorBase(NumSamples, RandomGenerator){
	maxF = MaxF;

	tmpGenotypeFrequencies = new double*[6];
	for(int i=0; i<6; ++i) tmpGenotypeFrequencies[i] = new double[3];
}

TMajorMinorEstimatorMLE::~TMajorMinorEstimatorMLE(){
	for(int i=0; i<6; ++i) delete[] tmpGenotypeFrequencies[i];
	delete[] tmpGenotypeFrequencies;
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
	estimateGenotypeFrequenciesEM(tmpGenotypeFrequencies[alleleicCombination]);
	return calculateLog10Likelihood(tmpGenotypeFrequencies[alleleicCombination]);
};

void TMajorMinorEstimatorMLE::findMLAllelicCombination(uint8_t** phred){
	//calculate L10L for each allelic combination
	fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[0]);
	L10L_perCombination[0] = estimateGenotypeFrequencies(phred, 0);
	L10L = L10L_perCombination[0];
	for(int i=1; i<6; ++i){
		fillLikelihoods(phred, genoMap.alleleicCombinationToGenotypes[i]);
		L10L_perCombination[i] = estimateGenotypeFrequencies(phred, i);
		if(L10L_perCombination[i] > L10L){
			L10L = L10L_perCombination[i];
		}
	}

	//pick combination
	chooseMLAllelicCombinationAmongThoseWithEqualLikelihood();

	genotypeFrequencies = tmpGenotypeFrequencies[allelicCombination];
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
TMajorMinor::TMajorMinor(TParameters & params, TLog* Logfile){
	logfile = Logfile;
	vcfOpened = false;

	//initialize random generator
	//TODO: do the random generator initialization in the task switcher?
	logfile->listFlush("Initializing random generator ...");
	if(params.parameterExists("fixedSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator=new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
};

void TMajorMinor::openVCF(std::string filenameTag, TGlfMultiReader & glfReader, bool usePhredLikelihoods){
	if(vcfOpened) closeVCF();

	//open vcf file
	filenameTag += ".vcf.gz";
	vcf.open(filenameTag.c_str());
	vcfOpened = true;

	//write info
	//TODO: create VCF class to harmonize code across different uses. Also include code in Tiger and other
	vcf << "##fileformat=VCFv4.3\n";
	vcf << "##source=ATLAS_GLF_Caller\n";
	glfReader.writeVCFHeader(vcf, usePhredLikelihoods);
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
		MMEstimator = new TMajorMinorEstimatorSkotte(glfReader.numSamples(), randomGenerator);
	} else if(method == "MLE"){
		double maxF = params.getParameterDoubleWithDefault("maxF", 0.00001);
		logfile->list("Will estimate major / minor alleles using the MLE method with maxF " + toString(maxF) + ".");
		MMEstimator = new TMajorMinorEstimatorMLE(glfReader.numSamples(), randomGenerator, maxF);
	} else throw "Unknown MajorMinor method '" + method + "'!";

	//output
	bool usePhredLikelihoods = params.parameterExists("phredLik");
	if(usePhredLikelihoods)
		logfile->list("Will write likelihoods as integers in phred format (PL tag in VCF).");
	else
		logfile->list("Will write log10(likelihoods) as float (GL tag in VCF).");

	//think about filters
	int minSamplesWithData = params.getParameterIntWithDefault("minSamplesWithData", 0);
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
			MMEstimator->estimateMajorMinor(glfReader.data);

			//filter on variant quality
			if(MMEstimator->variantQuality >= minVariantQuality){

				//write to VCF
				glfReader.writeSiteToVCF(vcf, MMEstimator->variantQuality, genoMap.genotypeMap[MMEstimator->major][MMEstimator->major], genoMap.genotypeMap[MMEstimator->major][MMEstimator->minor], genoMap.genotypeMap[MMEstimator->minor][MMEstimator->minor], randomGenerator, usePhredLikelihoods);
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
	closeVCF();
};



//---------------------------------------------------
// Skotte
//---------------------------------------------------

