/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

#include "GenotypeTypes.h"
#include "debugtools.h"

namespace PopulationTools{

//---------------------------------------------------
// TMajorMinorEstimatorBase
//---------------------------------------------------
TMajorMinorEstimatorBase::TMajorMinorEstimatorBase(TRandomGenerator* RandomGenerator){
	bestAllelicCombination = genometools::AllelicCombination::NN;
	minor = genometools::Base::N;
	major = genometools::Base::N;
	L10L = 0.0;
	randomGenerator = RandomGenerator;
	variantQuality = 0;
};

void TMajorMinorEstimatorBase::useAllAlleleicCombinations(){
	usedAllelicCombinations.clear();
	for(AllelicCombination ac = AllelicCombination::min; ac < AllelicCombination::max; ++ac){
		usedAllelicCombinations.push_back(ac);
	}
};

void TMajorMinorEstimatorBase::useAllelicCombinationsThatContain(const Base & base){
	usedAllelicCombinations.clear();

	//get the three alleleic combinations that contain base
	for(Base b = Base::min; b < Base::max; ++b){
		if(b != base){
			usedAllelicCombinations.push_back(allelicCombination(base, b));
		}
	}
};

void TMajorMinorEstimatorBase::chooseBestAllelicCombinationAmongThoseWithEqualScores(){
	//identify max L10L
	L10L = L10L_perCombination.max();

	//select MLE combination
	std::vector<AllelicCombination> best_combinations;
	for(AllelicCombination ac = AllelicCombination::min; ac < AllelicCombination::max; ++ac){
		if(L10L_perCombination[ac] == L10L)
			best_combinations.push_back(ac);
	}
	bestAllelicCombination = best_combinations[randomGenerator->sample(best_combinations.size())];
};

void TMajorMinorEstimatorBase::findMLAllelicCombination(const TMultiGLFData &){
	throw "Function TMajorMinorEstimatorBase::findMLAllelicCombination(TGlfMultiReader & glfReader, TGlfConverter & glfConverter) not implemented for base class!";
};

void  TMajorMinorEstimatorBase::_estimateMajorMinor(const TMultiGLFData & data){
	findMLAllelicCombination(data);

	//which one is major?
	if(genotypeFrequencies.alleleFrequency() < 0.5){
		major = first(bestAllelicCombination);
		minor = second(bestAllelicCombination);
	} else {
		major = second(bestAllelicCombination);
		minor = first(bestAllelicCombination);

		//also flip frequencies
		genotypeFrequencies.flip();
	}

	//calculate variant quality
	const auto refHom = genometools::genotype(major, major);
	coretools::Log10Probability LL_fixed_glfPhred = 0.0;
	for(uint32_t i=0; i<data.size(); ++i){
		if(data[i].hasData()){
			if(data[i].isHaploid())
				LL_fixed_glfPhred += (coretools::Log10Probability) data[i][major];
			else
				LL_fixed_glfPhred += (coretools::Log10Probability) data[i][refHom];
		}
	}

	//set variant quality
	if(LL_fixed_glfPhred > L10L){
		//can happen if MLE estimation is not perfect
		//difference is usually super smnall but enough to cause Log210Probability to complain it it out of range.
		variantQuality = coretools::Log10Probability(0.0);
	} else {
		variantQuality = LL_fixed_glfPhred - L10L;
	}
};

void  TMajorMinorEstimatorBase::estimateMajorMinor(const TMultiGLFData & data){
	//use all allelic combinations
	useAllAlleleicCombinations();

	//now estimate
	_estimateMajorMinor(data);
};

void  TMajorMinorEstimatorBase::estimateMajorMinor(const TMultiGLFData & data, const Base & base){
	//use all allelic combinations
	useAllelicCombinationsThatContain(base);

	//now estimate
	_estimateMajorMinor(data);
};

//---------------------------------------------------
// TMajorMinorEstimatorSkotte
//---------------------------------------------------
TMajorMinorEstimatorSkotte::TMajorMinorEstimatorSkotte(TRandomGenerator* RandomGenerator, double EpsilonF):TMajorMinorEstimatorBase(RandomGenerator){
	using BG = genometools::BiallelicGenotype;
	 epsilonF = EpsilonF;

	//diploid
	priorGenotypeFrequencies[BG::homoFirst] = 0.25;
	priorGenotypeFrequencies[BG::het] = 0.50;
	priorGenotypeFrequencies[BG::homoSecond] = 0.25;

	//haploid
	priorGenotypeFrequencies[BG::haploidFirst] = 0.50;
	priorGenotypeFrequencies[BG::haploidSecond] = 0.50;
};

void TMajorMinorEstimatorSkotte::findMLAllelicCombination(const TMultiGLFData & data){
	//calculate L10L for each allelic combination used
	for(auto& ac : usedAllelicCombinations){
		data.fill(genotypeLikelihoods, ac);
		L10L_perCombination[ac] = priorGenotypeFrequencies.calculateLog10Likelihood(genotypeLikelihoods, genotypeLikelihoods.size());
		L10L = L10L_perCombination[ac];
	}

	//pick combination with highest likelihood
	chooseBestAllelicCombinationAmongThoseWithEqualScores();

	//now estimate genotype frequencies at MLE allelic combination
	data.fill(genotypeLikelihoods, bestAllelicCombination);
	genotypeFrequencies.estimate(genotypeLikelihoods, genotypeLikelihoods.size(), epsilonF);

	//calculate likelihood again with better genotype frequencies
	L10L_perCombination[bestAllelicCombination] = genotypeFrequencies.calculateLog10Likelihood(genotypeLikelihoods, genotypeLikelihoods.size());
	L10L = L10L_perCombination[bestAllelicCombination];
};

//---------------------------------------------------
// TMajorMinorEstimatorMLE
//---------------------------------------------------
TMajorMinorEstimatorMLE::TMajorMinorEstimatorMLE(TRandomGenerator* RandomGenerator, double EpsilonF):TMajorMinorEstimatorBase(RandomGenerator){
	epsilonF = EpsilonF;
};

coretools::Log10Probability TMajorMinorEstimatorMLE::estimateGenotypeFrequencies(const TMultiGLFData & data, const AllelicCombination & thisAlleleicCombination){
	using genometools::index;
	data.fill(genotypeLikelihoods, thisAlleleicCombination);
	tmpGenotypeFrequencies[index(thisAlleleicCombination)].estimate(genotypeLikelihoods, genotypeLikelihoods.size(), epsilonF);
	return tmpGenotypeFrequencies[index(thisAlleleicCombination)].calculateLog10Likelihood(genotypeLikelihoods, genotypeLikelihoods.size());
};

void TMajorMinorEstimatorMLE::findMLAllelicCombination(const TMultiGLFData & data){
	//calculate L10L for each allelic combination
	for(auto& ac : usedAllelicCombinations){
		L10L_perCombination[ac] = estimateGenotypeFrequencies(data, ac);
	}

	//pick combination
	chooseBestAllelicCombinationAmongThoseWithEqualScores();
	genotypeFrequencies.set(tmpGenotypeFrequencies[index(bestAllelicCombination)]);
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
TMajorMinor::TMajorMinor(TLog* Logfile, TParameters &, TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	hasReference = false;
	randomGenerator = RandomGenerator;
	minSamplesWithData = 1;
	minVariantQuality = 0;
};

void TMajorMinor::estimateMajorMinor(TParameters & params){
	//open GLF files
	GLF::TGlfMultiReader glfReader(params, logfile);
	glfReader.setAllActive();

	//add reference, if provided
	if(params.parameterExists("fasta")){
		logfile->list("Will only identify the most likely alternative allele (argument: fasta)");
		std::string fastaFile = params.getParameter<std::string>("fasta");
		logfile->list("Reading reference sequence from '" + fastaFile + "'");
		glfReader.addReference(fastaFile);
		hasReference = true;
	}

	//estimation method
	std::string method = params.getParameterWithDefault<std::string>("method", "MLE");
	TMajorMinorEstimatorBase* MMEstimator;
	double maxF = params.getParameterWithDefault("maxF", 0.0000001);
	if(method == "Skotte"){
		logfile->list("Will estimate major / minor alleles using the Skotte method with maxF ", maxF, ". (parameters method and maxF)");
		MMEstimator = new TMajorMinorEstimatorSkotte(randomGenerator, maxF);
	} else if(method == "MLE"){
		logfile->list("Will estimate major / minor alleles using the MLE method with maxF ", maxF, ". (parameters method and maxF)");
		MMEstimator = new TMajorMinorEstimatorMLE(randomGenerator, maxF);
	} else throw "Unknown MajorMinor method '" + method + "'!";

	bool usePhredLikelihoods = params.parameterExists("phredLik");
	if(usePhredLikelihoods){
		logfile->list("Will write phred-scaled likelihoods. (parameter phredLik)");
	} else {
		logfile->list("Will write raw likelihoods. (use phredLik to phred-scale)");
	}

	//read filters
	if(params.parameterExists("printAll")){
		logfile->list("Will all sites and samples. (parameter printAll)");
		minSamplesWithData = 0;
		minVariantQuality = 0;
	} else {
		minSamplesWithData = params.getParameterWithDefault<uint32_t>("minSamplesWithData", 1);
		if(minSamplesWithData > 0){
			logfile->list("Will only print sites for which at least ", minSamplesWithData, " samples have data. (parameter minSamplesWithData)");
		}

		minVariantQuality = params.getParameterWithDefault<genometools::PhredIntProbability>("minVariantQual", genometools::PhredIntProbability::highest());
		if(minVariantQuality > genometools::PhredIntProbability::highest()){
			logfile->list("Will only print sites with variant quality >= ", minVariantQuality, " samples have data. (parameter minVariantQual)");
		}
	}

	if(minSamplesWithData > 0){
		glfReader.onlyJumpToPositionsWithData(true);
	} else {
		glfReader.onlyJumpToPositionsWithData(false);
	}

	//limit input
	long limitSites = params.getParameterWithDefault("limitSites", 0);
	if(limitSites > 0)
		logfile->list("Will stop at input position ", limitSites, ". (parameter 'limitSites')");
	if(limitSites < 0)
		throw "maxPos cannot be negative!";

	//filename tag
	std::string outname = params.getParameterWithDefault<std::string>("out", "ATLAS_majorMinor");
	logfile->list("Will write output files with tag '" + outname + "'. (parameter 'out')");

	//get sample names. Make sure they are unique in the vcf
	std::vector<std::string> sampleNames;
	if(params.parameterExists("sampleNames")){
		logfile->startIndent("Using the following sample names (parameter 'sampleNames'):");
		params.fillParameterIntoContainer("sampleNames", sampleNames, ',');

		if(sampleNames.size() != glfReader.numActiveSamples()){
			throw "Number of provided same names does not match number of GLF files!";
		}

		//report
		auto glfNames = glfReader.namesOfActiveFiles();
		for(uint32_t i = 0; i < glfReader.numActiveSamples(); ++i){
			logfile->list(glfNames[i], " -> ", sampleNames[i]);
		}
		logfile->endIndent();
	} else {
		logfile->list("Will deduce sample names from GLF file names. (use 'sampleNames' to provide alternative names)");
		sampleNames = glfReader.sampleNamesOfActiveFiles();

		//if there are duplicates, add suffix
		bool foundDuplicates = false;
		for(auto& s : sampleNames){
			uint32_t c = 0;
			for(auto& t : sampleNames){
				if(t == s){
					++c;
					if(c > 1){
						t = t + "." + coretools::str::toString(c);

						 //report
						if(!foundDuplicates){
							logfile->startIndent("Duplicate samples will be rename as follows:");
							foundDuplicates = true;
						}
						logfile->list(s, " -> ", t);
					}
				}
			}
		}
		if(foundDuplicates){
			logfile->endIndent();
		}
	}

	//open vcf file
	GLF::TGlfMultiReaderVcf vcf(outname + ".vcf.gz", "ATLAS_GLF_Caller", sampleNames, usePhredLikelihoods, randomGenerator);

	//vars
	logfile->startIndent("Parsing through glf files:");
	coretools::TTimer timer;
    long counter = 0;

	while(glfReader.readNext()){
		++counter;
		//filter on missingness
		if(glfReader.numActiveSamplesWithData() >= minSamplesWithData){
			//do major minor
			if(hasReference){
				Base ref = glfReader.refBase();
				MMEstimator->estimateMajorMinor(glfReader.data, ref);

				//write to VCF
				if(MMEstimator->variantQuality >= minVariantQuality){
					if(MMEstimator->major == ref){
						vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality, glfReader.data, MMEstimator->major, MMEstimator->minor);
					} else {
						vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality, glfReader.data, MMEstimator->minor, MMEstimator->major);
					}
				}
			} else {
				MMEstimator->estimateMajorMinor(glfReader.data);

				//write to VCF
				if(MMEstimator->variantQuality >= minVariantQuality){
					vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality, glfReader.data, MMEstimator->major, MMEstimator->minor);
				}
			}
		} //end filter on missingness

		//report progress
		if(counter % 1000000 == 0){
			logfile->list("Parsed ", counter, " positions in ", timer.formattedTime(), " min.");
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
