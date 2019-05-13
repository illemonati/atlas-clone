/*
 * TAlleleFrequencyEstimator.cpp
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */
/*
#include "TAlleleFrequencyEstimator.h"

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleHardyWeinbergFreqEstimator                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleHardyWeinbergFreqEstimator::TAlleleHardyWeinbergFreqEstimator(){
	maxIter = 1000;
	f = 0.0;
};

double TAlleleHardyWeinbergFreqEstimator::estimate(TPopulationLikehoodStorage storage, double epsilonF){
	double pGenotype[3];
	pGenotype[0] = 1.0; pGenotype[1] = 1.0; pGenotype[2] = 1.0; //initialize all to equal
	double weights[3];

	//run EM
	int iter = 0;
	double epsilon = epsilonF + 1.0;
	f = 0.0;
	while(iter < maxIter && epsilon > epsilonF){
		double old_f = f;

		//calculate sums
		double sum_1 = 0.0; double sum_2 = 0.0;
		for(int i=0; i<storage.numSamples; i++){
			//calculate weights
			weights[0] = storage[i][0] * pGenotype[0];
			weights[1] = storage[i][1] * pGenotype[1];
			weights[2] = storage[i][2] * pGenotype[2];
			double sum = weights[0] + weights[1] + weights[2];

			//add to sums
			sum_1 += weights[1] / sum;
			sum_2 += weights[2] / sum;
		}

		//estimate f
		f = (sum_1 + 2.0 * sum_2) / (2.0 * storage.numSamples);

		//recaluclate pGenotype
		pGenotype[0] = (1.0 - f) * (1.0 - f);
		pGenotype[1] = 2.0 * f * (1.0 - f);
		pGenotype[2] = f * f;

		//calculate F
		epsilon = fabs(f - old_f);
	}

	//return estimate
	return f;
};


////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleFreqEstimator                                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleFreqEstimator::TAlleleFreqEstimator(TParameters & Parameters, TLog* Logfile){
	vcfRead = false;
	_numLoci = 0;
	logfile = Logfile;
};

void TAlleleFreqEstimator::estimateAlleleFreq(TParameters & Parameters){
	if(vcfRead)
		throw "VCF already read!";

	//create reader
	bool saveAlleleFrequencies = true;
	TPopulationLikelihoodReader reader(Parameters, logfile, saveAlleleFrequencies);
	reader.doEstimateGenotypeFrequencies();

	// open vcf file
	vcfFilename = Parameters.getParameterString("vcf");
	logfile->startIndent("Estimating allele population frequencies from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename, logfile);

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	TPopulationLikehoodStorage storage(samples.numSamples());

	//create allele freqeuncy estimators
	TAlleleHardyWeinbergFreqEstimator HWestimator;
	double epsF = Parameters.getParameterDoubleWithDefault("epsF", 0.0000001);
	reader.doEstimateGenotypeFrequencies();

	//output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterStringWithDefault("out", tmp) + "_alleleFreq.txt.gz";
	logfile->list("Will write allele frequencies to file '" + outputName + "'.");
	TOutputFileZipped out(outputName);

	//write header
	out.writeHeader({"chr", "pos", "ref", "alt", "freqAltHW", "freqRefRef", "freqRefAlt", "freqAltAlt", "freqAltGeno"});

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(storage, samples, logfile)){
    	//print SNP
 		reader.writePosition(out);

 		//write HW estimates
 		out << HWestimator.estimate(storage, epsF);

 		//write genotype frequency estimates
 		reader.genotypeFrequencies()->writeDiploidFrequencies(out);
 		reader.genotypeFrequencies()->writeHaploidFrequencies(out);

 		//write frequency estimate based on genotype estimates
 		out << reader.allelFrequency() << std::endl;

 		//update for next
 		++_numLoci;
     }

    //clean up
	vcfRead = true;
	out.close();

    //report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};
*/
