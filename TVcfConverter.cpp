/*
 * TVcfConverter.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: caduffm
 */
#include "TVcfConverter.h"

TVcfConverter::TVcfConverter(TParameters & Parameters, TLog* Logfile){
	logfile = Logfile;
}


void TVcfConverter::convertVcf(TParameters & parameters){
	std::string convertTo = parameters.getParameterString("convertTo");
	if (convertTo == "lfmm"){
		logfile->startIndent("Start converting VCF file to lfmm-format.");
		convertToLfmm(parameters);
	} else throw "Unknown file format '" + convertTo + "'!";
}


void TVcfConverter::convertToLfmm(TParameters & parameters){
	//read samples
	TPopulationSamples samples;
	if(parameters.parameterExists("samples"))
		samples.readSamples(parameters.getParameterString("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = parameters.getParameterString("vcf");
	logfile->startIndent("Reading genotypes from VCF file '" + vcfFilename + "':");
	TPopulationLikelihoodReader reader(parameters, logfile, false);
	reader.openVCF(vcfFilename, logfile);
	logfile->endIndent();

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//open output file
	std::string outname = parameters.getParameterStringWithDefault("out", "");
	if(outname == ""){
		// guess from filename
		outname = parameters.getParameterString("vcf");
		outname = extractBefore(outname, ".");
	}
	std::string filename = outname + ".lfmm";
	logfile->list("Will write the lfmm genotype matrix to file '" + outname + "'.");
	std::ofstream lfmmFile(filename.c_str());
	if(!lfmmFile)
		throw "Failed to open file '" + filename + "' for writing!";

	//write header // TODO: not sure whether lfmm format accepts a header, check!
	char sep = '\t';
	lfmmFile << "chr\tpos";
	for(int s=0; s<samples.numSamples(); s++)
		lfmmFile << "\t" << samples.getNameFromOrderedIndex(s);
	lfmmFile << "\n";

	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	uint8_t* curLocus = new uint8_t[samples.numSamples() * 3];
	bool* sampleIsMissing = new bool[samples.numSamples()];

	//run through VCF file
	logfile->startIndent("Parsing VCF file and reading genotypes:");
	while(reader.readDataFromVCF(curLocus, sampleIsMissing, samples, logfile)){
		//write chromosome and position
		lfmmFile << reader.chr() << sep << reader.position();

		//print genotype for every sample
		for(int s=0; s<samples.numSamples(); s++){
			if (sampleIsMissing[s])
				lfmmFile << "\t" << 9; // missing genotypes are encoded by 9 in lfmm
			else
				lfmmFile << "\t" << findMaxGenotype(&curLocus[3*s]);
		}
		lfmmFile << std::endl;
	}

	//clean up
	delete[] curLocus;
	delete[] sampleIsMissing;

	//report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	logfile->endIndent();
}

int TVcfConverter::findMaxGenotype(uint8_t* phred){
	int observedGenotype = 0;
	for (int genotype = 1; genotype < 3; genotype++){
		if (phred[genotype] < phred[observedGenotype])
			 observedGenotype = genotype;
	}
	return observedGenotype;
}

