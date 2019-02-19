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
	// open readers
	TPopulationSamples samples;
	TPopulationLikelihoodReader reader(parameters, logfile, false);
	prepareReadingVcf(parameters, samples, reader);

	std::string convertTo = parameters.getParameterString("convertTo");
	if (convertTo == "lfmm"){
		logfile->startIndent("Start converting VCF file to lfmm-format.");
		convertToLfmm(parameters, samples, reader);
	} else throw "Unknown file format '" + convertTo + "'!";
}

void TVcfConverter::prepareReadingVcf(TParameters & parameters, TPopulationSamples & samples, TPopulationLikelihoodReader & reader){
	//read samples
	if(parameters.parameterExists("samples"))
		samples.readSamples(parameters.getParameterString("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = parameters.getParameterString("vcf");
	logfile->startIndent("Reading genotypes from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename, logfile);
	logfile->endIndent();

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());
}

std::ofstream TVcfConverter::openOutputFile(TParameters & parameters, std::string fileExtension){
	//open output file
	std::string outname = parameters.getParameterStringWithDefault("out", "");
	if(outname == ""){
		// guess from filename
		outname = parameters.getParameterString("vcf");
		outname = extractBefore(outname, ".");
	}
	std::string filename = outname + "." + fileExtension;
	logfile->list("Will write the output to file '" + outname + "'.");
	std::ofstream file(filename.c_str());
	if(!file)
		throw "Failed to open file '" + filename + "' for writing!";
	return file;
}

void TVcfConverter::convertToLfmm(TParameters & parameters, TPopulationSamples & samples, TPopulationLikelihoodReader & reader){
	// open output file
	std::ofstream lfmmFile = openOutputFile(parameters, "lfmm");

	// write header
	char sep = '\t';
	lfmmFile << "chr\tpos";
	for(int s=0; s<samples.numSamples(); s++)
		lfmmFile << "\t" << samples.getNameFromOrderedIndex(s);
	lfmmFile << "\n";

	// initialize variables for vcf-file
	uint8_t* curLocus = new uint8_t[samples.numSamples() * 3];
	bool* sampleIsMissing = new bool[samples.numSamples()];

	// run through VCF file
	logfile->startIndent("Parsing VCF file and reading genotypes:");
	while(reader.readDataFromVCF(curLocus, sampleIsMissing, samples, logfile)){
		// write chromosome and position
		lfmmFile << reader.chr() << sep << reader.position();

		// print genotype for every sample
		for(int s=0; s<samples.numSamples(); s++){
			short geno = reader.genotype(s);
			if(geno == 3) lfmmFile << "\t" << 9;
			else lfmmFile << "\t" << geno;
		}
		lfmmFile << std::endl;
	}

	// clean up
	delete[] curLocus;
	delete[] sampleIsMissing;

	// report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	logfile->endIndent();
}
