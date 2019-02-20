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
		logfile->startIndent("Converting VCF file to lfmm-format...");
		convertToLfmm(parameters, samples, reader);
	} else throw "Unknown file format '" + convertTo + "'!";
}

void TVcfConverter::prepareReadingVcf(TParameters & parameters, TPopulationSamples & samples, TPopulationLikelihoodReader & reader){
	//read samples
	if(parameters.parameterExists("samples"))
		samples.readSamples(parameters.getParameterString("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = parameters.getParameterString("vcf");
	reader.openVCF(vcfFilename, logfile);
	logfile->endIndent();

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());
}

void TVcfConverter::openOutputFile(TParameters & parameters, std::string fileExtension, std::ofstream & file){
	//open output file
	std::string outname = parameters.getParameterStringWithDefault("out", "");
	if(outname == ""){
		// guess from filename
		outname = parameters.getParameterString("vcf");
		outname = extractBefore(outname, ".");
	}
	std::string filename = outname + "." + fileExtension;
	logfile->list("Will write the output to file '" + outname + "'.");
	file.open(filename.c_str());
	if(!file)
		throw "Failed to open file '" + filename + "' for writing!";
}

void TVcfConverter::convertToLfmm(TParameters & parameters, TPopulationSamples & samples, TPopulationLikelihoodReader & reader){
	// open output file
	std::ofstream lfmmFile;
	openOutputFile(parameters, "lfmm", lfmmFile);

	// initialize variables for vcf-file
	uint8_t* curLocus = new uint8_t[samples.numSamples() * 3];
	bool* sampleIsMissing = new bool[samples.numSamples()];
	std::vector<short*> genotypes;

	// run through VCF file
	logfile->startIndent("Parsing VCF file and reading genotypes:");
	while(reader.readDataFromVCF(curLocus, sampleIsMissing, samples, logfile)){
		// write header (chromosome and position) to output file
		lfmmFile << "\t" << reader.chr() << "/" << reader.position();

		// store genotype for every sample
		short* indivGenotypes = new short[samples.numSamples()];
		for(int s=0; s<samples.numSamples(); s++){
			short geno = reader.genotype(s);
			if(geno == 3) // replace missing genotype by 9 (lfmm format)
				indivGenotypes[s] = 9;
			else
				indivGenotypes[s] = geno;
		}
		genotypes.push_back(indivGenotypes);
	}
	lfmmFile << std::endl;

	// write to output file (lfmm: rows are individuals, columns are loci)
	for(int s=0; s<samples.numSamples(); s++){
		lfmmFile << samples.getNameFromOrderedIndex(s) << "\t";
		for (std::vector<short*>::iterator it = genotypes.begin(); it < genotypes.end(); it++){
			lfmmFile << (*it)[s] << "\t";
		}
		lfmmFile << std::endl;
	}

	// clean up
	delete[] curLocus;
	delete[] sampleIsMissing;
	for (std::vector<short*>::iterator it = genotypes.begin(); it < genotypes.end(); it++){
		delete[] *it;
	}

	// report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	logfile->endIndent();
}
