/*
 * TPolymorhhicWindowIdentifier.cpp
 *
 *  Created on: Feb 14, 2020
 *      Author: wegmannd
 */


#include <TPolymorhicWindowIdentifier.h>


TPolymorhicWindowIdentifier::TPolymorhicWindowIdentifier(TParameters & Parameters, TLog* Logfile){
	logfile = Logfile;
};


double TPolymorhicWindowIdentifier::_calcQualityPolymorphic(const TPopulationLikehoodWindow & window, const uint32_t & i){
	//calculate likelihood for polymorphic and monomorphic
};


void TPolymorhicWindowIdentifier::identifyPolymorphicWindows(TParameters & Parameters, TRandomGenerator* randomGenerator){
	//read samples
	TPopulationSamples samples;
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameterString("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = Parameters.getParameterString("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	TPopulationLikelihoodReaderWindow reader(Parameters, logfile, false);
	reader.openVCF(vcfFilename);
	logfile->endIndent();

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//output file
	std::string tmp = readBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterStringWithDefault("out", tmp) + "_polymorphicWindows.txt.gz";
	logfile->list("Will write polymoprhic state of windows to file '" + outputName + "'.");
	TOutputFileZipped out(outputName);

	//write header
	std::vector<std::string> header = {"chr", "start", "end"};
	samples.addOrderedSampleNamesToVector(header);
	out.writeHeader(header);

	//create likelihood window
	TPopulationLikehoodWindow window(0, samples.numSamples());

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(window, samples, glfConverter)){
    	//write window
 		reader.writeWindow(out);

 		//write polymoprhic state for each sample
 		for(uint32_t i = 0; i<samples.numSamples(); ++i){
 			if(window.individualHasMissingData(i)){
 				out << "NA";
 			} else {
 				uint32_t numPoly = 0;
 				for(uint32_t l = 0; l<window.numLoci(); ++l){
 					if(window[l][i].glfLikelihood_1 > window[l][i].glfLikelihood_0 && window[l][i].glfLikelihood_1 > window[l][i].glfLikelihood_2)
 						++numPoly;
 				}
 				if(numPoly > 0)
 					out << "1";
 				else
 					out << "0";
 			}
 		}

 		//end line
 		out << std::endl;
    }

    //report final status
	logfile->endIndent();
	reader.concludeFilters();
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};
