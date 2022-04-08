/*
 * TPolymorhhicWindowIdentifier.cpp
 *
 *  Created on: Feb 14, 2020
 *      Author: wegmannd
 */

#include "TPolymorhicWindowIdentifier.h"

#include <stdint.h>
#include <algorithm>
#include <exception>
#include <ostream>
#include <vector>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TFile.h"
#include "TPopulation.h"
#include "TPopulationLikelihoodLocus.h"
#include "TPopulationLikelihoods.h"
#include "TSampleLikelihoods.h"
#include "stringFunctions.h"
#include "strongTypes.h"

namespace PopulationTools{
using coretools::TParameters;
using coretools::TLog;
using coretools::TRandomGenerator;

TPolymorhicWindowIdentifier::TPolymorhicWindowIdentifier(TParameters &, TLog* Logfile){
	logfile = Logfile;
};

void TPolymorhicWindowIdentifier::identifyPolymorphicWindows(TParameters & Parameters, TRandomGenerator*){
	using BG = genometools::BiallelicGenotype;
	//read samples
	genometools::TPopulationSamples samples;
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameter<std::string>("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = Parameters.getParameter<std::string>("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
    genometools::TPopulationLikelihoodReaderWindow reader(Parameters, logfile, false);
	reader.openVCF(vcfFilename);
	logfile->endIndent();

	//Match samples
	if (samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	else
		samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//output file
	std::string tmp = coretools::str::readBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp) + "_polymorphicWindows.txt.gz";
	logfile->list("Will write polymoprhic state of windows to file '" + outputName + "'.");
	coretools::TOutputFile out(outputName);

	//write header
	std::vector<std::string> header = {"chr", "start", "end", "numWithData", "numMono", "numPoly"};
	samples.addSampleNamesToVector(header);
	out.writeHeader(header);

	//create likelihood window
    using TSampleLikelihoods = genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability>;
    genometools::TPopulationLikehoodWindow<TSampleLikelihoods> window(0, samples.numSamples());

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    uint64_t totalWindowsChecked = 0;
    uint64_t totalPolymorphicWindows = 0;
    while(reader.readDataFromVCF(window, samples)){
        int numWindowsWithData = 0;
        int numWindowsPoly = 0;
        std::vector<std::string> ind;

    	//write window
 		reader.writeWindow(out);

 		//write polymoprhic state for each sample
 		for(uint32_t i = 0; i<samples.numSamples(); ++i){
 			if(window.individualHasMissingData(i)){
 				ind.push_back("NA");
 			} else {
 				++numWindowsWithData;
 				uint32_t numPoly = 0;
 				for(uint32_t l = 0; l<window.numLoci(); ++l){
 					if(window[l][i][BG::het] < window[l][i][BG::homoFirst] && window[l][i][BG::het] < window[l][i][BG::homoSecond])
 						++numPoly;
 				}
 				if(numPoly > 0){
 					ind.push_back("1");
 					++numWindowsPoly;
 				} else {
 					ind.push_back("0");
 				}
 			}
 		}

 		//write data
 		out << numWindowsWithData << numWindowsWithData - numWindowsPoly << numWindowsPoly;
 		out.write(ind);
 		out << std::endl;

 		//update global counts
 		totalWindowsChecked += numWindowsWithData;
 		totalPolymorphicWindows += numWindowsPoly;
    }

    //report final status
	logfile->endIndent();
	reader.concludeFilters();
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";

	//print global stats
	logfile->conclude(totalPolymorphicWindows, " of ", totalWindowsChecked, " windows were found polymoprhic (", 100.0 * (double) totalPolymorphicWindows / (double) totalWindowsChecked, "%).");

	logfile->endIndent();
};

}; //end namespace
