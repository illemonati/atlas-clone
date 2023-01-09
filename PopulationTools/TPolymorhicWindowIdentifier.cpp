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

#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/strongTypes.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/TSampleLikelihoods.h"
#include "genometools/VCF/TPopulation.h"
#include "genometools/VCF/TPopulationLikelihoodLocus.h"
#include "genometools/VCF/TPopulationLikelihoods.h"

namespace PopulationTools{
using coretools::instances::parameters;
using coretools::instances::logfile;

void TPolymorhicWindowIdentifier::identifyPolymorphicWindows() {
	using BG = genometools::BiallelicGenotype;
	//read samples
	genometools::TPopulationSamples samples;
	if(parameters().parameterExists("samples"))
		samples.readSamples(parameters().getParameter<std::string>("samples"));

	//open VCF reader
	std::string vcfFilename = parameters().getParameter<std::string>("vcf");
	logfile().startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
    genometools::TPopulationLikelihoodReaderWindow reader(false);
	reader.openVCF(vcfFilename);
	logfile().endIndent();

	//Match samples
	if (samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	else
		samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//output file
	auto tmp = coretools::str::readBeforeLast(vcfFilename, ".vcf");
	std::string outputName = parameters().getParameterWithDefault("out", tmp) + "_polymorphicWindows.txt.gz";
	logfile().list("Will write polymoprhic state of windows to file '" + outputName + "'.");
	coretools::TOutputFile out(outputName);

	//write header
	std::vector<std::string> header = {"chr", "start", "end", "numWithData", "numMono", "numPoly"};
	samples.addSampleNamesToVector(header);
	out.writeHeader(header);

	//create likelihood window
    using TSampleLikelihoods = genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability>;
    genometools::TPopulationLikehoodWindow<TSampleLikelihoods> window(0, samples.numSamples());

    //run through VCF file
    logfile().startIndent("Parsing VCF file:");
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
 		out.writeln(numWindowsWithData, numWindowsWithData - numWindowsPoly, numWindowsPoly, ind);

 		//update global counts
 		totalWindowsChecked += numWindowsWithData;
 		totalPolymorphicWindows += numWindowsPoly;
    }

    //report final status
	logfile().endIndent();
	reader.concludeFilters();
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";

	//print global stats
	logfile().conclude(totalPolymorphicWindows, " of ", totalWindowsChecked, " windows were found polymoprhic (", 100.0 * (double) totalPolymorphicWindows / (double) totalWindowsChecked, "%).");

	logfile().endIndent();
};

}; //end namespace
