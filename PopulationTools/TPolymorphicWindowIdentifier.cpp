/*
 * TPolymorhhicWindowIdentifier.cpp
 *
 *  Created on: Feb 14, 2020
 *      Author: wegmannd
 */

#include "TPolymorphicWindowIdentifier.h"

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

void TPolymorphicWindowIdentifier::run() {
	using BG = genometools::BiallelicGenotype;
	//read samples
	genometools::TPopulationSamples samples;
	if(parameters().exists("samples"))
		samples.readSamples(parameters().get<std::string>("samples"));

	//open VCF reader
	std::string vcfFilename = parameters().get<std::string>("vcf");
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
	std::string outputName = parameters().get("out", tmp) + "_polymorphicWindows.txt.gz";
	logfile().list("Will write polymorphic state of windows to file '" + outputName + "'.");
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
    size_t totalWindowsChecked = 0;
    size_t totalPolymorphicWindows = 0;
    while(reader.readDataFromVCF(window, samples)){
        int numWindowsWithData = 0;
        int numWindowsPoly = 0;
        std::vector<std::string> ind;

    	//write window
 		reader.writeWindow(out);

 		//write polymoprhic state for each sample
 		for (size_t i = 0; i<samples.numSamples(); ++i){
 			if(window.individualHasMissingData(i)){
 				ind.push_back("NA");
 			} else {
 				++numWindowsWithData;
 				size_t numPoly = 0;
 				for(size_t l = 0; l<window.numLoci(); ++l){
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
		UERROR("No usable loci in VCF file '", vcfFilename, "'!");

	//print global stats
	logfile().conclude(totalPolymorphicWindows, " of ", totalWindowsChecked, " windows were found polymoprhic (", 100.0 * (double) totalPolymorphicWindows / (double) totalWindowsChecked, "%).");

	logfile().endIndent();
};

}; //end namespace
