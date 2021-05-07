/*
 * TVcfDiagnostics.cpp
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */
#include "TVcfDiagnostics.h"

#include <algorithm>
#include <typeinfo>
#include <sstream>
#include "gzstream.h"

namespace VCF{

//---------------------------------------------------
//TAnnotator
//---------------------------------------------------

VcfDiagnostics::VcfDiagnostics(TParameters & Params, TLog* Logfile){
    logfile = Logfile;
    randomGeneratorInitialized = false;
    randomGenerator = NULL;
    chr=-1;

    //open vcf file
    openVCF(Params.getParameter<std::string>("vcf"), vcfFile);

    //outputname
    outname = Params.getParameterWithDefault<std::string>("out", "");
    if(outname == ""){
        //guess from filename
        outname = vcfFile.filename;
        outname = extractBeforeLast(outname, ".");
        if(isZipped)
            //if zipped there is extra .gz
            outname = extractBeforeLast(outname, ".");

    }
    logfile->list("Writing output files with prefix '" + outname + "'.");
}

//open input stream
void VcfDiagnostics::openVCF(std::string filename, TVcfFile_base & vcfFile){
    //open vcf file
    if(filename.find(".gz") == std::string::npos){
        logfile->list("Reading vcf from file '" + filename + "'.");
        isZipped = false;
    } else {
        logfile->list("Reading vcf from gzipped file '" + filename + "'.");
        isZipped = true;
    }

    vcfFile.openStream(filename, isZipped);

    //enable parsers
    vcfFile.enablePositionParsing();
    vcfFile.enableVariantParsing();
    vcfFile.enableVariantQualityParsing();
    vcfFile.enableFormatParsing();
    vcfFile.enableSampleParsing();
    vcfFile.enableInfoParsing();
}

void VcfDiagnostics::initializeRandomGenerator(){
	if(!randomGeneratorInitialized){
		randomGenerator=new TRandomGenerator();
		logfile->list("Random generator initialized with seed " + toString(randomGenerator->usedSeed));
		randomGeneratorInitialized=true;
	}
}

void VcfDiagnostics::vcfToInvariantBed(){
	//open vcf file
	logfile->list("Writing sites that are invariant across individuals to bed:");

	//open output
	logfile->list("Writing files to '" + outname + ".bed.gz'");
	gz::ogzstream bedFile((outname + (std::string) ".bed.gz").c_str());

	//parse vcf file
	std::map<char, int> bases = {{'A', 1}, {'C', 1}, {'G', 1}, {'T', 1}};
	int counter = 0;
	int curStartRegion = -1;
	bool previousStartIsVariant = false;
	std::string curChr;
	bool updateStart = true;
	long lastPosition;
	while(vcfFile.next()){
		++counter;
		if(updateStart){
			curStartRegion = vcfFile.position();
			curChr = vcfFile.chr();
			updateStart = false;
		}
		if(bases.find(vcfFile.getRefAllele()[0]) == bases.end()) {//!= 'A' && vcfFile.getRefAllele() != 'C' && vcfFile.getRefAllele() != 'G' && vcfFile.getRefAllele() != 'T'
			continue; //ignore indels
		}

		if(curChr != vcfFile.chr()){
			//is previous site invariant? -> write to file
			if(!previousStartIsVariant){
				bedFile << curChr << "\t" << curStartRegion - 1 << "\t" << lastPosition << std::endl;
			} else {
			}
			//update start directly, don't just set to true
			curStartRegion = vcfFile.position();
			curChr = vcfFile.chr();
		}

		int indCounter = 0;
		char allele = vcfFile.getFirstAlleleOfSample(0)[0];
		for(int ind=0; ind<vcfFile.numSamples(); ++ind){
			++indCounter;
			if(vcfFile.getFirstAlleleOfSample(ind)[0] != allele || vcfFile.getSecondAlleleOfSample(ind)[0] != allele){
				//there was a variant site, is previous site invariant? -> write to file
				if(previousStartIsVariant == false && counter != 1){
					bedFile << vcfFile.chr() << "\t" << curStartRegion - 1 << "\t" << vcfFile.position() - 1 << std::endl;
				}
				updateStart = true;
				previousStartIsVariant = true;
				break;
			}
			else {
				previousStartIsVariant = false;
			}
		}


//		if(indCounter == vcfFile.numSamples()){
//			//bed is 0-based
//			bedFile << vcfFile.chr() << "\t" << vcfFile.position() - 1 << "\t" << vcfFile.position() << std::endl;
//		}

		//in case of chr change, this is needed to write last region of previous chr
		lastPosition = vcfFile.position();

	}

	//write last region to file
	if(!previousStartIsVariant){
		bedFile << vcfFile.chr() << "\t" << curStartRegion - 1 << "\t" << vcfFile.position() << std::endl;
	}
	bedFile.close();
}

int VcfDiagnostics::findLastPassedFilterIndex(int obsValue, std::vector<int> & filtersAscendingOrder){
	int lastPassedFilterIndex = 0;
	for(unsigned int i=0; i<filtersAscendingOrder.size(); ++i){
		if(obsValue < filtersAscendingOrder[i])
			break;
		else
			lastPassedFilterIndex = i;
	}
	return lastPassedFilterIndex;
}

void VcfDiagnostics::assessAllelicImbalance(TParameters & Params){
	//output
	logfile->list("Writing files to '" + outname + "_allelicDepth.txt'");

	//limit input?
	int maxDP = Params.getParameterWithDefault<int>("maxDepth", 100);
	logfile->list("Ignoring sites with depth larger than " + toString(maxDP) + ".");

	int inputLines = Params.getParameterWithDefault<int>("inputLines", -1);
	if(inputLines <= 0){
		logfile->list("Reading whole vcf.");
	} else
		logfile->list("Limiting input to " + toString(inputLines) + " lines.");

	//initialize tables
	logfile->startIndent("Initializing count tables:");
	std::string qualityString = Params.getParameterWithDefault<std::string>("qualities", "0,10,20,30,40,50");
	std::vector<int> qualities;
	fillContainerFromString(qualityString, qualities, ',');

	std::vector<TCountTable*> countTables;
	for(unsigned int i=0; i<qualities.size(); ++i){
		countTables.emplace_back(new TCountTable(maxDP, maxDP, outname + "_qual" + toString(qualities[i]) + "_allelicDepth.txt", logfile));
	}
	logfile->endIndent();

	//temp variables
	int counter = 0;
	int numHet = 0;

	logfile->startIndent("Parsing vcf file:");
	while(vcfFile.next()){
		++counter;

		//get allelic depth at hetero sites across all samples. One table per quality filter.
		for(int i=0; i < vcfFile.numSamples(); ++i){
			if(!vcfFile.sampleIsMissing(i) && vcfFile.sampleIsHeteroRefNonref(i)){
				++numHet;

				//which position in table does site correspond to?
				std::vector<std::string> tmp;
				std::string tag="AD";
				fillContainerFromString(vcfFile.getSampleContentAt(tag, i), tmp, ',');
				int numRef = convertString<int>(tmp[0]);
				int numAlt = convertString<int>(tmp[1]);
				if(numRef == 0 || numAlt == 0)
					throw "Call at position " + toString(vcfFile.position()) + " is heterozygous but reference or alternative allelic depth is 0!";
				if(vcfFile.depthAsIntNoCheckForMissingSample("DP", i) > maxDP){
					logfile->warning("DP is " + toString(vcfFile.depthAsIntNoCheckForMissingSample("DP", i)) + " at pos " + toString(vcfFile.position()) + ". This site will be ignored.");
					continue;
				}

				//add count to correct table
				int quality = convertString<int>(vcfFile.getSampleContentAt("GQ", i));
				int index = findLastPassedFilterIndex(quality, qualities);
				for(int i=0; i<(index+1); ++i){
					++(countTables.at(i))->table[numRef][numAlt];
				}
			}
		}
		if(verbose && counter % 1000000 == 0)
			logfile->list("read " + toString(counter) + " lines, " + toString(numHet) + " heterozygous sites were found.");
		if(inputLines != -1 && counter > inputLines){
			logfile->list("reached input limit!");
			break;
		}
	}

	logfile->endIndent("done!");

	std::string description = "ref_depth/alt_depth";
	std::string rowPrefix = "ref_";
	std::string colPrefix = "alt_";


	//write tables
	for(unsigned int i=0; i<countTables.size(); ++i){
		countTables[i]->writeTable(description, rowPrefix, colPrefix);
		delete countTables[i];
	}

	//clean up
//	for(int i=0; i<countTables.size(); ++i){
//		delete countTables.at(i);
//	}

}

/*
void VcfDiagnostics::filterAllelicImbalance(){
	//open vcf file
	if(verbose) std::cerr << " - Filtering sites with allelic imbalance:" << std::endl;

	double pValFilter = params->getParameterWithDefault("pValFilter", 0.001);

	openVCF(vcfFile);

	//enable parsers
	vcfFile.enablePositionParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();
	vcfFile.enableVariantParsing();
	vcfFile.enableWriting();

	//prepare output vcf
	prepareVcfOutput(vcfFile);
	vcfFile.writeHeaderVCF_4_0();

	//how to filter?
	int counter = 0;
	std::string tag="AD";
	int numRef; int numAlt; int tot;
	int numHet;
	std::vector<int> tmp;
	double logHalf = log(0.5);
	double cumul;
	int i;
	long numFiltered = 0;
	long numPassedFilter = 0;

	//open log file for sites
	std::string filename = params->getParameterWithDefault<std::string>("infoOut", "ImbalanceFilter_siteInfo.txt");
	std::ofstream sitesInfoOut(filename.c_str());
	if(!sitesInfoOut)
		throw "Failed to open file '" + filename + " for writing!";
	sitesInfoOut << "chr\tpos\tnumHet\tnumRef\tnumAlt\taltFrac\tpVal\tstatus\n";

	while(vcfFile.next()){
		++counter;

		//count allelic depth across all samples
		numRef = 0; numAlt = 0;
		numHet = 0;
		for(int i=0; i < vcfFile.numSamples(); ++i){
			if(!vcfFile.sampleIsMissing(i) && vcfFile.sampleIsHeteroRefNonref(i)){
				fillContainerFromString(vcfFile.getSampleContentAt(tag, i), tmp, ',');
				numRef += tmp[0]; numAlt += tmp[1];
				++numHet;
			}
		}

		//write info
		sitesInfoOut << vcfFile.chr() << "\t" << vcfFile.position() << "\t" << numHet;

		//calculate p-value for sampling
		if(numHet > 0){
			sitesInfoOut << "\t" << numRef << "\t" << numAlt << "\t" << (double) numAlt / (double)(numAlt+numRef);
			cumul = 0.0;
			tot = numAlt + numRef;
			if(numAlt < numRef){
				for(i = 0; i <= numAlt; ++i)
					cumul += exp(randomGenerator->binomCoeffLn(tot, i) + logHalf*tot);
			} else {
				for(i = numAlt; i <= tot; ++i)
					cumul += exp(randomGenerator->binomCoeffLn(tot, i) + logHalf*tot);
			}

			//now filter
			if(cumul < pValFilter){
				++numFiltered;
				vcfFile.written = true;
				sitesInfoOut << "\t" << cumul << "\tfiltered\n";
			} else {
				vcfFile.writeLine();
				++numPassedFilter;
				sitesInfoOut << "\t" << cumul << "\tpassed\n";
			}
		} else {
			++numPassedFilter;
			sitesInfoOut << "\t-\t-\t-\t-\tpassed\n";
		}

		if(verbose && counter % 1000 == 0)	std::cerr << "    - read " << counter << " lines, " << numPassedFilter << " passed, " << numFiltered << " were filtered out." << std::endl;
	}

	//clean up
	sitesInfoOut.close();
}

*/


void VcfDiagnostics::fixIntAsFloat(){
	//open vcf file
	logfile->list("Fixing integers that are printed as floats:");

	//open output file
	std::string filename = outname + (std::string) "_fixed.vcf.gz";
	gz::ogzstream out(filename.c_str());
	if(!out) throw "Failed to open outputfile '" + filename + "'!";
	vcfFile.setOutStream(out);
	vcfFile.writeHeaderVCF_4_0();

	//tmp vars
	int counter = 0;
	std::vector<int> vec;

	//parse VCF

	logfile->startIndent("Parsing vcf file:");
	while(vcfFile.next()){
		++counter;

		//fix GP field
		std::string gp = vcfFile.fieldContentAsString("GP", 0);
		fillContainerFromString(gp, vec, ',');
		gp = std::to_string(vec[0]) + ',' + std::to_string(vec[1]) + ',' + std::to_string(vec[2]);
		vcfFile.updateField("GP", gp, 0);

		//write output
		vcfFile.writeLine();

		//report progress
		if(verbose && counter % 1000000 == 0)
			logfile->list("read " + toString(counter) + " lines.");
	}
	logfile->endIndent();
};

}; //end namespace
