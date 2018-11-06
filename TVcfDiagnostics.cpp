/*
 * TVcfDiagnostics.cpp
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */
#include <algorithm>
#include <typeinfo>
#include <sstream>
#include "../gzstream.h"
#include "TVcfDiagnostics.h"

//---------------------------------------------------
//TAnnotator
//---------------------------------------------------
VcfDiagnostics::VcfDiagnostics(TParameters* Params, TLog* Logfile){
	params = Params;
	logfile = Logfile;
	verbose = params->parameterExists("verbose");
	randomGeneratorInitialized = false;
	randomGenerator = NULL;
	chr=-1;

	//open vcf file
	openVCF(vcfFile);

	//outputname
	outname = params->getParameterStringWithDefault("out", "");
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
void VcfDiagnostics::openVCF(TVcfFile_base & vcfFile){
	//open vcf file
	std::string filename = params->getParameterString("vcf");
	if(filename.find(".gz") == std::string::npos){
		logfile->list("Reading vcf from file '" + filename + "'.");
		vcfFile.openStream(filename, false);
		isZipped = false;
	} else {
		logfile->list("Reading vcf from gzipped file '" + filename + "'.");
		vcfFile.openStream(filename, true);
		isZipped = true;
	}
}

void VcfDiagnostics::initializeRandomGenerator(){
	if(!randomGeneratorInitialized){
		randomGenerator=new TRandomGenerator();
		logfile->list("Random generator initialized with seed " + toString(randomGenerator->usedSeed));
		randomGeneratorInitialized=true;
	}
}

int VcfDiagnostics::baseToNumber(char base, std::string & marker){
	if(base == 'A') return 0;
	else if(base == 'C') return 1;
	else if(base == 'G') return 2;
	else if (base == 'T') return 3;
	else throw "unknown base " + toString(base) + " at marker " + marker;
}

//----------------------
// countstable
//----------------------

void VcfDiagnostics::vcfToBeagle(){
	//enable parsers
	vcfFile.enablePositionParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableVariantParsing(); //has to come before sample parsing!!
	vcfFile.enableSampleParsing();
	vcfFile.enableInfoParsing();

	//open output files:
	gz::ogzstream beagleFile((outname + (std::string) ".beagle.gz").c_str());

	//other variables
	int vcfLines = 0;

	//header string
	std::string header = ""; header += "marker\tallele1\tallele2";
	for(int i=0; i<vcfFile.numSamples(); ++i){
		for(int r=0; r<3; ++r){
			header += "\t" + vcfFile.sampleName(i);
		}
	}
	header += "\n";

	beagleFile << header;

	while(vcfFile.next()){
		//report progress
		++vcfLines;
		std::cout << "vcfLines " << vcfLines << std::endl;
		if(verbose && fmod(vcfLines, 100000000) == 0)
			std::cerr << " - Progress: Lines read from vcfFile = " << vcfLines << std::endl;

		//write line
		std::string marker = ""; marker += vcfFile.chr();marker += "_"; marker += toString(vcfFile.position());
		beagleFile << marker;
		beagleFile << "\t" << baseToNumber(vcfFile.getRefAllele(), marker) << "\t" << baseToNumber(vcfFile.getFirstAltAllele(), marker);
		for(int i=0; i<vcfFile.numSamples(); ++i){
			if(vcfFile.sampleIsMissing(i)){
				beagleFile << "\t0.333\t0.333\t0.333";
			} else {
//				std::cout << "i " << i << " ll " << vcfFile.genotypeLikelihoods(i).AA << " phred " << vcfFile.genotypeLikelihoodsPhred(i).AA << std::endl;
				beagleFile << "\t" << vcfFile.genotypeLikelihoodsPhred(i).AA
					<< "\t" << vcfFile.genotypeLikelihoodsPhred(i).AB
					<< "\t" << vcfFile.genotypeLikelihoodsPhred(i).BB;
			}
		}
		beagleFile << "\n";
	}

	beagleFile.close();

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

void VcfDiagnostics::assessAllelicImbalance(){
	//enable parsers
	vcfFile.enablePositionParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();
	vcfFile.enableVariantParsing();

	//output
	logfile->list("Writing files to '" + outname + "_allelicDepth.txt'");

	//limit input?
	int maxDP = params->getParameterIntWithDefault("maxDepth", 100);
	logfile->list("Ignoring sites with depth larger than " + toString(maxDP) + ".");

	long inputLines = params->getParameterLongWithDefault("inputLines", -1);
	if(inputLines <= 0){
		logfile->list("Reading whole vcf.");
	} else
		logfile->list("Limiting input to " + toString(inputLines) + " lines.");

	//initialize tables
	logfile->startIndent("Initializing count tables:");
	std::string qualityString = params->getParameterStringWithDefault("qualities", "0,10,20,30,40,50");
	std::vector<int> qualities;
	fillVectorFromString(qualityString, qualities, ',');

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
				fillVectorFromString(vcfFile.getSampleContentAt(tag, i), tmp, ',');
				int numRef = stringToInt(tmp[0]);
				int numAlt = stringToInt(tmp[1]);
				if(numRef == 0 || numAlt == 0)
					throw "Call at position " + toString(vcfFile.position()) + " is heterozygous but reference or alternative allelic depth is 0!";
				if(vcfFile.depthAsIntNoCheckForMissingSample("DP", i) > maxDP){
					logfile->warning("DP is " + toString(vcfFile.depthAsIntNoCheckForMissingSample("DP", i)) + " at pos " + toString(vcfFile.position()) + ". This site will be ignored.");
					continue;
				}

				//add count to correct table
				int quality = stringToInt(vcfFile.getSampleContentAt("GQ", i));
				int index = findLastPassedFilterIndex(quality, qualities);
				if(index == 4) throw "stopped";
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

	double pValFilter = params->getParameterDoubleWithDefault("pValFilter", 0.001);

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
	std::string filename = params->getParameterStringWithDefault("infoOut", "ImbalanceFilter_siteInfo.txt");
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
				fillVectorFromString(vcfFile.getSampleContentAt(tag, i), tmp, ',');
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


