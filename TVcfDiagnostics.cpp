/*
 * TVcfDiagnostics.cpp
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */
#include <algorithm>
#include <typeinfo>
#include <sstream>
#include "gzstream.h"
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
}

//open input stream
void VcfDiagnostics::openVCF(TVcfFile_base & vcfFile){
	//open vcf file
	std::string filename = params->getParameterString("vcf");
	if(filename.find(".gz") == std::string::npos){
		logfile->list("Reading vcf from file '" + filename + "'.");
		vcfFile.openStream(filename, false);
	} else {
		logfile->list("Reading vcf from gzipped file '" + filename + "'.");
		vcfFile.openStream(filename, true);
	}
}

void VcfDiagnostics::initializeRandomGenerator(){
	if(!randomGeneratorInitialized){
		randomGenerator=new TRandomGenerator();
		if(verbose) std::cerr << " - Randomgenerator initialized with seed " << randomGenerator->usedSeed << std::endl;
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

void VcfDiagnostics::vcfToBeagle(){
	//open vcf file
	openVCF(vcfFile);

	//enable parsers
	vcfFile.enablePositionParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableVariantParsing(); //has to come before sample parsing!!
	vcfFile.enableSampleParsing();
	vcfFile.enableInfoParsing();

	//open output files:
	std::string outname = params->getParameterString("out");
	if(verbose) std::cerr << " - Writing genotype likelihoods in beagle format to '" << outname << ".beagle.gz" << std::endl;
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

void VcfDiagnostics::assessAllelicImbalance(){
	//open vcf file
	if(verbose) std::cerr << " - assessing allelic imbalance:" << std::endl;

	openVCF(vcfFile);

	//enable parsers
	vcfFile.enablePositionParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();
	vcfFile.enableVariantParsing();

	//prepare output
	std::string outname = params->getParameterString("outname") + "_allelicDepth.txt";
	if(verbose) std::cerr << "    - Writing files to '" << outname  << std::endl;
	int maxDP = params->getParameterIntWithDefault("maxDepth", 1000);

	//limit input?
	long inputLines = params->getParameterLongWithDefault("inputLines", -1);
	if(inputLines <= 0){
		if(verbose)
			std::cerr << "    - Reading whole vcf." << std::endl;
	} else
		std::cerr << "    - Limiting input to " << inputLines << " lines." << std::endl;

	//initialize table
	long** table = new long*[maxDP];
	for(int i=0; i<maxDP; ++i){
		table[i] = new long[maxDP];
	}

	for(int i=0; i<maxDP; ++i){
		for(int j=0; j<maxDP; ++j)
			table[i][j] = 0;
	}

	//temp variables
	int counter = 0;
	int numHet = 0;

	//open log file for sites
	std::ofstream outTable(outname.c_str());
	if(!outTable)
		throw "Failed to open file '" + outname + " for writing!";

	while(vcfFile.next()){
		++counter;

		//count allelic depth across all samples
		for(int i=0; i < vcfFile.numSamples(); ++i){
			if(!vcfFile.sampleIsMissing(i) && vcfFile.sampleIsHeteroRefNonref(i)){
				++numHet;
				std::vector<std::string> tmp;
				std::string tag="AD";
				fillVectorFromString(vcfFile.getSampleContentAt(tag, i), tmp, ',');
				int numRef = stringToInt(tmp[0]);
				int numAlt = stringToInt(tmp[1]);
				if(vcfFile.depthAsIntNoCheckForMissingSample("DP", i) > maxDP){
					std::cerr << "WARNING: DP is " + toString(vcfFile.depthAsIntNoCheckForMissingSample("DP", i)) + " at pos " + toString(vcfFile.position()) + ". This site will be ignored.";
					continue;
				}
				++table[numRef][numAlt];
			}
		}
		if(verbose && counter % 1000000 == 0)
			std::cerr << "    - read " << counter << " lines, " << numHet << " heterozygous sites were found." << std::endl;
		if(inputLines != -1 && counter > inputLines){
			std::cerr << "    - reached input limit!" << std::endl;
			break;
		}
	}

	std::cerr << "    - Done reading vcf!" << std::endl;;

	//write table
	//header
	outTable << "ref_depth/alt_depth";
	for(int i=0; i<maxDP; ++i){
		outTable << "\talt_" << i;
	}
	outTable << "\n";

	for(int i=0; i<maxDP; ++i){
		outTable << "ref_" << i;
		for(int j=0; j<maxDP; ++j)
			outTable << "\t" << table[i][j];
		outTable << "\n";
	}

	//clean up
	outTable.close();

	for(int i = 0; i < maxDP; ++i)
	    delete[] table[i];
	delete[] table;
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


