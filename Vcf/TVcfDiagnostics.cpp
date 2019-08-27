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
VcfDiagnostics::VcfDiagnostics(TParameters & Params, TLog* Logfile){
	logfile = Logfile;
	randomGeneratorInitialized = false;
	randomGenerator = NULL;
	chr=-1;

	//open vcf file
	openVCF(Params.getParameterString("vcf"), vcfFile);

	//outputname
	outname = Params.getParameterStringWithDefault("out", "");
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

int VcfDiagnostics::baseToNumber(char base, std::string & marker){
	if(base == 'A') return 0;
	else if(base == 'C') return 1;
	else if(base == 'G') return 2;
	else if (base == 'T') return 3;
	else throw "unknown base " + toString(base) + " at marker " + marker;
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
		if(bases.find(vcfFile.getRefAllele()) == bases.end()) {//!= 'A' && vcfFile.getRefAllele() != 'C' && vcfFile.getRefAllele() != 'G' && vcfFile.getRefAllele() != 'T'
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
		char allele = vcfFile.getFirstAlleleOfSample(0);
		for(int ind=0; ind<vcfFile.numSamples(); ++ind){
			++indCounter;
			if(vcfFile.getFirstAlleleOfSample(ind) != allele || vcfFile.getSecondAlleleOfSample(ind) != allele){
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
	int maxDP = Params.getParameterIntWithDefault("maxDepth", 100);
	logfile->list("Ignoring sites with depth larger than " + toString(maxDP) + ".");

	long inputLines = Params.getParameterLongWithDefault("inputLines", -1);
	if(inputLines <= 0){
		logfile->list("Reading whole vcf.");
	} else
		logfile->list("Limiting input to " + toString(inputLines) + " lines.");

	//initialize tables
	logfile->startIndent("Initializing count tables:");
	std::string qualityString = Params.getParameterStringWithDefault("qualities", "0,10,20,30,40,50");
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

/***************************************
 * 									   *
 * 			 Vcf to Beagle			   *
 * 									   *
 ***************************************/

void VcfDiagnostics::writeBeagleHeader(TOutputFileZipped & beagleFile){
	//header string
	std::vector <std::string> header {"marker", "allele1", "allele2"};
	for(int i=0; i<vcfFile.numSamples(); ++i){
		for(int r=0; r<3; ++r)
			header.push_back(vcfFile.sampleName(i));
	}
	beagleFile.writeHeader(header);
}

void VcfDiagnostics::printProgressFrequencyFiltering(const struct timeval & startTime, long lineCounter, long numAcceptedLoci){
	struct timeval end;
	gettimeofday(&end, nullptr);
	float runtime = (end.tv_sec  - startTime.tv_sec)/60.0;
	logfile->list("Parsing line " + toString(lineCounter) + ", retained " + toString(numAcceptedLoci) + " loci in " + toString(runtime) + " min");
}

void VcfDiagnostics::concludeFilters(TVcfFilters & vcfFilters, const struct timeval & startTime,
		long lineCounter, long notBialleleicCounter,
		long missingSNPCounter, long lowFreqSNPCounter,
		long lowVariantQualityCounter,
		long noPLCounter, long notOnChrCounter, long numAcceptedLoci){
	printProgressFrequencyFiltering(startTime, lineCounter, numAcceptedLoci);
	if(notBialleleicCounter > 0)
		logfile->conclude(toString(notBialleleicCounter) + " loci were not bi-allelic.");
	if(lowVariantQualityCounter > 0)
		logfile->conclude(toString(lowVariantQualityCounter) + " loci had variant quality < " + toString(vcfFilters.minVariantQuality) + ".");
	if(noPLCounter > 0)
		logfile->conclude(toString(noPLCounter) + " loci had no PL or GL field.");
	if(missingSNPCounter > 0)
		logfile->conclude(toString(missingSNPCounter) + " loci had < " + toString(vcfFilters.minNumSamplesWithData) + " samples with data.");
	if(notOnChrCounter > 0)
		logfile->conclude(toString(notOnChrCounter) + " loci were on other chromosomes than specified.");
	/*if(lowFreqSNPCounter > 0)
		logfile->conclude(toString(lowFreqSNPCounter) + " loci had MAF < " + toString(vcfFilters.freqFilter) + ".");*/
}

int VcfDiagnostics::filterOnDepth(int numSamples, double * data, TVcfFilters & vcfFilters){
	int numIndividualsWithData = 0;
	unsigned int index = 0;
	for(int s = 0; s < numSamples; ++s, index+=3){
		// depth filter: if a locus has < minDepth reads, flag locus as missing (set all genotype likelihoods = 0.333)
		if (vcfFile.sampleDepth(s) < vcfFilters.minDepth)
			vcfFile.setSampleMissing(s);
		else
			numIndividualsWithData++;
	}

	return numIndividualsWithData;
}

inline double VcfDiagnostics::phredToError(double phred){
	return pow(10.0, -phred/10.0);
};

bool VcfDiagnostics::readVcfAndWriteBeagle(int numSamples, TVcfFilters & vcfFilters,
		TOutputFileZipped & beagleFile,
		long & lineCounter, long & notBialleleicCounter,
		long & missingSNPCounter, long & lowFreqSNPCounter,
		long & lowVariantQualityCounter,
		long & noPLCounter, long & notOnChrCounter, long & numAcceptedLoci, struct timeval & startTime){

    double data[3*numSamples];

	//read next
	while(vcfFile.next()){ // new line in vcf-file (= new locus)
		++lineCounter;

		//print progress
		if(lineCounter % vcfFilters.progressFrequency == 0)
			printProgressFrequencyFiltering(startTime, lineCounter, numAcceptedLoci);

		// limit lines
		if(vcfFilters.limitLines > 0 && lineCounter > vcfFilters.limitLines){
			logfile->list("Reached limit of " + toString(vcfFilters.limitLines) + " lines.");
			return false;
		}

		// keep chromosomes
		if (!vcfFilters.chromosomesToKeep.empty() && // don't keep this chromosome
				std::find(vcfFilters.chromosomesToKeep.begin(), vcfFilters.chromosomesToKeep.end(), vcfFile.chr()) == vcfFilters.chromosomesToKeep.end()){
			notOnChrCounter ++;
			continue;
		}

        //skip sites with != 2 alleles
        if(vcfFile.getNumAlleles() != 2){
        	notBialleleicCounter++;
        	continue;
        }

        //skip sites with too low variant quality
        if(vcfFilters.minVariantQuality > 0 && (vcfFile.variantQualityIsMissing() || vcfFile.variantQuality() < vcfFilters.minVariantQuality)){
        	lowVariantQualityCounter++;
        	continue;
        }

		//check if GL or PL is given
        int numIndividualsWithData = 0;
		if(vcfFile.formatColExists("GL")){
			numIndividualsWithData = filterOnDepth(numSamples, data, vcfFilters);
			unsigned int index = 0;
			for(int s = 0; s < numSamples; ++s, index+=3){
				if (vcfFile.sampleIsMissing(s)){
					data[index] = 0.333; data[index + 1] = 0.333; data[index + 2] = 0.333;
				}
				else
					vcfFile.fillLog10GenotypeLikelihoods(s, data[index], data[index + 1], data[index + 2]);
			}
		} else if(vcfFile.formatColExists("PL")){
			numIndividualsWithData = filterOnDepth(numSamples, data, vcfFilters);
			uint8_t tmp[3];
			unsigned int index = 0;
			for(int s = 0; s < numSamples; ++s, index+=3){
				if (vcfFile.sampleIsMissing(s)){
					data[index] = 0.333; data[index + 1] = 0.333; data[index + 2] = 0.333;
				}
				else{
					vcfFile.fillPhredScore(s, tmp[0], tmp[1], tmp[2]);
					data[index] = phredToError(tmp[0]);
					data[index + 1] = phredToError(tmp[1]);
					data[index + 2] = phredToError(tmp[2]);
				}
			}
		} else {
			++noPLCounter;
			continue;
		}

		// missingness filter: if less than <minNumSamplesWithData> individuals per locus have are missing, remove locus
		if (numIndividualsWithData < vcfFilters.minNumSamplesWithData){
			missingSNPCounter++;
			continue;
		}

		/*//filter in MAF --> TODO: check with Dan if this is necessary
		if(vcfFilters.freqFilter > 0.0 || vcfFilters.estimateGenotypeFrequencies){
			genoFrequencies.estimate(data, glfConverter, vcfFilters.epsilonF);

			if(genoFrequencies.MAF < vcfFilters.freqFilter){
				lowFreqSNPCounter++;
				continue;
			}
		}*/

		//SNP is accepted!
		++numAcceptedLoci;
		writeLocusToBeagleFile(beagleFile, data, numSamples);
		return true;
    }

	//return false at end of file
	logfile->endIndent("Reached end of VCF file.");
	concludeFilters(vcfFilters, startTime,
			lineCounter, notBialleleicCounter,
			missingSNPCounter, lowFreqSNPCounter,
			lowVariantQualityCounter, noPLCounter, notOnChrCounter, numAcceptedLoci);

	return false;
}

void VcfDiagnostics::writeLocusToBeagleFile(TOutputFileZipped & beagleFile,
		double * data, int numSamples){
	//write line
	std::string marker = vcfFile.chr() + ":" + toString(vcfFile.position());
	beagleFile << marker; // marker
	beagleFile << baseToNumber(vcfFile.getRefAllele(), marker) << baseToNumber(vcfFile.getFirstAltAllele(), marker); // ref and alt allele

	for (int s = 0; s < numSamples * 3; s+=3){
		beagleFile << data[s] << data[s + 1] << data[s + 2];
	}

	beagleFile.endLine();
}

void VcfDiagnostics::vcfToBeagle(TParameters & Params){
	// read vcf filters
	TVcfFilters vcfFilters(Params, logfile);

	//open output files
	TOutputFileZipped beagleFile(outname + ".beagle.gz");
	writeBeagleHeader(beagleFile);

	// initialize counters (filtering)
	long lineCounter = 0;
	long notBialleleicCounter = 0;
	long missingSNPCounter = 0;
	long lowFreqSNPCounter = 0;
	long lowVariantQualityCounter = 0;
	long noPLCounter = 0;
	long notOnChrCounter = 0;
	long numAcceptedLoci = 0;

	//set time at beginning
	logfile->startIndent("Start parsing VCF-file...");
	struct timeval startTime;
	gettimeofday(&startTime, nullptr);

    // initialize storage (store locus, then write, then overwrite)
    int numSamples = vcfFile.numSamples();

	// parse vcf-file and write to beagle
	while(readVcfAndWriteBeagle(numSamples, vcfFilters, beagleFile,
			lineCounter, notBialleleicCounter,
			missingSNPCounter, lowFreqSNPCounter,
			lowVariantQualityCounter, noPLCounter, notOnChrCounter, numAcceptedLoci, startTime))
	{ continue; }


	// clean up
	beagleFile.close();
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
		fillVectorFromString(gp, vec, ',');
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

/***************************************
 * 									   *
 * 			 TVcfFiltering			   *
 * 									   *
 ***************************************/

TVcfFilters::TVcfFilters(TParameters & Params, TLog * logfile){
	//settings
	limitLines = 0;
	minDepth = 1;
	minNumSamplesWithData = 0;
	//freqFilter = 0.0;
	//epsilonF = 0.001; //F for EM algorithm to estimate allele frequencies
	minVariantQuality = 0;
	//estimateGenotypeFrequencies = false;

	initialize(Params, logfile);
}

void TVcfFilters::initialize(TParameters & Parameters, TLog * logfile){
	//read parsing parameters
	// do we limit the lines to read?
	limitLines = Parameters.getParameterLongWithDefault("limitLines", -1);
	if(limitLines > 0)
		logfile->list("Will limit analysis to the first " + toString(limitLines) + " lines of the VCF file.");

	// do we set a depth filter?
	minDepth = Parameters.getParameterIntWithDefault("minDepth", 1);
	if(minDepth < 1)
		throw "minDepth must be >= 1!";
	if(minDepth > 1)
		logfile->list("Will filter samples to a minimum depth of " + toString(minDepth) + ".");

	// do we set a missingness filter?
	minNumSamplesWithData = Parameters.getParameterIntWithDefault("minSamplesWithData", 1);
	if(minNumSamplesWithData < 1)
		throw "minNumSamplesWithData must be >= 1!";
	if(minNumSamplesWithData > 1)
		logfile->list("Will remove loci where less than " + toString(minNumSamplesWithData) + " samples have data.");

	// parameters to set a filter on the allele frequency?
	/*freqFilter = Parameters.getParameterDoubleWithDefault("minMAF", 0.0); // MAF = minor allele frequency
	if(freqFilter < 0.0 || freqFilter >= 0.5)
		throw "MAF filter must be within (0.0,0.5)!";
	if(freqFilter > 0.0){
		estimateGenotypeFrequencies = true;
		epsilonF = Parameters.getParameterDoubleWithDefault("epsF", 0.0000001);
		logfile->list("Will filter on an allele frequency of " + toString(freqFilter) + ".");
	} else {
		estimateGenotypeFrequencies = false;
	}*/

	//filter on variant quality?
	minVariantQuality = Parameters.getParameterIntWithDefault("minVariantQuality", 0);
	if(minVariantQuality < 0) throw "minVariantQuality must be >= 0!";
	if(minVariantQuality > 0){
		logfile->list("Will only keep sites with variant quality >= " + toString(minVariantQuality) + ".");
	}

	// filter for specific chromosomes?
	if(Parameters.parameterExists("keepChromosomes"))
		specifyChromosomesToKeep(logfile, Parameters);

	//set progress frequency
	progressFrequency = Parameters.getParameterIntWithDefault("reportFreq", 10000);
}


void TVcfFilters::specifyChromosomesToKeep(TLog* logfile, TParameters & Parameters){
	std::string argument = Parameters.getParameterString("keepChromosomes");
	if(stringContains(argument, ".txt")){ // specified as a file name
		logfile->startIndent("Reading chromosomes that should be kept from '" + argument + "'");
		std::ifstream keepChromosomesFile(argument.c_str());
		if(!keepChromosomesFile)
			throw "Failed to open file '" + argument + "'!";
		while(keepChromosomesFile.good() && !keepChromosomesFile.eof()){
			std::string line;
			std::getline(keepChromosomesFile, line);
			std::vector<std::string> vec;
			fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
			//skip empty lines
			if(vec.size() > 0)
				chromosomesToKeep.push_back(vec[0]);
		}
		keepChromosomesFile.close();
	}
	else { // specified as a vector on command line
		logfile->startIndent("Reading chromosomes from command line.");
		fillVectorFromString(Parameters.getParameterString("keepChromosomes"), chromosomesToKeep, ',');
	}

	// write to logfile
	logfile->startIndent("Will keep the following chromosomes in the beagle file:");
	for (std::vector<std::string>::iterator it = chromosomesToKeep.begin(); it < chromosomesToKeep.end(); it ++)
		logfile->list(*it);

	logfile->endIndent();
	logfile->endIndent();
}

