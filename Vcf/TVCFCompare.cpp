/*
 * TCompareVCF.cpp
 *
 *  Created on: May 7, 2019
 *      Author: wegmannd
 */

#include "TVCFCompare.h"

//--------------------------------------------------------------
// TGenotypeComparisonTable
//--------------------------------------------------------------
TGenotypeComparisonTable::TGenotypeComparisonTable(){
	//set size
	size = 15; //4 haploid + 10 diploid + 1 missing
	missingIndex = 14;
	firstDiploidIndex = 4;

	//allocate count table
	counts = new int*[size];
	for(int i=0; i<size; i++){
		counts[i] = new int[size];

		//set to zero
		for(int j=0; j<size; ++j){
			counts[i][j] = 0;
		}
	}
};

TGenotypeComparisonTable::~TGenotypeComparisonTable(){
	for(int i=0; i<size; i++){
		delete[] counts[i];
	}
	delete[] counts;
};

//add haploid genotypes
void TGenotypeComparisonTable::add(const Base b1, const Base b2){
	++counts[b1][b2];
};

void TGenotypeComparisonTable::addOtherMissing(const int sample, const Base b){
	if(sample == 0){
		++counts[b][missingIndex];
	} else {
		++counts[missingIndex][b];
	}
};


void TGenotypeComparisonTable::addFirstMissing(const Base b2){
	++counts[missingIndex][b2];
};

void TGenotypeComparisonTable::addSecondMissing(const Base b1){
	++counts[b1][missingIndex];
};

//add diploid genotypes
void TGenotypeComparisonTable::add(Genotype g1, Genotype g2){
	++counts[firstDiploidIndex + g1][firstDiploidIndex + g2];
};

void TGenotypeComparisonTable::addOtherMissing(const int sample, const Genotype g){
	if(sample == 0){
		++counts[firstDiploidIndex + g][missingIndex];
	} else {
		++counts[missingIndex][firstDiploidIndex + g];
	}
};

void TGenotypeComparisonTable::addFirstMissing(Genotype g2){
	++counts[missingIndex][firstDiploidIndex + g2];
};

void TGenotypeComparisonTable::addSecondMissing(Genotype g1){
	++counts[firstDiploidIndex + g1][missingIndex];
};


//add haploid / diploid combination of genotypes
void TGenotypeComparisonTable::add(const Genotype g1, const Base b2){
	++counts[firstDiploidIndex + g1][b2];
};

void TGenotypeComparisonTable::add(const Base b1, const Genotype g2){
	++counts[b1][firstDiploidIndex + g2];
};

//write
void TGenotypeComparisonTable::write(const std::string filename, TGenotypeMap & genoMap){
	//open output file
	TOutputFilePlain out(filename);

	//write header
	std::vector<std::string> header = {"vcf1/vcf2"};
	//haploid bases
	for(int i=0; i<4; i++){
		header.push_back(std::string(1, genoMap.getBaseAsChar(i)));
	}

	//diploid genotypes
	for(int i=0; i<10; i++){
		header.push_back(genoMap.getGenotypeString(i));
	}

	//missing
	header.push_back("N/NN");

	//write header
	out.writeHeader(header);

	//write counts
	for(int i=0; i<size; i++){
		//write row name
		out << header[i+1];

		for(int j=0; j<size; j++){
			out << counts[i][j];
		}

		out.endLine();
	};
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------
TVCFCompare::TVCFCompare(TLog* Logfile){
	logfile = Logfile;
};

int TVCFCompare::openVCF(std::string & filename, TVcfFileSingleLine & vcfFile, std::string & sampleName){
	//open vcf file
	if(filename.find(".gz") == std::string::npos){
		logfile->list("Reading sample '" + sampleName + "' from VCF file '" + filename + "'.");
		vcfFile.openStream(filename, false);
	} else {
		logfile->list("Reading sample '" + sampleName + "' from gzipped VCF file '" + filename + "'.");
		vcfFile.openStream(filename, true);
	}

	vcfFile.enablePositionParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();
	vcfFile.enableVariantParsing();

	//return sample number
	return vcfFile.sampleNumber(sampleName);
};

void TVCFCompare::addToOtherMissing(TGenotypeComparisonTable & counts, const int sample, const int & sampleInVCF,  TVcfFileSingleLine & vcfFile){
	if(!vcfFile.sampleIsMissing(sampleInVCF)){
		if(vcfFile.sampleIsDiploid(sampleInVCF)){
			counts.addOtherMissing(sample, vcfFile.sampleGenotype(sampleInVCF, genoMap));
		} else {
			counts.addOtherMissing(sample, vcfFile.getFirstAlleleOfSample(sampleInVCF, genoMap));
		}
	}
};

void TVCFCompare::moveVcfFile(TVcfFileSingleLine & vcfFile, std::vector<std::string> & parsedChromosomesVcf){
	vcfFile.next();
	if(!vcfFile.eof && vcfFile.chr() != parsedChromosomesVcf.back()){
		parsedChromosomesVcf.push_back(vcfFile.chr());
	}
};

void TVCFCompare::compareVCFFiles(TParameters & parameters){
	//open vcf files
	logfile->startIndent("Open VCF files to compare:");
	std::vector<std::string> fileNames;
	parameters.fillParameterIntoVector("vcf", fileNames, ',');

	std::vector<std::string> sampleNames;
	parameters.fillParameterIntoVector("samples", sampleNames, ',');

	//currently only implemented for comparing two VCFs
	if(fileNames.size() != 2)
		throw "VCF comparison requires two VCF file names (not " + toString(fileNames.size()) + ")!";
	if(fileNames.size() != 2)
		throw "VCF comparison requires two sample names (not " + toString(sampleNames.size()) + ")!";

	//open VCF files
	TVcfFileSingleLine vcfFile0, vcfFile1;
	int sample0 = openVCF(fileNames[0], vcfFile0, sampleNames[0]);
	int sample1 = openVCF(fileNames[1], vcfFile1, sampleNames[1]);
	logfile->endIndent();

	//prepare count table
	TGenotypeComparisonTable counts;

	//now parse VCf files and compare calls
	logfile->startIndent("Parsing vcf file:");

	//read first for both VCF
	vcfFile0.next();
	vcfFile1.next();

	//save names of alreadyparsed chromosomes
	std::vector<std::string> parsedChromosomesVcf0, parsedChromosomesVcf1;
	parsedChromosomesVcf0.push_back(vcfFile0.chr());
	parsedChromosomesVcf1.push_back(vcfFile1.chr());

	//parse VCF
	uint32_t numLines = 0;
	struct timeval start;
	gettimeofday(&start, NULL);
	while(!vcfFile0.eof || !vcfFile1.eof){
		//is one end of file?
		if(vcfFile0.eof){
			addToOtherMissing(counts, 1, sample1, vcfFile1);
			moveVcfFile(vcfFile1, parsedChromosomesVcf1);
		} else if(vcfFile1.eof){
			addToOtherMissing(counts, 0, sample0, vcfFile0);
			moveVcfFile(vcfFile0, parsedChromosomesVcf0);
		}

		//are we on the same chromosome?
		else if(vcfFile0.chr() == vcfFile1.chr()){
			//same chromosome
			if(vcfFile0.position() == vcfFile1.position()){
				if(vcfFile0.sampleIsMissing(sample0)){
					//do not add comparisons where both are missing
					if(!vcfFile1.sampleIsMissing(sample1)){
						addToOtherMissing(counts, 1, sample1, vcfFile1);
					}
				} else {
					if(vcfFile1.sampleIsMissing(sample1)){
						addToOtherMissing(counts, 0, sample0, vcfFile0);
					} else {
						//both have calls
						if(vcfFile0.sampleIsDiploid(sample0)){
							if(vcfFile1.sampleIsDiploid(sample1)){
								counts.add(vcfFile0.sampleGenotype(sample0, genoMap), vcfFile1.sampleGenotype(sample1, genoMap));
							} else {
								counts.add(vcfFile0.sampleGenotype(sample0, genoMap), vcfFile1.getFirstAlleleOfSample(sample1, genoMap));
							}
						} else {
							if(vcfFile1.sampleIsDiploid(sample1)){
								counts.add(vcfFile0.getFirstAlleleOfSample(sample0, genoMap), vcfFile1.sampleGenotype(sample1, genoMap));
							} else {
								counts.add(vcfFile0.getFirstAlleleOfSample(sample0, genoMap), vcfFile1.getFirstAlleleOfSample(sample1, genoMap));
							}
						}
					}
				}

				//advance both
				moveVcfFile(vcfFile0, parsedChromosomesVcf0);
				moveVcfFile(vcfFile1, parsedChromosomesVcf1);
			} else {
				if(vcfFile0.position() < vcfFile1.position()){
					//position is missing in vcfFile1
					addToOtherMissing(counts, 0, sample0, vcfFile0);
					moveVcfFile(vcfFile0, parsedChromosomesVcf0);
				} else {
					//position is missing in vcfFile0
					addToOtherMissing(counts, 1, sample1, vcfFile1);
					moveVcfFile(vcfFile1, parsedChromosomesVcf1);
				}
			}
		} else {
			//we are on different chromosomes
			//has VCFFile1 already parsed chromosome of vcfFile0?
			if(std::find(parsedChromosomesVcf1.begin(), parsedChromosomesVcf1.end(), vcfFile0.chr()) != parsedChromosomesVcf1.end()){
				//vcfFile1 is on next chromosome, hence, position is missing in vcfFile1
				addToOtherMissing(counts, 0, sample0, vcfFile0);
				moveVcfFile(vcfFile0, parsedChromosomesVcf0);
			} else if(std::find(parsedChromosomesVcf0.begin(), parsedChromosomesVcf0.end(), vcfFile1.chr()) != parsedChromosomesVcf0.end()){
				//vcfFile0 is on next chromosome, hence, position is missing in vcfFile0
				addToOtherMissing(counts, 1, sample1, vcfFile1);
				moveVcfFile(vcfFile1, parsedChromosomesVcf1);
			} else {
				throw "Chromosomes differ between the two VCF files!";
			}
		}

		//report progress
		++numLines;
		if(numLines % 10000 == 0){
			struct timeval end;
			gettimeofday(&end, NULL);
			float runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(numLines) + " lines in " + toString(runtime) + " min.");
		}
	}
	struct timeval end;
	gettimeofday(&end, NULL);
	float runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(numLines) + " lines in " + toString(runtime) + " min.");
	logfile->list("Reached end of files.");
	logfile->endIndent();

	//write output file
	std::string out = parameters.getParameterString("out", false);
	if(out.empty()){
		//guess from filename
		//get base name of first VCF file
		out = fileNames[0];
		out = extractBeforeLast(out, '.');
		if(fileNames[0].find(".gz") != std::string::npos){
			//if zipped there is extra .gz
			out = extractBeforeLast(out, ".");
		}

		//get base name of first VCF file
		std::string tmp = fileNames[1];
		tmp = extractBeforeLast(tmp, '.');
		if(fileNames[1].find(".gz") != std::string::npos){
			//if zipped there is extra .gz
			tmp = extractBeforeLast(tmp, ".");
		}

		out += "_" + tmp;
	}

	out += "_CallComparison.txt";

	//writing output file
	logfile->listFlush("Writing counts to file '" + out + "' ...");
	counts.write(out, genoMap);
	logfile->done();
};



























