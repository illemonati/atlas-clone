/*
 * TCompareVCF.cpp
 *
 *  Created on: May 7, 2019
 *      Author: wegmannd
 */

#include "TVcfCompare.h"

#include <stddef.h>
#include <sys/time.h>
#include <algorithm>
#include <exception>
#include <memory>

#include "TFile.h"
#include "stringFunctions.h"

namespace VCF{

using genometools::Base;
using genometools::Genotype;
using coretools::TLog;
using coretools::TParameters;
using coretools::str::toString;

//--------------------------------------------------------------
// TGenotypeComparisonTable
//--------------------------------------------------------------
TGenotypeComparisonTable::TGenotypeComparisonTable(){
	//set to zero
	for(int i=0; i<size; i++){
		for(int j=0; j<size; ++j){
			counts[i][j] = 0;
		}
	}
};

//add haploid genotypes
void TGenotypeComparisonTable::add(const Base b1, const Base b2){
	++counts[index(b1)][index(b2)];
};

void TGenotypeComparisonTable::addOtherMissing(const int sample, const Base b){
	if(sample == 0){
		++counts[index(b)][missingIndex];
	} else {
		++counts[missingIndex][index(b)];
	}
};


void TGenotypeComparisonTable::addFirstMissing(const Base b2){
	++counts[missingIndex][index(b2)];
};

void TGenotypeComparisonTable::addSecondMissing(const Base b1){
	++counts[index(b1)][missingIndex];
};

//add diploid genotypes
void TGenotypeComparisonTable::add(Genotype g1, Genotype g2){
	++counts[firstDiploidIndex + index(g1)][firstDiploidIndex + index(g2)];
};

void TGenotypeComparisonTable::addOtherMissing(const int sample, const Genotype g){
	if(sample == 0){
		++counts[firstDiploidIndex + index(g)][missingIndex];
	} else {
		++counts[missingIndex][firstDiploidIndex + index(g)];
	}
};

void TGenotypeComparisonTable::addFirstMissing(Genotype g2){
	++counts[missingIndex][firstDiploidIndex + index(g2)];
};

void TGenotypeComparisonTable::addSecondMissing(Genotype g1){
	++counts[firstDiploidIndex + index(g1)][missingIndex];
};


//add haploid / diploid combination of genotypes
void TGenotypeComparisonTable::add(const Genotype g1, const Base b2){
	++counts[firstDiploidIndex + index(g1)][index(b2)];
};

void TGenotypeComparisonTable::add(const Base b1, const Genotype g2){
	++counts[index(b1)][firstDiploidIndex + index(g2)];
};

//write
void TGenotypeComparisonTable::write(const std::string filename){
	//open output file
	coretools::TOutputFile out(filename);

	//write header
	std::vector<std::string> header = {"vcf1/vcf2"};
	//haploid bases
	for(auto b = Base::min; b < Base::max; ++b){
		header.push_back(toString(b));
	}

	//diploid genotypes
	for(auto g = Genotype::min; g < Genotype::max; ++g){
		header.push_back(toString(g));
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
// TVCFComapreVCF
//--------------------------------------------------------------
TVcfComapreVCF::TVcfComapreVCF(){
	sampleIndex = 0;
	vcfFile = nullptr;
	vcfFileOpen = false;
	minDepth = 0;
	minQual = 0.0;
};

TVcfComapreVCF::TVcfComapreVCF(std::string & filename, std::string & sampleName, TLog* logfile){
	//open vcf file
	if(filename.find(".gz") == std::string::npos){
		logfile->list("Reading sample '" + sampleName + "' from VCF file '" + filename + "'.");
		vcfFile = new genometools::TVcfFileSingleLine(filename, false);
	} else {
		logfile->list("Reading sample '" + sampleName + "' from gzipped VCF file '" + filename + "'.");
		vcfFile = new genometools::TVcfFileSingleLine(filename, true);
	}
	vcfFileOpen = true;

	vcfFile->enablePositionParsing();
	vcfFile->enableFormatParsing();
	vcfFile->enableSampleParsing();
	vcfFile->enableVariantParsing();

	//add sample index
	sampleIndex = vcfFile->sampleNumber(sampleName);

	//move to first position
	vcfFile->next();

	if(vcfFile->eof)
		throw "Vcf file '" + filename + "' is empty!";

	//store first chr
	parsedChromosomes.push_back(vcfFile->chr());

	//set filters to zero
	minDepth = 0;
	minQual = 0.0;
};

TVcfComapreVCF::TVcfComapreVCF(TVcfComapreVCF&& other){
	sampleIndex = other.sampleIndex;
	other.sampleIndex = 0;

	vcfFile = other.vcfFile;
	other.vcfFile = nullptr;

	vcfFileOpen = other.vcfFileOpen;
	other.vcfFileOpen = false;

	minDepth = other.minDepth;
	other.minDepth = 0;

	minQual = other.minQual;
	other.minQual = 0.0;

	parsedChromosomes = other.parsedChromosomes;
	other.parsedChromosomes.clear();
};

TVcfComapreVCF& TVcfComapreVCF::operator=(TVcfComapreVCF&& other){
	 if(this != &other){
		delete vcfFile;
		sampleIndex = other.sampleIndex;
		other.sampleIndex = 0;

		vcfFile = other.vcfFile;
		other.vcfFile = nullptr;

		vcfFileOpen = other.vcfFileOpen;
		other.vcfFileOpen = false;

		minDepth = other.minDepth;
		other.minDepth = 0;

		minQual = other.minQual;
		other.minQual = 0.0;

		parsedChromosomes = other.parsedChromosomes;
		other.parsedChromosomes.clear();
	 }

	 return *this;
};

TVcfComapreVCF::~TVcfComapreVCF(){
	if(vcfFileOpen)
		delete vcfFile;
};

void TVcfComapreVCF::next(){
	vcfFile->next();
	if(!vcfFile->eof){
		//store chr if new
		if(vcfFile->chr() != parsedChromosomes.back())
			parsedChromosomes.push_back(vcfFile->chr());

		//filter
		if(!vcfFile->sampleIsMissing(sampleIndex) && minDepth > 0){
			if(vcfFile->sampleDepth(sampleIndex) < minDepth){
				vcfFile->setSampleMissing(sampleIndex);
			}
		}

		if(!vcfFile->sampleIsMissing(sampleIndex) && minQual > 0){
			if(vcfFile->sampleGenotypeQuality(sampleIndex) < minQual){
				vcfFile->setSampleMissing(sampleIndex);
			}
		}
	}
};

void TVcfComapreVCF::setFilters(const int MinDepth, const double MinQual){
	minDepth = MinDepth;
	minQual = MinQual;
};

bool TVcfComapreVCF::chrParsed(const std::string chr){
	return std::find(parsedChromosomes.begin(), parsedChromosomes.end(), chr) != parsedChromosomes.end();
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------
TVcfCompare::TVcfCompare(TLog* Logfile){
	logfile = Logfile;
};

void TVcfCompare::addToOtherMissing(TGenotypeComparisonTable & counts, const int sample){
	if(!vcfFiles[sample].isMissing()){
		if(vcfFiles[sample].isDiploid()){
			counts.addOtherMissing(sample, vcfFiles[sample].genotype());
		} else {
			counts.addOtherMissing(sample, vcfFiles[sample].base());
		}
	}
};

void TVcfCompare::compareVCFFiles(TParameters & parameters){
	//open vcf files
	logfile->startIndent("Open VCF files to compare:");
	std::vector<std::string> fileNames;
	parameters.fillParameterIntoContainer("vcf", fileNames, ',');

	std::vector<std::string> sampleNames;
	parameters.fillParameterIntoContainer("samples", sampleNames, ',');

	//currently only implemented for comparing two VCFs
	if(fileNames.size() != 2)
		throw "VCF comparison requires two VCF file names (not " + toString(fileNames.size()) + ")!";
	if(fileNames.size() != 2)
		throw "VCF comparison requires two sample names (not " + toString(sampleNames.size()) + ")!";
	if(fileNames.size() != sampleNames.size())
		throw "The number of vcf file names " + toString(fileNames.size()) + " does not the match the number of provided sample names " + toString(sampleNames.size()) + "!";

	//open VCF files
	for(size_t i=0; i<fileNames.size(); i++){
		vcfFiles.emplace_back(fileNames[i], sampleNames[i], logfile);
	}
	logfile->endIndent();

	//are filters in place?
	int minDepth = parameters.getParameterWithDefault<int>("minDepth", 0);
	if(minDepth > 0){
		logfile->list("Will consider genotypes with depth < " + toString(minDepth) + " as missing.");
	}
	double minQual = parameters.getParameterWithDefault("minQual", 0.0);
	if(minQual > 0){
		logfile->list("Will consider genotypes with quality < " + toString(minQual) + " as missing.");
	}

	//limitLines
	bool limitLines = false;
	long lineLimit = -1;
	if(parameters.parameterExists("limitLines")){
		limitLines = true;
		logfile->list("Will stop reading after ", limitLines, " lines.");
		lineLimit = parameters.getParameter<int>("limitLines");
	}

	//set filters in VCF files
	for(TVcfComapreVCF& it : vcfFiles){
		it.setFilters(minDepth, minQual);
	}

	//prepare count table
	TGenotypeComparisonTable counts;

	//now parse VCf files and compare calls
	logfile->startIndent("Parsing vcf file:");
	uint32_t numLines = 0;
	struct timeval start;
	gettimeofday(&start, NULL);
	while((!vcfFiles[0].eof() || !vcfFiles[1].eof()) && (!limitLines  || numLines < lineLimit)){
		//is one end of file?
		if(vcfFiles[0].eof()){
			addToOtherMissing(counts, 1);
			vcfFiles[1].next();
		} else if(vcfFiles[1].eof()){
			addToOtherMissing(counts, 0);
			vcfFiles[0].next();
		}

		//are we on the same chromosome?
		else if(vcfFiles[0].chr() == vcfFiles[1].chr()){
			//same chromosome
			if(vcfFiles[0].position() == vcfFiles[1].position()){
				if(vcfFiles[0].isMissing()){
					//do not add comparisons where both are missing
					if(!vcfFiles[1].isMissing()){
						addToOtherMissing(counts, 1);
					}
				} else {
					if(vcfFiles[1].isMissing()){
						addToOtherMissing(counts, 0);
					} else {
						//both have calls
						if(vcfFiles[0].isDiploid()){
							if(vcfFiles[1].isDiploid()){
								counts.add(vcfFiles[0].genotype(), vcfFiles[1].genotype());
							} else {
								counts.add(vcfFiles[0].genotype(), vcfFiles[1].base());
							}
						} else {
							if(vcfFiles[1].isDiploid()){
								counts.add(vcfFiles[0].base(), vcfFiles[1].genotype());
							} else {
								counts.add(vcfFiles[0].base(), vcfFiles[1].base());
							}
						}
					}
				}

				//advance both
				vcfFiles[0].next();
				vcfFiles[1].next();
			} else {
				if(vcfFiles[0].position() < vcfFiles[1].position()){
					//position is missing in vcfFile1
					addToOtherMissing(counts, 0);
					vcfFiles[0].next();
				} else {
					//position is missing in vcfFile0
					addToOtherMissing(counts, 1);
					vcfFiles[1].next();
				}
			}
		} else {
			//we are on different chromosomes
			//has VCFFile1 already parsed chromosome of vcfFile0?
			if(vcfFiles[1].chrParsed(vcfFiles[0].chr())){
				//vcfFile1 is on next chromosome, hence, position is missing in vcfFile1
				addToOtherMissing(counts, 0);
				vcfFiles[0].next();
			} else if(vcfFiles[0].chrParsed(vcfFiles[1].chr())){
				//vcfFile0 is on next chromosome, hence, position is missing in vcfFile0
				addToOtherMissing(counts, 1);
				vcfFiles[1].next();
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
	std::string out = parameters.getParameter<std::string>("out", false);
	if(out.empty()){
		//guess from filename
		//get base name of first VCF file
		out = fileNames[0];
		out = coretools::str::extractBeforeLast(out, '.');
		if(fileNames[0].find(".gz") != std::string::npos){
			//if zipped there is extra .gz
			out = coretools::str::extractBeforeLast(out, ".");
		}

		//get base name of first VCF file
		std::string tmp = fileNames[1];
		tmp = coretools::str::extractBeforeLast(tmp, '.');
		if(fileNames[1].find(".gz") != std::string::npos){
			//if zipped there is extra .gz
			tmp = coretools::str::extractBeforeLast(tmp, ".");
		}

		out += "_" + tmp;
	}

	out += "_CallComparison.txt";

	//writing output file
	logfile->listFlush("Writing counts to file '" + out + "' ...");
	counts.write(out);
	logfile->done();
};

}; //end namespace
