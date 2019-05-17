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
// TVCFComapreVCF
//--------------------------------------------------------------
TVCFComapreVCF::TVCFComapreVCF(){
	sampleIndex = 0;
	vcfFile = nullptr;
	vcfFileOpen = false;
	minDepth = 0;
	minQual = 0.0;
};

TVCFComapreVCF::TVCFComapreVCF(std::string & filename, std::string & sampleName, TLog* logfile){
	//open vcf file
	if(filename.find(".gz") == std::string::npos){
		logfile->list("Reading sample '" + sampleName + "' from VCF file '" + filename + "'.");
		vcfFile = new TVcfFileSingleLine(filename, false);
	} else {
		logfile->list("Reading sample '" + sampleName + "' from gzipped VCF file '" + filename + "'.");
		vcfFile = new TVcfFileSingleLine(filename, true);
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

TVCFComapreVCF::TVCFComapreVCF(TVCFComapreVCF&& other){
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

TVCFComapreVCF& TVCFComapreVCF::operator=(TVCFComapreVCF&& other){
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

TVCFComapreVCF::~TVCFComapreVCF(){
	if(vcfFileOpen)
		delete vcfFile;
};

void TVCFComapreVCF::next(){
	vcfFile->next();
	if(!vcfFile->eof){
		//store chr if new
		if(vcfFile->chr() != parsedChromosomes.back())
			parsedChromosomes.push_back(vcfFile->chr());

		//filter
		if(!vcfFile->sampleIsMissing(sampleIndex)){
			if(minDepth > 0){
				if(vcfFile->sampleDepth(sampleIndex) < minDepth){
					vcfFile->setSampleMissing(sampleIndex);
				}
			}

			if(minQual > 0){
				if(vcfFile->sampleGenotypeQuality(sampleIndex) < minQual){
					vcfFile->setSampleMissing(sampleIndex);
				}
			}
		}
	}
};

void TVCFComapreVCF::setFilters(const int MinDepth, const double MinQual){
	minDepth = MinDepth;
	minQual = MinQual;
};

bool TVCFComapreVCF::chrParsed(const std::string chr){
	return std::find(parsedChromosomes.begin(), parsedChromosomes.end(), chr) != parsedChromosomes.end();
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------
TVCFCompare::TVCFCompare(TLog* Logfile){
	logfile = Logfile;
};

void TVCFCompare::addToOtherMissing(TGenotypeComparisonTable & counts, const int sample){
	if(!vcfFiles[sample].isMissing()){
		if(vcfFiles[sample].isDiploid()){
			counts.addOtherMissing(sample, vcfFiles[sample].genotype(genoMap));
		} else {
			counts.addOtherMissing(sample, vcfFiles[sample].base(genoMap));
		}
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
	if(fileNames.size() != sampleNames.size())
		throw "The number of vcf file names " + toString(fileNames.size()) + " does not the match the number of provided sample names " + toString(sampleNames.size()) + "!";

	//open VCF files
	for(size_t i=0; i<fileNames.size(); i++){
		vcfFiles.emplace_back(fileNames[i], sampleNames[i], logfile);
	}
	logfile->endIndent();

	//are filters in place?
	int minDepth = parameters.getParameterIntWithDefault("minDepth", 0);
	if(minDepth > 0){
		logfile->list("Will consider genotypes with depth < " + toString(minDepth) + " as missing.");
	}
	double minQual = parameters.getParameterDoubleWithDefault("minQual", 0.0);
	if(minQual > 0){
		logfile->list("Will consider genotypes with quality < " + toString(minQual) + " as missing.");
	}

	//limitLines
	bool limitLines = false;
	long lineLimit = -1;
	if(parameters.parameterExists("limitLines")){
		limitLines = true;
		logfile->list("Will stop reading after " + toString(limitLines) + " lines.");
		lineLimit = parameters.getParameterLong("limitLines");
	}

	//set filters in VCF files
	for(TVCFComapreVCF& it : vcfFiles){
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
								counts.add(vcfFiles[0].genotype(genoMap), vcfFiles[1].genotype(genoMap));
							} else {
								counts.add(vcfFiles[0].genotype(genoMap), vcfFiles[1].base(genoMap));
							}
						} else {
							if(vcfFiles[1].isDiploid()){
								counts.add(vcfFiles[0].base(genoMap), vcfFiles[1].genotype(genoMap));
							} else {
								counts.add(vcfFiles[0].base(genoMap), vcfFiles[1].base(genoMap));
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



























