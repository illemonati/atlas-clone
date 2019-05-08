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
	size = 11;
	missingIndex = 10;

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

void TGenotypeComparisonTable::add(Genotype g1, Genotype g2){
	++counts[g1][g2];
};

void TGenotypeComparisonTable::addFirstMissing(Genotype g2){
	++counts[missingIndex][g2];
};

void TGenotypeComparisonTable::addSecondMissing(Genotype g1){
	++counts[g1][missingIndex];
};

void TGenotypeComparisonTable::add(const char ind1_first, const char ind1_second, const char ind2_first, const char ind2_second){
	add(genoMap.getGenotype(ind1_first, ind1_second), genoMap.getGenotype(ind2_first, ind2_second));
};

void TGenotypeComparisonTable::addFirstMissing(const char ind2_first, const char ind2_second){
	addFirstMissing(genoMap.getGenotype(ind2_first, ind2_second));
};

void TGenotypeComparisonTable::addSecondMissing(const char ind1_first, const char ind1_second){
	addSecondMissing(genoMap.getGenotype(ind1_first, ind1_second));
};

void TGenotypeComparisonTable::add(TVcfFileSingleLine & vcfFile1, TVcfFileSingleLine & vcfFile2){
	add(genoMap.getGenotype(vcfFile1.getFirstAlleleOfSample(0), vcfFile1.getSecondAlleleOfSample(0)), genoMap.getGenotype(vcfFile2.getFirstAlleleOfSample(0), vcfFile2.getSecondAlleleOfSample(0)));

	vcfFile1.next();
	vcfFile2.next();
};

void TGenotypeComparisonTable::addFirstMissing(TVcfFileSingleLine & vcfFile2){
	addFirstMissing(genoMap.getGenotype(vcfFile2.getFirstAlleleOfSample(0), vcfFile2.getSecondAlleleOfSample(0)));
	vcfFile2.next();
};

void TGenotypeComparisonTable::addSecondMissing(TVcfFileSingleLine & vcfFile1){
	addSecondMissing(genoMap.getGenotype(vcfFile1.getFirstAlleleOfSample(0), vcfFile1.getSecondAlleleOfSample(0)));
	vcfFile1.next();
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------
TVCFCompare::TVCFCompare(TLog* Logfile){
	logfile = Logfile;
};

void TVCFCompare::openVCF(std::string & filename, TVcfFileSingleLine & vcfFile){
	//open vcf file
	if(filename.find(".gz") == std::string::npos){
		logfile->list("Reading vcf from file '" + filename + "'.");
		vcfFile.openStream(filename, false);
	} else {
		logfile->list("Reading vcf from gzipped file '" + filename + "'.");
		vcfFile.openStream(filename, true);
	}

	vcfFile.enablePositionParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();
	vcfFile.enableVariantParsing();
};

void TVCFCompare::compareVCFFiles(TParameters & parameters){
	//open vcf files
	logfile->startIndent("Open VCF files to compare:");
	std::vector<std::string> filenames;
	parameters.fillParameterIntoVector("vcf", filenames, ',');

	//currently only implemented for comparing two VCFs
	if(filenames.size() != 2)
		throw "VCF comparison requires two VCF filenames (not " + toString(filenames.size()) + ")!";

	//open VCF files
	TVcfFileSingleLine vcfFile1, vcfFile2;
	openVCF(filenames[0], vcfFile1);
	openVCF(filenames[1], vcfFile2);
	logfile->endIndent();

	//prepare count table
	TGenotypeComparisonTable counts;

	//now parse VCf files and compare calls
	logfile->startIndent("Parsing vcf file:");

	//read first for both VCF
	vcfFile1.next();
	vcfFile2.next();

	std::string curChr = vcfFile1.chr();

	while(!vcfFile1.eof || !vcfFile2.eof){
		//is one end of file?
		if(vcfFile1.eof){
			//first is missing
			counts.addFirstMissing(vcfFile2);

			//advance second
			vcfFile2.next();
		} else if(vcfFile2.eof){
			//second is missing
			counts.addSecondMissing(vcfFile1);

			//advance first
			vcfFile1.next();
		}

		//are we on the same chromosome?
		else if(vcfFile1.chr() != vcfFile2.chr()){
			if(vcfFile1.chr() < vcfFile2.chr()){
				//second is on next chromosome
				//hence, second is missing
				counts.addSecondMissing(vcfFile1);

				//advance first
				vcfFile1.next();
			} else {
				//first is on next chromosome
				//first is missing
				counts.addFirstMissing(vcfFile2);

				//advance second
				vcfFile2.next();
			}
		} else {
			//same chromosome
			if(vcfFile1.position() == vcfFile2.position()){
				if(vcfFile1.sampleIsMissing(0)){
					//do not add comparisons where both are missing
					if(!vcfFile2.sampleIsMissing(0)){
						counts.addFirstMissing(vcfFile2);
					} else {
						if(vcfFile2.sampleIsMissing(0)){
							counts.addSecondMissing(vcfFile1);
						} else {
							counts.add(vcfFile1, vcfFile2);
						}
					}
				}

				//advance both
				vcfFile1.next();
				vcfFile2.next();
			} else {
				if(vcfFile1.position() < vcfFile2.position()){
					//second is missing
					counts.addSecondMissing(vcfFile1);

					//advance first
					vcfFile1.next();
				} else {
					//first is missing
					counts.addFirstMissing(vcfFile2);

					//advance second
					vcfFile2.next();
				}
			}
		}
	}
	logfile->endIndent();

};
