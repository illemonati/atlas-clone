/*
 * TCompareVCF.cpp
 *
 *  Created on: May 7, 2019
 *      Author: wegmannd
 */

#include "TVCFCompare.h"


TGenotypeComparisonTable::TGenotypeComparisonTable(){
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


	//now parse VCf files and compare calls
	logfile->startIndent("Parsing vcf file:");

	//read first for both VCF
	vcfFile1.next();
	vcfFile2.next();

	while(!vcfFile1.eof || !vcfFile2.eof){

		//if pos is same

	}
	logfile->endIndent();

};
