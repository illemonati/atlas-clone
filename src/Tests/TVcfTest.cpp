/*
 * TVcfTest.cpp
 *
 *  Created on: Nov 27, 2018
 *      Author: vivian
 */

#include "TVcfTest.h"

TAtlasTest_invariantBed::TAtlasTest_invariantBed(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "invariantBed";

	//variables
	filenameTag = _testingPrefix + _name;
	vcfFileName = filenameTag + ".vcf";
	bedFileName = filenameTag + ".bed.gz";
	logfile->list("Will write test vcf to " + vcfFileName);
};

void TAtlasTest_invariantBed::writeLineVcf(std::string chr, std::string pos, std::string gt){
	vcf << chr << "\t" << pos << "\t.\tA\tT\t.\t.\tDP=1\tGT:DP:AD:GQ:PL\t" << gt << ":1:1,0:3:0,3,40\t0/0:1:1,0:3:0,3,40\n";
}

void TAtlasTest_invariantBed::writeTestVcf(){
	vcf.open(vcfFileName.c_str());
	if(!vcf) throw "Failed to open output file '" + vcfFileName + "'!";

	vcf << "##fileformat=VCFv4.2\n";
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n";
	writeLineVcf("1", "1", "0/1");
	writeLineVcf("1", "4", "0/0");
	writeLineVcf("1", "5", "0/0");
	writeLineVcf("1", "7", "0/0");
	writeLineVcf("1", "9", "0/1");
	writeLineVcf("1", "11", "0/0");
	writeLineVcf("1", "41", "0/0");
	writeLineVcf("3", "3", "0/0");
	writeLineVcf("3", "12", "0/0");
	writeLineVcf("3", "13", "0/1");
	writeLineVcf("4", "1", "0/0");

	vcf.close();
}

bool TAtlasTest_invariantBed::checkBedFile(){
	logfile->list("reading conserved regions from '" + bedFileName + "'.");
	gz::igzstream bed(bedFileName.c_str());
	if(!bed)
		throw "Failed to open file '" + bedFileName + "!";

	std::string line;

//	while(bed.good() && !bed.eof()){
//		std::getline(bed, line);
//		std::cout << line << std::endl;
//	}

	//1st line
	if(!std::getline(bed, line))
		throw "file is already finished";
	if(line != "1\t3\t8"){
		logfile->newLine();
		logfile->conclude("line 1 is " + line + " instead of 1\t3\t9");
		return false;
	}

	//2nd line
	std::string truth = "1\t10\t41";
	getline(bed, line);
	if(line != (truth)){
		logfile->newLine();
		logfile->conclude("line 2 is " + line + " instead of " + truth);
		return false;
	}

	//3rd line
	truth = "3\t2\t12";
	getline(bed, line);
	if(line != (truth)){
		logfile->newLine();
		logfile->conclude("line 3 is " + line + " instead of " + truth);
		return false;
	}

	//4th line
	truth = "4\t0\t1";
	getline(bed, line);
	if(line != (truth)){
		logfile->newLine();
		logfile->conclude("line 3 is " + line + " instead of " + truth);
		return false;
	}

	else
		return true;
}

bool TAtlasTest_invariantBed::run(){
	//create vcf file
	writeTestVcf();

	//write invariant bed
	_testParams.addParameter("vcf", vcfFileName);
	if(!runTGenomeFromInputfile("VCFToInvariantBed"))
		return false;

	//check bed file
	if(checkBedFile())
		return true;
	else
		return false;
}
