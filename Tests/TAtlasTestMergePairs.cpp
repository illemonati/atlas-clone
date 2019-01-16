/*
 * TAtlasTestMergePairs.cpp
 *
 *  Created on: Jan 16, 2019
 *      Author: linkv
 */

#include "TAtlasTestMergePairs.h"


TAtlasTest_mergePairs::TAtlasTest_mergePairs(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "pmdEmpiricSimulation";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + "_mergedReads.bam";
	readGroupName = "TestReadGroup";
	readLength = params.getParameterIntWithDefault("pileupTest_readLength", 100);
	chrLength = readLength * 5;
	phredError = params.getParameterIntWithDefault("pileupTest_qual", 50);
}


bool TAtlasTest_mergePairs::run(){
	//1) create a bam and fasta file with known pileup results
	//----------------------------------------------
	writeBAM();

	//2) Run ATLAS to create pileup
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("maxReadLength", toString(readLength));
	_testParams.addParameter("window", toString(2*readLength));

	if(!runTGenomeFromInputfile("mergeReads"))
//		return false;
		std::cout << "didnt work" << std::endl;

	//3) check if results are OK
	//--------------------------
	return checkMergedBAMFile();
}

void TAtlasTest_mergePairs::writeBAM(){
	//create a bam file with known merging results
	logfile->startIndent("Writing a test BAM file:");
	logfile->listFlush("Opening bam file '" + filenameTag + ".bam' for writing ...");

	//prepare header
	BamTools::SamHeader header("");
	header.Version = "1.4";
	header.GroupOrder = "none";
	header.SortOrder = "coordinate";
	header.ReadGroups.Add(readGroupName + "\tPU:UNKNOWN\tLB:UNKNOWN\tSM:Sim1\tCN:UNKNOWN\tPL:ILLUMINA");
	header.Sequences.Add(BamTools::SamSequence("Chr1", chrLength));
	header.Sequences.Add(BamTools::SamSequence("Chr2", chrLength));

	BamTools::RefVector references;
	references.push_back(BamTools::RefData("Chr1", chrLength));
	references.push_back(BamTools::RefData("Chr2", chrLength));

	//now open file
	BamTools::BamWriter bamWriter;
	if (!bamWriter.Open(filenameTag + ".bam", header, references))
		throw "Failed to open BAM file '" + filenameTag + ".bam" + "'!";
	logfile->done();

	//create alignments

	//1st mate
	logfile->listFlush("Writing reads to BAM ...");
	BamTools::BamAlignment bamAlignment;
	bamAlignment.AddTag("RG", "Z", readGroupName);
	bamAlignment.MapQuality = 50;
	bamAlignment.RefID = 0;
	bamAlignment.Position = 558;
	bamAlignment.InsertSize = 64;
	bamAlignment.MatePosition = 559;
	bamAlignment.Bin = 163;
	bamAlignment.Name = "ST-E00201:183:HFH3GALXX:1:2218:21349:71348_fwd";
	bamAlignment.QueryBases = "ACAGAGAAAGCAAGGGAATTCCAGAAAAATATCCACTTCTGCTTCACTGACTACACTTAAAGCC";
	bamAlignment.Qualities = "DAFEGAJGGHDIGHHHJGGGHHHHJGGGGGCGHHHGJGHJGDJGHHGJGJGJFGHGJGGGGHDH";
	bamAlignment.Length = bamAlignment.QueryBases.size();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);

	//2nd mate
	bamAlignment.Position = 559;
	bamAlignment.InsertSize = 64;
	bamAlignment.MatePosition = 558;
	bamAlignment.Bin = 83;
	bamAlignment.Name = "ST-E00201:183:HFH3GALXX:1:2218:21349:71348_rev";
	bamAlignment.QueryBases = "CAGAGAAAGCAAGGGAATTCCAGAAAAATATCCACTTCTGCTTCACTGACTACACTTAAAGCC";
	bamAlignment.Qualities = "HJHJHGGJDHGJHHHGGGJHGJHGGGGGGGJHGHHGJHIDHGJGHHIHHHGHGHHGGGDGA@=";
	bamAlignment.Length = bamAlignment.QueryBases.size();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);

	//close BAM file
	bamWriter.Close();
	logfile->done();

	//index BAM file
	logfile->listFlush("Creating index of BAM file '" + bamFileName + "' ...");
	BamTools::BamReader reader;
	if(!reader.Open(bamFileName))
		throw "Failed to open BAM file '" + bamFileName + "' for indexing!";

	reader.CreateIndex(BamTools::BamIndex::STANDARD);
	reader.Close();
	logfile->done();

	//done!
	logfile->endIndent();
};

bool TAtlasTest_mergePairs::checkMergedBAMFile(){
	//BamFile stuff
	BamTools::BamReader bamReader;
 	BamTools::BamAlignment bamAlignment;

	//open BAM file
	logfile->list("Reading data from BAM file '" + bamFileName + "'.");
	if (!bamReader.Open(bamFileName))
		throw "Failed to open BAM file '" + bamFileName + "'!";
	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + bamFileName + "'!";

	//read through BAM
	while (bamReader.GetNextAlignment(bamAlignment)){
		std::cout << bamAlignment.QueryBases << std::endl;
	}
	return false;
}
