/*
 * TAtlasTestRemoveSoftClips.cpp
 *
 *  Created on: Jul 26, 2019
 *      Author: linkv
 */

#include "TAtlasTestRemoveSoftClips.h"

//------------------------------------------
//TAtlasTest_RemoveSoftClips
//------------------------------------------

TAtlasTest_removeSoftClips::TAtlasTest_removeSoftClips(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "removeSoftClips";

	//variables
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
};

bool TAtlasTest_removeSoftClips::run(){
	//1) create a bam and fasta file with known pileup results
	//----------------------------------------------
	writeBAM();

	//2a) Run ATLAS to create pileup with bamTools
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("method", "bamtools");

	if(!runTGenomeFromInputfile("removeSoftClippedBases"))
		return false;

	//2b) check if results are OK
	//--------------------------
	if(!checkNewBAM())
		return false;

	return true;
};

void TAtlasTest_removeSoftClips::writeBAM(){
	//create a bam file with known pileup results
	logfile->startIndent("Writing a test BAM file:");
	logfile->listFlush("Opening bam file '" + bamFileName + "' for writing ...");

	//prepare header
	std::string RGName = "RG_1";
	int chrLength = 1000000;
	int readLength = 100;

	BamTools::SamHeader header("");
	header.Version = "1.4";
	header.GroupOrder = "none";
	header.SortOrder = "coordinate";
	header.ReadGroups.Add(RGName + "\tPU:UNKNOWN\tLB:UNKNOWN\tSM:Sim1\tCN:UNKNOWN\tPL:ILLUMINA");
	header.Sequences.Add(BamTools::SamSequence("Chr1", chrLength));

	BamTools::RefVector references;
	references.push_back(BamTools::RefData("Chr1", chrLength));

	//now open file
	BamTools::BamWriter bamWriter;
	if (!bamWriter.Open(bamFileName, header, references))
		throw "Failed to open BAM file '" + bamFileName + "'!";
	logfile->done();

	//create alignment
	logfile->listFlush("Writing reads to BAM ...");
	BamTools::BamAlignment bamAlignment;
	bamAlignment.Length = readLength;
	bamAlignment.AddTag("RG", "Z", RGName);
	bamAlignment.MapQuality = 50;
	bamAlignment.Name = "*";
	bamAlignment.RefID = 0;
	bamAlignment.Position = 0;
	bamAlignment.Qualities = std::string(bamAlignment.Length, 30);
	bamAlignment.QueryBases = "CC" + std::string(bamAlignment.Length - 4, 'A') + "CC";
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 2));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 4));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 2));
}

bool TAtlasTest_removeSoftClips::checkNewBAM(){

}
