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

	logfile->listFlush("Writing reads to BAM ...");
	BamTools::BamAlignment bamAlignment;

	//first alignment: soft clips on both sides
	bamAlignment.Name = "softClips_bothsides";
	bamAlignment.Length = readLength;
	bamAlignment.AddTag("RG", "Z", RGName);
	bamAlignment.MapQuality = 50;
	bamAlignment.RefID = 0;
	bamAlignment.Position = 0;
	bamAlignment.Qualities = std::string(bamAlignment.Length, (char) 63);
	bamAlignment.QueryBases = "CC" + std::string(bamAlignment.Length - 4, 'A') + "CC";
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 2));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 4));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 2));
	trueQueryBases.push_back(std::string(bamAlignment.Length - 4, 'A'));
	bamWriter.SaveAlignment(bamAlignment);

	//second alignment: soft clips on left side
	bamAlignment.Name = "softClips_left";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 2));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 2));
	bamAlignment.QueryBases = "CC" + std::string(bamAlignment.Length - 2, 'A');
	trueQueryBases.push_back(std::string(bamAlignment.Length - 2, 'A'));
	bamWriter.SaveAlignment(bamAlignment);

	//third alignment: soft clips on right side
	bamAlignment.Name = "softClips_right";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 2));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 2));
	bamAlignment.QueryBases = std::string(bamAlignment.Length - 2, 'A') + "CC";
	trueQueryBases.push_back(std::string(bamAlignment.Length - 2, 'A'));
	bamWriter.SaveAlignment(bamAlignment);

	//fourth alignment: only soft clips
	bamAlignment.Name = "softClips_all";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', bamAlignment.Length));
	bamAlignment.QueryBases = std::string(bamAlignment.Length , 'C');
	trueQueryBases.push_back("");
	bamWriter.SaveAlignment(bamAlignment);

	//fifth alignment: no soft clips
	bamAlignment.Name = "softClips_none";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'A'));
	bamWriter.SaveAlignment(bamAlignment);

	//close BAM file
	bamWriter.Close();
	logfile->done();

	//index BAM file
	logfile->listFlush("Creating index of BAM file '" + filenameTag + ".bam' ...");
	BamTools::BamReader reader;
	if(!reader.Open(filenameTag + ".bam"))
		throw "Failed to open BAM file '" + filenameTag + ".bam' for indexing!";
	reader.CreateIndex(BamTools::BamIndex::STANDARD);
	reader.Close();
	logfile->done();

	//done!
	logfile->endIndent();

}

bool TAtlasTest_removeSoftClips::checkNewBAM(){
	//BamFile stuff
	BamTools::BamReader bamReader;
 	BamTools::BamAlignment bamAlignment;

	//open BAM file
	logfile->list("Reading data from BAM file '" + bamFileName + "'.");
	if (!bamReader.Open(bamFileName))
		throw "Failed to open BAM file '" + bamFileName + "'!";
	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + filenameTag + ".bam'!";

	//read through BAM
	int counter = 0;
	while(bamReader.GetNextAlignment(bamAlignment)){
		if(!basicChecks(bamAlignment, counter))
			return false;
		if(bamAlignment.QueryBases != trueQueryBases.at(counter)){
			logfile->newLine();
			logfile->conclude("Read " + bamAlignment.Name + ", isRev = " + toString(bamAlignment.IsReverseStrand()) + ": query bases not same as true bases! Read " + bamAlignment.QueryBases + " but was expecting " + trueQueryBases[counter]);
			return false;
		}

		++counter;
	}
	if((unsigned) counter != trueQueryBases.size()){
		logfile->newLine();
		logfile->conclude("Incorrect number of alignments in produced BAM file");
		return false;
	}

	return true;
};

bool TAtlasTest_removeSoftClips::basicChecks(BamTools::BamAlignment bamAlignment, const int pairNumber){
	if(bamAlignment.QueryBases.size() != bamAlignment.Qualities.size()){
		logfile->newLine();
		logfile->conclude( "Read number " + toString(pairNumber) + ": query bases not same size as qualities!");
		return false;
	}

	return true;
}
