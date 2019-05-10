/*
 * TAtlasTestFilter.cpp
 *
 *  Created on: May 8, 2019
 *      Author: vivian
 */

#include "TAtlasTestFilter.h"


TAtlasTest_filter::TAtlasTest_filter(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "testFilter";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + "_filtered.bam";
	readGroupName = "TestReadGroup";
	readLength = params.getParameterIntWithDefault("mergingTest_readLength", 100);
	chrLength = readLength * 5;
	phredError = params.getParameterIntWithDefault("mergingTest_qual", 50);


	keepImproperPairs = params.parameterExists("filter_keepImproperPairs");
	keepUnmappedReads = params.parameterExists("filter_keepUnmappedReads");
	keepFailedQC = params.parameterExists("filter_keepFailedQC");
	keepSecondary = params.parameterExists("filter_keepSecondary");
	keepSupplementary = params.parameterExists("filter_keepSupplementary");
	keepDuplicates = params.parameterExists("filter_keepDuplicates");
	filterSoftClips = params.parameterExists("filter_filterSoftClips");
	filterOrphanedReads = params.parameterExists("keepOrphans");

//	readGroups.readGroupInUse(curReadGroupID)
//					&& useStrand[bamAlignment.IsReverseStrand()]
//					&& useMate[bamAlignment.IsSecondMate()];
}

bool TAtlasTest_filter::run(){
	//1) create a bam and fasta file with known pileup results
	//----------------------------------------------
	writeBAM();

	//2) Run ATLAS to create filtered BAM
	//-----------------------------
	_testParams.addParameter("bam", filenameTag + ".bam");

	if(keepImproperPairs){
		_testParams.addParameter("keepImproperPairs", "");
	}
	if(keepUnmappedReads){
		_testParams.addParameter("keepUnmappedReads", "");
	}
	if(keepFailedQC){
		_testParams.addParameter("keepFailedQC", "");
	}
	if(keepSecondary){
		_testParams.addParameter("keepSecondary", "");
	}
//	if(keepSupplementary){
//		_testParams.addParameter("keepSupplementary", "");
//	}
	if(keepDuplicates){
		_testParams.addParameter("keepDuplicates", "");
	}
	if(filterSoftClips){
		_testParams.addParameter("filterSoftClips", "");
	}
	if(!filterOrphanedReads){
		_testParams.addParameter("keepOrphans", "");
	}

	if(!runTGenomeFromInputfile("filter"))
		return false;

	//3) check if results are OK
	//--------------------------
	return checkfilteredBAMFile();
}

bool TAtlasTest_filter::checkfilteredBAMFile(){
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

		if(bamAlignment.Name != shouldKeep.at(counter)){
			logfile->newLine();
			logfile->conclude("Read " + bamAlignment.Name +  " but was expecting " + shouldKeep[counter]);
			return false;
		}

		++counter;
	}
	if((unsigned) counter != shouldKeep.size()){
		logfile->newLine();
		logfile->conclude("Incorrect number of alignments in merged BAM file");
		return false;
	}

	if(trueIgnoredReadMessages.size() > 0){
		//check ignored reads file
		std::string ignoredReadsFile = filenameTag + "_ignoredReads.txt.gz";
		logfile->listFlush("Reading ignored reads from '" + ignoredReadsFile + "...");
		gz::igzstream file(ignoredReadsFile.c_str());
		if(!file) throw "Failed to open file '" + ignoredReadsFile + "!";

		size_t lineNum = 0;
		std::vector<std::string> vec;

		//fill list of reads to omit
		while(file.good() && !file.eof()){
			std::string line;
			if(getline(file, line)){
				if(lineNum >= trueIgnoredReadMessages.size()){
					logfile->newLine();
					logfile->conclude("Too many entries in ignored messages file. Expected " + toString(trueIgnoredReadMessages.size()) + " but reading line " + toString(lineNum + 1));
					return false;
				}
				if(line != trueIgnoredReadMessages.at(lineNum)){
					logfile->newLine();
					logfile->conclude("Incorrect entry in ignored reads file on line " + toString(lineNum) + ". Was expecting '" + trueIgnoredReadMessages.at(lineNum) + "'  but read '" + line + "'");
					return false;
				}
				++lineNum;
			}
		}
		logfile->write("done! Read " + toString(lineNum) + " read names");

		if((unsigned) lineNum != trueIgnoredReadMessages.size()){
			logfile->newLine();
			logfile->conclude("Incorrect number of alignments in merged BAM file");
			return false;
		}
	}

	return true;
}

void TAtlasTest_filter::setToProperPairEtc(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsProperPair(true);
	bamAlignment.SetIsPaired(true);
	bamAlignment.SetIsMapped(true);
	bamAlignment.SetIsPrimaryAlignment(true);
	bamAlignment.SetIsDuplicate(false);
	bamAlignment.SetIsFailedQC(false);
	bamAlignment.SetIsMateMapped(true);
}

void TAtlasTest_filter::setToFwdMate(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsReverseStrand(false);
	bamAlignment.SetIsMateReverseStrand(true);
	bamAlignment.SetIsFirstMate(true);
	bamAlignment.SetIsSecondMate(false);
}

void TAtlasTest_filter::setToRevMate(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsReverseStrand(true);
	bamAlignment.SetIsMateReverseStrand(false);
	bamAlignment.SetIsFirstMate(false);
	bamAlignment.SetIsSecondMate(true);
}

void TAtlasTest_filter::writeBAM(){
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

	//--------------------------------------------------------
	//create alignments
	//--------------------------------------------------------

	// duplicate
	logfile->listFlush("Writing reads to BAM ...");
	BamTools::BamAlignment bamAlignment;
	setToProperPairEtc(bamAlignment);
	bamAlignment.SetIsDuplicate(true);
	setToFwdMate(bamAlignment);
	bamAlignment.AddTag("RG", "Z", readGroupName);
	bamAlignment.MapQuality = 50;
	bamAlignment.RefID = 0;
	bamAlignment.Position = 558;
	bamAlignment.InsertSize = 64;
	bamAlignment.MatePosition = 559;
	bamAlignment.Length = 64;
	bamAlignment.Name = "duplicate";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'C');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(50));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepDuplicates)
		shouldKeep.push_back("duplicate");

	// improper pair
	setToProperPairEtc(bamAlignment);
	bamAlignment.SetIsProperPair(false);
	bamAlignment.Name = "improperPair";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepImproperPairs)
		shouldKeep.push_back("improperPair");

	// unmapped
	setToProperPairEtc(bamAlignment);
	bamAlignment.SetIsMapped(false);
	bamAlignment.Name = "unmapped";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepUnmappedReads)
		shouldKeep.push_back("unmapped");

	// failedQC
	setToProperPairEtc(bamAlignment);
	bamAlignment.SetIsMapped(false);
	bamAlignment.Name = "failedQC";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepFailedQC)
		shouldKeep.push_back("failedQC");

	// secondary
	setToProperPairEtc(bamAlignment);
	bamAlignment.SetIsPrimaryAlignment(false);
	bamAlignment.Name = "secondaryAlignment";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepSecondary)
		shouldKeep.push_back("secondaryAlignment");

//	// supplementary
//	setToProperPairEtc(bamAlignment);
//	bamAlignment.SetIs;
//	bamAlignment.Name = "supplementaryAlignment";
//
//	bamWriter.SaveAlignment(bamAlignment);
//	if(keepSupplementary)
//		shouldKeep.push_back("supplementaryAlignment");

	// soft clips left
	setToProperPairEtc(bamAlignment);
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 3));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 3));
	bamAlignment.Name = "softClipsLeft";
	bamWriter.SaveAlignment(bamAlignment);
	if(!filterSoftClips)
		shouldKeep.push_back("softClipsLeft");

	// soft clips right
	setToProperPairEtc(bamAlignment);
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 3));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 3));
	bamAlignment.Name = "softClipsRight";
	bamWriter.SaveAlignment(bamAlignment);
	if(!filterSoftClips)
		shouldKeep.push_back("softClipsRight");


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
