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

	keepReadsLongerThanFragment = params.parameterExists("filter_keepReadsLongerThanFragment");
	keepImproperPairs = params.parameterExists("filter_keepImproperPairs");
	keepUnmappedReads = params.parameterExists("filter_keepUnmappedReads");
	keepFailedQC = params.parameterExists("filter_keepFailedQC");
	keepSecondary = params.parameterExists("filter_keepSecondary");
	keepSupplementary = params.parameterExists("filter_keepSupplementary");
	keepDuplicates = params.parameterExists("filter_keepDuplicates");
	filterSoftClips = params.parameterExists("filter_filterSoftClips");
	keepOrphanedReads = params.parameterExists("filter_keepOrphans");
	minMQ = 0;
	if(params.parameterExists("filter_minMQ")){
		minMQ = params.getParameterInt("filter_minMQ");
	}

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

	if(keepReadsLongerThanFragment){
		_testParams.addParameter("keepReadsLongerThanFragment", "");
	}
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
	if(keepOrphanedReads){
		_testParams.addParameter("keepOrphans", "");
	}

	_testParams.addParameter("minMQ", toString(minMQ));

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

void TAtlasTest_filter::setToSingleEnd(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsProperPair(false);
	bamAlignment.SetIsPaired(false);
	bamAlignment.SetIsMapped(true);
	bamAlignment.SetIsPrimaryAlignment(true);
	bamAlignment.SetIsDuplicate(false);
	bamAlignment.SetIsFailedQC(false);
	bamAlignment.SetIsMateMapped(false);
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
	header.Sequences.Add(BamTools::SamSequence("Chr3", chrLength));


	BamTools::RefVector references;
	references.push_back(BamTools::RefData("Chr1", chrLength));
	references.push_back(BamTools::RefData("Chr2", chrLength));
	references.push_back(BamTools::RefData("Chr3", chrLength));

	//now open file
	BamTools::BamWriter bamWriter;
	if (!bamWriter.Open(filenameTag + ".bam", header, references))
		throw "Failed to open BAM file '" + filenameTag + ".bam" + "'!";
	logfile->done();

	//--------------------------------------------------------
	//create alignments
	//--------------------------------------------------------

	// duplicate fwd
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
	else
		trueIgnoredReadMessages.push_back("Read duplicate, fwd : did not pass parser filters");

	// duplicate
	setToProperPairEtc(bamAlignment);
	setToRevMate(bamAlignment);
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
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepOrphanedReads)
		shouldKeep.push_back("duplicate");
	else
		trueIgnoredReadMessages.push_back("Read duplicate, rev : not a proper pair (orphan)");

	//------------------------------------------------
	// improper pair
	setToProperPairEtc(bamAlignment);
	setToFwdMate(bamAlignment);
	bamAlignment.SetIsProperPair(false);
	bamAlignment.Name = "improperPair";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepImproperPairs)
		shouldKeep.push_back("improperPair");
	else
		trueIgnoredReadMessages.push_back("Read improperPair, fwd : did not pass parser filters");

	// unmapped
	setToSingleEnd(bamAlignment);
	bamAlignment.SetIsMapped(false);
	bamAlignment.Name = "unmapped";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepUnmappedReads)
		shouldKeep.push_back("unmapped");
	else
		trueIgnoredReadMessages.push_back("Read unmapped, fwd : did not pass parser filters");

	// failedQC
	setToSingleEnd(bamAlignment);
	bamAlignment.SetIsMapped(false);
	bamAlignment.Name = "failedQC";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepFailedQC)
		shouldKeep.push_back("failedQC");
	else
		trueIgnoredReadMessages.push_back("Read failedQC, fwd : did not pass parser filters");

	// secondary
	setToSingleEnd(bamAlignment);
	bamAlignment.SetIsPrimaryAlignment(false);
	bamAlignment.Name = "secondaryAlignment";

	bamWriter.SaveAlignment(bamAlignment);
	if(keepSecondary)
		shouldKeep.push_back("secondaryAlignment");
	else
		trueIgnoredReadMessages.push_back("Read secondaryAlignment, fwd : did not pass parser filters");


//	// supplementary
//	setToSingleEnd(bamAlignment);
//	bamAlignment.SetIs;
//	bamAlignment.Name = "supplementaryAlignment";
//
//	bamWriter.SaveAlignment(bamAlignment);
//	if(keepSupplementary)
//		shouldKeep.push_back("supplementaryAlignment");

	// soft clips left
	setToSingleEnd(bamAlignment);
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 3));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 3));
	bamAlignment.Name = "softClipsLeft";
	bamWriter.SaveAlignment(bamAlignment);
	if(!filterSoftClips){
		shouldKeep.push_back("softClipsLeft");
	} else {
		trueIgnoredReadMessages.push_back("Read softClipsLeft, fwd : did not pass parser filters");
	}

	// soft clips right
	setToSingleEnd(bamAlignment);
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length - 3));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('S', 3));
	bamAlignment.Name = "softClipsRight";
	bamWriter.SaveAlignment(bamAlignment);
	if(!filterSoftClips)
		shouldKeep.push_back("softClipsRight");
	else {
		trueIgnoredReadMessages.push_back("Read softClipsRight, fwd : did not pass parser filters");
	}

	//1) basic overlap, rev read is completely set to zero
	//1st mate
	setToProperPairEtc(bamAlignment);
	setToFwdMate(bamAlignment);
	bamAlignment.AddTag("RG", "Z", readGroupName);
	bamAlignment.MapQuality = 50;
	bamAlignment.RefID = 0;
	bamAlignment.Position = 558;
	bamAlignment.InsertSize = 64;
	bamAlignment.MatePosition = 559;
	bamAlignment.Length = 64;
	bamAlignment.Name = "1st_pair";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'C');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(50));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("1st_pair");

	//2nd mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 559;
	bamAlignment.InsertSize = -64;
	bamAlignment.MatePosition = 558;
	bamAlignment.Length = 63;
	bamAlignment.Name = "1st_pair";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("1st_pair");

	//--------------------------------------------------------

	//2) No overlap
	//1st mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 565;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 625;
	bamAlignment.Length = 20;
	bamAlignment.Name = "2nd_pair_noOverlap";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("2nd_pair_noOverlap");

	//2nd mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 625;
	bamAlignment.InsertSize = -100;
	bamAlignment.MatePosition = 565;
	bamAlignment.Length = 20;
	bamAlignment.Name = "2nd_pair_noOverlap";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(50));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("2nd_pair_noOverlap");

	//--------------------------------------------------------
	//4) Not consecutive
	//Not consecutive 1st mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 662;
	bamAlignment.InsertSize = 105;
	bamAlignment.MatePosition = 767;
	bamAlignment.Length = 20;
	bamAlignment.Name = "4th_pair_notConsecutive";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'T');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("4th_pair_notConsecutive");

	//3) Wrong order
	//Wrong order 1st mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 665;
	bamAlignment.InsertSize = -100;
	bamAlignment.MatePosition = 765;
	bamAlignment.Length = 20;
	bamAlignment.Name = "3rd_pair_wrongOrder";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("3rd_pair_wrongOrder");

	//Wrong order 2nd mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 765;
	bamAlignment.MatePosition = 665;
	bamAlignment.Length = 100;
	bamAlignment.Name = "3rd_pair_wrongOrder";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(50));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("3rd_pair_wrongOrder");

//	trueIgnoredReadMessages.push_back("Blacklist: Forward read of pair with name 3rd_pair_wrongOrder because it was in the blacklist");

	//Not consecutive 2nd mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 767;
	bamAlignment.InsertSize = -105;
	bamAlignment.MatePosition = 662;
	bamAlignment.Length = 20;
	bamAlignment.Name = "4th_pair_notConsecutive";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'C');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("4th_pair_notConsecutive");


	//--------------------------------------------------------

	// 5) longer than insert size
	// first mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 768;
	bamAlignment.InsertSize = 20;
	bamAlignment.MatePosition = 770;
	bamAlignment.Length = 100;
	bamAlignment.Name = "5th_pair_longerThanInsert";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepReadsLongerThanFragment)
		shouldKeep.push_back("5th_pair_longerThanInsert");
	else
		trueIgnoredReadMessages.push_back("Read 5th_pair_longerThanInsert, fwd : longer than insert size (TLEN)");

	// second mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 770;
	bamAlignment.InsertSize = -20;
	bamAlignment.MatePosition = 768;
	bamAlignment.Length = 10;
	bamAlignment.Name = "5th_pair_longerThanInsert";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'T');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepOrphanedReads || keepReadsLongerThanFragment){
		shouldKeep.push_back("5th_pair_longerThanInsert");
	} else {
		trueIgnoredReadMessages.push_back("Read 5th_pair_longerThanInsert, rev : not a proper pair (orphan)");
	}

	//--------------------------------------------------------

	//6) mate too far away
	// first mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 771;
	bamAlignment.InsertSize = 3000;
	bamAlignment.MatePosition = 3751;
	bamAlignment.Length = 20;
	bamAlignment.Name = "6th_pair_mateTooFarAway";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepOrphanedReads){
		shouldKeep.push_back("6th_pair_mateTooFarAway");
	} else {
		trueIgnoredReadMessages.push_back("Read 6th_pair_mateTooFarAway, fwd : orphaned read: mate is farther away than 2000 bp");
	}

	//7) second too far away
	// first mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 772;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 872;
	bamAlignment.Length = 20;
	bamAlignment.Name = "7th_pair_secondTooFarAway";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(40));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("7th_pair_secondTooFarAway");

	//second too far away second mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 872;
	bamAlignment.InsertSize = -100;
	bamAlignment.MatePosition = 772;
	bamAlignment.Length = 20;
	bamAlignment.Name = "7th_pair_secondTooFarAway";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(40));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("7th_pair_secondTooFarAway");


	// 6 Mate too far away second mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 3751;
	bamAlignment.InsertSize = -3000;
	bamAlignment.MatePosition = 771;
	bamAlignment.Length = 20;
	bamAlignment.Name = "6th_pair_mateTooFarAway";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepOrphanedReads){
		shouldKeep.push_back("6th_pair_mateTooFarAway");
	} else {
		trueIgnoredReadMessages.push_back("Read 6th_pair_mateTooFarAway, rev : not a proper pair (orphan)");
	}

	//--------------------------------------------------------
	// 8) //deletion in overlap
	setToFwdMate(bamAlignment);
	bamAlignment.AddTag("RG", "Z", readGroupName);
	bamAlignment.MapQuality = 50;
	bamAlignment.Position = 3752;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 3772;
	bamAlignment.Length = 70;
	bamAlignment.Name = "8th_pair_indels";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', 10));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('I', 5));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', 10));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('D', 5));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', 45));
	bamAlignment.QueryBases = std::string(10, 'A') + std::string(5, 'C') + std::string(10, 'A') + std::string(45, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(50));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("8th_pair_indels");

	//2nd mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 3772;
	bamAlignment.InsertSize = -100;
	bamAlignment.MatePosition = 3752;
	bamAlignment.Length = 80 + 5;
	bamAlignment.Name = "8th_pair_indels";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', 10));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('I', 5));
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', 70));

	bamAlignment.QueryBases = std::string(10, 'A') + std::string(5, 'C') + std::string(70, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));

	bamWriter.SaveAlignment(bamAlignment);
	shouldKeep.push_back("8th_pair_indels");

	//--------------------------------------------------------
	//9) mate on different chr
	// first mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 4000;
	bamAlignment.InsertSize = -100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.Name = "9th_pair_mateOnDiffChr_first";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'C');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepOrphanedReads){
		shouldKeep.push_back("9th_pair_mateOnDiffChr_first");
	} else {
		trueIgnoredReadMessages.push_back("Read 9th_pair_mateOnDiffChr_first, fwd : mate on different chromosome");
	}


	//mate on diff chr second mate
	setToRevMate(bamAlignment);
	bamAlignment.RefID = 1;
	bamAlignment.Position = 20;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.Name = "9th_pair_mateOnDiffChr_second";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(keepOrphanedReads){
		shouldKeep.push_back("9th_pair_mateOnDiffChr_second");
	} else {
		trueIgnoredReadMessages.push_back("Read 9th_pair_mateOnDiffChr_second, rev : mate on different chromosome");
	}

	//--------------------------------------------------------

	//nothing, just on new chr
	setToSingleEnd(bamAlignment);
	bamAlignment.RefID = 2;
	bamAlignment.Position = 20;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.Name = "normal_mapping_qual";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(minMQ <= bamAlignment.MapQuality){
		shouldKeep.push_back("normal_mapping_qual");
	} else {
		trueIgnoredReadMessages.push_back("Read normal_mapping_qual, fwd : did not pass parser filters");
	}


	//small mapping quality
	setToSingleEnd(bamAlignment);
	bamAlignment.RefID = 2;
	bamAlignment.Position = 20;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.MapQuality = 10;
	bamAlignment.Name = "low_mapping_qual";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(minMQ <= bamAlignment.MapQuality){
		shouldKeep.push_back("low_mapping_qual");
	} else {
		trueIgnoredReadMessages.push_back("Read low_mapping_qual, rev : did not pass parser filters");
	}


	//--------------------------------------------------------


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
