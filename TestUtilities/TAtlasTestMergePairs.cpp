/*
 * TAtlasTestMergePairs.cpp
 *
 *  Created on: Jan 16, 2019
 *      Author: linkv
 */

#include "../TestUtilities/TAtlasTestMergePairs.h"


TAtlasTest_mergePairs::TAtlasTest_mergePairs():TAtlasTest(){
	_name = "testMerging";
	filenameTag = _testingPrefix + _name;

	chrLength = -1;
	readLength = -1;
	phredError = -1.0;
	filterOrphanedReads = false;
}

void TAtlasTest_mergePairs::setVariables(TParameters & params, TLog* Logfile, TTaskList* TaskList){
	logfile = Logfile;
	taskList = TaskList;
	bamFileName = filenameTag + "_mergedReads.bam";
	readGroupName = "ReadGroup";
	readLength = params.getParameterWithDefault<int>("mergingTest_readLength", 100);
	chrLength = readLength * 5;
	phredError = params.getParameterWithDefault<int>("mergingTest_qual", 50);
//	filterPairsDiffChr = params.parameterExists("filterPairsDiffChr");
	filterOrphanedReads = !params.parameterExists("keepOrphans");
};

bool TAtlasTest_mergePairs::run(TParameters & params, TLog* Logfile, TTaskList* TaskList){
	//1) Define variables
	setVariables(params, Logfile, TaskList);

	//2) create a bam and fasta file with known pileup results
	//----------------------------------------------
	writeBAM();

	//3) Run ATLAS to create BAM
	//-----------------------------
	_testParams.addParameter("bam", filenameTag + ".bam");
//	if(filterPairsDiffChr){
//		_testParams.addParameter("filterPairsDiffChr", "");
//	}
	if(!filterOrphanedReads){
		_testParams.addParameter("keepOrphans", "");
	}
	_testParams.addParameter("mergingMethod", "keepHighestQualBase");

	_testParams.addParameter("keepOriginalQuality", "");

	if(!runMain("mergeReads"))
		return false;

	//4) check if results are OK
	//--------------------------
	return checkMergedBAMFile();
};

void TAtlasTest_mergePairs::writePairedEndReads(BamTools::BamWriter & bamWriter){
	//1) basic overlap, rev read is completely set to zero
	//1st mate
	logfile->listFlush("Writing reads to BAM ...");
	BamTools::BamAlignment bamAlignment;
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
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'C'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(50)));
	trueIsProper.push_back(true);

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
	trueQueryBases.push_back(bamAlignment.QueryBases);
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(0)));
	trueIsProper.push_back(true);

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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'A'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
	trueIsProper.push_back(true);

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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'A'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(50)));
	trueIsProper.push_back(true);

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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'T'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
	trueIsProper.push_back(true);

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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'A'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
	trueIsProper.push_back(true);
//	trueIgnoredReadMessages.push_back("OrderError: Reverse read of pair with name 3rd_pair_wrongOrder is ignored because its forward mate has not been read");


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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'A'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(50)));
	trueIsProper.push_back(true);
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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'C'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
	trueIsProper.push_back(true);

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
	if(!filterOrphanedReads){
		trueQueryBases.push_back(std::string(bamAlignment.Length, 'T'));
		trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
		trueIsProper.push_back(false);
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
	if(!filterOrphanedReads){
		trueQueryBases.push_back(std::string(bamAlignment.Length, 'G'));
		trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
		trueIsProper.push_back(false);
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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'A'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(40)));
	trueIsProper.push_back(true);

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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'A'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(40)));
	trueIsProper.push_back(true);

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
	if(!filterOrphanedReads){
		trueQueryBases.push_back(std::string(bamAlignment.Length, 'G'));
		trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
		trueIsProper.push_back(false);
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
	trueQueryBases.push_back(bamAlignment.QueryBases);
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(50)));
	trueIsProper.push_back(true);

	//single end read in between, short
	setToSingleEndEtc(bamAlignment);
	bamAlignment.EditTag("RG", "Z", readGroupName + "_single");
	bamAlignment.Position = 3762;
	bamAlignment.Length = 60;
	bamAlignment.InsertSize = bamAlignment.Length;
	bamAlignment.Name = "Single_end_short";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));

	bamWriter.SaveAlignment(bamAlignment);
	trueQueryBases.push_back(bamAlignment.QueryBases);
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
	trueIsProper.push_back(false);

	//single end read in between, should be ignored
	setToSingleEndEtc(bamAlignment);
	bamAlignment.EditTag("RG", "Z", readGroupName + "_ignored");
	bamAlignment.Position = 3763;
	bamAlignment.Length = 200;
	bamAlignment.InsertSize = bamAlignment.Length;
	bamAlignment.Name = "Single_end_ignored";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));

	//single end read in between, long
	setToSingleEndEtc(bamAlignment);
	bamAlignment.EditTag("RG", "Z", readGroupName + "_single");
	bamAlignment.Position = 3763;
	bamAlignment.Length = 102;
	bamAlignment.InsertSize = bamAlignment.Length;
	bamAlignment.Name = "Single_end_long";
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'A');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));

	bamWriter.SaveAlignment(bamAlignment);
	trueQueryBases.push_back(bamAlignment.QueryBases);
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));

	trueIsProper.push_back(false);

	//2nd mate
	setToProperPairEtc(bamAlignment);
	setToRevMate(bamAlignment);
	bamAlignment.EditTag("RG", "Z", readGroupName);
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
	//overlap with deletion + overlap, lower qual + insertion + overlap, lower Qual + end not overlapping
	trueQueryBases.push_back(std::string(5, 'A')
						+ std::string(5, 'A')
						+ std::string(5, 'C')
						+ std::string(40, 'A')
						+ std::string(30, 'A'));
	trueQualities.push_back(std::string(5, qualMap.phredIntToQuality(30))
						+ std::string(5, qualMap.phredIntToQuality(0))
						+ std::string(5, qualMap.phredIntToQuality(30))
						+ std::string(40, qualMap.phredIntToQuality(0))
						+ std::string(30, qualMap.phredIntToQuality(30)));
	trueIsProper.push_back(true);

	//--------------------------------------------------------
	//9) mate on different chr
	// first mate
	setToFwdMate(bamAlignment);
	bamAlignment.Position = 4000;
	bamAlignment.InsertSize = -100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.Name = "9th_pair_mateOnDiffChr";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'C');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(!filterOrphanedReads){
		trueQueryBases.push_back(std::string(bamAlignment.Length, 'C'));
		trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
		trueIsProper.push_back(false);
	} else {
		trueIgnoredReadMessages.push_back("Read 9th_pair_mateOnDiffChr, fwd : mate on different chromosome");
	}


	//mate on diff chr second mate
	setToRevMate(bamAlignment);
	bamAlignment.RefID = 1;
	bamAlignment.Position = 20;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.Name = "9th_pair_mateOnDiffChr";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	if(!filterOrphanedReads){
		trueQueryBases.push_back(std::string(bamAlignment.Length, 'G'));
		trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(30)));
		trueIsProper.push_back(false);
	} else {
		trueIgnoredReadMessages.push_back("Read 9th_pair_mateOnDiffChr, rev : not a proper pair (orphan)");
	}



	//--------------------------------------------------------


	//alignment that is not a proper pair
	//deletions and insertions
};

void TAtlasTest_mergePairs::writeAlignments(BamTools::BamWriter & bamWriter){
	writePairedEndReads(bamWriter);
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
	header.ReadGroups.Add(readGroupName + "_single" + "\tPU:UNKNOWN\tLB:UNKNOWN\tSM:Sim1\tCN:UNKNOWN\tPL:ILLUMINA");
	header.ReadGroups.Add(readGroupName + "_ignored" + "\tPU:UNKNOWN\tLB:UNKNOWN\tSM:Sim1\tCN:UNKNOWN\tPL:ILLUMINA");
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

	//write reads
	writeAlignments(bamWriter);

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
};

void TAtlasTest_mergePairs::setToProperPairEtc(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsProperPair(true);
	bamAlignment.SetIsPaired(true);
	bamAlignment.SetIsMapped(true);
	bamAlignment.SetIsPrimaryAlignment(true);
	bamAlignment.SetIsDuplicate(false);
	bamAlignment.SetIsFailedQC(false);
	bamAlignment.SetIsMateMapped(true);
}

void TAtlasTest_mergePairs::setToSingleEndEtc(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsProperPair(false);
	bamAlignment.SetIsPaired(false);
	bamAlignment.SetIsMapped(true);
	bamAlignment.SetIsPrimaryAlignment(true);
	bamAlignment.SetIsDuplicate(false);
	bamAlignment.SetIsFailedQC(false);
	bamAlignment.SetIsMateMapped(false);
	bamAlignment.InsertSize = 0;
}

void TAtlasTest_mergePairs::setToFwdMate(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsReverseStrand(false);
	bamAlignment.SetIsMateReverseStrand(true);
	bamAlignment.SetIsFirstMate(true);
	bamAlignment.SetIsSecondMate(false);
};

void TAtlasTest_mergePairs::setToRevMate(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsReverseStrand(true);
	bamAlignment.SetIsMateReverseStrand(false);
	bamAlignment.SetIsFirstMate(false);
	bamAlignment.SetIsSecondMate(true);
};

bool TAtlasTest_mergePairs::basicChecks(BamTools::BamAlignment & bamAlignment, const int pairNumber){
	if(bamAlignment.QueryBases.size() != bamAlignment.Qualities.size()){
		logfile->newLine();
		logfile->conclude( "Read number " + toString(pairNumber) + ": query bases not same size as qualities!");
		return false;
	}
	if(bamAlignment.Qualities.size() > bamAlignment.InsertSize){
		logfile->newLine();
		logfile->conclude( "Read number " + toString(pairNumber) + ": longer than insert size!");
		return false;
	}

	return true;
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
		} if(bamAlignment.Qualities != trueQualities.at(counter)){
			logfile->newLine();
			logfile->conclude("Read " + bamAlignment.Name + ", isRev = " + toString(bamAlignment.IsReverseStrand()) + ": qualities not same as true qualities! Read " + bamAlignment.Qualities + " but was expecting " + trueQualities[counter]);
//			std::cout << "true qualities " << trueQualities.at(counter) << std::endl;
			return false;
		} if(bamAlignment.IsProperPair() != trueIsProper.at(counter)){
			logfile->newLine();
			logfile->conclude("Read " + bamAlignment.Name + ", isRev = " + toString(bamAlignment.IsReverseStrand()) + ": proper pair flag is " + toString(bamAlignment.IsProperPair()) + " but was expecting " + toString(trueIsProper.at(counter)));
			return false;
		}

		++counter;
	}
	if((unsigned) counter != trueQualities.size()){
		logfile->newLine();
		logfile->conclude("Incorrect number of alignments in merged BAM file");
		return false;
	}

	if(trueIgnoredReadMessages.size() > 0){
		//check ignored reads file
		std::string ignoredReadsFile = filenameTag + "_ignoredReads.txt.gz";
		logfile->listFlush("Reading ignored reads from '" + ignoredReadsFile + "'...");
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
		logfile->write("done!");

		if((unsigned) lineNum != trueIgnoredReadMessages.size()){
			logfile->newLine();
			logfile->conclude("Incorrect number of alignments in merged BAM file");
			return false;
		}
	}

	return true;
};


//----------------------------------
// TAtlasTest_mergeSplitPairs
//----------------------------------

TAtlasTest_mergeSplitPairs::TAtlasTest_mergeSplitPairs():TAtlasTest_mergePairs(){
	_name = "testSplitMerging";
	filenameTag = _testingPrefix + _name;
};

void TAtlasTest_mergeSplitPairs::setVariables(TParameters & params, TLog* Logfile, TTaskList* TaskList){
	logfile = Logfile;
	bamFileName = filenameTag + "_mergedReads.bam";
	readGroupName = "ReadGroup";
	readLength = params.getParameterWithDefault<int>("mergingTest_readLength", 100);
	chrLength = readLength * 5;
	phredError = params.getParameterWithDefault<int>("mergingTest_qual", 50);
	filterOrphanedReads = !params.parameterExists("keepOrphans");
};

void TAtlasTest_mergeSplitPairs::writeRGSpecsFile(){
	std::ofstream out;
	std::string outName = _name + "_RGSpecs.txt";
	out.open(outName.c_str());
	if(!out) throw "Failed to open file '" + outName + "'!";
	out << readGroupName << "\tpaired"  << std::endl;
	out << readGroupName + "_single" << "\tsingle" << "\t" << 102 << std::endl;
	out.close();

	outName = _name + "_ignoreThese.txt";
	out.open(outName.c_str());
	if(!out) throw "Failed to open file '" + outName + "'!";
	out << readGroupName + "_ignored" << std::endl;
	out.close();
};

bool TAtlasTest_mergeSplitPairs::run(TParameters & params, TLog* Logfile, TTaskList* TaskList){
	//1) define variables
	setVariables(params, Logfile, TaskList);

	//TODO: add function to test header
	//2) create a bam and fasta file with known pileup results
	//----------------------------------------------
	writeBAM();
	writeRGSpecsFile();

	//3) Run ATLAS to create BAM
	//-----------------------------
	_testParams.addParameter("bam", filenameTag + ".bam");
	_testParams.addParameter("readGroupSettings", _name + "_RGSpecs.txt");

//	if(filterPairsDiffChr){
//		_testParams.addParameter("filterPairsDiffChr", "");
//	}
	if(!filterOrphanedReads){
		_testParams.addParameter("keepOrphans", "");
	}
	_testParams.addParameter("mergingMethod", "keepHighestQualBase");
	_testParams.addParameter("keepOriginalQuality", "");
	_testParams.addParameter("ignoreReadGroups", _name + "_ignoreThese.txt");

	if(!runMain("splitMerge"))
		return false;

	//4) check if results are OK
	//--------------------------
	return checkMergedBAMFile();
}


