/*
 * TAtlasTestMergePairs.cpp
 *
 *  Created on: Jan 16, 2019
 *      Author: linkv
 */

#include "TAtlasTestMergePairs.h"


TAtlasTest_mergePairs::TAtlasTest_mergePairs(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "testMerging";
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
	_testParams.addParameter("bam", filenameTag + ".bam");
	_testParams.addParameter("maxReadLength", toString(readLength));
	_testParams.addParameter("window", toString(2*readLength));

	if(!runTGenomeFromInputfile("mergeReads"))
		return false;

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

	//--------------------------------------------------------
	//create alignments
	//--------------------------------------------------------

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
	trueQueryBases.push_back(std::string(bamAlignment.Length, 'N'));
	trueQualities.push_back(std::string(bamAlignment.Length, qualMap.phredIntToQuality(1)));

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
	//TODO: parser should write out message to ignoredReads file stating that read was too long

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
	trueIgnoredReadMessages.push_back("Blacklist: Reverse read of pair with name 5th_pair_longerThanInsert because it was in the blacklist\n");

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


	// Mate too far away second mate
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
	trueIgnoredReadMessages.push_back("DistanceError: Read with name 6th_pair_mateTooFarAway has a mate that is farther away than 2000 bp\n");

	//--------------------------------------------------------
	//8) mate on different chr
	// first mate
	setToRevMate(bamAlignment);
	bamAlignment.Position = 3751;
	bamAlignment.InsertSize = -100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.Name = "8th_pair_mateOnDiffChr";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);

	//mate on diff chr second mate
	setToRevMate(bamAlignment);
	bamAlignment.RefID = 1;
	bamAlignment.Position = 20;
	bamAlignment.InsertSize = 100;
	bamAlignment.MatePosition = 20;
	bamAlignment.Length = 20;
	bamAlignment.Name = "8th_pair_mateOnDiffChr";
	bamAlignment.QueryBases = std::string(bamAlignment.Length, 'G');
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(30));
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));

	bamWriter.SaveAlignment(bamAlignment);
	trueIgnoredReadMessages.push_back("DistanceError: Read with name 6th_pair_mateTooFarAway has a mate that is farther away than 2000 bp\n");

	//--------------------------------------------------------


	//alignment that is on other chromosome, with one mate that is on new chr

	//alignment that is not a proper pair


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

void TAtlasTest_mergePairs::setToFwdMate(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsReverseStrand(false);
	bamAlignment.SetIsMateReverseStrand(true);
	bamAlignment.SetIsFirstMate(true);
	bamAlignment.SetIsSecondMate(false);
}

void TAtlasTest_mergePairs::setToRevMate(BamTools::BamAlignment & bamAlignment){
	bamAlignment.SetIsReverseStrand(true);
	bamAlignment.SetIsMateReverseStrand(false);
	bamAlignment.SetIsFirstMate(false);
	bamAlignment.SetIsSecondMate(true);
}

bool TAtlasTest_mergePairs::basicChecks(BamTools::BamAlignment & bamAlignment, const int pairNumber){
	if(bamAlignment.QueryBases.size() != bamAlignment.Qualities.size()){
		logfile->newLine();
		logfile->conclude( "Read number " + toString(pairNumber) + ": query bases not same size as qualities!");
		return false;
	} if(bamAlignment.Qualities.size() > bamAlignment.InsertSize){
		logfile->newLine();
		logfile->conclude( "Read number " + toString(pairNumber) + ": longer than insert size!");
		return false;
	}

	return true;
}
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
			logfile->conclude("Read number " + toString(counter) + ": query bases not same as true bases!");
			return false;
		} if(bamAlignment.Qualities != trueQualities.at(counter)){
			logfile->newLine();
			logfile->conclude("Read number " + toString(counter) + ": qualities not same as true qualities!");
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
		logfile->listFlush("Reading ignored reads from '" + ignoredReadsFile + "...");
		gz::igzstream file(ignoredReadsFile.c_str());
		if(!file) throw "Failed to open file '" + ignoredReadsFile + "!";

		int lineNum = 0;
		std::vector<std::string> vec;

		//fill list of reads to omit
		while(file.good() && !file.eof()){
			std::string line;
			if(getline(file, line)){
				if(line != trueIgnoredReadMessages.at(lineNum)){
					logfile->newLine();
					logfile->conclude("Incorrect entry in ignored reads file on line " + toString(lineNum) + ". Was expecting '" + trueIgnoredReadMessages.at(lineNum) + "'  but read '" + line + "'");
	//				return false;
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
