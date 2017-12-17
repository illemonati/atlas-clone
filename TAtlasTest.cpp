/*
 * TAtlasTest.cpp
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */


#include "TAtlasTest.h"


//------------------------------------------
//TAtlasTest
//------------------------------------------
bool TAtlasTest::runTGenomeFromInputfile(){

}


//------------------------------------------
//TAtlasTest_pileup
//------------------------------------------
TAtlasTest_pileup::TAtlasTest_pileup(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "pileup";

	//variables
	readLength = params.getParameterIntWithDefault("pileupTest_readLength", 100);
	qual = params.getParameterIntWithDefault("pileupTest_qual", 50);
	params.fillParameterIntoVectorWithDefault("pileupTest_depths", depths, ',', "2,4,10,20,40");
	base = {'A', 'C', 'G', 'T'};
	chrLength = readLength * 5;
	filenameTag = "ATLAS_test_" + _name;
	bamFileName = filenameTag + ".bam";
	readGroupName = "TestReadGroup";
};

bool TAtlasTest_pileup::run(){
	//1) create a bam file with known pileup results
	//----------------------------------------------
	BamTools::BamWriter bamWriter;
	writeBAM(bamWriter);

	//2) Run ATLAS to create pileup
	//-----------------------------
	logfile->startIndent("Running task 'pileup':");
	_testParams.clear();
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("verbose", "");
	_testParams.addParameter("maxReadLength", toString(readLength));
	_testParams.addParameter("window", toString(2*readLength));

	atlasTaskSwitcher taskSwitcher(&_testParams, logfile);
	taskSwitcher.runTask("pileup");
	logfile->endIndent();

	//3) check if results are OK
	//--------------------------
	return checkPileupFile();
};

void TAtlasTest_pileup::writeBAM(BamTools::BamWriter & bamWriter){
	//create a bam file with known pileup results
	logfile->startIndent("Writing a test BAM file:");
	logfile->listFlush("Opening bam file '" + bamFileName + "' for writing ...");

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
	if (!bamWriter.Open(bamFileName, header, references))
		throw "Failed to open BAM file '" + bamFileName + "'!";
	logfile->done();

	//create alignment
	logfile->listFlush("Writing reads to BAM ...");
	BamTools::BamAlignment bamAlignment;
	bamAlignment.Length = readLength;
	bamAlignment.AddTag("RG", "Z", readGroupName);
	bamAlignment.MapQuality = 50;
	bamAlignment.Name = "*";
	bamAlignment.RefID = 0;
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));
	bamAlignment.Qualities = std::string(bamAlignment.Length, (char) qual + 33);
	bamAlignment.Position = 0;

	//homozygous
	int i;
	for(size_t d=0; d<depths.size(); ++d){
		bamAlignment.QueryBases = std::string(readLength, base[d % 4]);
		for(i=0; i<depths[d]; ++i)
			bamWriter.SaveAlignment(bamAlignment);
		bamAlignment.Position += readLength;
	}

	//heterozygous
	bamAlignment.RefID = 1;
	bamAlignment.Position = 0;
	int m;
	for(size_t d=0; d<depths.size(); ++d){
		//firts base
		bamAlignment.QueryBases = std::string(readLength, base[d % 4]);
		m = depths[d]/2;
		for(i=0; i<m; ++i)
			bamWriter.SaveAlignment(bamAlignment);

		//second base
		bamAlignment.QueryBases = std::string(readLength, base[(d+1) % 4]);
		for(i=m; i<depths[d]; ++i)
			bamWriter.SaveAlignment(bamAlignment);

		bamAlignment.Position += readLength;
	}

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

bool TAtlasTest_pileup::checkPileupFile(){
	logfile->startIndent("Checking pileup file:");

	//open pileup file
	std::string filename = filenameTag + "_pileup.txt";
	logfile->listFlush("Opening file '" + filename + "' for reading ...");
	std::ifstream in(filename.c_str());
	if(!in)
		throw "Failed to open file '" + filename + "'!";
	logfile->done();

	//skip header
	std::string tmp;
	getline(in, tmp);

	//some variables
	std::vector<std::string> line;
	unsigned int trueDepth;
	int numLines = 0;
	std::string chr = "Chr1";
	int truePos = 0;
	std::vector<size_t> baseCounts  = {0, 0, 0, 0};
	int firstBase;
	int depthFirstBase, depthSecondBase;

	//parse file line by line check contents
	logfile->listFlush("Parsing file ...");
	while(in.good() && !in.eof()){
		//read line into vector
		++numLines;
		++truePos;
		fillVectorFromLineWhiteSpaceSkipEmpty(in, line);

		//skip empty
		if(line.size() == 0) continue;

		//check columns
		if(line.size() != 14){
			logfile->conclude("Wrong number of columns in pileup file '" + filename + "' on line " + toString(numLines) + "!");
			return false;
		}


		//already on second chromosome?
		if(truePos > chrLength){
			truePos = 1;
			chr = "Chr2";
		}

		//check chr
		if(line[0] != chr){
			logfile->conclude("Wrong chromosome in pileup file '" + filename + "' on line " + toString(numLines) + "!");
			return false;
		}

		//check position
		if(stringToInt(line[1]) != truePos){
			logfile->conclude("Wrong position in pileup file '" + filename + "' on line " + toString(numLines) + ": expected " + toString(truePos) + ", found " + line[1] + "!");
			return false;
		}

		//check depth
		trueDepth = depths[(truePos-1) / readLength];
		if(stringToInt(line[2]) != trueDepth){
			logfile->conclude("Wrong depth in pileup file '" + filename + "' on line " + toString(numLines) + "!");
			return false;
		}

		//check bases
		baseCounts[0] = std::count(line[3].begin(), line[3].end(), 'A');
		baseCounts[1] = std::count(line[3].begin(), line[3].end(), 'C');
		baseCounts[2] = std::count(line[3].begin(), line[3].end(), 'G');
		baseCounts[3] = std::count(line[3].begin(), line[3].end(), 'T');

		firstBase = (int) ((truePos-1) / readLength) % 4;

		if(chr == "Chr1"){
			//is homozygous
			if(baseCounts[firstBase] != trueDepth){
				logfile->conclude("Wrong homozygous bases in pileup file '" + filename + "' on line " + toString(numLines) + "!");
				return false;
			}
		} else {
			//is heterozygous
			depthFirstBase = trueDepth / 2;
			depthSecondBase = trueDepth - depthFirstBase;
			if(baseCounts[firstBase] != depthFirstBase || baseCounts[(firstBase + 1) % 4] != depthSecondBase){
				logfile->conclude("Wrong heterozygous bases in pileup file '" + filename + "' on line " + toString(numLines) + "!");
				return false;
			}
		}
	}
	logfile->done();
	logfile->endIndent();

	return true;
}

