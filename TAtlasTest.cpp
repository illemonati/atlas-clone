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
TAtlasTest::TAtlasTest(TParameters & params, TLog* Logfile){
	logfile = Logfile;
	_name = "empty";
	_testingPrefix = params.getParameterStringWithDefault("prefix", "ATLAS_testing_");
	logfile->list("Will use prefix '" + _testingPrefix + "' for all testing files.");
};

bool TAtlasTest::runTGenomeFromInputfile(std::string task){
	logfile->startIndent("Running task '" + task + "':");
	_testParams.addParameter("task", task);
	_testParams.addParameter("verbose", "");

	//open task switcher and run task
	atlasTaskSwitcher taskSwitcher(&_testParams, logfile);
	bool returnVal = true;
	try{
		taskSwitcher.runTask();
	}
	catch (std::string & error){
		logfile->conclude(error);
		returnVal = false;
	}
	catch (const char* error){
		logfile->conclude(error);
		returnVal = false;
	}
	catch(std::exception & error){
		logfile->conclude(error.what());
		returnVal = false;
	}
	catch (...){
		logfile->conclude("unhandeled error!");
		returnVal = false;
	}
	logfile->endIndent();
	return returnVal;
};

//------------------------------------------
//TAtlasTest_pileup
//------------------------------------------
TAtlasTest_pileup::TAtlasTest_pileup(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "pileup";

	//variables
	readLength = params.getParameterIntWithDefault("pileupTest_readLength", 100);
	logfile->list("Will simulate reads of length " + toString(readLength) + ".");
	phredError = params.getParameterIntWithDefault("pileupTest_qual", 50);
	logfile->list("Will test with quality " + toString(phredError) + ".");
	params.fillParameterIntoVectorWithDefault("pileupTest_depths", depths, ',', "2,4,10,20,40");
	logfile->list("Will test the following depths: " + concatenateString(depths, ", ") + ".");
	chrLength = readLength * 5;
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
	readGroupName = "TestReadGroup";

	emissionTolerance = params.getParameterDoubleWithDefault("pileupTest_qual", 0.0001);
	logfile->list("Will allow for a relative error in emission probabilities up to " + toString(emissionTolerance) + ".");
};

bool TAtlasTest_pileup::run(){
	//1) create a bam file with known pileup results
	//----------------------------------------------
	writeBAM();

	//2) Run ATLAS to create pileup
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("maxReadLength", toString(readLength));
	_testParams.addParameter("window", toString(2*readLength));

	if(!runTGenomeFromInputfile("pileup"))
		return false;

	//3) check if results are OK
	//--------------------------
	return checkPileupFile();
};

void TAtlasTest_pileup::writeBAM(){
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
	BamTools::BamWriter bamWriter;
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
	bamAlignment.Qualities = std::string(bamAlignment.Length, qualMap.phredIntToQuality(phredError));
	bamAlignment.Position = 0;

	//homozygous
	int i;
	for(size_t d=0; d<depths.size(); ++d){
		bamAlignment.QueryBases = std::string(readLength, genoMap.getBaseAsChar(d % 4));
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
		bamAlignment.QueryBases = std::string(readLength, genoMap.getBaseAsChar(d % 4));
		m = depths[d]/2;
		for(i=0; i<m; ++i)
			bamWriter.SaveAlignment(bamAlignment);

		//second base
		bamAlignment.QueryBases = std::string(readLength, genoMap.getBaseAsChar((d+1) % 4));
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
	int firstBase, secondBase;
	int depthFirstBase, depthSecondBase;
	int b;
	double error = qualMap.phredIntToError(phredError);
	double relDiff;

	std::vector<double> emissionProbs(genoMap.numGenotypes, 0.0);

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
			logfile->newLine();
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
			logfile->newLine();
			logfile->conclude("Wrong chromosome in pileup file '" + filename + "' on line " + toString(numLines) + "!");
			return false;
		}

		//check position
		if(stringToInt(line[1]) != truePos){
			logfile->newLine();
			logfile->conclude("Wrong position in pileup file '" + filename + "' on line " + toString(numLines) + ": expected " + toString(truePos) + ", found " + line[1] + "!");
			return false;
		}

		//check depth
		trueDepth = depths[(truePos-1) / readLength];
		if(stringToInt(line[2]) != trueDepth){
			logfile->newLine();
			logfile->conclude("Wrong depth in pileup file '" + filename + "' on line " + toString(numLines) + "!");
			return false;
		}

		//check bases and emission probabilities
		for(b=0; b<4; ++b)
			baseCounts[b] = std::count(line[3].begin(), line[3].end(), genoMap.getBaseAsChar(b));

		firstBase = (int) ((truePos-1) / readLength) % 4;

		if(chr == "Chr1"){
			//homozgous case
			//bases
			if(baseCounts[firstBase] != trueDepth){
				logfile->newLine();
				logfile->conclude("Wrong homozygous bases in pileup file '" + filename + "' on line " + toString(numLines) + "!");
				return false;
			}

			//calc emission probs
			//set all to full error
			for(b=0; b<genoMap.numGenotypes; ++b)
				emissionProbs[b] = pow(error/3.0, trueDepth);

			//correct homozygous genotype
			emissionProbs[genoMap.getGenotype(firstBase,firstBase)] = pow(1.0-error, trueDepth);

			//all heterozygous that contain the correct base
			for(b=1; b<4; ++b)
				emissionProbs[genoMap.getGenotype(firstBase,(firstBase + b) % 4)] = pow(0.5 - error/3.0, trueDepth);
		} else {
			//heterozygous case
			//bases
			depthFirstBase = trueDepth / 2;
			depthSecondBase = trueDepth - depthFirstBase;
			secondBase = (firstBase + 1) % 4;
			if(baseCounts[firstBase] != depthFirstBase || baseCounts[secondBase] != depthSecondBase){
				logfile->newLine();
				logfile->conclude("Wrong heterozygous bases in pileup file '" + filename + "' on line " + toString(numLines) + "!");
				return false;
			}

			//calc emission probs
			//set all to full error
			for(b=0; b<genoMap.numGenotypes; ++b)
				emissionProbs[b] = pow(error/3.0, trueDepth);

			//all heterozygous with one correct base
			for(b=1; b<4; ++b){
				emissionProbs[genoMap.getGenotype(firstBase,(firstBase + b) % 4)] = pow(0.5 * (1.0 - error), depthFirstBase);
				emissionProbs[genoMap.getGenotype(firstBase,(firstBase + b) % 4)] *= pow(error / 3.0, depthSecondBase);

				emissionProbs[genoMap.getGenotype(secondBase,(secondBase + b) % 4)] = pow(0.5 * (1.0 - error), depthSecondBase);
				emissionProbs[genoMap.getGenotype(secondBase,(secondBase + b) % 4)] *= pow(error / 3.0, depthFirstBase);
			}

			//all homozygous that contain a correct base
			emissionProbs[genoMap.getGenotype(firstBase,firstBase)] = pow(0.5 * (1.0 - error), trueDepth);
			emissionProbs[genoMap.getGenotype(secondBase,secondBase)] = pow(0.5 * (1.0 - error), trueDepth);

			//correct heterozygous genotype
			emissionProbs[genoMap.getGenotype(firstBase,secondBase)] = pow(1.0-error, trueDepth);
		}

		//now check emission probabilities
		for(b=0; b<genoMap.numGenotypes; ++b){
			relDiff = (stringToDouble(line[b+4]) - emissionProbs[b]) / emissionProbs[b];
			if(relDiff > emissionTolerance){
				logfile->newLine();
				logfile->conclude("Wrong emission probability for genotype " + genoMap.getGenotypeString(b) + " in pileup file '" + filename + "' on line " + toString(numLines) + ": expected " + toString(emissionProbs[b]) + ", found " + line[b+4] + "!");
				return false;
			}
		}
	}
	logfile->done();
	logfile->endIndent();

	return true;
}

TAtlasTest_theta::TAtlasTest_theta(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "theta";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
	simTheta = 	params.getParameterDoubleWithDefault("thetaTest_theta", 0.001);

};

bool TAtlasTest_theta::run(){
	//1) Run ATLAS to simulate BAM file if not yet existant
	//-----------------------------
	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("chrLength", "50000000");
	_testParams.addParameter("ploidy", "2");
	_testParams.addParameter("depth", "2");
	_testParams.addParameter("theta", toString(simTheta));

	//only simulate BAM if it does not already exist
	std::string filenameBAM = filenameTag + ".bam";
	logfile->listFlush("Opening file '" + filenameBAM + "' for reading ...");
	std::ifstream bam(filenameBAM.c_str());
	if(!bam){
		if(!runTGenomeFromInputfile("simulate"))
			return false;
	} else
		logfile->flush("file already exists");

	logfile->newLine();

	//2) Run ATLAS to estimateTheta
	//-----------------------------
	//only simulate BAM if it does not already exist
	std::string filenameTheta = filenameTag + "_theta_estimates.txt";
	logfile->listFlush("Opening file '" + filenameTheta + "' for reading ...");
	std::ifstream in(filenameTheta.c_str());
	if(!in){
		_testParams.addParameter("bam", bamFileName);
		if(!runTGenomeFromInputfile("estimateTheta"))
			return false;
	} else
		logfile->flush("file already exists");
	logfile->newLine();

	//3) check if results are OK
	//--------------------------
	return checkThetaFile();
};

bool TAtlasTest_theta::checkThetaFile(){
	//open theta file
	std::string filename = filenameTag + "_theta_estimates.txt";
	logfile->listFlush("Opening file '" + filename + "' for reading ...");
	std::ifstream in(filename.c_str());
	if(!in)
		throw "Failed to open file '" + filename + "'!";
	logfile->done();

	//some variables
	std::string tmp;
	std::vector<std::string> line;
	float numLines = 0.0;
	float sum = 0.0;
	float mean;

	//skip header
	getline(in, tmp);

	while(in.good() && !in.eof()){
		//read line into vector
		++numLines;
		fillVectorFromLineWhiteSpaceSkipEmpty(in, line);

		//skip empty
		if(line.size() == 0) continue;

		//check columns
		if(line.size() != 14){
			logfile->newLine();
			logfile->conclude("Wrong number of columns in theta file '" + filename + "' on line " + toString(numLines) + "!");
			return false;
		}
		//add new theta value to sum
		sum += stringToDouble(line[10]);

	}
	mean = sum / numLines;
	if(abs(mean - simTheta) > (simTheta / 100.0)){
		logfile->newLine();
		logfile->conclude("Theta was NOT estimated within range! mean estimated theta = " + toString(mean) + " with simulated theta = " + toString(simTheta));
		return false;
	} else {
		logfile->newLine();
		logfile->conclude("Theta was estimated within range! mean estimated theta = " + toString(mean) + " with simulated theta = " + toString(simTheta));
		return true;
	}
}




