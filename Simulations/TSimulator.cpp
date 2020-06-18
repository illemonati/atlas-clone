/*
 * TSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#include "TSimulator.h"

namespace Simulations{

//---------------------------------------------------
//TSimulator
//---------------------------------------------------
TSimulator::TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	randomGenerator = RandomGenerator;

	//set basic things to empty
	refInitialized = false;
	writeTrueGenotypes = false;
	writeVariantInvariantBedFiles = false;
	sampleSize = 0;
	seqDepth = 0;
	averageReadLength = 0;
	maxReadLength = 0;
	referenceDivergence = 0.0;
};

void TSimulator::initializeCommonSettings(TParameters & params){
	//depth
	float depth = params.getParameterDoubleWithDefault("depth", 10.0);
	logfile->list("Will simulate to an average depth of " + toString(depth) + ".");
	setDepth(depth);

	//base frequencies
	std::vector<float> freq;
	std::string tmp = params.getParameterStringWithDefault("baseFreq", "0.25,0.25,0.25,0.25");
	fillVectorFromString(tmp, freq, ',');
	if(freq.size() != 4) throw "baseFreq vector must have size = 4!";
	setBaseFreq(freq);

	//reference divergence
	referenceDivergence = params.getParameterDoubleWithDefault("refDiv", 0.01);
	logfile->list("Will simulate data with reference divergence = " + toString(referenceDivergence) + ".");

	//fill cumul table for reference divergence
	cumulRef[0] = 1.0 - referenceDivergence;
	cumulRef[1] = cumulRef[0] + referenceDivergence / 3.0;
	cumulRef[2] = cumulRef[1] + referenceDivergence / 3.0;
	cumulRef[3] = 1.0;

	//read groups
	initializeReadSimulator(params);

	//chromosomes
	initializeChromosomes(params, logfile);

	//extra output on sites
	writeTrueGenotypes = params.parameterExists("writeTrueGenotypes");
	writeVariantInvariantBedFiles = params.parameterExists("writeVariantBED");

	//output name
	outname = params.getParameterStringWithDefault("out", "ATLAS_simulations");
	logfile->list("Will write output files with tag '" + outname + "'.");

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	referenceObj.initialize(filename, logfile);
};

//--------------------------------------------------------------
//Function to initialize read groups
//--------------------------------------------------------------
void TSimulator::saveToMap(std::string & name, std::string args, std::map<std::string, std::string> & map, std::string & filename){
	//save, but check if name already exists!
	if(map.find(name) != map.end())
			throw "Duplicated read group name '" + name + "'in file '" + filename + "'!";
	map.emplace(name, args);
};

void TSimulator::initializeReadLengthDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & readLengthMap){
	logfile->startIndent("Reading read length distribution:");
	std::string s = params.getParameterStringWithDefault("readLength", "single:fixed(100)");

	//check if it is a file
	size_t pos = s.find("(");
	if(pos == std::string::npos){
		//is a file
		logfile->listFlush("Reading distribution(s) from file '" + s + "' ...");
		std::ifstream file(s.c_str());
		if(!file)
			throw "Failed to open read-length file '" + s + "!\nEither provide a valid read length distribution, or a valid file name listing this distribution for each read group.";

		//variables
		int lineNum = 0;
		std::string line;
		std::vector<std::string> vec;

		//now parse file
		while(file.good() && !file.eof()){
			++lineNum;
			//skip empty lines or those that start with //
			std::getline(file, line);
			line = extractBefore(line, "//");
			trimString(line);
			if(!line.empty()){
				fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
				if(vec.size() != 2)
					throw "Found " + toString(vec.size()) + " instead of 2 columns in '" + s + "' on line " + toString(lineNum) + "!\n Expect 1) read group name and 2) function string.";

				//check format
				if(vec[1].find("single:") == std::string::npos && vec[1].find("paired:") == std::string::npos)
					throw "Failed to understand read length distribution '" + vec[1] + "': keyword single or paired missing!";

				//save to map
				saveToMap(vec[0], vec[1], readLengthMap, s);
			}
		}
		logfile->done();
		logfile->conclude("Read distributions for " + toString(readLengthMap.size()) + " read groups.");
		perReadGroup = true;
	} else {
		//is a function on the command line
		logfile->list("Will use '" + s + " for all read groups.");
		readLengthMap.emplace("-", s);
		perReadGroup = false;
	}
	logfile->endIndent();
};

void TSimulator::initializeQualityDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & qualityDistMap){
	logfile->startIndent("Reading quality distribution:");
	std::string s = params.getParameterStringWithDefault("qualityDist", "normal(30,10)[0,93]");

	//check if it is a specific function
	size_t pos = s.find("(");
	if(pos == std::string::npos){
		//is a file
		logfile->listFlush("Reading distribution(s) from file '" + s + "' ...");
		std::ifstream file(s.c_str());
		if(!file)
			throw "Failed to open quality distribution file '" + s + "!\nEither provide a valid quality distribution, or a valid file name listing this distribution for each read group.";

		//variables
		int lineNum = 0;
		std::string line;
		std::vector<std::string> vec;

		//now parse file
		while(file.good() && !file.eof()){
			++lineNum;
			//skip empty lines or those that start with //
			std::getline(file, line);
			line = extractBefore(line, "//");
			trimString(line);
			if(!line.empty()){
				fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
				if(vec.size() != 2)
					throw "Found " + toString(vec.size()) + " instead of 2 columns in '" + s + "' on line " + toString(lineNum) + "!\n Expect 1) read group name and 2) function string.";

				//save to map
				saveToMap(vec[0], vec[1], qualityDistMap, s);
			}
		}
		logfile->done();
		logfile->conclude("Read distributions for " + toString(qualityDistMap.size()) + " read groups.");
		perReadGroup = true;
	} else {
		//is a function on the command line
		logfile->list("Will use '" + s + "' for all read groups.");
		qualityDistMap.emplace("-", s);
		perReadGroup = false;
	}
	logfile->endIndent();
}

void TSimulator::initializeQualityTransformations(TParameters & params, bool & perReadGroup, std::map<std::string, TSimulatorQualityTransformParameters > & qualTransformMap){
	//initialize quality transformation
	//Currently we allow for five options:
	//  1) recal transformation initialized from the command line (one for all read groups)//
	//  2) read-group specific recal transformation provided via a recal file	//
	//  3) BQSR transformation initialized from the command line (one for all read groups)
	//  4) read-group specific BQSR transformation from a file (NOT IMPLEMENTED!)
	//  5) no quality transformation

	//map has format: < readGroup, < type, args > >

	logfile->startIndent("Reading quality transformation:");
	std::vector<std::string> string_vec;

	//tmp vars
	std::string::size_type pos;
	std::string line;
	std::vector<std::string> vec;
	std::string rgName;

	//Recal
	//*****
	if(params.parameterExists("recal")){
		std::string recalString = params.getParameterString("recal");
		pos = recalString.find_first_of('[');

		//Option 1: recal from numbers: a single string valid for all read groups.
		//---------------------------------------------------------------------
		if(pos != std::string::npos){

			//save to map
			logfile->list("Will use '" + recalString + "' for all read groups.");
			qualTransformMap["-"] = TSimulatorQualityTransformParameters("recal", recalString, recalString);
			perReadGroup = false;
		}

		//Option 2: a recal input file
		//----------------------------
		else {
			//check if file exists
			logfile->listFlush("Reading transformation(s) from file '" + recalString + "' ...");
			std::ifstream recalFile(recalString.c_str());
			if(!recalFile)
				throw "Failed to open recal file '" + recalString + "!\nEither provide a valid file name or a string of the format model[quality parameters;position parameters; context parameters].";

			//now parse file
			int lineNum = 0;
			std::string tmpString;
			while(recalFile.good() && !recalFile.eof()){
				++lineNum;
				//skip empty lines or those that start with //
				std::getline(recalFile, line);
				line = extractBefore(line, "//");
				trimString(line);
				if(!line.empty()){
					fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);

					//remove header
					if(vec[0] == "readGroup") continue;

					//make sure all columns are given
					if(vec.size() != 6)
						throw "Found " + toString(vec.size()) + " instead of 6 columns in '" + recalString + "' on line " + toString(lineNum) + "!\n expect read group name, mate, model, quality parameters, position parameters and context parameters!";

					//create model string: model[qual;pos;context]
					std::string modelString = vec[2] + "[" + vec[3] + ";" + vec[4] + ";" + vec[5] + "]";

					//check if one mates were already added for this read group
					std::map<std::string, TSimulatorQualityTransformParameters >::iterator mapIt = qualTransformMap.find(vec[0]);
					if(mapIt == qualTransformMap.end()){
						//add new read group
						if(vec[1] == "first")
							qualTransformMap[vec[0]] = TSimulatorQualityTransformParameters("recal", modelString, "-");
						else if(vec[1] == "second")
							qualTransformMap[vec[0]] = TSimulatorQualityTransformParameters("recal", "-", modelString);
						else if(vec[1] == "both")
							qualTransformMap[vec[0]] = TSimulatorQualityTransformParameters("recal", modelString, modelString);
						else
							throw "Unknown mate '" + vec[1] + "' in '" + recalString + "' on line " + toString(lineNum) + "!";
					} else{
						//add to existing read group
						if(vec[1] == "first"){
							if(mapIt->second.parameters_firstMate != "-")
								throw "Duplicate entry for first mate of read group '" + vec[0] + "' in '" + recalString + "' on line " + toString(lineNum) + "!";
							mapIt->second.parameters_firstMate = modelString;
						} else if(vec[1] == "second"){
							if(mapIt->second.parameters_secondMate != "-")
								throw "Duplicate entry for second mate of read group '" + vec[0] + "' in '" + recalString + "' on line " + toString(lineNum) + "!";
							mapIt->second.parameters_secondMate = modelString;
						} else if(vec[1] == "both"){
							if(mapIt->second.parameters_firstMate != "-")
								throw "Duplicate entry for first mate of read group '" + vec[0] + "' in '" + recalString + "' on line " + toString(lineNum) + "!";
							else if(mapIt->second.parameters_secondMate != "-")
								throw "Duplicate entry for second mate of read group '" + vec[0] + "' in '" + recalString + "' on line " + toString(lineNum) + "!";
							mapIt->second.parameters_firstMate = modelString;
							mapIt->second.parameters_secondMate = modelString;
						} else
							throw "Unknown mate '" + vec[1] + "' in '" + recalString + "' on line " + toString(lineNum) + "!";
					}
				}
			}
			logfile->done();
			logfile->conclude("Read transformations for " + toString(qualTransformMap.size()) + " read groups.");
			perReadGroup = true;
		}
	}

	//BQSR
	//****
	else if(params.parameterExists("BQSRTransformation")){
		std::string BQSRString = params.getParameterString("BQSRTransformation");
		pos = BQSRString.find_first_of('[');

		//Option 3: BQSR from numbers: a single one valid for all read groups.
		//---------------------------------------------------------------------
		if(pos != std::string::npos){
			BQSRString.erase(0, pos+1);
			pos = BQSRString.find_first_of(']');
			if(pos == std::string::npos)
				throw "Failed to understand BQSR string: missing '['!\nEither provide a valid file name or the BQSR parameters as '[phi1,phi2,revIntercept]";

			BQSRString.erase(pos, 1);

			//save to map
			logfile->list("Will use '" + BQSRString + "' for all read groups.");
			qualTransformMap["-"] = TSimulatorQualityTransformParameters("BQSR", BQSRString, BQSRString);
			perReadGroup = false;
		} else throw "Failed to understand BQSR string: missing '['!\nEither provide a valid file name or the BQSR parameters as '[phi1,phi2,revIntercept]";
	}

	//No transformation
	//*****************
	else {
		//Option 5: no transformation
		//---------------------------------------------------------------------
		logfile->list("Will print original quality scores for all read groups.");
		qualTransformMap["-"] = TSimulatorQualityTransformParameters("none", "-", "-");
		perReadGroup = false;
	}
	logfile->endIndent();
};


void TSimulator::initializePMD(TParameters & params, bool & perReadGroup, std::map<std::string, std::pair<std::string, std::string> > & pmdMap){
	//map has format: < readGroup, < pmdCT, pmdGA > >

	logfile->startIndent("Reading PMD:");
	//Check if PMD is provided in a file
	if(params.parameterExists("pmdFile")){
		//read from file for each read group
		std::string filename = params.getParameterString("pmdFile");
		logfile->listFlush("Reading PMD from file '" + filename + "' ...");
		std::ifstream file(filename.c_str());
		if(!file) throw "Failed to open PMD file '" + filename + "'!";

		//variables
		int lineNum = 0;
		std::string line;
		std::vector<std::string> vec;
		std::vector<std::string>::iterator nameIt;

		//parse file that has structure: readGroup PMD(CT) PMD(GA)
		while(file.good() && !file.eof()){
			++lineNum;
			//skip empty lines or those that start with //
			std::getline(file, line);
			line = extractBefore(line, "//");
			trimString(line);
			if(!line.empty()){
				fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
				if(vec.size() != 3) throw "Found " + toString(vec.size()) + " instead of 3 columns in '" + filename + "' on line " + toString(lineNum) + "!";

				//save read group name, but check if name already exists!
				if(pmdMap.find(vec[0]) != pmdMap.end())
					throw "Duplicated read group name '" + vec[0] + "'in file '" + filename + "'!";

				//save to map
				pmdMap[vec[0]] = std::pair<std::string, std::string>(vec[1], vec[2]);
			}
		}

		//close file
		file.close();

		logfile->done();
		logfile->conclude("Read PMD for " + toString(pmdMap.size()) + " read groups.");
		perReadGroup = true;
	}

	//Read from command line
	else {
		if(params.parameterExists("pmd")){
			std::string pmdString = params.getParameterString("pmd");
			pmdMap["-"] = std::pair<std::string, std::string>(pmdString, pmdString);
		} else {
			std::string pmdCTString;
			std::string pmdGAString;

			if(params.parameterExists("pmdCT")){
				pmdCTString = params.getParameterString("pmdCT");
			} else pmdCTString = "none";

			if(params.parameterExists("pmdGA")){
				pmdGAString = params.getParameterString("pmdGA");
			} else pmdGAString = "none";

			//add to map
			pmdMap["-"] = std::pair<std::string, std::string>(pmdCTString, pmdGAString);
			perReadGroup = false;
		}
		logfile->list("Will use PMD as provided on the command line for all read groups.");
	}
	logfile->endIndent();
};

void TSimulator::initializeContamination(TParameters & params, bool & perReadGroup, std::map<std::string, double> & contaminationMap){
	logfile->startIndent("Reading contamination:");
	std::string s = params.getParameterStringWithDefault("contamination", "0.0");

	//check if it is a single number or a file
	if(stringIsProbablyANumber(s)){
		//is a numberon the command line
		double rate = convertString<double>(s);
		logfile->list("Will use a contamination rate of " + toString(rate) + " for all read groups.");
		contaminationMap.emplace("-", rate);
		perReadGroup = false;
	} else {
		//is a file
		logfile->listFlush("Reading contamination from file '" + s + "' ...");
		std::ifstream file(s.c_str());
		if(!file)
			throw "Failed to open quality distribution file '" + s + "!\nEither provide a valid quality distribution, or a valid file name listing this distribution for each read group.";

		//variables
		int lineNum = 0;
		std::string line;
		std::vector<std::string> vec;

		//now parse file
		while(file.good() && !file.eof()){
			++lineNum;
			//skip empty lines or those that start with //
			std::getline(file, line);
			line = extractBefore(line, "//");
			trimString(line);
			if(!line.empty()){
				fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
				if(vec.size() != 2)
					throw "Found " + toString(vec.size()) + " instead of 2 columns in '" + s + "' on line " + toString(lineNum) + "!\n Expect 1) read group name and 2) contamination rate.";

				//save to map
				if(contaminationMap.find(vec[0]) != contaminationMap.end())
						throw "Duplicated read group name '" + vec[0] + "'in file '" + s + "'!";
				double rate = convertString<double>(s);
				contaminationMap.emplace(vec[0], rate);
			}
		}
		logfile->done();
		logfile->conclude("Read distributions for " + toString(contaminationMap.size()) + " read groups.");
		perReadGroup = true;
	}
	logfile->endIndent();
}

void TSimulator::addToReadGroupVector(std::vector<std::string> & vec, const std::string & rg){
	//add read group if it does not exist yet
	if(std::find(vec.begin(), vec.end(), rg) == vec.end())
		vec.push_back(rg);
};

void TSimulator::initializeReadGroup(const std::string & readLengthString, std::string & readGroupName, int rgNumber, int maxPrintQual){
	//single or paired end? Is indicated at beginning of readLengthString!
	if(readLengthString.find("single:") == 0){
		readSimulators.push_back(new TSimulatorSingleEndRead(readGroupName, rgNumber, maxPrintQual, randomGenerator));
	} else if(readLengthString.find("paired:") == 0){
		readSimulators.push_back(new TSimulatorPairedEndReads(readGroupName, rgNumber, maxPrintQual, randomGenerator));
	} else
		throw "Unable to understand read length string '" + readLengthString + "'!";

	//add read Length distribution
	std::string readLengthDist = readAfterLast(readLengthString,':');
	readSimulators.back()->setReadLengthDistribution(readLengthDist, logfile);
};

void TSimulator::initializeReadSimulator(TParameters & params){
	// A) read length
	//---------------
	std::map<std::string, std::string> readLengthMap;
	bool readLengthPerReadGroup = false;
	initializeReadLengthDistribution(params, readLengthPerReadGroup, readLengthMap);

	//add read group names to list
	if(readLengthPerReadGroup){
		for(std::map<std::string, std::string>::iterator it=readLengthMap.begin(); it!=readLengthMap.end(); ++it)
			addToReadGroupVector(readGroupNames, it->first);
	}

	// B) initialize quality distribution
	//-----------------------------------
	std::map<std::string, std::string> qualityMap;
	bool qualityPerReadGroup = false;
	initializeQualityDistribution(params, qualityPerReadGroup, qualityMap);

	//add read group names to list
	if(qualityPerReadGroup){
		for(std::map<std::string, std::string>::iterator it=qualityMap.begin(); it!=qualityMap.end(); ++it)
			addToReadGroupVector(readGroupNames, it->first);
	}

	// C) initialize quality transformation
	//-------------------------------------
	std::map<std::string, TSimulatorQualityTransformParameters > qualTransformMap;
	bool qualTransformPerReadGroup = false;
	initializeQualityTransformations(params, qualTransformPerReadGroup, qualTransformMap);

	//add read group names to list
	if(qualTransformPerReadGroup){
		for(std::map<std::string, TSimulatorQualityTransformParameters >::iterator it=qualTransformMap.begin(); it!=qualTransformMap.end(); ++it)
			addToReadGroupVector(readGroupNames, it->first);
	}

	// D) initialize PMD
	//------------------
	std::map<std::string, std::pair<std::string, std::string> > pmdMap;
	bool pmdPerReadGroup = false;
	initializePMD(params, pmdPerReadGroup, pmdMap);

	//add read group names to list
	if(pmdPerReadGroup){
		for(std::map<std::string, std::pair<std::string, std::string> >::iterator it=pmdMap.begin(); it!=pmdMap.end(); ++it)
			addToReadGroupVector(readGroupNames, it->first);
	}

	// E) initialize contamination
	//----------------------------
	std::map<std::string, double> contaminationMap;
	bool contaminationPerReadGroup = false;
	initializeContamination(params, contaminationPerReadGroup, contaminationMap);

	//add read group names to list
	if(contaminationPerReadGroup){
		for(std::map<std::string, double>::iterator it=contaminationMap.begin(); it!=contaminationMap.end(); ++it)
			addToReadGroupVector(readGroupNames, it->first);
	}

	// F) other things
	//----------------
	int maxPrintQual = params.getParameterIntWithDefault("maxPrintQual", 93);
	logfile->list("Will print quality scores up to " + toString(maxPrintQual) + ".");
	logfile->endIndent();

	//now check for read groups: which ones do we simulate?
	//-----------------------------------------------------
	//Option 1: at least one file was given specifying multiple read groups
	if(readGroupNames.size() > 0){
		//create read groups as specified in the files
		logfile->startIndent("Initializing " + toString(readGroupNames.size()) + " read groups:");

		//now initialize
		int rgNumber = 1;
		for(std::vector<std::string>::iterator it=readGroupNames.begin(); it!=readGroupNames.end(); ++it, ++rgNumber){
			logfile->startIndent("Initializing readgroup '" + *it + "':");

			// A) read length
			if(readLengthPerReadGroup){
				std::map<std::string, std::string>::iterator rlIt = readLengthMap.find(*it);
				if(rlIt == readLengthMap.end())
					throw "Read length distribution not specified for read group '" + *it + "'!";

				initializeReadGroup(rlIt->second, *it, rgNumber, maxPrintQual);
			} else
				initializeReadGroup(readLengthMap.begin()->second, *it, rgNumber, maxPrintQual);

			// B) quality dist
			if(qualityPerReadGroup){
				std::map<std::string, std::string>::iterator qIt = qualityMap.find(*it);
				if(qIt == qualityMap.end())
					throw "Read quality distribution not specified for read group '" + *it + "'!";
				readSimulators.back()->setReadLengthDistribution(qIt->second, logfile);
			} else
				readSimulators.back()->setQualityDistribution(qualityMap.begin()->second);

			// C) quality transformation
			if(qualTransformPerReadGroup){
				std::map<std::string, TSimulatorQualityTransformParameters >::iterator qtIt = qualTransformMap.find(*it);
				if(qtIt == qualTransformMap.end()){
					//initialize without transformation
					TSimulatorQualityTransformParameters tp("none", "-", "-");
					readSimulators.back()->setQualityTransformation(tp, logfile);
				} else
					readSimulators.back()->setQualityTransformation(qtIt->second, logfile);
			} else{
				if(readSimulators.back()->type == "single"){
					TSimulatorQualityTransformParameters tp(qualTransformMap.begin()->second.type, qualTransformMap.begin()->second.parameters_firstMate, "-");
					readSimulators.back()->setQualityTransformation(tp, logfile);
				} else
					readSimulators.back()->setQualityTransformation(qualTransformMap.begin()->second, logfile);

			}

			// D) PMD
			if(pmdPerReadGroup){
				std::map<std::string, std::pair<std::string, std::string> >::iterator pmdIt = pmdMap.find(*it);
				if(pmdIt == pmdMap.end()){
					//initialize without transformation
					std::string type="none";
					readSimulators.back()->setPMD(type, type);
				} else
					readSimulators.back()->setPMD(pmdIt->second.first, pmdIt->second.second);
			} else
				readSimulators.back()->setPMD(pmdMap.begin()->second.first, pmdMap.begin()->second.second);

			// E) contamination
			if(contaminationPerReadGroup){
				std::map<std::string, double>::iterator contaminationIt = contaminationMap.find(*it);
				if(contaminationIt != contaminationMap.end())
					readSimulators.back()->setContamination(contaminationIt->second, &referenceObj);
			} else
				readSimulators.back()->setContamination(contaminationMap.begin()->second, &referenceObj);

			//check and print
			readSimulators.back()->printDetails(logfile);
			logfile->endIndent();
		}
		logfile->endIndent();
	}
	
	//Option 2: everything provided on command line
	else {
		//If everything was provided on the command line, allow for replicate read groups
		int numRG = params.getParameterIntWithDefault("numReadGroups", 1);
		std::string name;
		logfile->startIndent("Initializing " + toString(numRG) + " identical read group(s):");

		//now initialize
		for(int i=0; i<numRG; ++i){
			name = "SimReadGroup" + toString(i+1);
			readGroupNames.push_back(name);
			logfile->startIndent("Initializing read group '" + name + "':");
			initializeReadGroup(readLengthMap.begin()->second, name, i+1, maxPrintQual);
			readSimulators.back()->setQualityDistribution(qualityMap.begin()->second);
			readSimulators.back()->setQualityTransformation(qualTransformMap.begin()->second, logfile);
			readSimulators.back()->setPMD(pmdMap.begin()->second.first, pmdMap.begin()->second.second);
			readSimulators.back()->setContamination(contaminationMap.begin()->second, &referenceObj);

			//check and print
			readSimulators.back()->printDetails(logfile);
			logfile->endIndent();
		}
	}

	//initialize read group frequencies frequencies
	initializeReadGroupFrequencies(params);
}

void TSimulator::initializeReadGroupFrequencies(TParameters & params){
	cumulSimGroupFrequenies.reserve(readSimulators.size());
	simGroupFrequencies.reserve(readSimulators.size());
	if(params.parameterExists("readGroupFreq")){
		//read frequencies
		std::vector<std::string> vec;
		params.fillParameterIntoVector("readGroupFreq", vec, true);
		std::vector<double> freq;
		repeatIndexes(vec, freq);
		if(freq.size() != readSimulators.size())
			throw "Provided read group frequencies do not match number of read groups!";

		//normalize and print
		double sum = 0;
		for(size_t i=1; i<readSimulators.size(); ++i)
			sum += freq[i];

		logfile->startIndent("Will simulate read groups with the following frequencies:");
		for(size_t i=1; i<readSimulators.size(); ++i){
			simGroupFrequencies[i] = freq[i] / sum;
			logfile->list(toString(simGroupFrequencies[i]) + " " + readSimulators[i]->name());
		}
		logfile->endIndent();

		//fill cumulative
		cumulSimGroupFrequenies[0] = simGroupFrequencies[0];
		for(size_t i=1; i<readSimulators.size(); ++i)
			cumulSimGroupFrequenies[i] = cumulSimGroupFrequenies[i-1] + simGroupFrequencies[i];
		cumulSimGroupFrequenies[readSimulators.size() - 1] = 1.0; //ensure last entry is 1.0
	} else{
		//equal frequencies
		logfile->list("Will simulate reads equally distributed among read groups.");
		for(size_t i=0; i<readSimulators.size(); ++i){
			simGroupFrequencies[i] = (double) 1.0 / (double) readSimulators.size();
			cumulSimGroupFrequenies[i] = (double) (i+1) / (double) readSimulators.size();
		}
	}

	//precalculate some stuff
	averageReadLength = 0;
	maxReadLength = 0;
	int i=0;

	for(TSimulatorSingleEndRead* readSimsIt : readSimulators){
		averageReadLength += simGroupFrequencies[i] * readSimsIt->meanReadLength();
		if(readSimsIt->maxReadLength() > maxReadLength)
			maxReadLength = readSimsIt->maxReadLength();
		i++;
	}
}

//--------------------------------------------------------------
//Initialize chromosomes, depth and base frequencies
//--------------------------------------------------------------
void TSimulator::initializeChromosomes(TParameters & params, TLog* logfile){
	std::vector<std::string> string_vec;
	std::vector<long> chrLength;
	params.fillParameterIntoVectorWithDefault("chrLength", string_vec, ',', "1000000");
	repeatIndexes(string_vec, chrLength);
	std::vector<int> ploidy;
	if(params.parameterExists("ploidy")){
		params.fillParameterIntoVector("ploidy", string_vec, ',');
		repeatIndexes(string_vec, ploidy);
	} else {
		for(size_t i=0; i<chrLength.size(); ++i)
			ploidy.push_back(2);
	}
	if(ploidy.size() != chrLength.size())
		throw "List of chromosome lengths and ploidies differ in length!";
	std::vector<bool> haploid;
	for(std::vector<int>::iterator it=ploidy.begin(); it!=ploidy.end(); ++it){
		if(*it == 1) haploid.push_back(true);
		else if(*it == 2) haploid.push_back(false);
		else throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!";
	}

	if(chrLength.size() < 1)
		throw "Issue understanding length of chromosomes!";
	if(chrLength.size() == 1){
		int numChr = params.getParameterIntWithDefault("numChr", 1);
		std::string text = "Will simulate " + toString(numChr) ;
		if(haploid[0]) text += " haploid";
		else text += " diploid";
		text += " chromosome(s) of length " + toString(chrLength[0]) + " each.";
		logfile->list(text);
		initializeChromosomes(numChr, chrLength[0], haploid[0]);
	} else {
		logfile->startIndent("Will simulate " + toString(chrLength.size()) + " chromosome(s) of the following length:");
		std::vector<bool>::iterator hIt=haploid.begin();
		std::string text;
		for(std::vector<long>::iterator it=chrLength.begin(); it!=chrLength.end(); ++it, ++hIt){
			text = toString(*it) + " (";
			if(*hIt) text += "haploid)";
			else text += "diploid)";
			logfile->list(text);
		}
		initializeChromosomes(chrLength, haploid);
		logfile->endIndent();
	}
}

void TSimulator::initializeChromosomes(int numChr, long chrLength, bool haploid){
	chromosomes.clear();
	for(int i=0; i<numChr; ++i){
		chromosomes.push_back(TSimulatorChromosome("chr" + toString(i+1), i, chrLength, haploid));
	}
}

void TSimulator::initializeChromosomes(std::vector<long> & chrLength, std::vector<bool> haploid){
	chromosomes.clear();
	for(unsigned int i=0; i<chrLength.size(); ++i){
		chromosomes.push_back(TSimulatorChromosome("chr" + toString(i+1), i, chrLength[i], haploid[i]));
	}
}

void TSimulator::setDepth(float depth){
	seqDepth = depth;
}

void TSimulator::setBaseFreq(std::vector<float> & freq){
	float sum = 0.0;
	for(int i=0; i<4; ++i){
		baseFreq[i] = freq[i];
		sum += baseFreq[i];
	}
	for(int i=0; i<4; ++i){
		baseFreq[i] /= sum;
	}
	cumulBaseFreq[0] = baseFreq[0];
	cumulBaseFreq[1] = cumulBaseFreq[0] + baseFreq[1];
	cumulBaseFreq[2] = cumulBaseFreq[1] + baseFreq[2];
	cumulBaseFreq[3] = 1.0;

	logfile->list("Simulating with base frequencies A:" + toString(baseFreq[0]) + " C:" + toString(baseFreq[1])+ " G:" + toString(baseFreq[1])+ " T:" + toString(baseFreq[1]));
}

//--------------------------------------------------------------
//Run simulations
//--------------------------------------------------------------
void TSimulator::simulateReadsFromHaplotypes(std::vector<TSimulatorChromosome>::iterator & thisChr, Base** haplotypes, TSimulatorBamFile & bamFile, std::string extraProgressText){
	//Initialize probabilities to simulate reads
	uint64_t numReads;
	if(averageReadLength == 0) numReads = 0;
	else numReads = thisChr->length * seqDepth / averageReadLength;

	uint64_t chrLengthForStart = thisChr->length - maxReadLength + 1;
	double probReadPerSite = 1.0 / (double) chrLengthForStart;
	uint64_t numReadsSimulated = 0;
	uint32_t numReadsHere;

	//prepare bam alignment
	for(std::vector<TSimulatorSingleEndRead*>::iterator readSimsIt = readSimulators.begin(); readSimsIt!=readSimulators.end(); ++readSimsIt)
		(*readSimsIt)->setRefId(thisChr->refID);

	//initialize progress reporting
	int prog;
	int oldProg = 0;
	std::string progressString = "Simulating about " + toString(numReads) + " reads" + extraProgressText + " ...";

	logfile->listFlush(progressString);

	//now simulate
	for(uint64_t l=0; l<chrLengthForStart; ++l){
		//write unwritten alignments
		for(TSimulatorSingleEndRead* rs : readSimulators)
			rs->writeUnwrittenAlignments(l, bamFile);

		//draw random number to get number of reads starting at this position
		numReadsHere = randomGenerator->getBiomialRand(probReadPerSite, numReads);
		//now simulate
		if(numReadsHere > 0){
			numReadsSimulated += numReadsHere;
			for(uint32_t r=0; r<numReadsHere; ++r){
				readSimulators[randomGenerator->pickOne(readSimulators.size(), cumulSimGroupFrequenies.data())]->simulate(haplotypes[randomGenerator->pickOne(2)], l, bamFile);
			}

			//report progress
			prog = 100.0 * (float) numReadsSimulated / (float) numReads;
			if(prog > oldProg){
				oldProg = prog;
				logfile->listOverFlush(progressString + "(" + toString(prog) + "%)");
			}
		}
	}
	//write unwritten alignments
	for(TSimulatorSingleEndRead* rs : readSimulators)
		rs->writeUnwrittenAlignments(thisChr->length, bamFile);

	logfile->overList(progressString + " done!  ");
	logfile->conclude("Simulated a total of " + toString(numReadsSimulated) + " reads.");
};

void TSimulator::runSimulations(){
	//open bam files
	TSimulatorBamFiles bamFiles(sampleSize, outname, readGroupNames, chromosomes, logfile);

	//prepare haplotypes and
	TSimulatorHaplotypes haplotypes(sampleSize);

	//open files to store extra info on sites
	if(writeTrueGenotypes){
		//open file for true genotypes
		std::string filename = outname + "_trueGenotypes.vcf.gz";
		haplotypes.openTrueGenotypeVCF(filename);
	}

	TSimulatorVariantInvariantBedFiles bedFiles;
	if(writeVariantInvariantBedFiles)
		bedFiles.open(outname);

	//simulate sequences
	int refId = 0;
	for(std::vector<TSimulatorChromosome>::iterator chrIt = chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId){
		logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

		//update reference storage and update haplotype lengths
		referenceObj.setChr(chrIt->name, chrIt->length);
		haplotypes.setLength(chrIt->length);

		//simulate genotypes
		logfile->listFlush("Simulating genotypes ...");
		if(chrIt->haploid)
			simulateHaplotypesHaploid(haplotypes, *chrIt, referenceObj.getPointerToRef());
		else
			simulateHaplotypesDiploid(haplotypes, *chrIt, referenceObj.getPointerToRef());
		logfile->done();

		//write true genotypes
		if(writeTrueGenotypes){
			logfile->listFlush("Writing true genotypes ...");
			haplotypes.writeTrueGenotypes(*chrIt, referenceObj.getPointerToRef(), genoMap);
			logfile->done();
		}

		//write BED files
		if(writeVariantInvariantBedFiles)
			bedFiles.write(haplotypes, *chrIt);

		//now simulate and write reads
		logfile->startIndent("Simulating reads:");
		for(int i=0; i<sampleSize; ++i)
			simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesOfIndividual(i), bamFiles[i], " for individual " + toString(i+1));
		logfile->endIndent();

		//end of chromosome
		logfile->endIndent();
	}

	//close stuff
	bamFiles.close();
	logfile->endIndent();
	haplotypes.closeTrueGenotypeVCF();
	referenceObj.close();
};

//---------------------------------------------------------
//TSimulatorOneIndividual
//---------------------------------------------------------
TSimulatorOneIndividual::TSimulatorOneIndividual(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator):TSimulator(Logfile, RandomGenerator){
	logfile->startIndent("Reading parameters to simulate a single individual:");

	//first common stuff
	initializeCommonSettings(params);
	sampleSize = 1;

	//now theta
	std::vector<std::string> tmp;
	params.fillParameterIntoVectorWithDefault("theta", tmp, ',', "0.001");
	repeatIndexes(tmp, thetas);
	if(thetas.size() == 1){
		logfile->list("Will simulate a single individual with theta = " + toString(thetas[0]) + ".");
		for(unsigned int i=1; i<chromosomes.size(); ++i)
			thetas.push_back(thetas[0]);
	} else {
		logfile->list("Will simulate a single individual with chromosome specific thetas " + concatenateString(thetas, ", "));
	}

	//one theta per chromosome
	if(thetas.size() != chromosomes.size())
		throw "Number of theta values provided does not match number of chromosomes to simulate!";

	//done
	logfile->endIndent();
};

TSimulatorOneIndividual::~TSimulatorOneIndividual(){
	thetas.clear();
};

void TSimulatorOneIndividual::simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	//fill mutation table
	mutTable.fill(baseFreq, thetas[chromosome.refID]);

	for(uint64_t l=0; l<chromosome.length; ++l){
		haplotypes(0,0,l) = static_cast<Base>(randomGenerator->pickOne(4, cumulBaseFreq));
		haplotypes(0,1,l) = static_cast<Base>(randomGenerator->pickOne(4, mutTable[haplotypes(0,0,l)]));

		//decide on reference sequence
		if(haplotypes(0,0,l) == haplotypes(0,1,l))
			ref[l] = static_cast<Base> ((haplotypes(0,0,l) + randomGenerator->pickOne(4, cumulRef)) % 4);
		else
			ref[l] = static_cast<Base> (haplotypes(0,randomGenerator->pickOne(2),l));
	}
};

void TSimulatorOneIndividual::simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	//fill mutation table
	mutTable.fill(baseFreq, thetas[chromosome.refID]);

	//now simulate genotypes
	for(uint64_t l=0; l<chromosome.length; ++l){
		haplotypes(0,0,l) = static_cast<Base> (randomGenerator->pickOne(4, cumulBaseFreq));
		haplotypes(0,1,l) = haplotypes(0,0,l);

		//decide on ref
		ref[l] = static_cast<Base> ((haplotypes(0,0,l) + randomGenerator->pickOne(4, cumulRef)) % 4);
	}
};

//---------------------------------------------------------
//TSimulatorPairOfIndividuals
//---------------------------------------------------------
TSimulatorPairOfIndividuals::TSimulatorPairOfIndividuals(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator):TSimulator(Logfile, RandomGenerator){
	logfile->startIndent("Reading parameters to simulate two individuals with a specific genetic distance:");

	//first common stuff
	initializeCommonSettings(params);
	sampleSize = 2;

	//Initialize phis
	std::vector<std::string> tmp;
	params.fillParameterIntoVector("phi", tmp, ',');
	repeatIndexes(tmp, phis);

	if(phis.size() != 9)
		throw "Wrong number of phi! Required are nine values for genotype combinations 00/00, 00/01, 01/00, 00/11, 01/01, 01/02, 00/12, 01/22, 01/23";

	//normalize phis
	double sum = 0.0;
	for(std::vector<double>::iterator it=phis.begin(); it!=phis.end(); ++it)
		sum += *it;
	if(sum != 1.0){
		logfile->list("Normalizing phi to sum to one (currently summing to " + toString(sum) + ").");
		for(std::vector<double>::iterator it=phis.begin(); it!=phis.end(); ++it)
			*it /= sum;
	}
	logfile->list("Used phi are: " + concatenateString(phis, ", "));

	//initializes tables
	fillTables();

	//done
	logfile->endIndent();
};


void TSimulatorPairOfIndividuals::fillTables(){
	//file cumulative frequencies of cases (phis)
	double sum = 0.0;
	int genoCase = 0;
	for(std::vector<double>::iterator it = phis.begin(); it != phis.end(); ++it, ++genoCase){
		sum += *it;
		cumulGenoCaseFrequencies[genoCase] = sum;
	}
	cumulGenoCaseFrequencies[genoCase] = 1.0;
	if(fabs(sum - 1.0) > 0.0000000001)
		throw "Phis do not sum to 1.0! They sum to " + toString(sum) + ".";

	//prepare genotype frequency tables for each case
	cumulGenoCombinationFreq = new double*[9];
	genoTrans = new Base**[9];

	//some variables
	int a,b,c,d;
	int index;
	double* cumul = new double[24];

	//case 0: aa/aa
	//-----------------------------------------
	cumulGenoCombinationFreq[0] = new double[4];
	numGenotypeCombinations[0] = 4;
	genoTrans[0] = new Base*[4];
	sum = 0.0;
	for(a=0; a<4; ++a){
		sum += baseFreq[a];
		cumulGenoCombinationFreq[0][a] = sum;
		genoTrans[0][a] = new Base[4];
		genoTrans[0][a][0] = static_cast<Base>(a);
		genoTrans[0][a][1] = static_cast<Base>(a);
		genoTrans[0][a][2] = static_cast<Base>(a);
		genoTrans[0][a][3] = static_cast<Base>(a);
	}

	//cases 1 to 3: aa/ab, ab/aa, aa/bb
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	index = 0;
	sum = 0.0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				sum += baseFreq[a] * baseFreq[b];
				cumul[index] = sum;
				++index;
			}
		}
	}
	//normalize
	for(index=0; index<12; ++index)
		cumul[index] /= sum;

	//now initialize
	for(int ca = 1; ca<4; ++ca){
		cumulGenoCombinationFreq[ca] = new double[12];
		numGenotypeCombinations[ca] = 12;
		genoTrans[ca] = new Base*[12];
		for(index=0; index<12; ++index){
			genoTrans[ca][index] = new Base[4];
			cumulGenoCombinationFreq[ca][index] = cumul[index];
		}
	}

	//assign genotype translations
	index = 0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				genoTrans[1][index][0] = static_cast<Base>(a); genoTrans[1][index][1] = static_cast<Base>(a); genoTrans[1][index][2] = static_cast<Base>(a); genoTrans[1][index][3] = static_cast<Base>(b);
				genoTrans[2][index][0] = static_cast<Base>(a); genoTrans[2][index][1] = static_cast<Base>(b); genoTrans[2][index][2] = static_cast<Base>(a); genoTrans[2][index][3] = static_cast<Base>(a);
				genoTrans[3][index][0] = static_cast<Base>(a); genoTrans[3][index][1] = static_cast<Base>(a); genoTrans[3][index][2] = static_cast<Base>(b); genoTrans[3][index][3] = static_cast<Base>(b);
				++index;
			}
		}
	}

	//cases 4: ab/ab
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	index = 0;
	sum = 0.0;
	genoTrans[4] = new Base*[6];
	for(a=0; a<3; ++a){
		for(b=a+1; b<4; ++b){
			sum += baseFreq[a] * baseFreq[b];
			cumul[index] = sum;
			genoTrans[4][index] = new Base[4];
			genoTrans[4][index][0] = static_cast<Base>(a); genoTrans[4][index][1] = static_cast<Base>(b); genoTrans[4][index][2] = static_cast<Base>(a); genoTrans[4][index][3] = static_cast<Base>(b);
			++index;
		}
	}
	//normalize
	for(index=0; index<6; ++index)
		cumul[index] /= sum;

	//now initialize
	cumulGenoCombinationFreq[4] = new double[6];
	numGenotypeCombinations[4] = 6;
	for(index=0; index<6; ++index){
		cumulGenoCombinationFreq[4][index] = cumul[index];

	}

	//case 5: ab/ac
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	index = 0;
	sum = 0.0;
	genoTrans[5] = new Base*[24];
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						sum += baseFreq[a] * baseFreq[b] * baseFreq[c];
						cumul[index] = sum;
						genoTrans[5][index] = new Base[4];
						genoTrans[5][index][0] = static_cast<Base>(a); genoTrans[5][index][1] = static_cast<Base>(b); genoTrans[5][index][2] = static_cast<Base>(a); genoTrans[5][index][3] = static_cast<Base>(c);
						++index;
					}
				}
			}
		}
	}
	//normalize
	for(index=0; index<24; ++index)
		cumul[index] /= sum;

	//now initialize
	cumulGenoCombinationFreq[5] = new double[24];
	numGenotypeCombinations[5] = 24;
	for(index=0; index<24; ++index){
		cumulGenoCombinationFreq[5][index] = cumul[index];
	}

	//cases 6 and 7: aa/bc, ab/cc
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	index = 0;
	sum = 0.0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						sum += baseFreq[a] * baseFreq[b] * baseFreq[c];
						cumul[index] = sum;
						++index;
					}
				}
			}
		}
	}
	//normalize
	for(index=0; index<24; ++index)
		cumul[index] /= sum;

	//now initialize
	for(int ca = 6; ca<8; ++ca){
		cumulGenoCombinationFreq[ca] = new double[24];
		numGenotypeCombinations[ca] = 24;
		genoTrans[ca] = new Base*[24];
		for(index=0; index<24; ++index){
			cumulGenoCombinationFreq[ca][index] = cumul[index];
			genoTrans[ca][index] = new Base[4];
		}
	}

	//assign genotype translations
	index = 0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						genoTrans[5][index][0] = static_cast<Base>(a); genoTrans[5][index][1] = static_cast<Base>(b); genoTrans[5][index][2] = static_cast<Base>(a); genoTrans[5][index][3] = static_cast<Base>(c);
						genoTrans[6][index][0] = static_cast<Base>(a); genoTrans[6][index][1] = static_cast<Base>(a); genoTrans[6][index][2] = static_cast<Base>(b); genoTrans[6][index][3] = static_cast<Base>(c);
						genoTrans[7][index][0] = static_cast<Base>(a); genoTrans[7][index][1] = static_cast<Base>(b); genoTrans[7][index][2] = static_cast<Base>(c); genoTrans[7][index][3] = static_cast<Base>(c);
						++index;
					}
				}
			}
		}
	}

	//case 8: ab/cd
	//-----------------------------------------
	cumulGenoCombinationFreq[8] = new double[24];
	numGenotypeCombinations[8] = 24;
	genoTrans[8] = new Base*[24];
	index = 0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						for(d=0; d<4; ++d){
							if(d!=a && d!=b && d!=c){
								cumulGenoCombinationFreq[8][index] = (double) (index+1.0) / 24.0;
								genoTrans[8][index] = new Base[4];
								genoTrans[8][index][0] = static_cast<Base>(a);
								genoTrans[8][index][1] = static_cast<Base>(b);
								genoTrans[8][index][2] = static_cast<Base>(c);
								genoTrans[8][index][3] = static_cast<Base>(d);
								++index;
							}
						}
					}
				}
			}
		}
	}

	//prepare haplotype order table to randomly pick
	//----------------------------------------------
	orderLookup = new short*[4];
	orderLookup[0] = new short[4];
	orderLookup[0][0] = 0; orderLookup[0][1] = 1; orderLookup[0][2] = 2; orderLookup[0][3] = 3;
	orderLookup[1] = new short[4];
	orderLookup[1][0] = 0; orderLookup[1][1] = 1; orderLookup[1][2] = 3; orderLookup[1][3] = 2;
	orderLookup[2] = new short[4];
	orderLookup[2][0] = 1; orderLookup[2][1] = 0; orderLookup[2][2] = 2; orderLookup[2][3] = 3;
	orderLookup[3] = new short[4];
	orderLookup[3][0] = 1; orderLookup[3][1] = 0; orderLookup[3][2] = 3; orderLookup[3][3] = 2;

	//clean up
	delete[] cumul;

	//set as initialized
	tablesInitialized = true;
};

void TSimulatorPairOfIndividuals::deleteTables(){
	if(tablesInitialized){
		for(int i=0; i<9; ++i){
			delete[] cumulGenoCombinationFreq[i];
			for(int j=0; j<numGenotypeCombinations[i]; ++j)
				delete[] genoTrans[i][j];
			delete[] genoTrans[i];
		}
		delete[] cumulGenoCombinationFreq;
		delete[] genoTrans;

		//and lookup
		for(int i=0; i<4; ++i)
			delete[] orderLookup[i];
		delete[] orderLookup;

		tablesInitialized = false;
	}
};

void TSimulatorPairOfIndividuals::simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	//first run diploid
	simulateHaplotypesDiploid(haplotypes, chromosome, ref);

	//now set homozygous
	for(uint64_t l=0; l<chromosome.length; ++l){
		//assign to haplotypes
		haplotypes(0,1,l) = haplotypes(0,0,l);
		haplotypes(1,1,l) = haplotypes(1,0,l);
	}
};

void TSimulatorPairOfIndividuals::simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	//run across loci
	for(uint64_t l=0; l<chromosome.length; ++l){
		//pick a case
		int c = randomGenerator->pickOne(9, cumulGenoCaseFrequencies);

		//pick genotypes
		int g = randomGenerator->pickOne(numGenotypeCombinations[c], cumulGenoCombinationFreq[c]);

		//pick order
		int o = randomGenerator->pickOne(4);

		//assign to haplotypes
		haplotypes(0,0,l) = genoTrans[c][g][orderLookup[o][0]];
		haplotypes(0,1,l) = genoTrans[c][g][orderLookup[o][1]];
		haplotypes(1,0,l) = genoTrans[c][g][orderLookup[o][2]];
		haplotypes(1,1,l) = genoTrans[c][g][orderLookup[o][3]];

		//simulate reference
		if(c == 0){
			ref[l] = static_cast<Base>((genoTrans[c][g][0] + randomGenerator->pickOne(4, cumulRef)) % 4);
		} else {
			int r = randomGenerator->pickOne(4);
			ref[l] = genoTrans[c][g][r];
		}
	}
};

//---------------------------------------------------------
//TSimulatorSFS
//---------------------------------------------------------
TSimulatorSFS::TSimulatorSFS(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator):TSimulator(Logfile, RandomGenerator){
	logfile->startIndent("Reading parameters to simulate a population sample given an SFS:");

	//first common stuff
	initializeCommonSettings(params);

	//sample size
	sampleSize = params.getParameterIntWithDefault("sampleSize", 10);

	//read SFS
	logfile->startIndent("Initializing SFS:");
	if(params.parameterExists("sfs")){
		logfile->startIndent("Reading SFS from files:");

		std::vector<std::string> tmp;
		std::vector<std::string> sfsFileNames;
		params.fillParameterIntoVector("sfs", tmp, ',');
		repeatIndexes(tmp, sfsFileNames);

		//if a single SFS is given: use it for all chromosomes
		if(sfsFileNames.size() == 1){
			for(size_t i=1; i<chromosomes.size(); ++i)
				sfsFileNames.emplace_back(sfsFileNames[0]);
		}

		//check if numbe rof chromosomes given matches number of chromosomes
		if(sfsFileNames.size() != chromosomes.size())
			throw "Number of SFS files does not match number of chromosomes!";

		//initialize SFS from files
		bool folded = params.parameterExists("folded");
		initializeSFS(sfsFileNames, folded);
	} else if(params.parameterExists("theta")){
		//parse theta from command line
		std::vector<std::string> tmp;
		params.fillParameterIntoVector("theta", tmp, ',');
		std::vector<double> thetas;
		repeatIndexes(tmp, thetas);
		if(thetas.size() == 1){
			logfile->list("Will simulate from SFS with theta = " + toString(thetas[0]) + ".");
			for(unsigned int i=1; i<chromosomes.size(); ++i)
				thetas.push_back(thetas[0]);
		} else {
			logfile->list("Will simulate data from chromosome specific SFS with thetas " + concatenateString(thetas, ", "));
		}
		initializeSFS(thetas);
	} else throw "Either argument sfs or theta must be provided to simulate population samples!";


	//fill mutation table
	mutTable.fill(baseFreq);

	//done
	logfile->endIndent();
};

TSimulatorSFS::~TSimulatorSFS(){
	//deleting SFS
	for(std::vector<SFS*>::iterator it=sfs.begin(); it!=sfs.end(); ++it)
		delete *it;
};

void TSimulatorSFS::initializeSFS(std::vector<double> & thetas){
	if(thetas.size() != chromosomes.size())
		throw "Number of theta values does not match number of chromosomes!";

	//generate SFS for each chromosome
	logfile->listFlush("Initializing SFS ...");
	int chr = 1;
	std::string filename;
	std::vector<double>::iterator thetaIt = thetas.begin();
	for(std::vector<TSimulatorChromosome>::iterator chrIt = chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++thetaIt, ++chr){
		sfs.push_back(new SFS((2 - chrIt->haploid) * sampleSize, (float) *thetaIt));

		//save true SFS
		filename = outname + "_trueSFS_chr" + toString(chr) + ".txt";
		(*sfs.rbegin())->writeToFile(filename);
	}
	logfile->done();
	logfile->conclude("True SFS written to '" + outname + "_trueSFS_chr*.txt'.");
};

void TSimulatorSFS::initializeSFS(std::vector<std::string> & sfsFileNames, bool folded){
	if(sfsFileNames.size() != chromosomes.size())
		throw "Number of SFS files does not match number of chromosomes!";

	//read the SFS of each chromosome from the corresponding file
	std::vector<SFS*> sfs;
	std::vector<std::string>::iterator it = sfsFileNames.begin();
	int nChr;
	for(std::vector<TSimulatorChromosome>::iterator chrIt = chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++it){
		logfile->listFlush("Reading the sfs of chromosome '" + chrIt->name + "' from file '" + *it + "' ...");
		if(folded) sfs.push_back(new SFSfolded(*it));
		else sfs.push_back(new SFS(*it));
		logfile->done();

		nChr = (2-chrIt->haploid) * sampleSize;
		if((*sfs.rbegin())->numChromosomes != nChr)
			throw "SFS does not match sample size! It contains data for " + toString((*sfs.rbegin())->numChromosomes) + " instead of " + toString(nChr) + " chromosomes.";
	}
};

static inline int is_odd(int x){ return x % 2 != 0; }

void TSimulatorSFS::simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	//now simulate haplotypes
	for(uint64_t l=0; l<chromosome.length; ++l){
		//pick alleles
		Base ancestral = static_cast<Base>(randomGenerator->pickOne(4, cumulBaseFreq));
		Base derived = static_cast<Base>(randomGenerator->pickOne(4, mutTable[ancestral]));

		//pick derived allele frequency
		int alleleCount = sfs[chromosome.refID]->getRandomAlleleCount(randomGenerator);

		//pick haplotypes that are derived
		int numNeeded = alleleCount;
		for(int i=0; i<sampleSize; ++i){
			if(randomGenerator->getRand() < (double) numNeeded / (double) (sampleSize - i)){
				haplotypes(i,0,l) = derived;
				--numNeeded;
				if(numNeeded == 0){
					for(int j=i+1; j<sampleSize; ++j){
						haplotypes(i,0,l) = ancestral;
						haplotypes(i,1,l) = ancestral;
					}
					break;
				}
			} else
				haplotypes(i,0,l) = ancestral;

			//make homozygous
			haplotypes(i,1,l) = haplotypes(i,0,l);
		}

		//decide on reference sequence
		if(alleleCount > 0){
			if(randomGenerator->getRand() < (double) alleleCount / (double) sampleSize)
				ref[l] = derived;
			else
				ref[l] = ancestral;
		} else
			ref[l] = static_cast<Base>((ancestral + randomGenerator->pickOne(4, cumulRef)) % 4);
	}

};

void TSimulatorSFS::simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	int numHaplotypes = 2 * sampleSize;
	for(uint64_t l=0; l<chromosome.length; ++l){
		//pick alleles
		Base ancestral = static_cast<Base>(randomGenerator->pickOne(4, cumulBaseFreq));
		Base derived = static_cast<Base>(randomGenerator->pickOne(4, mutTable[ancestral]));

		//pick derived allele frequency
		int alleleCount = sfs[chromosome.refID]->getRandomAlleleCount(randomGenerator);
		//oo << alleleCount << "\n";

		//pick haplotypes that are derived
		if(alleleCount == 0){
			for(int i=0; i<sampleSize; ++i){
				haplotypes(i,0,l) = ancestral;
				haplotypes(i,1,l) = ancestral;
			}

			//decide on reference sequence
			ref[l] = static_cast<Base>((ancestral + randomGenerator->pickOne(4, cumulRef)) % 4);
		} else {
			int numNeeded = alleleCount;
			for(int i=0; i<numHaplotypes; ++i){
				double prob = (double) numNeeded / (double) (numHaplotypes - i);
				if(randomGenerator->getRand() < prob){
					haplotypes(i / 2, is_odd(i), l) = derived;
					--numNeeded;
					if(numNeeded == 0){
						for(int j=i+1; j<numHaplotypes; ++j)
							haplotypes(j / 2, is_odd(j), l) = ancestral;
						break;
					}
				} else
					haplotypes(i / 2, is_odd(i), l) = ancestral;
			}

			//decide on reference sequence
			if(randomGenerator->getRand() < (double) alleleCount / (double) numHaplotypes)
				ref[l] = derived;
			else
				ref[l] = ancestral;
		}
	}
};

//---------------------------------------------------------
//TSimulatorHardyWeinberg
//---------------------------------------------------------
TSimulatorHardyWeinberg::TSimulatorHardyWeinberg(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator):TSimulator(Logfile, RandomGenerator){
	logfile->startIndent("Reading parameters to simulate a population sample under Hardy-Weinberg equilibrium:");

	//first common stuff
	initializeCommonSettings(params);

	//sample size
	sampleSize = params.getParameterIntWithDefault("sampleSize", 10);

	//parameters of beta distribution
	fracPoly = params.getParameterDoubleWithDefault("fracPoly", 0.1);
	logfile->list("Will simulate " + toString(fracPoly) + " of all sites as polymorphic.");
	alpha = params.getParameterDoubleWithDefault("alpha", 0.5);
	if(alpha <= 0.0) throw "Alpha must be > 0!";
	beta = params.getParameterDoubleWithDefault("beta", 0.5);
	if(beta <= 0.0) throw "Beta must be > 0!";
	logfile->list("Polymoprhic sites will have derived allele frequencies f~Beta(" + toString(alpha) + ", " + toString(beta) + ").");
	F = params.getParameterDoubleWithDefault("F", 0.0);
	if(F > 0.0) logfile->list("Will use an inbreeding coefficient of " + toString(F) + ".");
	if(F < 0.0 || F > 1.0) throw "Inbreeding coefficient F must be within [0,1]!";

	//write true allele freq?
	writeTrueAlleleFreq = false;
	if(params.parameterExists("writeTrueAlleleFreq")){
		std::string alleleFreqFile = outname + "_trueAlleleFreq.txt.gz";
		logfile->list("Will write true allele frequencies to file '" + alleleFreqFile + "'.");
		trueFreqFile.open(alleleFreqFile);
		trueFreqFile.writeHeader({"Chr", "Pos", "Ancestral", "Derived", "derivedFreq", "MAF"});
		writeTrueAlleleFreq = true;
	}

	//fill mutation table
	mutTable.fill(baseFreq);

	//done
	logfile->endIndent();
};

void TSimulatorHardyWeinberg::fillCumulGenoProb(const double & f){
	double oneMinus_f = 1.0 - f;
	cumulGenoProb[0] = F * oneMinus_f + (1.0 - F) * oneMinus_f * oneMinus_f;
	cumulGenoProb[1] = cumulGenoProb[0] + (1.0 - F) * 2.0 * f * oneMinus_f;
	cumulGenoProb[2] = 1.0;
};

void TSimulatorHardyWeinberg::simulateSite(TSimulatorHardyWeinbergSite & site, const std::string & chr, const uint64_t & pos, Base* & ref){
	//simulate bases
	site.reference = static_cast<Base>(randomGenerator->pickOne(4, cumulBaseFreq));
	site.alternative = static_cast<Base>(randomGenerator->pickOne(4, mutTable[site.reference]));

	//is the site polymorphic?
	if(randomGenerator->getRand() < fracPoly){
		site.isPolymorphic = true;
		site.f = randomGenerator->getBetaRandom(alpha, beta);

		//reference is a random sample from pop with frequency f: flip if ref is alt!
		if(randomGenerator->getRand() < site.f){
			Base tmp = site.reference;
			site.reference = site.alternative;
			site.alternative = tmp;
			site.f = 1 - site.f;
		}
	} else {
		site.isPolymorphic = false;

		//is reference diverged?
		if(randomGenerator->getRand() < referenceDivergence){
			site.f = 1.0;
		} else {
			site.f = 0.0;
		}
	}

	//store reference
	ref[pos] = site.reference;

	//write true frequency: pos is 1 based!
	if(writeTrueAlleleFreq){
		trueFreqFile << chr << pos+1 << toBase[site.reference] << toBase[site.alternative] << site.f;
		if(site.f < 0.5){
			trueFreqFile << site.f << std::endl;
		} else {
			trueFreqFile << 1.0 - site.f << std::endl;
		}
	}
};

void TSimulatorHardyWeinberg::fillhaplotypesMonomoprhic(TSimulatorHaplotypes & haplotypes, const uint64_t & locus, TSimulatorHardyWeinbergSite & site){
	if(site.f == 0.0){
		for(int i=0; i<sampleSize; ++i){
			haplotypes(i,0,locus) = site.reference;
			haplotypes(i,1,locus) = site.reference;
		}
	} else {
		for(int i=0; i<sampleSize; ++i){
			haplotypes(i,0,locus) = site.alternative;
			haplotypes(i,1,locus) = site.alternative;
		}
	}
};

void TSimulatorHardyWeinberg::simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	//storage
	TSimulatorHardyWeinbergSite site;

	//now simulate haplotypes
	for(uint64_t l=0; l<chromosome.length; ++l){
		//simulate site
		simulateSite(site, chromosome.name, l, ref);

		//polymoprhic or not?
		if(site.isPolymorphic){
			//simulate genotypes
			for(int i=0; i<sampleSize; ++i){
				if(randomGenerator->getRand() < site.f){
					haplotypes(i,0,l) = site.alternative;
					haplotypes(i,1,l) = site.alternative;
				} else {
					haplotypes(i,0,l) = site.reference;
					haplotypes(i,1,l) = site.reference;
				}
			}
		} else {
			fillhaplotypesMonomoprhic(haplotypes, l, site);
		}
	}
};

void TSimulatorHardyWeinberg::simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){
	//storage
	TSimulatorHardyWeinbergSite site;

	//now simulate haplotypes
	for(uint64_t l=0; l<chromosome.length; ++l){
		//simulate site
		simulateSite(site, chromosome.name, l, ref);

		//polymoprhic or not?
		if(site.isPolymorphic){
			fillCumulGenoProb(site.f);

			//simulate genotypes
			for(int i=0; i<sampleSize; ++i){
				int geno = randomGenerator->pickOne(3, cumulGenoProb);
				if(geno == 0){
					haplotypes(i,0,l) = site.reference;
					haplotypes(i,1,l) = site.reference;
				} else if(geno == 1){
					if(randomGenerator->getRand() < 0.5){
						haplotypes(i,0,l) = site.reference;
						haplotypes(i,1,l) = site.alternative;
					} else {
						haplotypes(i,0,l) = site.alternative;
						haplotypes(i,1,l) = site.reference;
					}
				} else {
					haplotypes(i,0,l) = site.alternative;
					haplotypes(i,1,l) = site.alternative;
				}
			}
		} else {
			fillhaplotypesMonomoprhic(haplotypes, l, site);
		}
	}
};

//--------------------------------------------------------------------
//Functions to simulate pooled data
//--------------------------------------------------------------------
//TODO: Need to switch to haplotype model

/*
void TSimulator::simulatePooledData(int sampleSize, SFS & sfs, std::string outname){
	//open BAM file
	openBamFile(outname + ".bam");

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	openFastaFile(filename);

	//prepare variables
	float* altFreq = NULL;
	long numReads;
	long chrLengthForStart;
	double probReadPerSite;
	int numReadsHere;
	long numReadsSimulated;
	initializeQualToErrorTable();

	//open frequency file
	filename = outname + "_frequencies.txt";
	std::ofstream freqFile(filename.c_str());

	//simulate sequences
	int refId = 0;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId){
		logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

		//simulate reference and alternative sequence
		simulateReferenceAndAlternativeSequenceCurChromosome();

		//simulate alternative frequencies (and write to file)
		logfile->listFlush("Simulating alternative allele frequencies ...");
		delete[] altFreq;
		altFreq = new float[chrIt->length];
		for(int l=0; l<chrIt->length; ++l){
			altFreq[l] = sfs.getRandomFrequency(randomGenerator);
			freqFile << chrIt->name << "\t" << l+1 << altFreq[l] << "\n";
		}
		logfile->done();

		//simulating reads
		numReads = chrIt->length * seqDepth / readLength;
		chrLengthForStart = chrIt->length - readLength;
		probReadPerSite = 1.0 / (double) chrLengthForStart;
		numReadsSimulated = 0;
		bamAlignment.RefID = refId;
		int prog;
		int oldProg = 0;
		std::string progressString = "Simulating about " + toString(numReads) + " reads ...";
		logfile->listFlush(progressString);
		for(long l=0; l<chrLengthForStart; ++l){
			//draw random number to get number of reads starting at this position
			numReadsHere = randomGenerator->getBiomialRand(probReadPerSite, numReads);

			//now simulate
			if(numReadsHere > 0){
				simulateReads(numReadsHere, l, altFreq);
				numReadsSimulated += numReadsHere;

				//report progress
				prog = 100.0 * (float) numReadsSimulated / (float) numReads;
				if(prog > oldProg){
					oldProg = prog;
					logfile->listOverFlush(progressString + "(" + toString(prog) + "%)");
				}
			}
		}
		logfile->overList(progressString + " done!  ");
		logfile->conclude("Simulated a total of " + toString(numReadsSimulated) + " reads.");
		logfile->endIndent();
	}

	//close stuff
	closeBamFile();
	closeFastaFile();
	freqFile.close();

	//clear memory
	delete[] altFreq;
}
*/

}; //end namespace
