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
	_logfile = Logfile;
	_randomGenerator = RandomGenerator;

	//set basic things to empty
	_refInitialized = false;
	_writeTrueGenotypes = false;
	_writeVariantInvariantBedFiles = false;
	_sampleSize = 0;
	_seqDepth = 0;
	_averageReadLength = 0;
	_maxReadLength = 0;
	_referenceDivergence = 0.0;
};

void TSimulator::_initializeCommonSettings(TParameters & params){
	//depth
	float depth = params.getParameterDoubleWithDefault("depth", 10.0);
	_logfile->list("Will simulate to an average depth of " + toString(depth) + ".");
	setDepth(depth);

	//base frequencies
	std::vector<float> freq;
	std::string tmp = params.getParameterStringWithDefault("baseFreq", "0.25,0.25,0.25,0.25");
	fillVectorFromString(tmp, freq, ',');
	if(freq.size() != 4) throw "baseFreq vector must have size = 4!";
	setBaseFreq(freq);

	//reference divergence
	_referenceDivergence = params.getParameterDoubleWithDefault("refDiv", 0.01);
	_logfile->list("Will simulate data with reference divergence = " + toString(_referenceDivergence) + ".");

	//fill cumul table for reference divergence
	_cumulRef[0] = 1.0 - _referenceDivergence;
	_cumulRef[1] = _cumulRef[0] + _referenceDivergence / 3.0;
	_cumulRef[2] = _cumulRef[1] + _referenceDivergence / 3.0;
	_cumulRef[3] = 1.0;

	//read groups
	_initializeReadSimulator(params);

	//chromosomes
	initializeChromosomes(params, _logfile);

	//extra output on sites
	_writeTrueGenotypes = params.parameterExists("writeTrueGenotypes");
	_writeVariantInvariantBedFiles = params.parameterExists("writeVariantBED");

	//output name
	_outname = params.getParameterStringWithDefault("out", "ATLAS_simulations");
	_logfile->list("Will write output files with tag '" + _outname + "'.");

	//open FASTA file for reference sequences
	std::string filename = _outname + ".fasta";
	_referenceObj.initialize(filename, _logfile);
};

//--------------------------------------------------------------
//Function to initialize read groups
//--------------------------------------------------------------
std::vector<std::string>  TSimulator::_readSimInfoPerReadGroup(const std::string & Filename, const std::string & Column, const std::string & Name){
	_logfile->listFlush("Reading " + Name + "(s) from file '" + Filename + "' ...");

	TInputFile in(Filename, {"ReadGroup", Column}, "\t", "//");
	std::vector<std::string> vec;
	std::vector<bool> found(_readGroups.size(), false);

	//return map
	std::vector<std::string> ret(_readGroups.size());

	//now parse file
	while(in.read(vec)){
		//find read group
		uint32_t rg = _readGroups.getId(vec[0]);
		found[rg] = true;
		ret[rg] = vec[1];
	}
	_logfile->done();
	_logfile->conclude("Read " + Name + "s for " + toString(in.lineNumber()) + " read groups.");

	//check if there was data for each read group
	for(size_t i = 0; i < found.size(); ++i){
		if(!found[i]){
			throw "No " + Name + " given for read group '" + _readGroups.getName(i) + "' in file '" + Filename + "'!";
		}
	}

	//return
	return ret;
};

void TSimulator::_initializeReadGroup(const std::string & readLengthString, const BAM::TReadGroup & ReadGroup){
	//single or paired end? Is indicated at beginning of readLengthString!
	if(readLengthString.find("single:") == 0){
		_readSimulators.push_back(new TSimulatorSingleEndRead(ReadGroup, _randomGenerator, _genoMap));
	} else if(readLengthString.find("paired:") == 0){
		_readSimulators.push_back(new TSimulatorPairedEndReads(ReadGroup, _randomGenerator, _genoMap));
	} else
		throw "Unable to understand string '" + readLengthString + "'!";

	//add read Length distribution
	std::string readLengthDist = readAfterLast(readLengthString,':');
	_readSimulators.back()->setReadLengthDistribution(readLengthDist, _logfile);
};

void TSimulator::_initializeReadGroupsFromReadLengthDistribution(TParameters & params,
																 const std::string & ParameterName,
																 const std::string & DefaultValue,
																 const std::string & Name){
	_logfile->startIndent("Parsing read length distribution (parameter '" + ParameterName + "'):");
	std::string s = params.getParameterStringWithDefault(ParameterName, DefaultValue);

	_readSimulators.clear();

	//We allow for two options:
	//  1) initialized from the command line (one for all read groups)
	//  2) read-group specific as given in a file

	//check if it is a file (should not contain a ':')
	size_t pos = s.find(":");
	if(pos != std::string::npos){
		//Option 1: a single read length distribution for all
		//---------------------------------------------------------------------
		_logfile->list("Will use '" + s + " for all read groups.");

		//create read groups
		for(auto& r : _readGroups){
			_initializeReadGroup(s, r);
		}
	} else {
		//Option 2: read group specific, given in a file
		//---------------------------------------------------------------------
		std::vector<std::string> dist = _readSimInfoPerReadGroup(s, ParameterName, Name);

		for(uint32_t r = 0; r < dist.size(); ++r){
			_initializeReadGroup(dist[r], _readGroups[r]);
		}
	}
	_logfile->endIndent();
};

void TSimulator::_initializeDistribution(TParameters & params,
										 const std::string & ParameterName,
										 const std::string & DefaultValue,
										 const std::string & Name,
										 void (TSimulatorSingleEndRead::*function)(std::string string)){
	_logfile->startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");
	std::string s = params.getParameterStringWithDefault(ParameterName, DefaultValue);

	//We allow for two options:
	//  1) initialized from the command line (one for all read groups)
	//  2) read-group specific as given in a file

	//check if it is a file (should not contain a ':')
	size_t pos = s.find(":");
	if(pos != std::string::npos){
		//Option 1: a single read distribution for all
		//---------------------------------------------------------------------
		_logfile->list("Will use '" + s + " for all read groups.");

		//create read groups
		for(auto& r : _readSimulators){
			r->*(function)(s);
		}
	} else {
		//Option 2: read group specific, given in a file
		//---------------------------------------------------------------------
		std::vector<std::string> dist = _readSimInfoPerReadGroup(s, ParameterName, Name);

		for(uint32_t r = 0; r < _readSimulators.size().size(); ++r){
			_readSimulators[r]->*(function)(dist[r]);
		}
	}
	_logfile->endIndent();
};

void TSimulator::_initializePMD(TParameters & params,
		 	 	 	 	 	 	const std::string & ParameterName,
		                        const std::string & Name){

	_logfile->startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");

	if(params.parameterExists(ParameterName)){
		std::string pmdString = params.getParameterString(ParameterName);
		_PMD.initialize(pmdString, _readGroups, _logfile);

		//add PMD to simulators
		for(uint32_t r = 0; r < _readSimulators.size(); ++r){
			_readSimulators[r]->setPMD(&_PMD[r]);
		}
	} else {
		_logfile->list("Not simulating any PMD.");
	}
};

void TSimulator::_initializeQualityTransformations(TParameters & params,
	 	 										   const std::string & ParameterName,
												   const std::string & Name){

	_logfile->startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");

	if(params.parameterExists(ParameterName)){
		std::string recalString = params.getParameterString("recal");
		_recal.createModels(recalString, &_readGroups, &_readGroupMap, _logfile);

		//add recal to simulators
		for(uint32_t r = 0; r < _readSimulators.size(); ++r){
			_readSimulators[r]->setPMD(&_PMD[r]);
		}
	} else {
		_logfile->list("Not simulating any quality transformation.");
	}

	_logfile->endIndent();
};

void TSimulator::_initializeContamination(TParameters & params, bool & perReadGroup, std::map<std::string, double> & contaminationMap){
	_logfile->startIndent("Reading contamination:");
	std::string s = params.getParameterStringWithDefault("contamination", "0.0");

	//check if it is a single number or a file
	if(stringIsProbablyANumber(s)){
		//is a numberon the command line
		double rate = convertString<double>(s);
		_logfile->list("Will use a contamination rate of " + toString(rate) + " for all read groups.");
		contaminationMap.emplace("-", rate);
		perReadGroup = false;
	} else {
		//is a file
		_logfile->listFlush("Reading contamination from file '" + s + "' ...");
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
				fillVectorFromStringWhiteSpace(line, vec, true);
				if(vec.size() != 2)
					throw "Found " + toString(vec.size()) + " instead of 2 columns in '" + s + "' on line " + toString(lineNum) + "!\n Expect 1) read group name and 2) contamination rate.";

				//save to map
				if(contaminationMap.find(vec[0]) != contaminationMap.end())
						throw "Duplicated read group name '" + vec[0] + "'in file '" + s + "'!";
				double rate = convertString<double>(s);
				contaminationMap.emplace(vec[0], rate);
			}
		}
		_logfile->done();
		_logfile->conclude("Read distributions for " + toString(contaminationMap.size()) + " read groups.");
		perReadGroup = true;
	}
	_logfile->endIndent();
}

void TSimulator::_addToReadGroupVector(std::vector<std::string> & vec, const std::string & rg){
	//add read group if it does not exist yet
	if(std::find(vec.begin(), vec.end(), rg) == vec.end())
		vec.push_back(rg);
};

void TSimulator::_addReadGroupsIfFile(const std::string & ParameterName, TParameters & Parameters, BAM::TReadGroups & ReadGroups){
	//check if parameter is given
	if(Parameters.parameterExists(ParameterName)){
		std::string s = Parameters.getParameterString(ParameterName);

		//check if string s provides a definition (contains a ':') or is a file (does not contain a ':')
		if(!stringContains(s, ":")){
			//is probably a file -> try to open it
			if(std::filesystem::exists(s)){
				TInputFile in(s, {"ReadGroup"}, "\t", "//");
				std::vector<std::string> tmp;
				while(in.read(tmp)){
					//add all non-existing elemets
					ReadGroups.add(tmp[0]);
				}
			}
		}
	}
};


void TSimulator::_initializeReadSimulator(TParameters & params){
	// For which read groups?
	// Check for each parameter if it is given per read group (a file) or common to all
	_readGroups.clear();
	_addReadGroupsIfFile("readLength", params, _readGroups);
	_addReadGroupsIfFile("qualityDist", params, _readGroups);
	_addReadGroupsIfFile("pmd", params, _readGroups);
	_addReadGroupsIfFile("recal", params, _readGroups);
	_addReadGroupsIfFile("readGroupFreq", params, _readGroups);

	//any read groups specified?
	if(_readGroups.empty()){
		int numRG = params.getParameterIntWithDefault("numReadGroups", 1);
		for(int i=0; i < numRG; ++i){
			_readGroups.add("SimReadGroup" + toString(i+1));
		}

		//report
		if(numRG == 1){
			_logfile->startIndent("Initializing one read group (parameter 'numreadGroups'):");
		} else if(numRg > 1){
			_logfile->startIndent("Initializing " + toString(numRG) + " identical read groups (parameter 'numreadGroups'):");
		} else {
			throw "numReadGroups must be at least 1!";
		}

	} else {
		_logfile->startIndent("Initializing " + toString(_readGroups.size()) + " individual read group(s):");
	}

	// A) read length: used to create simulator read groups
	//---------------
	_initializeReadGroupsFromReadLengthDistribution(params, "readLength", "single:fixed(100)", "read length");

	// B) initialize quality distribution
	//-----------------------------------
	_initializeDistribution(params, "baseQuality", "normal(30,10)[0,93]", "base quality distribution", &TSimulatorSingleEndRead::setQualityDistribution);

	// C) initialize mapping quality distribution
	//-----------------------------------
	_initializeDistribution(params, "mappingQuality", "normal(60,10)[1,255]", "mapping quality distribution", &TSimulatorSingleEndRead::setMappingQualityDistribution);

	// D) initialize PMD
	//------------------
	_initializePMD(params, "pmd", "post-mortem damage");

	// E) initialize quality transformation
	//-------------------------------------
	std::map<std::string, TSimulatorQualityTransformParameters > qualTransformMap;
	bool qualTransformPerReadGroup = false;
	_initializeQualityTransformations(params, qualTransformPerReadGroup, qualTransformMap);

	//add read group names to list
	if(qualTransformPerReadGroup){
		for(std::map<std::string, TSimulatorQualityTransformParameters >::iterator it=qualTransformMap.begin(); it!=qualTransformMap.end(); ++it)
			_addToReadGroupVector(_readGroupNames, it->first);
	}


	// E) initialize contamination
	//----------------------------
	std::map<std::string, double> contaminationMap;
	bool contaminationPerReadGroup = false;
	_initializeContamination(params, contaminationPerReadGroup, contaminationMap);

	//add read group names to list
	if(contaminationPerReadGroup){
		for(std::map<std::string, double>::iterator it=contaminationMap.begin(); it!=contaminationMap.end(); ++it)
			_addToReadGroupVector(_readGroupNames, it->first);
	}

	// F) other things
	//----------------
	int maxPrintQual = params.getParameterIntWithDefault("maxPrintQual", 93);
	_logfile->list("Will print quality scores up to " + toString(maxPrintQual) + ".");
	_logfile->endIndent();

	//now check for read groups: which ones do we simulate?
	//-----------------------------------------------------
	//Option 1: at least one file was given specifying multiple read groups
	if(_readGroupNames.size() > 0){
		//create read groups as specified in the files
		_logfile->startIndent("Initializing " + toString(_readGroupNames.size()) + " read groups:");

		//now initialize
		int rgNumber = 0;
		for(std::vector<std::string>::iterator it=_readGroupNames.begin(); it!=_readGroupNames.end(); ++it, ++rgNumber){
			_logfile->startIndent("Initializing readgroup '" + *it + "':");

			// A) read length
			if(readLengthPerReadGroup){
				std::map<std::string, std::string>::iterator rlIt = readLengthMap.find(*it);
				if(rlIt == readLengthMap.end())
					throw "Read length distribution not specified for read group '" + *it + "'!";

				_initializeReadGroup(rlIt->second, *it, rgNumber, maxPrintQual);
			} else
				_initializeReadGroup(readLengthMap.begin()->second, *it, rgNumber, maxPrintQual);

			// B) quality dist
			if(qualityPerReadGroup){
				std::map<std::string, std::string>::iterator qIt = qualityMap.find(*it);
				if(qIt == qualityMap.end())
					throw "Read quality distribution not specified for read group '" + *it + "'!";
				_readSimulators.back()->setReadLengthDistribution(qIt->second, _logfile);
			} else
				_readSimulators.back()->setQualityDistribution(qualityMap.begin()->second);

			// C) quality transformation
			if(qualTransformPerReadGroup){
				std::map<std::string, TSimulatorQualityTransformParameters >::iterator qtIt = qualTransformMap.find(*it);
				if(qtIt == qualTransformMap.end()){
					//initialize without transformation
					TSimulatorQualityTransformParameters tp("none", "-", "-");
					_readSimulators.back()->setQualityTransformation(tp, _logfile);
				} else
					_readSimulators.back()->setQualityTransformation(qtIt->second, _logfile);
			} else{
				if(_readSimulators.back()->type() == "single"){
					TSimulatorQualityTransformParameters tp(qualTransformMap.begin()->second.type, qualTransformMap.begin()->second.parameters_firstMate, "-");
					_readSimulators.back()->setQualityTransformation(tp, _logfile);
				} else
					_readSimulators.back()->setQualityTransformation(qualTransformMap.begin()->second, _logfile);

			}

			// D) PMD
			if(pmdPerReadGroup){
				std::map<std::string, std::pair<std::string, std::string> >::iterator pmdIt = pmdMap.find(*it);
				if(pmdIt == pmdMap.end()){
					//initialize without transformation
					std::string type="none";
					_readSimulators.back()->setPMD(type, type);
				} else
					_readSimulators.back()->setPMD(pmdIt->second.first, pmdIt->second.second);
			} else
				_readSimulators.back()->setPMD(pmdMap.begin()->second.first, pmdMap.begin()->second.second);

			// E) contamination
			if(contaminationPerReadGroup){
				std::map<std::string, double>::iterator contaminationIt = contaminationMap.find(*it);
				if(contaminationIt != contaminationMap.end())
					_readSimulators.back()->setContamination(contaminationIt->second, &_referenceObj);
			} else
				_readSimulators.back()->setContamination(contaminationMap.begin()->second, &_referenceObj);

			//check and print
			_readSimulators.back()->printDetails(_logfile);
			_logfile->endIndent();
		}
		_logfile->endIndent();
	}
	
	//Option 2: everything provided on command line
	else {
		//If everything was provided on the command line, allow for replicate read groups
		int numRG = params.getParameterIntWithDefault("numReadGroups", 1);
		std::string name;
		_logfile->startIndent("Initializing " + toString(numRG) + " identical read group(s):");

		//now initialize
		for(int i=0; i<numRG; ++i){
			name = "SimReadGroup" + toString(i+1);
			_readGroupNames.push_back(name);
			_logfile->startIndent("Initializing read group '" + name + "':");
			_initializeReadGroup(readLengthMap.begin()->second, name, i, maxPrintQual);
			_readSimulators.back()->setQualityDistribution(qualityMap.begin()->second);
			_readSimulators.back()->setQualityTransformation(qualTransformMap.begin()->second, _logfile);
			_readSimulators.back()->setPMD(pmdMap.begin()->second.first, pmdMap.begin()->second.second);
			_readSimulators.back()->setContamination(contaminationMap.begin()->second, &_referenceObj);

			//check and print
			_readSimulators.back()->printDetails(_logfile);
			_logfile->endIndent();
		}
		_logfile->endIndent();
	}

	//initialize read group frequencies frequencies
	_initializeReadGroupFrequencies(params);
};

void TSimulator::_initializeReadGroupFrequencies(TParameters & params){
	_cumulSimGroupFrequenies.reserve(_readSimulators.size());
	_simGroupFrequencies.reserve(_readSimulators.size());
	if(params.parameterExists("readGroupFreq")){
		//read frequencies
		std::vector<std::string> vec;
		params.fillParameterIntoVector("readGroupFreq", vec, true);
		std::vector<double> freq;
		repeatIndexes(vec, freq);
		if(freq.size() != _readSimulators.size())
			throw "Provided read group frequencies do not match number of read groups!";

		//normalize and print
		double sum = 0;
		for(size_t i=1; i<_readSimulators.size(); ++i)
			sum += freq[i];

		_logfile->startIndent("Will simulate read groups with the following frequencies:");
		for(size_t i=1; i<_readSimulators.size(); ++i){
			_simGroupFrequencies[i] = freq[i] / sum;
			_logfile->list(toString(_simGroupFrequencies[i]) + " " + _readSimulators[i]->name());
		}
		_logfile->endIndent();

		//fill cumulative
		_cumulSimGroupFrequenies[0] = _simGroupFrequencies[0];
		for(size_t i=1; i<_readSimulators.size(); ++i)
			_cumulSimGroupFrequenies[i] = _cumulSimGroupFrequenies[i-1] + _simGroupFrequencies[i];
		_cumulSimGroupFrequenies[_readSimulators.size() - 1] = 1.0; //ensure last entry is 1.0
	} else{
		//equal frequencies
		_logfile->list("Will simulate reads equally distributed among read groups.");
		for(size_t i=0; i<_readSimulators.size(); ++i){
			_simGroupFrequencies[i] = (double) 1.0 / (double) _readSimulators.size();
			_cumulSimGroupFrequenies[i] = (double) (i+1) / (double) _readSimulators.size();
		}
	}

	//precalculate some stuff
	_averageReadLength = 0;
	_maxReadLength = 0;
	int i=0;

	for(TSimulatorSingleEndRead* readSimsIt : _readSimulators){
		_averageReadLength += _simGroupFrequencies[i] * readSimsIt->meanReadLength();
		if(readSimsIt->maxReadLength() > _maxReadLength)
			_maxReadLength = readSimsIt->maxReadLength();
		i++;
	}
}

//--------------------------------------------------------------
//Initialize chromosomes, depth and base frequencies
//--------------------------------------------------------------
void TSimulator::initializeChromosomes(TParameters & params, TLog* logfile){
	std::vector<std::string> string_vec;
	std::vector<uint32_t> chrLength;
	params.fillParameterIntoVectorWithDefault("chrLength", string_vec, ',', "1000000");
	repeatIndexes(string_vec, chrLength);
	std::vector<uint8_t> ploidy;
	if(params.parameterExists("ploidy")){
		params.fillParameterIntoVector("ploidy", string_vec, ',');
		repeatIndexes(string_vec, ploidy);
	} else {
		for(size_t i=0; i<chrLength.size(); ++i)
			ploidy.push_back(2);
	}
	if(ploidy.size() != chrLength.size())
		throw "List of chromosome lengths and ploidies differ in length!";
	for(auto& p : ploidy){
		if(p != 1 && p!=2){
			throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!";
		}
	}

	if(chrLength.size() < 1)
		throw "Issue understanding length of chromosomes!";
	if(chrLength.size() == 1){
		int numChr = params.getParameterIntWithDefault("numChr", 1);
		std::string text = "Will simulate " + toString(numChr) ;
		if(ploidy[0] == 1) text += " haploid";
		else text += " diploid";
		text += " chromosome(s) of length " + toString(chrLength[0]) + " each.";
		logfile->list(text);
		initializeChromosomes(numChr, chrLength[0], ploidy[0]);
	} else {
		logfile->startIndent("Will simulate " + toString(chrLength.size()) + " chromosome(s) of the following length:");
		std::vector<uint8_t>::iterator hIt=ploidy.begin();
		std::string text;
		for(std::vector<uint32_t>::iterator it=chrLength.begin(); it!=chrLength.end(); ++it, ++hIt){
			text = toString(*it) + " (";
			if(*hIt == 1) text += "haploid)";
			else text += "diploid)";
			logfile->list(text);
		}
		initializeChromosomes(chrLength, ploidy);
		logfile->endIndent();
	}
};

void TSimulator::initializeChromosomes(const uint32_t & numChr, const uint32_t & chrLength, const uint8_t & ploidy){
	_chromosomes.clear();
	for(uint32_t i=0; i<numChr; ++i){
		_chromosomes.appendChromosome("chr" + toString(i+1), chrLength, ploidy);
	}
};

void TSimulator::initializeChromosomes(std::vector<uint32_t> & chrLength, std::vector<uint8_t> haploid){
	_chromosomes.clear();
	for(size_t i=0; i<chrLength.size(); ++i){
		_chromosomes.appendChromosome("chr" + toString(i+1), chrLength[i], haploid[i]);
	}
};

void TSimulator::setDepth(float depth){
	_seqDepth = depth;
};

void TSimulator::setBaseFreq(std::vector<float> & freq){
	float sum = 0.0;
	for(int i=0; i<4; ++i){
		_baseFreq[i] = freq[i];
		sum += _baseFreq[i];
	}
	for(int i=0; i<4; ++i){
		_baseFreq[i] /= sum;
	}
	_cumulBaseFreq[0] = _baseFreq[0];
	_cumulBaseFreq[1] = _cumulBaseFreq[0] + _baseFreq[1];
	_cumulBaseFreq[2] = _cumulBaseFreq[1] + _baseFreq[2];
	_cumulBaseFreq[3] = 1.0;

	_logfile->list("Simulating with base frequencies A:" + toString(_baseFreq[0]) + " C:" + toString(_baseFreq[1])+ " G:" + toString(_baseFreq[1])+ " T:" + toString(_baseFreq[1]));
};

//--------------------------------------------------------------
//Run simulations
//--------------------------------------------------------------
void TSimulator::_simulateReadsFromHaplotypes(const BAM::TChromosome & thisChr, Base** haplotypes, TSimulatorBamFile & bamFile, std::string extraProgressText){
	//Initialize probabilities to simulate reads
	uint64_t numReads;
	if(_averageReadLength == 0) numReads = 0;
	else numReads = thisChr.length * _seqDepth / _averageReadLength;

	uint64_t chrLengthForStart = thisChr.length - _maxReadLength + 1;
	double probReadPerSite = 1.0 / (double) chrLengthForStart;
	uint64_t numReadsSimulated = 0;
	uint32_t numReadsHere;

	//initialize progress reporting
	int prog;
	int oldProg = 0;
	std::string progressString = "Simulating about " + toString(numReads) + " reads" + extraProgressText + " ...";

	_logfile->listFlush(progressString);

	//now simulate
	for(uint32_t l=0; l<chrLengthForStart; ++l){
		//write unwritten alignments
		for(TSimulatorSingleEndRead* rs : _readSimulators)
			rs->writeUnwrittenAlignments(l, bamFile);

		//draw random number to get number of reads starting at this position
		numReadsHere = _randomGenerator->getBinomialRand(probReadPerSite, numReads);
		//now simulate
		if(numReadsHere > 0){
			numReadsSimulated += numReadsHere;
			for(uint32_t r=0; r<numReadsHere; ++r){
				int rg = _randomGenerator->pickOne(_readSimulators.size(), _cumulSimGroupFrequenies.data());
				_readSimulators[rg]->simulate(haplotypes[_randomGenerator->sample(2)], thisChr.refID(), l, bamFile);
			}

			//report progress
			prog = 100.0 * (float) numReadsSimulated / (float) numReads;
			if(prog > oldProg){
				oldProg = prog;
				_logfile->listOverFlush(progressString + "(" + toString(prog) + "%)");
			}
		}
	}
	//write unwritten alignments
	for(TSimulatorSingleEndRead* rs : _readSimulators)
		rs->writeUnwrittenAlignments(thisChr.length, bamFile);

	_logfile->overList(progressString + " done!  ");
	_logfile->conclude("Simulated a total of " + toString(numReadsSimulated) + " reads.");
};

void TSimulator::runSimulations(){
	//open bam files
	TSimulatorBamFiles bamFiles(_sampleSize, _outname, _readGroupNames, _chromosomes, _logfile, _genoMap, _qualMap);

	//prepare haplotypes and
	TSimulatorHaplotypes haplotypes(_sampleSize);

	//open files to store extra info on sites
	if(_writeTrueGenotypes){
		//open file for true genotypes
		std::string filename = _outname + "_trueGenotypes.vcf.gz";
		haplotypes.openTrueGenotypeVCF(filename);
	}

	TSimulatorVariantInvariantBedFiles bedFiles;
	if(_writeVariantInvariantBedFiles)
		bedFiles.open(_outname);

	//simulate sequences
	for(auto& chr : _chromosomes){
		_logfile->startIndent("Simulating chromosome " + chr.name + ":");

		//update reference storage and update haplotype lengths
		_referenceObj.setChr(chr.name, chr.length);
		haplotypes.setLength(chr.length);

		//simulate genotypes
		_logfile->listFlush("Simulating genotypes ...");
		if(chr.ploidy == 1)
			_simulateHaplotypesHaploid(haplotypes, chr, _referenceObj);
		else
			_simulateHaplotypesDiploid(haplotypes, chr, _referenceObj);
		_logfile->done();

		//write true genotypes
		if(_writeTrueGenotypes){
			_logfile->listFlush("Writing true genotypes ...");
			haplotypes.writeTrueGenotypes(chr.name, _referenceObj.getPointerToRef(), _genoMap);
			_logfile->done();
		}

		//write BED files
		if(_writeVariantInvariantBedFiles)
			bedFiles.write(haplotypes, chr.name);

		//now simulate and write reads
		_logfile->startIndent("Simulating reads:");
		for(int i=0; i<_sampleSize; ++i)
			_simulateReadsFromHaplotypes(chr, haplotypes.getHaplotypesOfIndividual(i), bamFiles[i], " for individual " + toString(i+1));
		_logfile->endIndent();

		//end of chromosome
		_logfile->endIndent();
	}

	//close stuff
	bamFiles.close();
	_logfile->endIndent();
	haplotypes.closeTrueGenotypeVCF();
	_referenceObj.close();
};

//---------------------------------------------------------
//TSimulatorOneIndividual
//---------------------------------------------------------
TSimulatorOneIndividual::TSimulatorOneIndividual(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator):TSimulator(Logfile, RandomGenerator){
	_logfile->startIndent("Reading parameters to simulate a single individual:");

	//first common stuff
	_initializeCommonSettings(params);
	_sampleSize = 1;

	//now theta
	std::vector<std::string> tmp;
	params.fillParameterIntoVectorWithDefault("theta", tmp, ',', "0.001");
	repeatIndexes(tmp, thetas);
	if(thetas.size() == 1){
		_logfile->list("Will simulate a single individual with theta = " + toString(thetas[0]) + ".");
		for(unsigned int i=1; i<_chromosomes.size(); ++i)
			thetas.push_back(thetas[0]);
	} else {
		_logfile->list("Will simulate a single individual with chromosome specific thetas " + concatenateString(thetas, ", "));
	}

	//one theta per chromosome
	if(thetas.size() != _chromosomes.size())
		throw "Number of theta values provided does not match number of chromosomes to simulate!";

	//done
	_logfile->endIndent();
};

TSimulatorOneIndividual::~TSimulatorOneIndividual(){
	thetas.clear();
};

void TSimulatorOneIndividual::_simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	//fill mutation table
	mutTable.fill(_baseFreq, thetas[chromosome.refID()]);

	for(uint64_t l=0; l<chromosome.length; ++l){
		haplotypes(0,0,l) = static_cast<Base>(_randomGenerator->pickOne(4, _cumulBaseFreq));
		haplotypes(0,1,l) = static_cast<Base>(_randomGenerator->pickOne(4, mutTable[haplotypes(0,0,l)]));

		//decide on reference sequence
		if(haplotypes(0,0,l) == haplotypes(0,1,l))
			ref[l] = static_cast<Base> ((haplotypes(0,0,l) + _randomGenerator->pickOne(4, _cumulRef)) % 4);
		else
			ref[l] = static_cast<Base> (haplotypes(0,_randomGenerator->sample(2),l));
	}
};

void TSimulatorOneIndividual::_simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	//fill mutation table
	mutTable.fill(_baseFreq, thetas[chromosome.refID()]);

	//now simulate genotypes
	for(uint64_t l=0; l<chromosome.length; ++l){
		haplotypes(0,0,l) = static_cast<Base> (_randomGenerator->pickOne(4, _cumulBaseFreq));
		haplotypes(0,1,l) = haplotypes(0,0,l);

		//decide on ref
		ref[l] = static_cast<Base> ((haplotypes(0,0,l) + _randomGenerator->pickOne(4, _cumulRef)) % 4);
	}
};

//---------------------------------------------------------
//TSimulatorPairOfIndividuals
//---------------------------------------------------------
TSimulatorPairOfIndividuals::TSimulatorPairOfIndividuals(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator):TSimulator(Logfile, RandomGenerator){
	_logfile->startIndent("Reading parameters to simulate two individuals with a specific genetic distance:");

	//first common stuff
	_initializeCommonSettings(params);
	_sampleSize = 2;

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
		_logfile->list("Normalizing phi to sum to one (currently summing to " + toString(sum) + ").");
		for(std::vector<double>::iterator it=phis.begin(); it!=phis.end(); ++it)
			*it /= sum;
	}
	_logfile->list("Used phi are: " + concatenateString(phis, ", "));

	//initializes tables
	fillTables();

	//done
	_logfile->endIndent();
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
		sum += _baseFreq[a];
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
				sum += _baseFreq[a] * _baseFreq[b];
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
			sum += _baseFreq[a] * _baseFreq[b];
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
						sum += _baseFreq[a] * _baseFreq[b] * _baseFreq[c];
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
						sum += _baseFreq[a] * _baseFreq[b] * _baseFreq[c];
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

void TSimulatorPairOfIndividuals::_simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	//first run diploid
	_simulateHaplotypesDiploid(haplotypes, chromosome, ref);

	//now set homozygous
	for(uint64_t l=0; l<chromosome.length; ++l){
		//assign to haplotypes
		haplotypes(0,1,l) = haplotypes(0,0,l);
		haplotypes(1,1,l) = haplotypes(1,0,l);
	}
};

void TSimulatorPairOfIndividuals::_simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	//run across loci
	for(uint64_t l=0; l<chromosome.length; ++l){
		//pick a case
		int c = _randomGenerator->pickOne(9, cumulGenoCaseFrequencies);

		//pick genotypes
		int g = _randomGenerator->pickOne(numGenotypeCombinations[c], cumulGenoCombinationFreq[c]);

		//pick order
		int o = _randomGenerator->sample(4);

		//assign to haplotypes
		haplotypes(0,0,l) = genoTrans[c][g][orderLookup[o][0]];
		haplotypes(0,1,l) = genoTrans[c][g][orderLookup[o][1]];
		haplotypes(1,0,l) = genoTrans[c][g][orderLookup[o][2]];
		haplotypes(1,1,l) = genoTrans[c][g][orderLookup[o][3]];

		//simulate reference
		if(c == 0){
			ref[l] = static_cast<Base>((genoTrans[c][g][0] + _randomGenerator->pickOne(4, _cumulRef)) % 4);
		} else {
			int r = _randomGenerator->sample(4);
			ref[l] = genoTrans[c][g][r];
		}
	}
};

//---------------------------------------------------------
//TSimulatorSFS
//---------------------------------------------------------
TSimulatorSFS::TSimulatorSFS(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator):TSimulator(Logfile, RandomGenerator){
	_logfile->startIndent("Reading parameters to simulate a population sample given an SFS:");

	//first common stuff
	_initializeCommonSettings(params);

	//sample size
	_sampleSize = params.getParameterIntWithDefault("sampleSize", 10);

	//read SFS
	_logfile->startIndent("Initializing SFS:");
	if(params.parameterExists("sfs")){
		_logfile->startIndent("Reading SFS from files:");

		std::vector<std::string> tmp;
		std::vector<std::string> sfsFileNames;
		params.fillParameterIntoVector("sfs", tmp, ',');
		repeatIndexes(tmp, sfsFileNames);

		//if a single SFS is given: use it for all chromosomes
		if(sfsFileNames.size() == 1){
			for(size_t i=1; i<_chromosomes.size(); ++i)
				sfsFileNames.emplace_back(sfsFileNames[0]);
		}

		//check if numbe rof chromosomes given matches number of chromosomes
		if(sfsFileNames.size() != _chromosomes.size())
			throw "Number of SFS files does not match number of chromosomes!";

		//initialize SFS from files
		bool folded = params.parameterExists("folded");
		_initializeSFS(sfsFileNames, folded);
	} else if(params.parameterExists("theta")){
		//parse theta from command line
		std::vector<std::string> tmp;
		params.fillParameterIntoVector("theta", tmp, ',');
		std::vector<double> thetas;
		repeatIndexes(tmp, thetas);
		if(thetas.size() == 1){
			_logfile->list("Will simulate from SFS with theta = " + toString(thetas[0]) + ".");
			for(unsigned int i=1; i<_chromosomes.size(); ++i)
				thetas.push_back(thetas[0]);
		} else {
			_logfile->list("Will simulate data from chromosome specific SFS with thetas " + concatenateString(thetas, ", "));
		}
		_initializeSFS(thetas);
	} else throw "Either argument sfs or theta must be provided to simulate population samples!";


	//fill mutation table
	mutTable.fill(_baseFreq);

	//done
	_logfile->endIndent();
};

TSimulatorSFS::~TSimulatorSFS(){
	//deleting SFS
	for(std::vector<SFS*>::iterator it=sfs.begin(); it!=sfs.end(); ++it)
		delete *it;
};

void TSimulatorSFS::_initializeSFS(std::vector<double> & thetas){
	if(thetas.size() != _chromosomes.size())
		throw "Number of theta values does not match number of chromosomes!";

	//generate SFS for each chromosome
	_logfile->listFlush("Initializing SFS ...");
	std::string filename;
	std::vector<double>::iterator thetaIt = thetas.begin();
	for(auto& chr : _chromosomes){
		sfs.push_back(new SFS(chr.ploidy * _sampleSize, (float) *thetaIt));

		//save true SFS
		filename = _outname + "_trueSFS_chr" + toString(chr.refID()) + ".txt";
		(*sfs.rbegin())->writeToFile(filename);
		++thetaIt;
	}
	_logfile->done();
	_logfile->conclude("True SFS written to '" + _outname + "_trueSFS_chr*.txt'.");
};

void TSimulatorSFS::_initializeSFS(std::vector<std::string> & sfsFileNames, bool folded){
	if(sfsFileNames.size() != _chromosomes.size())
		throw "Number of SFS files does not match number of chromosomes!";

	//read the SFS of each chromosome from the corresponding file
	std::vector<SFS*> sfs;
	std::vector<std::string>::iterator it = sfsFileNames.begin();
	for(auto& chr : _chromosomes){
		_logfile->listFlush("Reading the sfs of chromosome '" + chr.name + "' from file '" + *it + "' ...");
		if(folded) sfs.push_back(new SFSfolded(*it));
		else sfs.push_back(new SFS(*it));
		_logfile->done();

		uint32_t nChr = chr.ploidy * _sampleSize;
		if(sfs.back()->numChromosomes != nChr){
			throw "SFS does not match sample size! It contains data for " + toString((*sfs.rbegin())->numChromosomes) + " instead of " + toString(nChr) + " chromosomes.";
		}
		++it;
	}
};

static inline int is_odd(int x){ return x % 2 != 0; }

void TSimulatorSFS::_simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	//now simulate haplotypes
	for(uint32_t l=0; l<chromosome.length; ++l){
		//pick alleles
		Base ancestral = static_cast<Base>(_randomGenerator->pickOne(4, _cumulBaseFreq));
		Base derived = static_cast<Base>(_randomGenerator->pickOne(4, mutTable[ancestral]));

		//pick derived allele frequency
		int alleleCount = sfs[chromosome.refID()]->getRandomAlleleCount(_randomGenerator);

		//pick haplotypes that are derived
		int numNeeded = alleleCount;
		for(int i=0; i<_sampleSize; ++i){
			if(_randomGenerator->getRand() < (double) numNeeded / (double) (_sampleSize - i)){
				haplotypes(i,0,l) = derived;
				--numNeeded;
				if(numNeeded == 0){
					for(int j=i+1; j<_sampleSize; ++j){
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
			if(_randomGenerator->getRand() < (double) alleleCount / (double) _sampleSize)
				ref[l] = derived;
			else
				ref[l] = ancestral;
		} else
			ref[l] = static_cast<Base>((ancestral + _randomGenerator->pickOne(4, _cumulRef)) % 4);
	}

};

void TSimulatorSFS::_simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	int numHaplotypes = 2 * _sampleSize;
	for(uint64_t l=0; l<chromosome.length; ++l){
		//pick alleles
		Base ancestral = static_cast<Base>(_randomGenerator->pickOne(4, _cumulBaseFreq));
		Base derived = static_cast<Base>(_randomGenerator->pickOne(4, mutTable[ancestral]));

		//pick derived allele frequency
		int alleleCount = sfs[chromosome.refID()]->getRandomAlleleCount(_randomGenerator);
		//oo << alleleCount << "\n";

		//pick haplotypes that are derived
		if(alleleCount == 0){
			for(int i=0; i<_sampleSize; ++i){
				haplotypes(i,0,l) = ancestral;
				haplotypes(i,1,l) = ancestral;
			}

			//decide on reference sequence
			ref[l] = static_cast<Base>((ancestral + _randomGenerator->pickOne(4, _cumulRef)) % 4);
		} else {
			int numNeeded = alleleCount;
			for(int i=0; i<numHaplotypes; ++i){
				double prob = (double) numNeeded / (double) (numHaplotypes - i);
				if(_randomGenerator->getRand() < prob){
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
			if(_randomGenerator->getRand() < (double) alleleCount / (double) numHaplotypes)
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
	_logfile->startIndent("Reading parameters to simulate a population sample under Hardy-Weinberg equilibrium:");

	//first common stuff
	_initializeCommonSettings(params);

	//sample size
	_sampleSize = params.getParameterIntWithDefault("sampleSize", 10);

	//parameters of beta distribution
	fracPoly = params.getParameterDoubleWithDefault("fracPoly", 0.1);
	_logfile->list("Will simulate " + toString(fracPoly) + " of all sites as polymorphic. (parameter fracPoly)");
	alpha = params.getParameterDoubleWithDefault("alpha", 0.5);
	if(alpha <= 0.0) throw "Alpha must be > 0!";
	beta = params.getParameterDoubleWithDefault("beta", 0.5);
	if(beta <= 0.0) throw "Beta must be > 0!";
	_logfile->list("Polymoprhic sites will have derived allele frequencies f~Beta(" + toString(alpha) + ", " + toString(beta) + "). (parameters alpha, beta)");
	F = params.getParameterDoubleWithDefault("F", 0.0);
	if(F == 0.0){
		_logfile->list("Will assume no inbreeding. (parameter F=0)");
	} else {
		_logfile->list("Will use an inbreeding coefficient of " + toString(F) + ". (parameter F)");
		if(F < 0.0 || F > 1.0) throw "Inbreeding coefficient F must be within [0,1]!";
	}

	//write true allele freq?
	writeTrueAlleleFreq = false;
	if(params.parameterExists("writeTrueAlleleFreq")){
		std::string alleleFreqFile = _outname + "_trueAlleleFreq.txt.gz";
		_logfile->list("Will write true allele frequencies to file '" + alleleFreqFile + "'.");
		trueFreqFile.open(alleleFreqFile);
		trueFreqFile.writeHeader({"Chr", "Pos", "Ancestral", "Derived", "derivedFreq", "MAF"});
		writeTrueAlleleFreq = true;
	}

	//fill mutation table
	mutTable.fill(_baseFreq);

	//done
	_logfile->endIndent();
};

void TSimulatorHardyWeinberg::_fillCumulGenoProb(const double & f){
	double oneMinus_f = 1.0 - f;
	cumulGenoProb[0] = F * oneMinus_f + (1.0 - F) * oneMinus_f * oneMinus_f;
	cumulGenoProb[1] = cumulGenoProb[0] + (1.0 - F) * 2.0 * f * oneMinus_f;
	cumulGenoProb[2] = 1.0;
};

void TSimulatorHardyWeinberg::_simulateSite(TSimulatorHardyWeinbergSite & site, const std::string & chr, const uint64_t & pos, TSimulatorReference & ref){
	//simulate bases
	site.reference = static_cast<Base>(_randomGenerator->pickOne(4, _cumulBaseFreq));
	site.alternative = static_cast<Base>(_randomGenerator->pickOne(4, mutTable[site.reference]));

	//is the site polymorphic?
	if(_randomGenerator->getRand() < fracPoly){
		site.isPolymorphic = true;
		site.f = _randomGenerator->getBetaRandom(alpha, beta);

		//reference is a random sample from pop with frequency f: flip if ref is alt!
		if(_randomGenerator->getRand() < site.f){
			Base tmp = site.reference;
			site.reference = site.alternative;
			site.alternative = tmp;
			site.f = 1 - site.f;
		}
	} else {
		site.isPolymorphic = false;

		//is reference diverged?
		if(_randomGenerator->getRand() < _referenceDivergence){
			site.f = 1.0;
		} else {
			site.f = 0.0;
		}
	}

	//store reference
	ref[pos] = site.reference;

	//write true frequency: pos is 1 based!
	if(writeTrueAlleleFreq){
		trueFreqFile << chr << pos+1 << _toBase[site.reference] << _toBase[site.alternative] << site.f;
		if(site.f < 0.5){
			trueFreqFile << site.f << std::endl;
		} else {
			trueFreqFile << 1.0 - site.f << std::endl;
		}
	}
};

void TSimulatorHardyWeinberg::_fillhaplotypesMonomoprhic(TSimulatorHaplotypes & haplotypes, const uint64_t & locus, TSimulatorHardyWeinbergSite & site){
	if(site.f == 0.0){
		for(int i=0; i<_sampleSize; ++i){
			haplotypes(i,0,locus) = site.reference;
			haplotypes(i,1,locus) = site.reference;
		}
	} else {
		for(int i=0; i<_sampleSize; ++i){
			haplotypes(i,0,locus) = site.alternative;
			haplotypes(i,1,locus) = site.alternative;
		}
	}
};

void TSimulatorHardyWeinberg::_simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	//storage
	TSimulatorHardyWeinbergSite site;

	//now simulate haplotypes
	for(uint64_t l=0; l<chromosome.length; ++l){
		//simulate site
		_simulateSite(site, chromosome.name, l, ref);

		//polymoprhic or not?
		if(site.isPolymorphic){
			//simulate genotypes
			for(int i=0; i<_sampleSize; ++i){
				if(_randomGenerator->getRand() < site.f){
					haplotypes(i,0,l) = site.alternative;
					haplotypes(i,1,l) = site.alternative;
				} else {
					haplotypes(i,0,l) = site.reference;
					haplotypes(i,1,l) = site.reference;
				}
			}
		} else {
			_fillhaplotypesMonomoprhic(haplotypes, l, site);
		}
	}
};

void TSimulatorHardyWeinberg::_simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){
	//storage
	TSimulatorHardyWeinbergSite site;

	//now simulate haplotypes
	for(uint64_t l=0; l<chromosome.length; ++l){
		//simulate site
		_simulateSite(site, chromosome.name, l, ref);

		//polymoprhic or not?
		if(site.isPolymorphic){
			_fillCumulGenoProb(site.f);

			//simulate genotypes
			for(int i=0; i<_sampleSize; ++i){
				int geno = _randomGenerator->pickOne(3, cumulGenoProb);
				if(geno == 0){
					haplotypes(i,0,l) = site.reference;
					haplotypes(i,1,l) = site.reference;
				} else if(geno == 1){
					if(_randomGenerator->getRand() < 0.5){
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
			_fillhaplotypesMonomoprhic(haplotypes, l, site);
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
