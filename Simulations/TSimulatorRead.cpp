/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"

namespace Simulations{

//----------------------------------
//TSimulatorSingleEndRead
//----------------------------------
TSimulatorSingleEndRead::TSimulatorSingleEndRead(std::string readGroupName, int readGroupID, int MaxPrintPhredInt, TRandomGenerator* RandomGenerator, GenotypeLikelihoods::TGenotypeMap & GenoMap):_genoMap(GenoMap){
	_type = "single-end";

	//set variables
	_randomGenerator = RandomGenerator;
	_maxPrintPhredInt = MaxPrintPhredInt;

	_readLengthDist = NULL;
	_readLengthInitialized = false;

	_qualityDist = NULL;
	_qualityDistInitialized = false;
	_qualDistType = "";

	_qualityTransform = NULL;
	_qualityTransformInitialized = false;

	_hasPMD = false;
	_isInitialized = false;

	_isContaminated = false;
	_contaminationRate = 0.0;
	contaminationSource = NULL;

	//initialize bamAlignment
	_name = readGroupName;
	_readGroupID = readGroupID;
	_readNamePrefix = "ATL:0:A:1:" + toString(_readGroupID) + ":"; //"<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	_readXPos = 1;
	_readYPos = 1;
	_alignment.setReadGroup(_readGroupID);
	_alignment.setMappingQuality(50);

	bases = NULL;
	phredIntQualities = NULL;
};

TSimulatorSingleEndRead::~TSimulatorSingleEndRead(){
	if(_readLengthInitialized){
		delete _readLengthDist;
		delete[] bases;
		delete[] phredIntQualities;
	}
	if(_qualityDistInitialized)
		delete _qualityDist;
	if(_qualityTransformInitialized)
		delete _qualityTransform;
};

bool TSimulatorSingleEndRead::checkInitialization(){
	_isInitialized = _readLengthInitialized && _qualityDistInitialized && _qualityTransformInitialized;
	return _isInitialized;
}

void TSimulatorSingleEndRead::setReadLengthDistribution(std::string s, TLog* logfile){
	size_t pos = s.find("(");
	std::string tmp;

	if(pos == std::string::npos)
		throw "Unable to understand read length distribution '" + s + "'!";

	//initialize appropriate function
	std::string type = s.substr(0, pos);
	s.erase(0, pos);

	if(type == "gamma")
		_readLengthDist = new TSimulatorReadLengthGamma(s, _randomGenerator, logfile);
	else if(type == "gammaMode"){
		_readLengthDist = new TSimulatorReadLengthGammaMode(s, _randomGenerator, logfile);

	}
	else if(type == "fixed")
		_readLengthDist = new TReadLengthDistribution(s, _randomGenerator);
	else throw "Unknown read length distribution '" + type + "'!";

	//initialize bases and qualities
	bases = new Base[_readLengthDist->max()];
	phredIntQualities = new int[_readLengthDist->max()];

	//update status
	_readLengthInitialized = true;
	checkInitialization();
};

void TSimulatorSingleEndRead::setQualityDistribution(std::string s){
	size_t pos = s.find("(");
	std::string tmp;

	if(pos == std::string::npos)
		throw "Unable to understand read length distribution '" + s + "'!";

	//initialize appropriate function
	std::string type = s.substr(0, pos);
	s.erase(0, pos);
	if(type == "fixed")
		_qualityDist = new TSimulatorQualityDist(s);
	else if(type == "normal")
		_qualityDist = new TSimulatorQualityDistNormal(s, _randomGenerator);
	else if(type == "binned")
		_qualityDist = new TSimulatorQualityDistBinned(s, _randomGenerator);
	else throw "Unknown read quality distribution '" + type + "'!";

	_qualDistType = type;

	_qualityDistInitialized = true;
	checkInitialization();
};

void TSimulatorSingleEndRead::setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile){
	if(!_qualityDistInitialized)
		throw "Can not initialize quality transformation in TSimulatorRead: quality distribution not initialized!";

	if(parameters.parameters_secondMate != "-" && parameters.parameters_secondMate != parameters.parameters_firstMate)
		throw "Quality transformation for second mate provided, but read group is single end!";

	if(parameters.type == "none")
		_qualityTransform = new TSimulatorQualityTransformation(_qualityDist, _randomGenerator);
	else {
		if(parameters.parameters_firstMate == "-")
			throw "Quality transformation for first mate not provided!";

		if(parameters.type == "recal"){
			_qualityTransform = new TSimulatorQualityTransformationRecal(parameters.parameters_firstMate, _readLengthDist->max(), _qualityDist, _randomGenerator);
		} else
			throw "Unknown quality transformation type '" + parameters.type + "'!";
	}

	_qualityTransformInitialized = true;
	checkInitialization();
};

void TSimulatorSingleEndRead::setPMD(const std::string & pmdStringCT, const std::string & pmdStringGA){
	_pmdObject.initializeFunction(pmdStringCT, GenotypeLikelihoods::pmdCT);
	_pmdObject.initializeFunction(pmdStringGA, GenotypeLikelihoods::pmdGA);

	_hasPMD = _pmdObject.hasDamage();
	checkInitialization();
};

void TSimulatorSingleEndRead::_applyPMD(Base* _bases, const TReadLength & readLength, const bool isReverseStrand){
	if(_hasPMD){
		if(isReverseStrand){
			for(int p=0; p<readLength.read; ++p){
				if(_bases[p] == C ){
					if(_randomGenerator->getRand() < _pmdObject.getProbThreePrime(readLength.fragment - readLength.read + p))
						_bases[p] = T;
				} else if(_bases[p] == G){
					if(_randomGenerator->getRand() < _pmdObject.getProbFivePrime(readLength.read - p - 1))
						_bases[p] = A;
				}
			}
		} else {
			for(int p=0; p<readLength.read; ++p){
				if(_bases[p] == C ){
					if(_randomGenerator->getRand() < _pmdObject.getProbFivePrime(p))
						_bases[p] = T;
				} else if(_bases[p] == G){
					if(_randomGenerator->getRand() < _pmdObject.getProbThreePrime(readLength.fragment - p - 1))
						_bases[p] = A;
				}
			}
		}
	}
};

void TSimulatorSingleEndRead::setContamination(double rate, TSimulatorReference* source){
	_isContaminated = true;
	_contaminationRate = rate;
	contaminationSource = source;

	//check
	if(_contaminationRate < 0.0)
		throw "Contamination rate must be >= 0.0!";
	if(_contaminationRate > 1.0)
		throw "Contamination rate must be <= 0.0!";
	if(_contaminationRate == 0.0)
		_isContaminated = false;
};

std::string TSimulatorSingleEndRead::_getNextReadName(){
	++_readXPos;
	if(_readXPos == 65536){
		++_readYPos;
		_readXPos = 1;
	}
	return _readNamePrefix + toString(_readXPos) + ":" + toString(_readYPos);
};

void TSimulatorSingleEndRead::_simulateBasesQualities(BAM::TAlignment & alignment,
													  Base* haplotype,
													  const uint64_t pos,
													  const TReadLength & readLength,
													  const bool readIsContaminated,
													  const bool isReverse,
													  TSimulatorQualityTransformation* qualityTransform){
	//copy bases
	if(_isContaminated && _randomGenerator->getRand() < _contaminationRate)
		memcpy(bases, contaminationSource->getPointerToRef() + pos, readLength.read);
	else
		memcpy(bases, haplotype + pos, readLength.read);

	//apply PMD
	_applyPMD(bases, readLength, isReverse);

	//simulate qualities and errors
	qualityTransform->simulateQualitiesAndErrors(bases, phredIntQualities, readLength.read, isReverse);

	//copy bases and qualities
	std::string queryBases, qualities;
	for(int p=0; p<readLength.read; ++p){
		queryBases += _genoMap.baseToChar[bases[p]];
		qualities += (char) (std::min(phredIntQualities[p], _maxPrintPhredInt) + 33);
	}

	//fill alignment
	//TODO: verify that it should be negative in case of single end!
	alignment.setSequenceQualitiesOnlyMatches(queryBases, qualities);
	alignment.setIsReverseStrand(isReverse);
	if(isReverse){
		alignment.setInsertSize(-readLength.fragment);
	} else {
		alignment.setInsertSize(readLength.fragment);
	}
};

void TSimulatorSingleEndRead::simulate(Base* haplotype, const uint32_t & refID, const uint32_t & pos, TSimulatorBamFile & bamFile){
	//pick a fragment and read length, strand and contamination
	TReadLength readLength = _readLengthDist->sample();
	bool isReverse = _randomGenerator->getRand() < 0.5;
	bool readIsContaminated = _isContaminated && _randomGenerator->getRand() < _contaminationRate;

	//simulated bases and qualities
	_simulateBasesQualities(_alignment, haplotype, pos, readLength, readIsContaminated, isReverse, _qualityTransform);

	//write bam alignment
	_alignment.move(refID, pos);
	_alignment.setName(_getNextReadName());
	bamFile.saveAlignment(_alignment);
};

void TSimulatorSingleEndRead::printDetails(TLog* logfile){
	logfile->list(_type + ".");
	if(_readLengthInitialized)
		_readLengthDist->printDetails(logfile);
	else throw "Read length distribution not initialized!";

	if(_qualityDistInitialized)
		_qualityDist->printDetails(logfile);
	else throw "Read quality distribution not initialized!";

	if(_qualityTransformInitialized)
		_qualityTransform->printDetails(logfile);
	else throw "Quality transformation not initialized!";

	if(_hasPMD){
		logfile->list("C->T PMD as " + _pmdObject.getFunctionString(GenotypeLikelihoods::pmdCT));
		logfile->list("G->A PMD as " + _pmdObject.getFunctionString(GenotypeLikelihoods::pmdGA));
	} else {
		logfile->list("No PMD.");
	}

	if(_isContaminated)
		logfile->list("Contaminated with rate " + toString(_contaminationRate) + ".");
	else
		logfile->list("Read group is not contaminated.");
};

//----------------------------------
// TSimulatorPairedEndReads
//----------------------------------
TSimulatorPairedEndReads::TSimulatorPairedEndReads(std::string readGroupName, int readGroupNumber, int MaxPrintQual, TRandomGenerator* RandomGenerator, GenotypeLikelihoods::TGenotypeMap & GenoMap)
:TSimulatorSingleEndRead(readGroupName, readGroupNumber, MaxPrintQual, RandomGenerator, GenoMap){
	_type = "paired-end";
	qualityTransform_secondMate = NULL;

	//set SAM flags
	_flags.setIsPaired(true);
	_flags.setIsProperPair(true);
	_flags.setIsRead1(true);
	_flags.setMateIsReverseStrand(true);

	//set SAM flags of second mate
	_mateFlags.setIsPaired(true);
	_mateFlags.setIsProperPair(true);
	_mateFlags.setIsRead2(true);
	_mateFlags.setIsReverseStrand(true);
};

TSimulatorPairedEndReads::~TSimulatorPairedEndReads(){
	if(_qualityTransformInitialized)
		delete qualityTransform_secondMate;

	for(auto& it : bamAlignmentSecondMates_idle)
		delete it;
};

void TSimulatorPairedEndReads::setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile){
	if(!_qualityDistInitialized)
		throw "Can not initialize quality transformation in TSimulatorRead: quality distribution not initialized!";

	if(parameters.type == "none"){
		_qualityTransform = new TSimulatorQualityTransformation(_qualityDist, _randomGenerator);
		qualityTransform_secondMate = new TSimulatorQualityTransformation(_qualityDist, _randomGenerator);
	} else {
		if(parameters.parameters_firstMate == "-")
			throw "Quality transformation for first mate not provided!";
		if(parameters.parameters_secondMate == "-")
			throw "Quality transformation for second mate not provided!";

		if(parameters.type == "recal"){
			_qualityTransform = new TSimulatorQualityTransformationRecal(parameters.parameters_firstMate, _readLengthDist->max(), _qualityDist, _randomGenerator);
			qualityTransform_secondMate = new TSimulatorQualityTransformationRecal(parameters.parameters_secondMate, _readLengthDist->max(), _qualityDist, _randomGenerator);
		} else
			throw "Unknown quality transformation type '" + parameters.type + "'!";
	}

	_qualityTransformInitialized = true;
	checkInitialization();
};

void TSimulatorPairedEndReads::simulate(Base* haplotype, const uint32_t refID, const uint32_t & pos, TSimulatorBamFile & bamFile){
	//pick a fragment, read length and contamination
	TReadLength readLength = _readLengthDist->sample();
	bool readIsContaminated = _isContaminated && _randomGenerator->getRand() < _contaminationRate;

	// Fill FIRST mate
	//------------------
	//simulated bases and qualities
	_simulateBasesQualities(_alignment, haplotype, pos, readLength, readIsContaminated, false, _qualityTransform);

	//write bam alignment
	_alignment.move(refID, pos);
	_alignment.setName(_getNextReadName());
	bamFile.saveAlignment(_alignment);

	// Fill SECOND mate
	//------------------
	//identify position
	uint32_t matePos = pos + readLength.diff();

	//create alignment
	BAM::TAlignment* secondMate;
	if(!bamAlignmentSecondMates_idle.empty()){
		//reuse
		auto it = bamAlignmentSecondMates_idle.begin();
		secondMate = *it;
		bamAlignmentSecondMates_idle.erase(it);
	} else {
		//create new
		secondMate = new BAM::TAlignment;
		secondMate->setReadGroup(_readGroupID);
		secondMate->setMappingQuality(50);
		secondMate->setIsReverseStrand(true);
	}

	//simulated bases and qualities
	_simulateBasesQualities(*secondMate, haplotype, matePos, readLength, readIsContaminated, true, qualityTransform_secondMate);

	//fill bam alignment
	secondMate->move(refID, matePos);
	secondMate->setName(_alignment.name());

	//write if it starts at same position as first, and keep for writing later otherwiese
	if(matePos == pos){
		bamFile.saveAlignment(*secondMate);
		bamAlignmentSecondMates_idle.push_back(secondMate);
	} else {
		bamAlignmentSecondMates.push_back(secondMate);
	}
};

void TSimulatorPairedEndReads::writeUnwrittenAlignments(const long & pos, TSimulatorBamFile & bamFile){
	auto it = bamAlignmentSecondMates.begin();
	while(it != bamAlignmentSecondMates.end()){
		if((*it)->position() <= pos){
			bamFile.saveAlignment(*(*it));
			bamAlignmentSecondMates_idle.push_back(*it);
			it = bamAlignmentSecondMates.erase(it);
		} else
			++it;
	}
};

}; // end namespace
