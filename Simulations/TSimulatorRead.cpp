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
TSimulatorSingleEndRead::TSimulatorSingleEndRead(const BAM::TReadGroup & ReadGroup, TRandomGenerator* RandomGenerator) :
	_readGroup(ReadGroup) {
	_type = "single-end";

	//set variables
	_randomGenerator = RandomGenerator;

	_pmdObject = NULL;
	_hasPMD = false;
	_recalObject = NULL;
	_hasRecal = false;
	_isInitialized = false;

	_isContaminated = false;
	_contaminationRate = 0.0;
	contaminationSource = NULL;

	//initialize bamAlignment
	_readNamePrefix = "ATL:0:A:1:" + coretools::str::toString(_readGroup.id) + ":"; //"<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	_readXPos = 1;
	_readYPos = 1;
	_alignment.setReadGroup(_readGroup.id);
	_alignment.setMappingQuality(50);

};

bool TSimulatorSingleEndRead::checkInitialization(){
	_isInitialized = _readLengthDist && _qualityDist && _mappingQualityDist;
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
		_readLengthDist = std::make_unique<TSimulatorReadLengthGamma>(s, _randomGenerator, logfile);
	else if(type == "gammaMode"){
		_readLengthDist = std::make_unique<TSimulatorReadLengthGammaMode>(s, _randomGenerator, logfile);
	}
	else if(type == "fixed")
		_readLengthDist = std::make_unique<TReadLengthDistribution>(s, _randomGenerator);
	else throw "Unknown read length distribution '" + type + "'!";

	//update status
	checkInitialization();
};

void TSimulatorSingleEndRead::_initializeQualityDistribution(std::string s, std::unique_ptr<TSimulatorQualityDist> & pointer){
	size_t pos = s.find("(");
	std::string tmp;

	if(pos == std::string::npos)
		throw "Unable to understand distribution '" + s + "'!";

	//initialize appropriate function
	std::string type = s.substr(0, pos);
	s.erase(0, pos);
	if(type == "fixed")
		pointer = std::make_unique<TSimulatorQualityDist>(s, _randomGenerator);
	else if(type == "normal")
		pointer = std::make_unique<TSimulatorQualityDistNormal>(s, _randomGenerator);
	else if(type == "binned")
		pointer = std::make_unique<TSimulatorQualityDistBinned>(s, _randomGenerator);
	else if(type == "freq")
		pointer = std::make_unique<TSimulatorQualityDistFreq>(s, _randomGenerator);
	else throw "Unknown read quality distribution '" + type + "'!";

	checkInitialization();
};

void TSimulatorSingleEndRead::setQualityDistribution(std::string s){
	_initializeQualityDistribution(s, _qualityDist);
};

void TSimulatorSingleEndRead::setMappingQualityDistribution(std::string s){
	_initializeQualityDistribution(s, _mappingQualityDist);
};

void TSimulatorSingleEndRead::setQualityTransformation(GenotypeLikelihoods::TSequencingErrorModelsOneReadGroup const* RecalObject){
	_recalObject = RecalObject;
	_hasRecal = _recalObject->recalibrates();
};

void TSimulatorSingleEndRead::setPMD(GenotypeLikelihoods::TPMDType const* PmdObject){
	_pmdObject = PmdObject;
	_hasPMD = _pmdObject->hasDamage();
};

void TSimulatorSingleEndRead::_applyPMD(std::vector<Base>& bases, const TReadLength & readLength, const bool isReverseStrand){
	if(_hasPMD){
		if(!isReverseStrand){
			for(uint16_t p = 0; p < readLength.read; ++p){
				_pmdObject->simulatePMD(bases[p], p, readLength.fragment - p - 1, isReverseStrand, *_randomGenerator);
			}
		} else {
			for(uint16_t p = 0; p < readLength.read; ++p){
				_pmdObject->simulatePMD(bases[p], readLength.read - p - 1, readLength.fragment - readLength.read + p, isReverseStrand, *_randomGenerator);
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
	return coretools::str::toString(_readNamePrefix, _readXPos, ":", _readYPos);
};

/*
void TSimulatorSingleEndRead::_simulateBasesQualities(BAM::TAlignment & alignment,
													  std::vector<Base> haplotype,
													  const uint64_t pos,
													  const TReadLength & readLength,
													  const bool readIsContaminated,
													  TSimulatorQualityTransformation* qualityTransform){

	//set read length
	if(alignment.isReverseStrand()){
		alignment.setInsertSize(-readLength.fragment);
	} else {
		alignment.setInsertSize(readLength.fragment);
	}

	//simulate true bases
	std::vector<Base> bases;
	if(readIsContaminated){
		std::vector<Base>::iterator start = contaminationSource->reference().begin() + pos;
		bases.assign(start, start + readLength.read);
	} else {
		std::vector<Base>::iterator start = haplotype.begin() + pos;
		bases.assign(start, start + readLength.read);
	}

	_cigar.clear();
	_cigar.add('M', readLength.read);

	//simulate true qualities
	std::vector<PhredIntProbability> phredIntQualities;
	_qualityDist->sample(phredIntQualities);
	_alignment.setSequenceQualities(_cigar, bases, phredIntQualities);

	//apply PMD

	I'm here!!!!!!
	Write functions to pass full alignments to PMD to simulate PMD
	Then pass full alignment to recal to simulate bases according to the error rates


	_applyPMD(_bases, readLength, isReverse);

	//simulate qualities and errors
	qualityTransform->simulateQualitiesAndErrors(_bases, _phredIntQualities, readLength.read, isReverse);

	//copy bases and qualities
	std::string queryBases, qualities;
	for(int p=0; p<readLength.read; ++p){
		queryBases += _genoMap.baseToChar[_bases[p]];

		CHECK!!! What do we have? Qual or phredInt? Use qualMap!

		qualities += (char) _phredIntQualities[p] + 33);
	}

	//--------------------------
	 *
	//adjust qualities for writing
	_qualMap.adjustQualitiesForWriting(qualities);

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

*/

void TSimulatorSingleEndRead::simulate(std::vector<Base>, uint32_t refID, uint32_t pos, TSimulatorBamFile & bamFile){
	//prepare alignment
	_alignment.move(refID, pos);
	_alignment.setName(_getNextReadName());
	_alignment.setIsReverseStrand(_randomGenerator->getRand() < 0.5);

	//pick a fragment and read length, strand and contamination
	//TReadLength readLength = _readLengthDist->sample();
	//bool readIsContaminated = _isContaminated && _randomGenerator->getRand() < _contaminationRate;

	//simulated bases and qualities
	//_simulateBasesQualities(_alignment, haplotype, pos, readLength, readIsContaminated, _qualityTransform);

	//write bam alignment
	bamFile.saveAlignment(_alignment);
};

void TSimulatorSingleEndRead::printDetails(TLog* logfile){
	logfile->list(_type + ".");
	if(_readLengthDist)
		_readLengthDist->printDetails(logfile);
	else throw "Read length distribution not initialized!";

	if(_mappingQualityDist)
		_mappingQualityDist->printDetails(logfile, "Mapping quality");
	else throw "Mapping quality distribution not initialized!";

	if(_qualityDist)
		_qualityDist->printDetails(logfile, "Base quality");
	else throw "Read quality distribution not initialized!";

	if(_hasRecal){
		//TODO: add recal string output
		logfile->list("Recal: ");
	}

	if(_hasPMD){
		logfile->list("PMD: " + _pmdObject->functionString());
	} else {
		logfile->list("No PMD.");
	}

	if(_isContaminated)
		logfile->list("Contaminated with rate ", _contaminationRate, ".");
	else
		logfile->list("Read group is not contaminated.");
};

//----------------------------------
// TSimulatorPairedEndReads
//----------------------------------
TSimulatorPairedEndReads::TSimulatorPairedEndReads(const BAM::TReadGroup& ReadGroup, coretools::TRandomGenerator* RandomGenerator) : TSimulatorSingleEndRead(ReadGroup, RandomGenerator){
	_type = "paired-end";

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
	for(auto& it : bamAlignmentSecondMates_idle)
		delete it;
};

void TSimulatorPairedEndReads::simulate(std::vector<Base> haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile & bamFile){
	//pick a fragment, read length and contamination
	TReadLength readLength = _readLengthDist->sample();
	bool readIsContaminated = _isContaminated && _randomGenerator->getRand() < _contaminationRate;

	// Fill FIRST mate
	//------------------
	//simulated bases and qualities
	//_simulateBasesQualities(_alignment, haplotype, pos, readLength, readIsContaminated, false, _qualityTransform);

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
		secondMate->setReadGroup(_readGroup.id);
		secondMate->setMappingQuality(50);
		secondMate->setIsReverseStrand(true);
	}

	//simulated bases and qualities
	//_simulateBasesQualities(*secondMate, haplotype, matePos, readLength, readIsContaminated, true, qualityTransform_secondMate);

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

void TSimulatorPairedEndReads::writeUnwrittenAlignments(long pos, TSimulatorBamFile & bamFile){
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
