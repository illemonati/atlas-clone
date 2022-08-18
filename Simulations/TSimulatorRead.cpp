/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"
#include <algorithm>
#include <memory>

#include "PhredProbabilityTypes.h"
#include "TLog.h"
#include "TPostMortemDamage.h"
#include "TRandomGenerator.h"
#include "TSequencedBase.h"
#include "SequencingError/TModel.h"
#include "TSimulatorAuxiliaryTools.h"
#include "stringFunctions.h"

namespace Simulations {
using genometools::Base;
using genometools::PhredIntProbability;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using BAM::RGInfo::TReadGroupInfoEntry;
using BAM::RGInfo::InfoType;


//------------------------------------------------
// TSimulatorRead
//------------------------------------------------
TSimulatorRead::TSimulatorRead(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)
	: _readGroup(ReadGroup),
	  _readGroupInfo(RGInfo),
	  _fragmentLengthDistr(_readGroupInfo[InfoType::fragmentLengthDistr]),
	  _qualityDist(_readGroupInfo[InfoType::baseQualityDistr]),
	  _mappingQualityDist(_readGroupInfo[InfoType::mappingQualityDistr])
	{

	// initialize bamAlignment
	_alignment.setReadGroup(_readGroup.id);

	//readNamePrefix: "<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	_readNamePrefix = "ATL:0:A:1:" + coretools::str::toString(_readGroup.id) + ":";

	//soft clip
	std::string sc = _readGroupInfo[InfoType::softClipDistr];
	if(!sc.empty()){
		//check if one or two values are given
		if(sc.find(':') == std::string::npos){
			//one distribution for both
			if(sc != "-" && sc != "fixed(0)"){
				_softClipDist3 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc);
				_softClipDist5 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc);
			} else {
				std::string sc3 = coretools::str::extractBefore(sc, ":");
				sc.erase(0,1);
				_softClipDist3 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc3);
				_softClipDist5 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc);
			}
		}
	}
}

double TSimulatorRead::_calcMeanReadLength(const uint16_t maxLen) const {
	// if fragments are always shorter than _numcycles, return mean fragment length
	if(_fragmentLengthDistr.max() < maxLen){
		return _fragmentLengthDistr.mean();
	}

	// else: take into account that read length is always <= _numCycles
	double m = 0.0;
	double cumul = 0.0;
	for(uint16_t i = 1; i <= maxLen; ++i){
		double f = _fragmentLengthDistr.density(i);
		m += f * (double) i;
		cumul += f;
	}

	//remaining are all of lenth _numCycles
	m += (1. - cumul) * maxLen;
	return m;
}

std::string TSimulatorRead::_getNextReadName() {
	++_readXPos;
	if (_readXPos == 65536) {
		++_readYPos;
		_readXPos = 1;
	}
	return coretools::str::toString(_readNamePrefix, _readXPos, ":", _readYPos);
}

void TSimulatorRead::_simulateAlignmentDetails(uint32_t refID, uint32_t pos){
	_alignment.move(refID, pos);
	_alignment.setName(_getNextReadName());

	//simulate mapping quality
	_alignment.setMappingQuality(_mappingQualityDist.sample().get());
}

bool TSimulatorRead::_simulateContamination(){
	return _contaminationRate > 0. && randomGenerator().getRand() < _contaminationRate;
}

void TSimulatorRead::_addSoftclippedBases(std::vector<Base> & bases, const std::unique_ptr<TCategoricalDistribution<uint16_t>> & softClippedDist){
	if(softClippedDist){
		auto len = softClippedDist->sample();
		if(len > 0){
			for (size_t i = 0; i < len; i++){
				bases.push_back(static_cast<Base>(randomGenerator().getRand<uint8_t>(0,4)));
			}
			_cigar.add('S', len);
		}
	}
}

void TSimulatorRead::_simulateBasesQualities(BAM::TAlignment & alignment,
						      const std::vector<Base>& haplotype,
						      const uint64_t pos,
							  const uint16_t fragmentLength,
						      const uint16_t readLength,
						      bool readIsContaminated){

	//prepare vector of bases
	std::vector<Base> bases;

	// set read length
	if (alignment.isReverseStrand()) {
		alignment.setInsertSize(-fragmentLength);
		_addSoftclippedBases(bases, _softClipDist3);
	} else {
		alignment.setInsertSize(fragmentLength);
		_addSoftclippedBases(bases, _softClipDist5);
	}

	// simulate true bases
	const auto start = readIsContaminated ? _contaminationSource->reference().cbegin() + pos : haplotype.cbegin() + pos;
	bases.insert(bases.end(), start, start + std::min(fragmentLength, readLength));

	if (alignment.isReverseStrand()) {
		_addSoftclippedBases(bases, _softClipDist5);
	} else {
		_addSoftclippedBases(bases, _softClipDist3);
	}
	
	// simulate true qualities
	std::vector<genometools::PhredIntProbability> phredIntQualities(bases.size());
	_qualityDist.sample(phredIntQualities);

	_alignment.setSequenceQualities(_cigar, bases, phredIntQualities);

	for (auto & b : _alignment) {
		if (_pmd && _pmd->hasDamage()) _pmd->simulate(b);

		const auto sm = b.isSecondMate();
		_recal[sm]->simulate(b);
	}
}

void TSimulatorRead::setRecal(
	GenotypeLikelihoods::SequencingError::TModel const *Recal1, GenotypeLikelihoods::SequencingError::TModel const *Recal2) {
	_recal[0] = Recal1;
	_recal[1] = Recal2;
}

void TSimulatorRead::setPMD(GenotypeLikelihoods::TPMDType const *Pmd) {
	_pmd = Pmd;
}

void TSimulatorRead::setContamination(double rate, TSimulatorReference *source) {
	_contaminationRate  = rate;
	_contaminationSource = source;

	// check
	if (_contaminationRate < 0.0) throw "Contamination rate must be >= 0.0!";
	if (_contaminationRate > 1.0) throw "Contamination rate must be <= 0.0!";
}

void TSimulatorRead::printDetails(double frequency) {
	//TODO: complete with all information
	/*
	logfile().startIndent("Read group '", _readGroup.name_ID, "':");
	logfile().list("Type: ", type(), ".");
	logfile().list("Frequency: ", frequency, ".");

	_fragmentLengthDist.printDetails();
	_mappingQualityDist.printDetails("Mapping quality");
	_qualityDist.printDetails("Base quality");

	if ((_recal[0] && _recal[0]->recalibrates())) {
		// TODO: add recal string output
		logfile().list("Recal First Mate: " + _recal[0]->getCovariateDefinition());
	}

	if ((_recal[1] && _recal[1]->recalibrates())) {
		// TODO: add recal string output
		logfile().list("Recal Second Mate: " + _recal[0]->getCovariateDefinition());
	}

	if (_pmd && _pmd->hasDamage()) {
		logfile().list("PMD: " + _pmd->functionString());
	} else {
		logfile().list("No PMD.");
	}

	if (_contaminationRate > 0.)
		logfile().list("Contaminated with rate ", _contaminationRate, ".");
	else
		logfile().list("Read group is not contaminated.");
	*/
}

//----------------------------------
// TSimulatorSingleEndRead
//----------------------------------
TSimulatorSingleEndRead::TSimulatorSingleEndRead(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)
	: TSimulatorRead(ReadGroup, RGInfo){

	//num cycles
	coretools::str::convertString< coretools::StrictlyPositive<uint16_t> >(_readGroupInfo[InfoType::numCycles],
			BAM::RGInfo::infos[InfoType::numCycles].description + " must be within [1,65535].", _numCycles);
}

double TSimulatorSingleEndRead::meanReadLength() const {
	return _calcMeanReadLength(_numCycles);
}

void TSimulatorSingleEndRead::simulate(const std::vector<Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile) {
	// prepare alignment
	_simulateAlignmentDetails(refID, pos);
	_alignment.setIsReverseStrand(randomGenerator().getRand() < 0.5);

	// simulate read length
	uint16_t fragmentLength = _fragmentLengthDistr.sample();

	// simulated bases and qualities
	_simulateBasesQualities(_alignment, haplotype, pos, fragmentLength, _numCycles, _simulateContamination());

	// write bam alignment
	bamFile.saveAlignment(_alignment);
}

//----------------------------------
// TSimulatorPairedEndReads
//----------------------------------
TSimulatorPairedEndReads::TSimulatorPairedEndReads(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)
	: TSimulatorRead(ReadGroup, RGInfo){

	//num cycles
	if(coretools::str::stringContains(_readGroupInfo[InfoType::numCycles], ',')){
		//two values: one for first and one for second mate
		coretools::str::convertString< coretools::StrictlyPositive<uint16_t> >(coretools::str::readBefore(_readGroupInfo[InfoType::numCycles], ','),
				BAM::RGInfo::infos[InfoType::numCycles].description + " must be within [1,65535].", _numCycles[0]);
		coretools::str::convertString< coretools::StrictlyPositive<uint16_t> >(coretools::str::readAfter(_readGroupInfo[InfoType::numCycles], ','),
				BAM::RGInfo::infos[InfoType::numCycles].description + " must be within [1,65535].", _numCycles[1]);
	} else {
		//one value to be used for both mates
		coretools::str::convertString< coretools::StrictlyPositive<uint16_t> >(_readGroupInfo[InfoType::numCycles],
				BAM::RGInfo::infos[InfoType::numCycles].description + " must be within [1,65535].", _numCycles[0]);
		_numCycles[1] = _numCycles[0];
	}

	// set SAM flags
	_flags.setIsPaired(true);
	_flags.setIsProperPair(true);
	_flags.setIsRead1(true);
	_flags.setMateIsReverseStrand(false);

	// set SAM flags of second mate
	_mateFlags.setIsPaired(true);
	_mateFlags.setIsProperPair(true);
	_mateFlags.setIsRead2(true);
	_mateFlags.setIsReverseStrand(true);
}

double TSimulatorPairedEndReads::meanReadLength() const {
	return _calcMeanReadLength(_numCycles[0] + _numCycles[1]);
}

void TSimulatorPairedEndReads::simulate(const std::vector<Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile) {
	// pick a fragment
	uint16_t fragmentLength = _fragmentLengthDistr.sample();
	bool readIsContaminated = _simulateContamination();

	// Fill FIRST mate
	//------------------
	// prepare alignment
	_simulateAlignmentDetails(refID, pos);

	// simulated bases and qualities
	_simulateBasesQualities(_alignment, haplotype, pos, fragmentLength, _numCycles[0], readIsContaminated);

	// write bam alignment
	bamFile.saveAlignment(_alignment);

	// Fill SECOND mate
	//------------------
	// identify position
	uint32_t matePos;
	if(fragmentLength <= _numCycles[1]){
		matePos = pos;
	} else {
		matePos = pos + (uint32_t) fragmentLength - (uint32_t) _numCycles[1];
	}

	// create new alignment
	BAM::TAlignment secondMate;
	secondMate.setReadGroup(_readGroup.id);
	secondMate.move(refID, matePos);
	secondMate.setName(_alignment.name());
	secondMate.setMappingQuality(_alignment.mappingQuality());
	secondMate.setSamFlags(_mateFlags);

	// simulated bases and qualities
	_simulateBasesQualities(secondMate, haplotype, matePos, fragmentLength, _numCycles[1], readIsContaminated);

	// write if it starts at same position as first, and keep for writing later otherwiese
	if (matePos == pos) {
		bamFile.saveAlignment(secondMate);
	} else {
		bamAlignmentSecondMates.push_back(secondMate);
	}
}

void TSimulatorPairedEndReads::writeUnwrittenAlignments(long , TSimulatorBamFile &bamFile) {
	for (auto & a: bamAlignmentSecondMates) {
		bamFile.saveAlignment(a);
	}
	bamAlignmentSecondMates.clear();
}

} // namespace Simulations
