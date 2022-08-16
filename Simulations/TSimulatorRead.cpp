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
	  _fragmentLengthDist(_readGroupInfo.get(InfoType::fragmentLengthDistr)),
	  _qualityDist(_readGroupInfo.get(InfoType::baseQualityDistr)),
	  _mappingQualityDist(_readGroupInfo.get(InfoType::mappingQualityDistr))
	{

	// initialize bamAlignment
	_alignment.setReadGroup(_readGroup.id);

	//readNamePrefix: "<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	_readNamePrefix = "ATL:0:A:1:" + coretools::str::toString(_readGroup.id) + ":";

	//soft clip
	std::string sc = _readGroupInfo.get(InfoType::softClipDistr);
	if(!sc.empty()){
		//check if one or two values are given
		if(sc.find(':') == std::string::npos){
			//one distribution for both
			if(sc != "-" && sc != "fixed(0)"){
				_softClipDist3 = std::make_unique(sc);
				_softClipDist5 = std::make_unique(sc);
			} else {
				std::string sc3 = coretools::str::extractBefore(sc, ":");
				sc.erase(0,1);
				_softClipDist3 = std::make_unique(sc3);
				_softClipDist5 = std::make_unique(sc);
			}
		}
	}
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
						      const TReadAndFragmentLength & readLength,
						      bool readIsContaminated/*,
										     TSimulatorQualityTransformation* qualityTransform*/){

	//prepare vector of bases
	std::vector<Base> bases;

	// set read length
	if (alignment.isReverseStrand()) {
		alignment.setInsertSize(-readLength.fragment);
		_addSoftclippedBases(bases, _softClipDist3);
	} else {
		alignment.setInsertSize(readLength.fragment);
		_addSoftclippedBases(bases, _softClipDist5);
	}

	// simulate true bases
	const auto start = readIsContaminated ? _contaminationSource->reference().cbegin() + pos : haplotype.cbegin() + pos;
	bases.insert(bases.end(), start, start + readLength.read);

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

std::string TSimulatorRead::_getNextReadName() {
	++_readXPos;
	if (_readXPos == 65536) {
		++_readYPos;
		_readXPos = 1;
	}
	return coretools::str::toString(_readNamePrefix, _readXPos, ":", _readYPos);
}


//----------------------------------
// TSimulatorSingleEndRead
//----------------------------------
TSimulatorSingleEndRead::TSimulatorSingleEndRead(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)
	: TSimulatorRead(ReadGroup, RGInfo){

	//num cycles
	_numCycles = coretools::str::convertString< coretools::StrictlyPositive<uint16_t> >(_readGroupInfo.get(InfoType::numCycles),
			BAM::RGInfo::infoType2Description(InfoType::numCycles) + " must be within [1,65535].", _numCycles);
}



void TSimulatorSingleEndRead::simulate(const std::vector<Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile) {
	// prepare alignment
	_alignment.move(refID, pos);
	_alignment.setName(_getNextReadName());
	_alignment.setIsReverseStrand(randomGenerator().getRand() < 0.5);

	// pick a fragment and read length, strand and contamination
	TReadAndFragmentLength readLength = _fragmentLengthDist.sample();
	bool readIsContaminated = _contaminationRate > 0. && randomGenerator().getRand() < _contaminationRate;

	// simulated bases and qualities
	_simulateBasesQualities(_alignment, haplotype, pos, readLength, readIsContaminated/*, _qualityTransform*/);

	//set mapping quality
	_alignment.setMappingQuality(static_cast<uint8_t>(_mappingQualityDist.sample()));

	// write bam alignment
	bamFile.saveAlignment(_alignment);
}

void TSimulatorSingleEndRead::printDetails(double frequency) {
	//TODO: complete with all information
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
}

//----------------------------------
// TSimulatorPairedEndReads
//----------------------------------
TSimulatorPairedEndReads::TSimulatorPairedEndReads(
		const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)

		const BAM::TReadGroup &, const uint16_t NumCyclesFirst, const uint16_t NumCyclesSecond) : TSimulatorSingleEndRead(ReadGroup, NumCclesFirst){
	_numCyclesSecond = NumCyclesSecond;
	// set SAM flags
	_flags.setIsPaired(true);
	_flags.setIsProperPair(true);
	_flags.setIsRead1(true);
	_flags.setMateIsReverseStrand(true);

	// set SAM flags of second mate
	_mateFlags.setIsPaired(true);
	_mateFlags.setIsProperPair(true);
	_mateFlags.setIsRead2(true);
	_mateFlags.setIsReverseStrand(true);
}

void TSimulatorPairedEndReads::simulate(
		const std::vector<Base>& /*haplotype*/, uint32_t refID, uint32_t pos,
					TSimulatorBamFile &bamFile) {
	// pick a fragment, read length and contamination
	TReadAndFragmentLength readLength  = _fragmentLengthDist.sample();
	//bool readIsContaminated = (_contaminationRate > 0.) && randomGenerator().getRand() < _contaminationRate;

	// Fill FIRST mate
	//------------------
	// simulated bases and qualities
	//_simulateBasesQualities(_alignment, haplotype, pos, readLength, readIsContaminated, false, _qualityTransform);

	// write bam alignment
	_alignment.move(refID, pos);
	_alignment.setName(_getNextReadName());
	bamFile.saveAlignment(_alignment);

	// Fill SECOND mate
	//------------------
	// identify position
	uint32_t matePos = pos + readLength.diff();

	// create alignment

	BAM::TAlignment secondMate;
	// create new
	secondMate.setReadGroup(_readGroup.id);
	secondMate.setMappingQuality(50);
	secondMate.setIsReverseStrand(true);

	// simulated bases and qualities
	//_simulateBasesQualities(*secondMate, haplotype, matePos, readLength, readIsContaminated, true,
	//qualityTransform_secondMate);

	// fill bam alignment
	secondMate.move(refID, matePos);
	secondMate.setName(_alignment.name());

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
