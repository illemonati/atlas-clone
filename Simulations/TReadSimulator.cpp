/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include <TReadSimulator.h>
#include <algorithm>
#include <memory>

#include "genometools/PhredProbabilityTypes.h"
#include "TPostMortemDamage.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSequencedBase.h"
#include "SequencingError/TModel.h"
#include "TSimulatorAuxiliaryTools.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/fromString.h"

namespace Simulations {
using genometools::Base;
using genometools::PhredIntProbability;
using genometools::TGenomePosition;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using BAM::RGInfo::TReadGroupInfoEntry;
using BAM::RGInfo::InfoType;

//------------------------------------------------
// TSimulatorRead
//------------------------------------------------

TReadSimulator::TReadSimulator(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)
	: _readGroup(ReadGroup), _recalModels(RGInfo) {

	// initialize bamAlignment
	_alignment.setReadGroup(_readGroup.id());

	//readNamePrefix: "<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	_readNamePrefix = "ATL:0:A:1:" + coretools::str::toString(_readGroup.id()) + ":";

	//initialize distributions
	_initDistribution(_fragmentLengthDistr, RGInfo, InfoType::fragmentLength);
	_initDistribution(_mappingQualityDist, RGInfo, InfoType::mappingQuality);
	_initDistribution(_qualityDist, RGInfo, InfoType::baseQuality);

	//soft clip
	logfile().listFlush(BAM::RGInfo::infos[InfoType::softClipping].description, ": ");
	if(RGInfo.has(InfoType::softClipping)){
		std::string sc = RGInfo.getString(InfoType::softClipping);
		logfile().write(sc);

		if (!sc.empty() && sc != "-") {
			// check if one or two values are given
			if (sc.find(':') == std::string::npos) {
				// one distribution for both
				_softClipDist3 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc);
				_softClipDist5 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc);
			} else {
				std::string sc3 = coretools::str::extractBefore(sc, ":");
				sc.erase(0, 1);
				_softClipDist3 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc3);
				_softClipDist5 = std::make_unique<TCategoricalDistribution<uint16_t>>(sc);
			}
		}
	} else {
		logfile().write("none");
	}
}

double TReadSimulator::_calcMeanReadLength(const uint16_t maxLen) const {
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

std::string TReadSimulator::_getNextReadName() {
	++_readXPos;
	if (_readXPos == 65536) {
		++_readYPos;
		_readXPos = 1;
	}
	return coretools::str::toString(_readNamePrefix, _readXPos, ":", _readYPos);
}

void TReadSimulator::_simulateAlignmentDetails(const TGenomePosition & Position){	;
	_alignment.move(Position);
	_alignment.setName(_getNextReadName());

	//simulate mapping quality
	_alignment.setMappingQuality(_mappingQualityDist.sample());
}

bool TReadSimulator::_simulateContamination(){
	return _contaminationRate > 0. && randomGenerator().getRand() < _contaminationRate;
}

void TReadSimulator::_addSoftclippedBases(std::vector<Base> & Bases, const std::unique_ptr<TCategoricalDistribution<uint16_t>> & SoftClippedDist, BAM::TCigar & Cigar){
	if(SoftClippedDist){
		auto len = SoftClippedDist->sample();
		if(len > 0){
			for (size_t i = 0; i < len; i++){
				Bases.push_back(static_cast<Base>(randomGenerator().getRand<uint8_t>(0,4)));
			}
			Cigar.add('S', len);
		}
	}
}

void TReadSimulator::_simulateBasesQualities(BAM::TAlignment & alignment,
						      const std::vector<Base>& haplotype,
							  const uint16_t fragmentLength,
						      const uint16_t readLength,
						      bool readIsContaminated){

	//prepare vector of bases
	std::vector<Base> bases;
	BAM::TCigar cigar;

	// set read length
	if (alignment.isReverseStrand()) {
		alignment.setInsertSize(-fragmentLength);
		_addSoftclippedBases(bases, _softClipDist3, cigar);
	} else {
		alignment.setInsertSize(fragmentLength);
		_addSoftclippedBases(bases, _softClipDist5, cigar);
	}

	// simulate true bases
	const auto start = readIsContaminated ? _contaminationSource->reference().cbegin() + alignment.position() : haplotype.cbegin() + alignment.position();
	auto len = std::min(fragmentLength, readLength);
	bases.insert(bases.end(), start, start + len);
	cigar.add('M', len);

	if (alignment.isReverseStrand()) {
		_addSoftclippedBases(bases, _softClipDist5, cigar);
	} else {
		_addSoftclippedBases(bases, _softClipDist3, cigar);
	}
	
	// simulate true qualities
	std::vector<genometools::PhredIntProbability> phredIntQualities(bases.size());
	_qualityDist.sample(phredIntQualities);

	alignment.setSequenceQualities(cigar, bases, phredIntQualities);

	//add PMD
	if (_pmd && _pmd->hasDamage()){
		for (auto & b : alignment) {
			_pmd->simulate(b);
		}
	}

	//add recal
	_recalModels.simulate(alignment);
}

void TReadSimulator::setPMD(GenotypeLikelihoods::TPMDType const *Pmd) {
	_pmd = Pmd;
}

void TReadSimulator::setContamination(double rate, TSimulatorReference *source) {
	_contaminationRate  = rate;
	_contaminationSource = source;

	// check
	if (_contaminationRate < 0.0) UERROR("Contamination rate must be >= 0.0!");
	if (_contaminationRate > 1.0) UERROR("Contamination rate must be <= 0.0!");
}

//----------------------------------
// TSimulatorSingleEndRead
//----------------------------------
TReadSimulatorSingleEnd::TReadSimulatorSingleEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)
	: TReadSimulator(ReadGroup, RGInfo) {

	_alignment.setSamFlags(_flags);

	//num cycles
	logfile().list(BAM::RGInfo::infos[InfoType::cycles].description, ": ", RGInfo.getString(InfoType::cycles));
	std::string error = "For single-end read groups, " + BAM::RGInfo::infos[InfoType::cycles].description + " must be a single integer within [1,65535].";
	auto& json = RGInfo[InfoType::cycles];

	std::string err = "Unable to understand " + BAM::RGInfo::infos[InfoType::cycles].description + ": ";
	std::string errRange = err + "expect a single integer within [1,65535].";

	if(json.is_number()){
		if(json.get<int>() < 1 || json.get<int>() > 65535){
			UERROR(errRange);
		}
		_numCycles = json.get<int>();
	} else if(json.is_array() && json.size() != 1){
		UERROR(errRange);
	} else if(json.is_string()){
		coretools::str::fromString(json.get<std::string>(), _numCycles, err);
	} else if(json.is_array() && json.size() == 1){
		coretools::str::fromString(json[0].get<std::string>(), _numCycles, err);
	} else {
		UERROR(errRange);
	}
}

double TReadSimulatorSingleEnd::meanReadLength() const {
	return _calcMeanReadLength(_numCycles);
}

void TReadSimulatorSingleEnd::simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype, BAM::TOutputBamFile & BamFile) {
	// pick a fragment
	const auto fragmentLength = _fragmentLengthDistr.sample();

	// prepare alignment
	_simulateAlignmentDetails(Position );
	_alignment.setIsReverseStrand(randomGenerator().getRand() < 0.5);

	// simulated bases and qualities
	_simulateBasesQualities(_alignment, Haplotype, fragmentLength, _numCycles, _simulateContamination());

	// write bam alignment
	BamFile.writeAlignment(_alignment);
}

/**
void TReadSimulatorSingleEnd::simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype, Simulations::TFASTQWriter &FASTQFile ) {
	// pick a fragment
	const auto fragmentLength = _fragmentLengthDistr.sample();

	// prepare alignment
	_simulateAlignmentDetails(Position);
	_alignment.setIsReverseStrand(randomGenerator().getRand() < 0.5);

	// simulated bases and qualities
	_simulateBasesQualities(_alignment, Haplotype, fragmentLength, _numCycles, _simulateContamination());

	// write bam alignment
	FASTQFile.writeAlignment(_alignment);
}
*/


//----------------------------------
// TSimulatorPairedEndReads
//----------------------------------
TReadSimulatorPairedEnd::TReadSimulatorPairedEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo)
	: TReadSimulator(ReadGroup, RGInfo){
	//num cycles
	logfile().list(BAM::RGInfo::infos[InfoType::cycles].description, ": ", RGInfo[InfoType::cycles]);
	auto& json = RGInfo[InfoType::cycles];

	std::string err = "Unable to understand " + BAM::RGInfo::infos[InfoType::cycles].description + ": ";
	std::string errRange = err + "expect one or two integers within [1,65535].";

	//TODO: probably need some json parsing functions to simplify this!
	if(json.is_array()){
		if(json.size() == 1){
			if(json[0].get<int>() < 1 || json[0].get<int>() > 65535){
				UERROR(errRange);
			}
			_numCycles[0] = json[0].get<int>();
			_numCycles[1] = _numCycles[0];
		} else if(json.size() == 2){
			if(json[0].get<int>() < 1 || json[1].get<int>() > 65535 || json[1].get<int>() < 1 || json[1].get<int>() > 65535){
				UERROR(errRange);
			}
			_numCycles[0] = json[0].get<int>();
			_numCycles[1] = json[1].get<int>();
		} else {
			UERROR(errRange);
		}
	} else if(json.is_number()){
		if(json[0].get<int>() < 1 || json[0].get<int>() > 65535){
			UERROR(errRange);
		}
		_numCycles[0] = json[0].get<int>();
		_numCycles[1] = _numCycles[0];
	} else if(json.is_string()){
		std::string ss = json.get<std::string>();
		if(coretools::str::stringContains(ss, ',')){
			//two values: one for first and one for second mate
			coretools::str::fromString(coretools::str::readBefore(ss, ','), _numCycles[0], err);
			coretools::str::fromString(coretools::str::readAfter(ss, ','), _numCycles[1], err);
		} else {
			//one value to be used for both mates
			coretools::str::fromString(ss, _numCycles[0], err);
			_numCycles[1] = _numCycles[0];
		}
	} else {
		UERROR(errRange);
	}

	// set SAM flags
	_flags.setIsPaired(true);
	_flags.setIsProperPair(true);
	_flags.setIsRead1(true);
	_flags.setMateIsReverseStrand(false);
	_alignment.setSamFlags(_flags);

	// set SAM flags of second mate
	_mateFlags.setIsPaired(true);
	_mateFlags.setIsProperPair(true);
	_mateFlags.setIsRead2(true);
	_mateFlags.setIsReverseStrand(true);
	_secondMate.setSamFlags(_mateFlags);
}

double TReadSimulatorPairedEnd::meanReadLength() const {
	return _calcMeanReadLength(_numCycles[0] + _numCycles[1]);
}

void TReadSimulatorPairedEnd::simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype, BAM::TOutputBamFile & BamFile) {
	// pick a fragment
	const auto fragmentLength     = _fragmentLengthDistr.sample();
	const auto readIsContaminated = _simulateContamination();

	// Fill FIRST mate

	// prepare alignment
	_simulateAlignmentDetails(Position);
	_simulateBasesQualities(_alignment, Haplotype, fragmentLength, _numCycles[0], readIsContaminated);

	// Fill SECOND mate

	// identify position
	_secondMate.move(_alignment);
	if(fragmentLength > _numCycles[1]){
		_secondMate += (uint32_t) fragmentLength - (uint32_t) _numCycles[1];
	}

	// create new alignment
	_secondMate.setReadGroup(_readGroup.id());
	_secondMate.setName(_alignment.name());
	_secondMate.setMappingQuality(_alignment.mappingQuality());

	// simulated bases and qualities
	_simulateBasesQualities(_secondMate, Haplotype, fragmentLength, _numCycles[1], readIsContaminated);

	// WRITE ALIGNMENTS
	//-----------------
	//set mate positions
	_alignment.setMateGenomicPosition(_secondMate);
	_secondMate.setMateGenomicPosition(_alignment);

	BamFile.writeAlignment(_alignment);

	// write mate if it starts at same position as first, and keep for writing later otherwise
	if (_secondMate == _alignment) {
		BamFile.writeAlignment(_secondMate);
	} else {
		BamFile.writeAlignmentLater(_secondMate);
	}
}

} // namespace Simulations
