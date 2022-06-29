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
#include "TSequencingErrorModel.h"
#include "TSimulatorAuxiliaryTools.h"
#include "stringFunctions.h"

namespace Simulations {
using genometools::Base;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;

//----------------------------------
// TSimulatorSingleEndRead
//----------------------------------
TSimulatorSingleEndRead::TSimulatorSingleEndRead(const BAM::TReadGroup &ReadGroup)
	: _readGroup(ReadGroup),
	  _readNamePrefix("ATL:0:A:1:" + coretools::str::toString(_readGroup.id) + ":") {
	//readNamePrefix: "<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	// initialize bamAlignment
	_alignment.setReadGroup(_readGroup.id);
}

bool TSimulatorSingleEndRead::checkInitialization() {
	return _readLengthDist && _qualityDist && _mappingQualityDist;
}

void TSimulatorSingleEndRead::setReadLengthDistribution(std::string s) {
	const auto pos = s.find("(");
	std::string tmp;

	if (pos == std::string::npos) throw "Unable to understand read length distribution '" + s + "'!";

	// initialize appropriate function
	const auto type = s.substr(0, pos);
	s.erase(0, pos);

	if (type == "gamma")
		_readLengthDist = std::make_unique<TSimulatorReadLengthGamma>(s);
	else if (type == "gammaMode") {
		_readLengthDist = std::make_unique<TSimulatorReadLengthGammaMode>(s);
	} else if (type == "fixed")
		_readLengthDist = std::make_unique<TReadLengthDistribution>(s);
	else
		throw "Unknown read length distribution '" + type + "'!";
}

std::unique_ptr<TSimulatorQualityDist> TSimulatorSingleEndRead::_initializeQualityDistribution(std::string s) {
	const auto pos = s.find("(");
	std::string tmp;
	if (pos == std::string::npos) throw "Unable to understand distribution '" + s + "'!";

	// initialize appropriate function
	const auto type = s.substr(0, pos);
	s.erase(0, pos);
	if (type == "fixed")
		return std::make_unique<TSimulatorQualityDistFixed>(s);
	else if (type == "normal")
		return std::make_unique<TSimulatorQualityDistNormal>(s);
	else if (type == "binned")
		return std::make_unique<TSimulatorQualityDistBinned>(s);
	else if (type == "freq")
		return std::make_unique<TSimulatorQualityDistFreq>(s);
	else
		throw "Unknown read quality distribution '" + type + "'!";
}

std::unique_ptr<TSimulatorSoftClipDist> TSimulatorSingleEndRead::_initializeSoftClipDistribution(std::string s) {
	const auto pos = s.find("(");
	std::string tmp;
	if (pos == std::string::npos) throw "Unable to understand distribution '" + s + "'!";

	// initialize appropriate function
	const auto type = s.substr(0, pos);
	s.erase(0, pos);
	if (type == "fixed")
		return std::make_unique<TSimulatorSoftClipDistFixed>(s);
	else if (type == "binned")
		return std::make_unique<TSimulatorSoftClipDistBinned>(s);
	else if (type == "freq")
		return std::make_unique<TSimulatorSoftClipDistFreq>(s);
	/*else if (type == "poisson")
		return std::make_unique<TSimulatorSoftClipDistPoisson>(s);*/
	else
		throw "Unknown read quality distribution '" + type + "'!";
}

void TSimulatorSingleEndRead::setQualityDistribution(std::string s) {
	_qualityDist = _initializeQualityDistribution(s);
}

void TSimulatorSingleEndRead::setMappingQualityDistribution(std::string s) {
	_mappingQualityDist = _initializeQualityDistribution(s);
}

void TSimulatorSingleEndRead::setSoftClipDistribution(std::string s) {
	_softClipDist = _initializeSoftClipDistribution(s);
}

void TSimulatorSingleEndRead::setRecal(
	GenotypeLikelihoods::SequencingError::TModel const *Recal1, GenotypeLikelihoods::SequencingError::TModel const *Recal2) {
	_recal[0] = Recal1;
	_recal[1] = Recal2;
}

void TSimulatorSingleEndRead::setPMD(GenotypeLikelihoods::TPMDType const *Pmd) {
	_pmd = Pmd;
}

void TSimulatorSingleEndRead::setContamination(double rate, TSimulatorReference *source) {
	_contaminationRate  = rate;
	_contaminationSource = source;

	// check
	if (_contaminationRate < 0.0) throw "Contamination rate must be >= 0.0!";
	if (_contaminationRate > 1.0) throw "Contamination rate must be <= 0.0!";
}

std::string TSimulatorSingleEndRead::_getNextReadName() {
	++_readXPos;
	if (_readXPos == 65536) {
		++_readYPos;
		_readXPos = 1;
	}
	return coretools::str::toString(_readNamePrefix, _readXPos, ":", _readYPos);
}

void TSimulatorSingleEndRead::_simulateBasesQualities(BAM::TAlignment & alignment,
						      const std::vector<Base>& haplotype,
						      const uint64_t pos,
						      const TReadLength & readLength,
						      bool readIsContaminated/*,
										     TSimulatorQualityTransformation* qualityTransform*/){

	// set read length
	if (alignment.isReverseStrand()) {
		alignment.setInsertSize(-readLength.fragment);
	} else {
		alignment.setInsertSize(readLength.fragment);
	}

	// simulate true bases
	uint16_t sClippedLength5 = static_cast<uint16_t>(_softClipDist->sample());
	uint16_t sClippedLength3 = static_cast<uint16_t>(_softClipDist->sample());
	if ((sClippedLength5 + sClippedLength3) >= readLength.read)
		throw "Number of softclipped bases either equal or exceed read length. Either increase read length or decrease softclipped length.";
	uint16_t mappedLength = readLength.read - sClippedLength5 - sClippedLength3;
	std::vector<Base> bases;

	const auto start = readIsContaminated ? _contaminationSource->reference().cbegin() + pos : haplotype.cbegin() + pos;

	for (size_t i = 1; i <= sClippedLength5; i++) bases.push_back(static_cast<Base>(randomGenerator().getRand<uint8_t>(0,4)));
	bases.insert(bases.end(), start, start + mappedLength);
	for (size_t i = 1; i <= sClippedLength3; i++) bases.push_back(static_cast<Base>(randomGenerator().getRand<uint8_t>(0,4)));

	_cigar.clear();
	if (sClippedLength5 > 0) _cigar.add('S', sClippedLength5);
	_cigar.add('M', mappedLength);
	if (sClippedLength3 > 0) _cigar.add('S', sClippedLength3);
	

	// simulate true qualities
	std::vector<genometools::PhredIntProbability> phredIntQualities(bases.size());
	_qualityDist->sample(phredIntQualities);

	_alignment.setSequenceQualities(_cigar, bases, phredIntQualities);

	for (auto & b : _alignment) {
		if (_pmd && _pmd->hasDamage()) _pmd->simulate(b);

		const auto sm = b.isSecondMate();
		_recal[sm]->simulate(b);
	}
}

void TSimulatorSingleEndRead::simulate(const std::vector<Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile) {
	// prepare alignment
	_alignment.move(refID, pos);
	_alignment.setName(_getNextReadName());
	_alignment.setIsReverseStrand(randomGenerator().getRand() < 0.5);

	// pick a fragment and read length, strand and contamination
	TReadLength readLength = _readLengthDist->sample();
	bool readIsContaminated = _contaminationRate > 0. && randomGenerator().getRand() < _contaminationRate;

	// simulated bases and qualities
	_simulateBasesQualities(_alignment, haplotype, pos, readLength, readIsContaminated/*, _qualityTransform*/);

	//set mapping quality
	_alignment.setMappingQuality(static_cast<uint8_t>(_mappingQualityDist->sample()));

	// write bam alignment
	bamFile.saveAlignment(_alignment);
}

void TSimulatorSingleEndRead::printDetails() {
	logfile().list(type() + ".");
	if (_readLengthDist)
		_readLengthDist->printDetails();
	else
		throw "Read length distribution not initialized!";

	if (_mappingQualityDist)
		_mappingQualityDist->printDetails("Mapping quality");
	else
		throw "Mapping quality distribution not initialized!";

	if (_qualityDist)
		_qualityDist->printDetails("Base quality");
	else
		throw "Read quality distribution not initialized!";

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
TSimulatorPairedEndReads::TSimulatorPairedEndReads(const BAM::TReadGroup &ReadGroup) : TSimulatorSingleEndRead(ReadGroup){
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

	void TSimulatorPairedEndReads::simulate(const std::vector<Base>& /*haplotype*/, uint32_t refID, uint32_t pos,
					TSimulatorBamFile &bamFile) {
	// pick a fragment, read length and contamination
	TReadLength readLength  = _readLengthDist->sample();
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
