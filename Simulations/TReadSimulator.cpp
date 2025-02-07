/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include <TReadSimulator.h>

#include "PMD/TModel.h"
#include "TOutputBamFile.h"
#include "TSimulatorReference.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"

namespace Simulations {
using BAM::RGInfo::InfoType;
using BAM::RGInfo::TReadGroupInfoEntry;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using coretools::instances::parameters;
using coretools::probdist::TCategoricalDistribution;
using coretools::P;
using genometools::Base;
using genometools::TGenomePosition;

//------------------------------------------------
// TSimulatorRead
//------------------------------------------------
TReadSimulator::TReadSimulator(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal)
		: _readGroup(&ReadGroup), _pmd(&Pmd), _recal(Recal) {

	// initialize bamAlignment
	_alignment.setReadGroup(_readGroup->id);

	//readNamePrefix: "<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	_readNamePrefix = "ATL:0:A:1:" + coretools::str::toString(_readGroup->id) + ":";

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
				_softClipDist3 = std::make_unique<TCategoricalDistribution<size_t>>(sc);
				_softClipDist5 = std::make_unique<TCategoricalDistribution<size_t>>(sc);
			} else {
				std::string sc3 = coretools::str::extractBefore(sc, ":");
				sc.erase(0, 1);
				_softClipDist3 = std::make_unique<TCategoricalDistribution<size_t>>(sc3);
				_softClipDist5 = std::make_unique<TCategoricalDistribution<size_t>>(sc);
			}
		}
	} else {
		logfile().write("none");
	}

	//duplication rate

	if(RGInfo.has(InfoType::duplicationRate)){
		coretools::str::fromString(RGInfo.getString(InfoType::duplicationRate), _duplicationRate, "duplication rate is not within [0,1]!");
		logfile().list(BAM::RGInfo::infos[InfoType::duplicationRate].description, ": ", _duplicationRate);
		if(_duplicationRate > 0.5) UERROR("Duplication rate must be within [0.0, 0.5]!");
		_duplicationRateAmongSimulated = P(_duplicationRate / (_duplicationRate.complement()));
	}

	if(parameters().exists("baseN")){
		_baseN = parameters().get<coretools::Probability>("baseN");
		logfile().list("Will simulate reads with base = N probability = ", _baseN, ". (parameter 'baseN')");
	} else {
		_baseN = P(0.001);
		logfile().list("Will simulate reads with base = N probability = ", _baseN, ". (set with 'baseN')");
	}
}

double TReadSimulator::_calcMeanReadLength(size_t maxLen) const {
	// if fragments are always shorter than _numcycles, return mean fragment length
	if(_fragmentLengthDistr.max() < maxLen){
		return _fragmentLengthDistr.mean();
	}

	// else: take into account that read length is always <= _numCycles
	double m = 0.0;
	double cumul = 0.0;
	for (size_t i = 1; i <= maxLen; ++i){
		const double f = _fragmentLengthDistr.density(i);
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

void TReadSimulator::_addSoftclippedBases(std::vector<Base> & Bases,const size_t &softClipLength, BAM::TCigar & Cigar){
	if(softClipLength > 0){
		for (size_t i = 0; i < softClipLength; i++){
			Bases.push_back(static_cast<Base>(randomGenerator().getRand<uint8_t>(0,4)));
		}
		Cigar.add('S', softClipLength);
	}
}

void TReadSimulator::_simulateBasesQualities(BAM::TAlignment &alignment, const std::vector<Base> &haplotype,
											 size_t fragmentLength, size_t readLength, bool readIsContaminated) {
	//prepare vector of bases
	std::vector<Base> bases;
	BAM::TCigar cigar;

	//sample softclip lengths
	const auto softClipLength3 = _softClipDist3 ? _softClipDist3->sample() : 0;
	const auto softClipLength5 = _softClipDist5 ? _softClipDist5->sample() : 0;

	// set read length
	if (alignment.isReverseStrand()) {
		alignment.setInsertSize(-fragmentLength);
		_addSoftclippedBases(bases, softClipLength3, cigar);
	} else {
		alignment.setInsertSize(fragmentLength);
		_addSoftclippedBases(bases, softClipLength5, cigar);
	}

	// simulate true bases
	const auto start = readIsContaminated ? _contaminationSource->reference().cbegin() + alignment.position() : haplotype.cbegin() + alignment.position();
	auto len = std::min(fragmentLength, readLength);
	assert(len > 0);
	bases.insert(bases.end(), start, start + len);
	cigar.add('M', len);

	if (alignment.isReverseStrand()) {
		_addSoftclippedBases(bases, softClipLength5, cigar);
	} else {
		_addSoftclippedBases(bases, softClipLength3, cigar);
	}

	// make some bases N
	for (auto& b: bases) {
		if (randomGenerator().getRand() < _baseN) b = Base::N;
	}
	
	// simulate true qualities
	std::vector<coretools::PhredInt> phredIntQualities(bases.size());
	_qualityDist.sample(phredIntQualities);

	alignment.setSequenceQualities(cigar, bases, phredIntQualities);

	_pmd->simulate(alignment);
	_recal[alignment.mate()]->simulate(alignment);
}

void TReadSimulator::setPMD(GenotypeLikelihoods::PMD::TModel const *Pmd) {
	_pmd = Pmd;
}

void TReadSimulator::setContamination(double rate, TSimulatorReference *source) {
	_contaminationRate  = rate;
	_contaminationSource = source;

	// check
	if (_contaminationRate < 0.0) UERROR("Contamination rate must be >= 0.0!");
	if (_contaminationRate > 1.0) UERROR("Contamination rate must be <= 0.0!");
}

size_t TReadSimulator::simulate(const TGenomePosition &Position, const std::vector<Base> &Haplotype,
							  BAM::TOutputBamFile &BamFile) {
	// Do not simulate fraction of reads that will be duplicates
	if (_duplicationRate == 0.0) {
		_simulate(Position, Haplotype);
		_writeSimulatedAlignments(BamFile);
		return 1;
	} else if (randomGenerator().getRand() > _duplicationRate) {
		_simulate(Position, Haplotype);
		_writeSimulatedAlignments(BamFile);

		if (randomGenerator().getRand() < _duplicationRateAmongSimulated) {
			_writeSimulatedAlignments(BamFile);
			return 2;
		}
		return 1;
	}
	return 0;
}

//----------------------------------
// TSimulatorSingleEndRead
//----------------------------------
	TReadSimulatorSingleEnd::TReadSimulatorSingleEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal)
		: TReadSimulator(ReadGroup, RGInfo, Pmd, Recal) {

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

void TReadSimulatorSingleEnd::_simulate(const TGenomePosition &Position, const std::vector<Base> &Haplotype) {
	// pick a fragment
	const auto fragmentLength = _fragmentLengthDistr.sample();

	// prepare alignment
	_simulateAlignmentDetails(Position);
	_alignment.setIsReverseStrand(randomGenerator().getRand() < 0.5);

	// simulated bases and qualities
	_simulateBasesQualities(_alignment, Haplotype, fragmentLength, _numCycles, _simulateContamination());
}

void TReadSimulatorSingleEnd::_writeSimulatedAlignments(BAM::TOutputBamFile & BamFile){
	assert(_alignment.size() > 0);
	BamFile.writeAlignment(_alignment);
}

//----------------------------------
// TSimulatorPairedEndReads
//----------------------------------
	TReadSimulatorPairedEnd::TReadSimulatorPairedEnd(const BAM::TReadGroup & ReadGroup, const TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal)
		: TReadSimulator(ReadGroup, RGInfo, Pmd, Recal){
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

	// set initial flags
	_alignment.setIsPaired(true);
	_alignment.setIsProperPair(true);
	_alignment.setIsRead1(true);
	_alignment.setIsReverseStrand(false);

	_secondMate.setIsPaired(true);
	_secondMate.setIsProperPair(true);
	_secondMate.setIsRead2(true);
	_secondMate.setIsReverseStrand(true);
}

double TReadSimulatorPairedEnd::meanReadLength() const {
	return _calcMeanReadLength(_numCycles[0] + _numCycles[1]);
}

void TReadSimulatorPairedEnd::_writeSimulatedAlignments(BAM::TOutputBamFile & BamFile){
	BamFile.writeAlignment(_alignment);

	// write mate if it starts at same position as first, and keep for writing later otherwise
	if (_secondMate == _alignment) {
		BamFile.writeAlignment(_secondMate);
	} else {
		BamFile.writeAlignmentLater(_secondMate);
	}
}

void TReadSimulatorPairedEnd::_simulate(const TGenomePosition & Position, const std::vector<Base> & Haplotype) {
	// pick a fragment
	const auto fragmentLength     = _fragmentLengthDistr.sample();
	const auto readIsContaminated = _simulateContamination();

	// first mate
	_simulateAlignmentDetails(Position);
	_alignment.setIsReverseStrand(randomGenerator().getRand() < 0.5);

	_simulateBasesQualities(_alignment, Haplotype, fragmentLength, _numCycles[0], readIsContaminated);

	// second mate
	// identify position
	_secondMate.move(_alignment);
	if(fragmentLength > _numCycles[1]){
		_secondMate += (size_t) fragmentLength - (size_t) _numCycles[1];
	}

	// create new alignment
	_secondMate.setReadGroup(_readGroup->id);
	_secondMate.setName(_alignment.name());
	_secondMate.setMappingQuality(_alignment.mappingQuality());
	_secondMate.setIsReverseStrand(!_alignment.isReverseStrand());
	assert(_alignment.isReverseStrand() != _secondMate.isReverseStrand());

	// simulated bases and qualities
	_simulateBasesQualities(_secondMate, Haplotype, fragmentLength, _numCycles[1], readIsContaminated);

	// WRITE ALIGNMENTS
	_alignment.setMateGenomicPosition(_secondMate);
	_secondMate.setMateGenomicPosition(_alignment);
}

} // namespace Simulations
