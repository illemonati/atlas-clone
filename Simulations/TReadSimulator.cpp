/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include <TReadSimulator.h>

#include "PMD/TModel.h"
#include "TOutputBamFile.h"
#include "TReadGroupInfo.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/TwoBases.h"

namespace Simulations {
using BAM::RGInfo::InfoType;
using BAM::RGInfo::TReadGroupInfoEntry;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using coretools::instances::parameters;
using coretools::probdist::TCategoricalDistribution;
using coretools::P;
using coretools::user_assert;
using genometools::Base;
using genometools::TwoBase;
using genometools::TGenomePosition;

namespace impl {
std::pair<size_t, size_t> refDiff(const std::vector<TwoBase> &Haplotype, coretools::TView<genometools::Base> Reference, size_t Pos, size_t Len) {
	std::pair<size_t, size_t> Ns{};
	for (size_t i = 0; i < Len; ++i) {
		const auto pi = Pos + i;
		const auto fi = first(Haplotype[pi]);
		const auto se = second(Haplotype[pi]);
		if (fi == se) continue;

		const auto re = Reference[pi];

		if (re == Base::N) continue;
		if ((fi != re) && (se != re)) continue;

		if (fi != re) ++Ns.first;
		else ++Ns.second;
	}
	return Ns;
}

//initialization functions
template <typename Distr>
void initDistribution(Distr & Dist, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const BAM::RGInfo::InfoType & Info){
	coretools::instances::logfile().list(coretools::str::capitalizeFirst(BAM::RGInfo::infos[Info].description), ": ", RGInfo.getString(Info));
	Dist.set(RGInfo.getString(Info));
}
} // namespace impl

//------------------------------------------------
// TSimulatorRead
//------------------------------------------------
TReadSimulator::TReadSimulator(const BAM::TReadGroup &ReadGroup, const TReadGroupInfoEntry &RGInfo,
							   const GenotypeLikelihoods::PMD::TModel &Pmd,
							   const GenotypeLikelihoods::SequencingError::RGModels &Recal)
	: _readGroup(&ReadGroup), _pmd(&Pmd), _recal(Recal) {

	// initialize bamAlignment

	//readNamePrefix: "<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	_readNamePrefix = "ATL:0:A:1:" + coretools::str::toString(_readGroup->id) + ":";

	//initialize distributions
	impl::initDistribution(_fragmentLengthDistr, RGInfo, InfoType::fragmentLength);
	impl::initDistribution(_mappingQualityDist, RGInfo, InfoType::mappingQuality);
	impl::initDistribution(_qualityDist, RGInfo, InfoType::baseQuality);

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
		user_assert(_duplicationRate <= 0.5, "Duplication rate must be within [0.0, 0.5]!");
		_duplicationRateAmongSimulated = P(_duplicationRate / (_duplicationRate.complement()));
	}

	if(parameters().exists("baseN")){
		_baseN = parameters().get<coretools::Probability>("baseN");
		logfile().list("Will simulate reads with base = N probability = ", _baseN, ". (parameter 'baseN')");
	} else {
		_baseN = P(0.001);
		logfile().list("Will simulate reads with base = N probability = ", _baseN, ". (set with 'baseN')");
	}

	if(parameters().exists("refBias")){
		_refBias = parameters().get<coretools::Probability>("refBias");
		logfile().list("Will simulate reads with a reference bias of ", _refBias, ". (parameter 'refBias')");
	} else {
		_refBias = P(0.5);
		logfile().list("Will simulate reads with a reference bias of ", _refBias, ". (set with 'refBias')");
	}
}

TReadSimulator::~TReadSimulator() {
	logfile().list(_refCount[0], " Reference difference in written haplotypes.");
	logfile().list(_refCount[1], " Reference difference in unwritten haplotypes.");
	logfile().list("Reference bias: ", double(_refCount[1])/(_refCount[0]+_refCount[1]));
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

bool TReadSimulator::_simulateContamination(){
	return _contaminationRate > 0. && randomGenerator().getRand() < _contaminationRate;
}

void TReadSimulator::_addSoftclippedBases(std::vector<Base> & Bases,const size_t &SoftClipLength, BAM::TCigar & Cigar) {
	if(SoftClipLength > 0){
		for (size_t i = 0; i < SoftClipLength; i++){
			Bases.push_back(static_cast<Base>(randomGenerator().getRand<uint8_t>(0,4)));
		}
		Cigar.add('S', SoftClipLength);
	}
}


void TReadSimulator::_simulateBasesQualities(BAM::TAlignment &Alignment, const std::vector<TwoBase> &Haplotype,
											 bool FirstHaplo, size_t ReadLength,
											 bool ReadIsContaminated) {

	assert(ReadLength > 0);
	user_assert(!ReadIsContaminated, "contaminated Reads not implemented yet!");
	// prepare vector of bases
	std::vector<Base> bases;
	BAM::TCigar cigar;

	//sample softclip lengths
	const auto softClipLength3 = _softClipDist3 ? _softClipDist3->sample() : 0;
	const auto softClipLength5 = _softClipDist5 ? _softClipDist5->sample() : 0;

	// set read length
	if (Alignment.isReverseStrand()) {
		_addSoftclippedBases(bases, softClipLength3, cigar);
	} else {
		_addSoftclippedBases(bases, softClipLength5, cigar);
	}

	// simulate true bases

	auto picker = [FirstHaplo](TwoBase tb) {
		return FirstHaplo ? first(tb) : second(tb);
	};
	for (size_t i = 0; i < ReadLength; ++i) {
		bases.push_back(picker(Haplotype[i + Alignment.position()]));
	}

	cigar.add('M', ReadLength);

	if (Alignment.isReverseStrand()) {
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

	Alignment.setSequenceQualities(cigar, bases, phredIntQualities);

	_pmd->simulate(Alignment);
	_recal[Alignment.mate()]->simulate(Alignment);
}

size_t TReadSimulator::simulate(const TGenomePosition &Position, const std::vector<TwoBase> &Haplotype,
								coretools::TView<genometools::Base> Reference, BAM::TOutputBamFile &BamFile) {
	// Do not simulate fraction of reads that will be duplicates
	if (_duplicationRate == 0.0) {
		if (!_simulate(Position, Haplotype, Reference)) return 0;
		_writeSimulatedAlignments(BamFile);
		return 1;
	} else if (randomGenerator().getRand() > _duplicationRate) {
		if (!_simulate(Position, Haplotype, Reference)) return 0;
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
		user_assert(json.get<int>() >= 1 && json.get<int>() <= 65535, errRange);

		_numCycles = json.get<int>();
	} else if(json.is_array() && json.size() != 1){
		throw coretools::TUserError(errRange);
	} else if(json.is_string()){
		coretools::str::fromString(json.get<std::string>(), _numCycles, err);
	} else if(json.is_array() && json.size() == 1){
		coretools::str::fromString(json[0].get<std::string>(), _numCycles, err);
	} else {
		throw coretools::TUserError(errRange);
	}

	_alignment.setReadGroup(_readGroup->id);
}

double TReadSimulatorSingleEnd::meanReadLength() const {
	return _calcMeanReadLength(_numCycles);
}

bool TReadSimulatorSingleEnd::_simulate(const TGenomePosition &Position, const std::vector<TwoBase> &Haplotype,
										coretools::TView<genometools::Base> Reference) {
	const auto fragmentLength = _fragmentLengthDistr.sample();
	const auto readLength     = std::min(fragmentLength, _numCycles);

	const bool firstHaplo = randomGenerator().getRand() < 0.5;
	const auto refDiff    = impl::refDiff(Haplotype, Reference, Position.position(), readLength);

	// Discard alignment?
	if (refDiff.first < refDiff.second) {
		// more ref Difference in second mate -> bias towards first
		const auto oRatio = firstHaplo ? _refBias.oddsRatio() : _refBias.complement().oddsRatio();
		if (oRatio < 1 && randomGenerator().getRand() > oRatio) return false;
	} else if (refDiff.first > refDiff.second) {
		// more ref Difference in first mate -> bias towards second
		const auto oRatio = firstHaplo ? _refBias.complement().oddsRatio() : _refBias.oddsRatio();
		if (oRatio < 1 && randomGenerator().getRand() > oRatio) return false;
	}

	// Statistics
	if (firstHaplo) {
		_refCount.front() += refDiff.first;
		_refCount.back() += refDiff.second;
	} else {
		_refCount.front() += refDiff.second;
		_refCount.back() += refDiff.first;
	}

	_alignment.move(Position);
	_alignment.setName(_getNextReadName());
	_alignment.setMappingQuality(_mappingQualityDist.sample());

	if (randomGenerator().getRand() < 0.5) {
		_alignment.setIsReverseStrand(true);
		_alignment.setInsertSize(-fragmentLength);
	} else {
		_alignment.setIsReverseStrand(false);
		_alignment.setInsertSize(fragmentLength);
	}

	_simulateBasesQualities(_alignment, Haplotype, firstHaplo, readLength, _simulateContamination());
	return true;
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
		using BAM::RGInfo::toString;

	//num cycles
	logfile().list(BAM::RGInfo::infos[InfoType::cycles].description, ": ", toString(RGInfo[InfoType::cycles]));
	auto& json = RGInfo[InfoType::cycles];

	std::string err = "Unable to understand " + BAM::RGInfo::infos[InfoType::cycles].description + ": ";
	std::string errRange = err + "expect one or two integers within [1,65535].";

	//TODO: probably need some json parsing functions to simplify this!
	if(json.is_array()){
		if(json.size() == 1){
			if(json[0].get<int>() < 1 || json[0].get<int>() > 65535){
				throw coretools::TUserError(errRange);
			}
			_numCycles[0] = json[0].get<int>();
			_numCycles[1] = _numCycles[0];
		} else if(json.size() == 2){
			if(json[0].get<int>() < 1 || json[1].get<int>() > 65535 || json[1].get<int>() < 1 || json[1].get<int>() > 65535){
				throw coretools::TUserError(errRange);
			}
			_numCycles[0] = json[0].get<int>();
			_numCycles[1] = json[1].get<int>();
		} else {
			throw coretools::TUserError(errRange);
		}
	} else if(json.is_number()){
		if(json.get<int>() < 1 || json.get<int>() > 65535){
			throw coretools::TUserError(errRange);
		}
		_numCycles[0] = json.get<int>();
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
		throw coretools::TUserError(errRange);
	}

	// set Alignment properties
	_fwdStrand.setIsPaired(true);
	_fwdStrand.setIsProperPair(true);
	_fwdStrand.setIsReverseStrand(false);
	_fwdStrand.setReadGroup(_readGroup->id);

	_revStrand.setIsPaired(true);
	_revStrand.setIsProperPair(true);
	_revStrand.setIsReverseStrand(true);
	_revStrand.setReadGroup(_readGroup->id);
}

double TReadSimulatorPairedEnd::meanReadLength() const {
	return _calcMeanReadLength(_numCycles[0] + _numCycles[1]);
}

void TReadSimulatorPairedEnd::_writeSimulatedAlignments(BAM::TOutputBamFile & BamFile){
	BamFile.writeAlignment(_fwdStrand);

	// write mate if it starts at same position as first, and keep for writing later otherwise
	if (_revStrand.from() == _fwdStrand.from()) {
		BamFile.writeAlignment(_revStrand);
	} else {
		BamFile.writeAlignmentLater(_revStrand);
	}
}

bool TReadSimulatorPairedEnd::_simulate(const TGenomePosition &Position, const std::vector<TwoBase> &Haplotype,
										coretools::TView<genometools::Base> Reference) {
	const auto fragmentLength     = _fragmentLengthDistr.sample();
	const auto readLength1        = std::min(fragmentLength, _numCycles.front());
	const auto readLength2        = std::min(fragmentLength, _numCycles.back());
	const auto readIsContaminated = _simulateContamination();

	const bool firstHaplo = randomGenerator().getRand() < 0.5;

	_fwdStrand.move(Position);
	_revStrand.move(_fwdStrand.from());
	if(fragmentLength > _numCycles.back()){
		_revStrand.advanceOnRef(fragmentLength - _numCycles.back());
	}

	const auto refDiff1 = impl::refDiff(Haplotype, Reference, _fwdStrand.position(), readLength1);
	const auto refDiff2 = impl::refDiff(Haplotype, Reference, _revStrand.position(), readLength2);

	// Discard alignment?
	if (refDiff1.first + refDiff2.first < refDiff1.second + refDiff2.second) {
		// more ref Difference in second mate -> bias towards first
		const auto oRatio = firstHaplo ? _refBias.oddsRatio() : _refBias.complement().oddsRatio();
		if (oRatio < 1 && randomGenerator().getRand() > oRatio) return false;
	} else if (refDiff1.first + refDiff2.first > refDiff1.second + refDiff2.second) {
		// more ref Difference in first mate -> bias towards second
		const auto oRatio = firstHaplo ? _refBias.complement().oddsRatio() : _refBias.oddsRatio();
		if (oRatio < 1 && randomGenerator().getRand() > oRatio) return false;
	}

	// Statistics
	if (firstHaplo) {
		_refCount.front() += refDiff1.first + refDiff2.first;
		_refCount.back() += refDiff1.second + refDiff2.second;
	} else {
		_refCount.front() += refDiff1.second + refDiff2.second;
		_refCount.back() += refDiff1.first + refDiff2.first;
	}

	_fwdStrand.setName(_getNextReadName());
	_fwdStrand.setMappingQuality(_mappingQualityDist.sample());

	_revStrand.setName(_fwdStrand.name());
	_revStrand.setMappingQuality(_fwdStrand.mappingQuality());

	_fwdStrand.setInsertSize(fragmentLength);
	_revStrand.setInsertSize(-fragmentLength);

	if (randomGenerator().getRand() < 0.5) {
		_fwdStrand.setIsSecondMate(true);
		_revStrand.setIsSecondMate(false);
	} else {
		_fwdStrand.setIsSecondMate(false);
		_revStrand.setIsSecondMate(true);
	}

	assert(_fwdStrand.isReverseStrand() != _revStrand.isReverseStrand());
	assert(_fwdStrand.isSecondMate() != _revStrand.isSecondMate());

	_simulateBasesQualities(_fwdStrand, Haplotype, firstHaplo, readLength1, readIsContaminated);
	_simulateBasesQualities(_revStrand, Haplotype, firstHaplo, readLength2, readIsContaminated);

	_fwdStrand.setMateGenomicPosition(_revStrand.from());
	_revStrand.setMateGenomicPosition(_fwdStrand.from());
	return true;
}

} // namespace Simulations
