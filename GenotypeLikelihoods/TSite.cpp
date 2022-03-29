/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"
#include "GenotypeTypes.h"
#include <numeric>

namespace GenotypeLikelihoods {

//-------------------------------------------------------
// TSite
//-------------------------------------------------------

void TSite::clear() noexcept {
	_bases.clear();
	refBase  = genometools::Base::N;
	genotype = genometools::Genotype::NN;
};

void TSite::add(const BAM::TSequencedBase &base) { _bases.push_back(base); };
void TSite::add(BAM::TSequencedBase &&base) { _bases.push_back(base); };

void TSite::addToBaseFrequencies(TBaseData &frequencies) const noexcept {
	if (!empty()) {
		const double weight = 1.0 / _bases.size();
		for (auto &b : _bases) { frequencies[b.base] += weight; }
	}
};

void TSite::downsample(uint32_t maxDepth, const coretools::TSubsamplePicker &picker) {
	// only subsample if depth > maxDepth
	if (_bases.size() > maxDepth) {
		// select subsample
		const auto &subsample = picker.pick(_bases.size(), maxDepth);

		// copy bases to new vector
		std::vector<BAM::TSequencedBase> newBases;
		newBases.reserve(subsample.size());
		for (auto i : subsample) { newBases.push_back(_bases[i]); }

		// swap vectors
		_bases = std::move(newBases);
	}
};

std::string TSite::getBases() const {
	return std::accumulate(_bases.cbegin(), _bases.cend(), "",
			       [](auto tot, auto b) { return tot + genometools::base2char(b.base); });
};

std::string TSite::getQualities() const {
	return std::accumulate(_bases.cbegin(), _bases.cend(), "", [](auto tot, auto b) {
		return tot + static_cast<char>(genometools::BaseQuality(b.recalibratedQualityAsPhredInt)); });
}

uint32_t TSite::refDepth() const {
	if (refBase == genometools::Base::N) return 0;

	return std::count_if(_bases.cbegin(), _bases.cend(), [this](auto b) {return b.base == refBase;});
};

void TSite::countAlleles(TBaseCounts &alleleCounts) const {
	alleleCounts.reset();
	for (const auto &b : _bases) { alleleCounts.add(b.base); }
};

void TSite::countMates(std::array<int, 2> &mateCounts) const {
	mateCounts.fill(0);
	for (const auto &b : _bases) { ++mateCounts[b.isSecondMate()]; }
};

void TSite::countFwdRev(std::array<int, 2> &frCounts) const {
	frCounts.fill(0);
	for (const auto &b : _bases) { ++frCounts[b.isReverseStrand()]; }
};

}; // namespace GenotypeLikelihoods
