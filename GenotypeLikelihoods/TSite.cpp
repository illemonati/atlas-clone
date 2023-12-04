/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"

#include <algorithm>
#include <numeric>
#include <utility>

#include "coretools/Containers/TStrongArray.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"

namespace GenotypeLikelihoods {

//-------------------------------------------------------
// TSite
//-------------------------------------------------------

void TSite::clear() noexcept {
	_bases.clear();
	refBase  = genometools::Base::N;
	genotype = genometools::Genotype::NN;
}

void TSite::add(const BAM::TSequencedBase &base) { _bases.push_back(base); }

TBaseData TSite::baseFrequencies() const noexcept {
	TBaseData bd{};
	if (!empty()) {
		const double weight = 1.0 / _bases.size();
		for (auto &b : _bases) { bd[b.base] += weight; }
	}
	return bd;
}

void TSite::downsample(size_t maxDepth, const coretools::TSubsamplePicker &picker) {
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
}

std::string TSite::getBases() const {
	if (empty()) return "-";
	return std::accumulate(_bases.cbegin(), _bases.cend(), std::string(""),
			       [](auto tot, auto b) { return tot + genometools::base2char(b.base); });
}

std::string TSite::getQualities() const {
	if (empty()) return "-";
	return std::accumulate(_bases.cbegin(), _bases.cend(), std::string(""), [](auto tot, auto b) {
		return tot + static_cast<char>(genometools::BaseQuality(b.recalibratedQualityAsPhredInt)); });
}

size_t TSite::refDepth() const {
	if (refBase == genometools::Base::N) return 0;

	return std::count_if(_bases.cbegin(), _bases.cend(), [this](auto b) {return b.base == refBase;});
}

TBaseCounts TSite::countAlleles() const {
	TBaseCounts alleleCounts{};
	for (const auto &b : _bases) { ++alleleCounts[b.base]; }
	return alleleCounts;
}

coretools::TStrongArray<size_t, BAM::Mate> TSite::countMates() const {
	coretools::TStrongArray<size_t, BAM::Mate> mateCounts{};
	for (const auto &b : _bases) { ++mateCounts[b.mate()]; }
	return mateCounts;
}

std::array<int, 2> TSite::countFwdRev() const {
	std::array<int, 2> frCounts{};
	for (const auto &b : _bases) { ++frCounts[b.isReverseStrand()]; }
	return frCounts;
}

}; // namespace GenotypeLikelihoods
