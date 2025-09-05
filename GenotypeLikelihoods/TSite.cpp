/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"

#include "TSequencedData.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "coretools/algorithms.h"
#include "genometools/Genotypes/Containers.h"
#include "genometools/Genotypes/Base.h"

namespace GenotypeLikelihoods {
using genometools::TBaseData;
using genometools::TBaseCounts;

//-------------------------------------------------------
// TSite
//-------------------------------------------------------

void TSite::clear() noexcept {
	_data.clear();
	refBase = genometools::Base::N;
	altBase = genometools::Base::N;
}

TBaseData TSite::baseFrequencies() const noexcept {
	TBaseData bd{};
	if (!empty()) {
		const double weight = 1.0 / _data.size();
		for (auto &b : _data) { bd[b.base] += weight; }
	}
	return bd;
}

void TSite::downsample(coretools::Probability p) {
	using coretools::instances::randomGenerator;
	const auto iMax  = _data.size() - 1;
	const auto pComp = p.complement();
	for (int i = iMax; i >= 0; --i) {
		if (randomGenerator().getRand() < pComp) {
			std::swap(_data[i], _data.back());
			_data.pop_back();
		}
	}
}


void TSite::shuffle() {
	using coretools::instances::randomGenerator;
	randomGenerator().shuffle(_data);
}

void TSite::limitDepth(size_t UpToDepth) {
	if (UpToDepth < _data.size()) _data.resize(UpToDepth);
}

std::string TSite::getBases() const {
	if (empty()) return "-";
	return std::accumulate(_data.cbegin(), _data.cend(), std::string(""),
			       [](auto tot, auto b) { return tot + genometools::base2char(b.base); });
}

std::vector<genometools::Base> TSite::sampleBases() const {
	std::vector<std::vector<genometools::Base>> b_bam;
	for (auto data: _data) {
		if (data.bamID >= b_bam.size()) b_bam.resize(data.bamID + 1);
		b_bam[data.bamID].push_back(data.base);
	}
	std::vector<genometools::Base> ret;
	for (const auto& s: b_bam) {
		if (s.empty()) continue;
		ret.push_back(s[coretools::instances::randomGenerator().sample(s.size())]);
	}
	return ret;
}

std::string TSite::getQualities() const {
	if (empty()) return "-";
	return std::accumulate(_data.cbegin(), _data.cend(), std::string(""), [](auto tot, auto b) {
		return tot + coretools::toChar(b.recalQuality); });
}

size_t TSite::refDepth() const {
	if (refBase == genometools::Base::N) return 0;

	return std::count_if(_data.cbegin(), _data.cend(), [this](auto b) {return b.base == refBase;});
}

TBaseCounts TSite::countAlleles() const {
	TBaseCounts alleleCounts{};
	for (const auto &b : _data) { ++alleleCounts[b.base]; }
	return alleleCounts;
}

coretools::PhredInt TSite::rmsMappingQual() const {
	if (_data.empty()) return {};

	double sqSum = 0;
	for (const auto &d : _data) {
		const auto mQ = d.mappingQuality.get();
		sqSum += mQ*mQ;
	}
	return coretools::PhredInt(static_cast<coretools::PhredInt>(sqrt(sqSum/_data.size())));
}

coretools::TStrongArray<size_t, BAM::Mate> TSite::countMates() const {
	coretools::TStrongArray<size_t, BAM::Mate> mateCounts{};
	for (const auto &b : _data) { ++mateCounts[b.mate()]; }
	return mateCounts;
}

std::array<int, 2> TSite::countFwdRev() const {
	std::array<int, 2> frCounts{};
	for (const auto &b : _data) { ++frCounts[b.get<BAM::Flags::ReversedStrand>()]; }
	return frCounts;
}

}; // namespace GenotypeLikelihoods
