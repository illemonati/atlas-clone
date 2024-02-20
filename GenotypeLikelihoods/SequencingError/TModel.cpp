/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "SequencingError/TModel.h"
#include "TAlignment.h"
#include "coretools/Main/TRandomGenerator.h"

namespace GenotypeLikelihoods {
namespace SequencingError {
using coretools::Probability;
using coretools::P;
using coretools::instances::randomGenerator;
using genometools::Base;

//*********************************************************
// TModelNoRecal
//*********************************************************
genometools::PhredIntProbability TNoRecal::phredInt(const BAM::TSequencedBase &data) const noexcept {
	if (data == Base::N) return genometools::PhredIntProbability::highest();
	return data.originalQuality;
}

TBaseLikelihoods TNoRecal::P_dij(const BAM::TSequencedBase &data) const noexcept {
	if (data == Base::N) { return TBaseLikelihoods{P(1.)}; }
	const auto eps = static_cast<Probability>(data.originalQuality);
	TBaseLikelihoods baseLikelihoods{P(1. / 3) * eps};
	baseLikelihoods[data.base] = eps.complement();
	return baseLikelihoods;
}

void TNoRecal::simulate(BAM::TAlignment &aln) const noexcept {
	for (auto &data : aln) {
		if (data.base == Base::N) continue;

		const auto e = static_cast<Probability>(data.originalQuality);
		if (randomGenerator().getRand() < e) {
			const int i = randomGenerator().getRand(1, 4); // 3 bases to choose from
			data.base   = Base((coretools::index(data.base) + i) % 4);
		}
	}
}

void TNoRecal::recalibrate(BAM::TAlignment &aln) const noexcept {
	for (auto &b : aln) b.recalQuality = b.originalQuality;
}

//*********************************************************
// TModelRecal
//*********************************************************

genometools::PhredIntProbability TWithRecal::phredInt(const BAM::TSequencedBase &data) const noexcept {
	if (data == Base::N) return genometools::PhredIntProbability::highest();
	return genometools::PhredIntProbability(_epsilon.calcErrorRate(data));
}

TBaseLikelihoods TWithRecal::P_dij(const BAM::TSequencedBase &data) const noexcept {
	if (data == Base::N) { return TBaseLikelihoods{P(1.)}; }
	const auto e = _epsilon.calcErrorRate(data);
	const auto l = data.base;
	TBaseLikelihoods baseLikelihoods;
	for (auto k = Base::min; k < Base::max; ++k) baseLikelihoods[k] = e * _rho[k][l];
	baseLikelihoods[l] = e.complement();
	return baseLikelihoods;
}

void TWithRecal::simulate(BAM::TAlignment &aln) const noexcept {
	for (auto &data : aln) {
		if (data.base != Base::N && randomGenerator().getRand() < _epsilon.calcErrorRate(data)) {
			constexpr coretools::TStrongArray<std::array<Base, 3>, Base> lss(
				{std::array<Base, 3>{Base::C, Base::G, Base::T},
				 {Base::A, Base::G, Base::T},
				 {Base::A, Base::C, Base::T},
				 {Base::A, Base::C, Base::G}});

			const auto k = data.base;
			const auto ls = lss[k];

			const double r = randomGenerator().getRand();
			if (r < _rho[k][ls[0]]) {
				data.base = ls[0];
			} else if (r < _rho[k][ls[0]] + _rho[k][ls[1]]) {
				data.base = ls[1];
			} else {
				data.base = ls[2];
			}
		}
	}
}

void TWithRecal::recalibrate(BAM::TAlignment &aln) const noexcept {
	for (auto &b : aln) b.recalQuality = phredInt(b);
}

BAM::RGInfo::TInfo TWithRecal::info() const {
	BAM::RGInfo::TInfo info;
	info        = _epsilon.info();
	info["rho"] = _rho.info();
	return info;
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
