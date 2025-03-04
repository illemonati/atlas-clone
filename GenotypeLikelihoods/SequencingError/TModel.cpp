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
using genometools::TBaseLikelihoods;

//*********************************************************
// TModelNoRecal
//*********************************************************
coretools::PhredInt TNoRecal::phredInt(const BAM::TSequencedData &data) const noexcept {
	if (data.base == Base::N) return coretools::PhredInt::highest();
	return data.originalQuality;
}

TBaseLikelihoods TNoRecal::P_dij(const BAM::TSequencedData &data) const noexcept {
	if (data.base == Base::N) { return TBaseLikelihoods{P(1.)}; }
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

coretools::PhredInt TWithRecal::phredInt(const BAM::TSequencedData &data) const noexcept {
	if (data.base == Base::N) return coretools::PhredInt::highest();
	return coretools::PhredInt(_epsilon.calcErrorRate(data));
}

TBaseLikelihoods TWithRecal::P_dij(const BAM::TSequencedData &data) const noexcept {
	if (data.base == Base::N) { return TBaseLikelihoods{P(1.)}; }
	const auto e = _epsilon.calcErrorRate(data);
	const auto l = data.base;
	TBaseLikelihoods baseLikelihoods;
	for (auto k = Base::min; k < Base::max; ++k) baseLikelihoods[k] = e * _rho.prob(k, data);
	baseLikelihoods[l] = e.complement();
	return baseLikelihoods;
}

void TWithRecal::simulate(BAM::TAlignment &aln) const noexcept {
	for (auto &d : aln) {
		if (d.get<BAM::Flags::SoftClipped>() || d.base == Base::N) continue;

		if (randomGenerator().getRand() < _epsilon.calcErrorRate(d)) {
			constexpr coretools::TStrongArray<std::array<Base, 3>, Base> lss(
				{std::array<Base, 3>{Base::C, Base::G, Base::T},
				 {Base::A, Base::G, Base::T},
				 {Base::A, Base::C, Base::T},
				 {Base::A, Base::C, Base::G}});

			const auto ls = lss[d.base];

			const double r = randomGenerator().getRand();
			double rhoSum  = _rho.prob(d, ls[0]);
			if (r < rhoSum) {
				d.base = ls[0];
			} else if (r < rhoSum + _rho.prob(d,ls[1])) {
				d.base = ls[1];
			} else {
				d.base = ls[2];
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
