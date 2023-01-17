/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "SequencingError/TModel.h"

#include <algorithm>
#include <array>
#include <math.h>
#include <memory>
#include <ostream>
#include <stdlib.h>
#include <tuple>
#include <type_traits>
#include <utility>

#include "RecalEstimatorTools.h"
#include "SequencingError/TCovariate.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "coretools/Types/weakTypes.h"
#include "coretools/devtools.h"
#include "coretools/enum.h"

namespace GenotypeLikelihoods {
namespace SequencingError {
using coretools::Probability;
using coretools::instances::randomGenerator;
using genometools::Base;
using namespace std::literals;
using namespace coretools::str;

//*********************************************************
// TModelNoRecal
//*********************************************************
Probability TModelNoRecal::errorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return Probability(1.0);
	return (Probability)base.originalQuality_phredInt;
}

genometools::PhredIntProbability TModelNoRecal::phredInt(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return genometools::PhredIntProbability::highest();
	return base.originalQuality_phredInt;
}

TBaseLikelihoods TModelNoRecal::baseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) { return TBaseLikelihoods{1.}; }
	const auto eps = static_cast<Probability>(base.originalQuality_phredInt);
	TBaseLikelihoods baseLikelihoods{(1. / 3) * eps};
	baseLikelihoods[base.base] = eps.complement();
	return baseLikelihoods;
}

void TModelNoRecal::simulate(BAM::TSequencedBase &base) const noexcept {
	if (base.base == Base::N) return;

	const auto e = static_cast<Probability>(base.originalQuality_phredInt);
	if (randomGenerator().getRand() < e) {
		const int i = randomGenerator().getRand(1, 4); // 3 bases to choose from
		base.base   = Base((coretools::index(base.base) + i) % 4);
	}
}

//*********************************************************
// TModelRecal
//*********************************************************

Probability TModelRecal::errorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return Probability::highest();
	return _epsilon.calcErrorRate(base);
}

genometools::PhredIntProbability TModelRecal::phredInt(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return genometools::PhredIntProbability::highest();
	return genometools::PhredIntProbability(_epsilon.calcErrorRate(base));
}

TBaseLikelihoods TModelRecal::baseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) { return TBaseLikelihoods{1.}; }
	const auto e = _epsilon.calcErrorRate(base);
	const auto l = base.base;
	TBaseLikelihoods baseLikelihoods;
	for (auto k = Base::min; k < Base::max; ++k) baseLikelihoods[k] = e * _rho[k][l];
	baseLikelihoods[l] = e.complement();
	return baseLikelihoods;
}

void TModelRecal::simulate(BAM::TSequencedBase &base) const noexcept {
	if (base.base == Base::N) return;

	const auto e = _epsilon.calcErrorRate(base);
	if (randomGenerator().getRand() < e) {
		const auto k = base.base;
		constexpr coretools::TStrongArray<std::array<Base, 3>, Base> lss(
			{std::array<Base, 3>{Base::C, Base::G, Base::T},
			 {Base::A, Base::G, Base::T},
			 {Base::A, Base::C, Base::T},
			 {Base::A, Base::C, Base::G}});
		const auto ls = lss[k];

		double r = randomGenerator().getRand();
		if (r < _rho[k][ls[0]]) {
			base.base = ls[0];
			return;
		}
		r -= _rho[k][ls[0]];
		if ((r < _rho[k][ls[1]])) {
			base.base = ls[1];
			return;
		}
		// else
		base.base = ls[2];
	}
}

BAM::RGInfo::TInfo TModelRecal::info() const {
	BAM::RGInfo::TInfo info;
	info        = _epsilon.info();
	info["rho"] = _rho.info();
	return info;
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
