/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "SequencingError/TModel.h"

#include <array>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include "RecalEstimatorTools.h"
#include "coretools/Main/TError.h"
#include "TGenotypeData.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSequencedBase.h"
#include "SequencingError/TCovariate.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Strings/fromString.h"
#include "coretools/devtools.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/weakTypes.h"
#include "coretools/enum.h"

namespace GenotypeLikelihoods {
namespace SequencingError {
using coretools::Probability;
using coretools::instances::randomGenerator;
using genometools::Base;
using namespace std::literals;
using namespace coretools::str;

//*********************************************************
// TRho
//*********************************************************

TRho::TRho(std::string_view Def) {
	using coretools::str::toString;
	//"default" implies default rho

	if (Def == "default" || Def == "-") {
		return;
	}

	TSplitter spl(Def, ';');
	size_t i = 0;
	for (auto s: spl) {
		if (i >= 4) UERROR("Too many rows given for rho, needed only 4!");

		std::array<double, 4> ar;
		TSplitter spl2(strip(s, "[]"), ',');
		size_t j = 0;

		for (auto ss : spl2) {
			if (j >= ar.size()) UERROR("Too many rho values given for row ", i, ", needed ", ar.size(), "!");

			if (strip(ss) == "-") {
				ar[j] = 0.;
			} else {
				coretools::str::fromString<true>(strip(ss), ar[j]);
			}
			++j;
		}
		if (j < ar.size()) UERROR("Too few(", j, ") rho values given, needed ", ar.size(), "!");

		ar[i]         = 0.;
		_rho[Base(i)] = TBaseProbabilities::normalize(ar);
		++i;
	}
	if (i < 4) UERROR("Too few rows given for rho, needed 4, not ", i, "!");
}

std::string TRho::definition() const noexcept {
	using coretools::str::toString;
	return "[[-,"s + toString(_rho[Base::A][Base::C]) + ',' + toString(_rho[Base::A][Base::G]) + ',' + toString(_rho[Base::A][Base::T]) + "];["
		+ toString(_rho[Base::C][Base::A]) + ",-," + toString(_rho[Base::C][Base::G]) + ',' + toString(_rho[Base::C][Base::T]) + "];["
		+ toString(_rho[Base::G][Base::A]) + ',' + toString(_rho[Base::G][Base::C]) + ",-," + toString(_rho[Base::G][Base::T]) + "];["
		+ toString(_rho[Base::T][Base::A]) + ',' + toString(_rho[Base::T][Base::C]) + ',' + toString(_rho[Base::T][Base::G]) + ",-]]";
}

BAM::RGInfo::TInfo TRho::info() const {
	return {
		{_rho[Base::A][Base::A].get(), _rho[Base::A][Base::C].get(), _rho[Base::A][Base::G].get(), _rho[Base::A][Base::T].get()},
		{_rho[Base::C][Base::A].get(), _rho[Base::C][Base::C].get(), _rho[Base::C][Base::G].get(), _rho[Base::C][Base::T].get()},
		{_rho[Base::G][Base::A].get(), _rho[Base::G][Base::C].get(), _rho[Base::G][Base::G].get(), _rho[Base::G][Base::T].get()},
		{_rho[Base::T][Base::A].get(), _rho[Base::T][Base::C].get(), _rho[Base::T][Base::G].get(), _rho[Base::T][Base::T].get()}
	};
}

void TRho::add(genometools::Base l, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept {
	for (auto k = Base::min; k < Base::max; ++k) {
		_rhoSum[k][l] += P_g_I_d*P_bbar_I_d[k];
	}
}

void TRho::estimate() noexcept {
	for (Base k = Base::min; k < Base::max; ++k) {
		_rhoSum[k][k] = 0.0;
		_rho[k] = TBaseProbabilities::normalize(_rhoSum[k]);
	}
	// reset
	_rhoSum.fill({});
}

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
	if (base == Base::N) {
		return TBaseLikelihoods{1.};
	}
	const auto eps = static_cast<Probability>(base.originalQuality_phredInt);
	TBaseLikelihoods baseLikelihoods{(1./3)*eps};
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

//-------------------------------------------------
// functions to calculate error rates
//-------------------------------------------------


constexpr Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return Probability(0.9999999999);
	if (eta < -23.03) return Probability(0.0000000001);

	return coretools::logistic(eta);
}

TModelRecal::TModelRecal(const BAM::RGInfo::TInfo &) : _epsilon("") {
	WINK();
}

Probability TModelRecal::errorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return Probability::highest();
	return _epsilon.calcErrorRate(base);
}

genometools::PhredIntProbability TModelRecal::phredInt(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return genometools::PhredIntProbability::highest();
	return genometools::PhredIntProbability(_epsilon.calcErrorRate(base));
}

TBaseLikelihoods TModelRecal::baseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) {
		return TBaseLikelihoods{1.};
	}
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
		constexpr coretools::TStrongArray<std::array<Base, 3>, Base> lss({std::array<Base, 3>{Base::C, Base::G, Base::T},
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
