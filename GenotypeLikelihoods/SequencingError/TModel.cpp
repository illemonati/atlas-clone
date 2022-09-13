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
#include "TError.h"
#include "TGenotypeData.h"
#include "TRandomGenerator.h"
#include "TSequencedBase.h"
#include "SequencingError/TCovariate.h"
#include "TStrongArray.h"
#include "devtools.h"
#include "mathFunctions.h"
#include "probability.h"
#include "stringFunctions.h"
#include "toString.h"
#include "weakTypes.h"
#include "enum.h"
#include "TError.h"

namespace GenotypeLikelihoods {
namespace SequencingError {
using coretools::Probability;
using coretools::index;
using coretools::instances::randomGenerator;
using genometools::Base;
using namespace std::literals;

//*********************************************************
// TRho
//*********************************************************

TRho::TRho(const std::string &Def) {
	using coretools::str::toString;
	using coretools::index;
	//"default" implies default rho

	if (Def == "default" || Def == "-") {
		return;
	}

	// otherwise: full matrix is provided
	std::vector<std::string> vec;
	std::string s = Def;
	coretools::str::fillContainerFromString(Def, vec, ';');
	if (vec.size() != 4) throw "Rho matrix has " + toString(vec.size()) + " instead of 4 rows!";

	// parse rows
	for (Base a = Base::min; a < Base::max; ++a) {
		std::string &row = vec[index(a)];
		coretools::str::trimString(row, "[]");
		std::vector<double> r;
		coretools::str::fillContainerFromString(row, r, ',');
		if (r.size() != 4)
			throw "Rho matrix has " + toString(r.size()) + " instead of 4 columns for row " + toString(index(a) + 1) + "!";

		r[index(a)] = 0.;
		_rho[a]     = r;
	}
}

void TRho::set(const BAM::RGInfo::TInfo & Def){
	if(Def.is_string() && Def.get() == "default"){
		return;
	}

	//check if definition is complete
	std::string expl = "Require four attributed 'A', 'C', 'G' and 'T', each specifying an array of four floating point numbers.";
	if(Def.size() != 4){
		UERROR("Unable to understand rho: ", expl);
	}

	//parse each FROM base
	for (Base from = Base::min; from < Base::max; ++from) {
		using genometools::toString;

		//check that attribute exists
		if(!Def.contains(toString(from))){
			UERROR("Unable to understand rho: missing attribute '", toString(from),"'! ", expl);
		}

		//ensure it is an array of length 4
		const BAM::RGInfo::TInfo& x = Def[toString(from)];
		if(!x.is_array() || x.size() != 4){
			UERROR("Unable to understand rho: attribute '", toString(from), "' does not specify an array of length 4! ", expl);
		}

		//parse row for each TO base
		for (Base to = Base::min; to < Base::max; ++to) {
			r[from][to] = x[index(to)];
			r[from][from] = 0.0;
		}
	}
}

std::string TRho::getDefinition() const noexcept {
	using coretools::str::toString;
	return "[[-,"s + toString(_rho[Base::A][Base::C]) + ',' + toString(_rho[Base::A][Base::G]) + ',' + toString(_rho[Base::A][Base::T]) + "];["
		+ toString(_rho[Base::C][Base::A]) + ",-," + toString(_rho[Base::C][Base::G]) + ',' + toString(_rho[Base::C][Base::T]) + "];["
		+ toString(_rho[Base::G][Base::A]) + ',' + toString(_rho[Base::G][Base::C]) + ",-," + toString(_rho[Base::G][Base::T]) + "];["
		+ toString(_rho[Base::T][Base::A]) + ',' + toString(_rho[Base::T][Base::C]) + ',' + toString(_rho[Base::T][Base::G]) + ",-]]";
}

BAM::RGInfo::TInfo TRho::getInfo() const noexcept{
	//TODO: find better conversion StribngArray -> JSON
	BAM::RGInfo::TInfo info;
	for (Base from = Base::min; from < Base::max; ++from) {
		std::vector<double> vec;
		vec.reserve(4);
		for (Base to = Base::min; to < Base::max; ++to) {
			vec.push_back(_rho[from][to]);
		}

		info[toString(from)] = vec;
	}
	return info;
}

void TRho::add(genometools::Base l, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept {
	for (auto k = Base::min; k < Base::max; ++k) {
		_rhoSum[k][l] += P_g_I_d*P_bbar_I_d[k];
	}
}

void TRho::estimate() noexcept {
	for (Base k = Base::min; k < Base::max; ++k) {
		_rhoSum[k][k] = 0.0;
		_rho[k] = _rhoSum[k];
	}
	// reset
	_rhoSum.fill({});
}

//*********************************************************
// TModelNoRecal
//*********************************************************
Probability TModelNoRecal::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return Probability(1.0);
	return (Probability)base.originalQuality_phredInt;
}

genometools::PhredIntProbability TModelNoRecal::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return genometools::PhredIntProbability::highest();
	return base.originalQuality_phredInt;
}

TBaseLikelihoods TModelNoRecal::getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
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


TModelRecal::TModelRecal(const std::string& RhoDef, const std::string &EpsilonDef): _rho(RhoDef), _epsilon(EpsilonDef) {}

TModelRecal::TModelRecal(BAM::RGInfo::TInfo & Def){
	//Extract RHO, if given
	if(Def.contains(TRho::name)){
		_rho.set(Def[TRho::name]);

		//remove rho from TInfo
		Def.erase(TRho::name);
	}

	//initialize covariates using provided definition
	_epsilon.initialize(Def);
}

BAM::RGInfo::TInfo TModelRecal::getInfo() const noexcept{
	BAM::RGInfo::TInfo info = _epsilon.getInfo();
	info.push_back(_rho.getInfo());
	return info;
}

//-------------------------------------------------
// functions to calculate error rates
//-------------------------------------------------

constexpr Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return Probability(0.9999999999);
	if (eta < -23.03) return Probability(0.0000000001);

	return coretools::logistic(eta);
}

Probability TModelRecal::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return Probability::highest();
	return _epsilon.calcErrorRate(base);
}

genometools::PhredIntProbability TModelRecal::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return genometools::PhredIntProbability::highest();
	return genometools::PhredIntProbability(_epsilon.calcErrorRate(base));
}

TBaseLikelihoods TModelRecal::getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
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

} // namespace SequencingError
} // namespace GenotypeLikelihoods
