/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Strings/toString.h"
#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "coretools/algorithms.h"
#include "coretools/Types/probability.h"
#include "stattools/MLEInference/TNelderMead.h"

#include <limits>
#include <numeric>
#include <random>
#include <armadillo>

namespace GenotypeLikelihoods {
using genometools::Base;
using genometools::Genotype;

namespace impl {
using coretools::TStrongArray;

TStrongArray<TGenotypeProbabilities, genometools::Base>	piTable(double mu, double theta_r, double theta_g) {
	using coretools::index;

	const arma::mat::fixed<4,4> l   = {{-2 - mu, 1, mu, 1}, {1, -2 - mu, 1, mu}, {mu, 1, -2 - mu, 1}, {1, mu, 1, -2 - mu}};
	const arma::mat::fixed<4,4> P_r = arma::expmat(theta_r*l);
	const arma::mat::fixed<4,4> P_g = arma::expmat(theta_g*l);

	coretools::TStrongArray<TGenotypeProbabilities, genometools::Base> pi;
	for (auto r = Base::min; r < Base::max; ++r) {
		TGenotypeData pi_r{};
		for (auto g = Genotype::min; g < Genotype::max; ++g) {
			const auto k = genometools::first(g);
			const auto l = genometools::second(g);
			const auto f = genometools::isHeterozygous(g) + 1; // homo = 1, het = 2
			for (auto R = Base::min; R < Base::max; ++R) {
				pi_r[g] += f*P_r(index(R),index(r))*P_g(index(k),index(R))*P_g(index(l),index(R));
			}
		}
		pi[r] = TGenotypeProbabilities::normalize(pi_r);
	}
	return pi;
}

double Q(double mu, double theta_r, double theta_g,
		 const coretools::TStrongArray<TGenotypeData, genometools::Base> &lkhSum) {
	try {
		const auto pi = impl::piTable(mu, theta_r, theta_g);
		double Q      = 0;
		for (auto r = Base::min; r < Base::max; ++r) {
			for (auto g = Genotype::min; g < Genotype::max; ++g) { Q += std::log(pi[r][g]) * lkhSum[r][g]; }
		}
		return Q;
	} catch (...) { return std::numeric_limits<double>::lowest(); }
}
} // namespace impl

TGenotypeLikelihoods THaploidDistribution::getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const {
	return TGenotypeLikelihoods({baseLikelihoods[Base::A], 0., 0., 0.,
			                     baseLikelihoods[Base::C], 0., 0.,
								 baseLikelihoods[Base::G], 0.,
								 baseLikelihoods[Base::T]});
}
coretools::Probability THaploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// first == second if Haploid
	return baseLikelihoods[genometools::first(genotype)];
}

double THaploidDistribution::normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base) {
	double sum = 0.;
	// only four
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g = genometools::genotype(b, b);
		likelihoods[g] *= _pi[b];
		sum += likelihoods[g];
	}
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g = genometools::genotype(b, b);
		_piSum[b] += likelihoods[g].scale(sum);
	}
	return sum;
}

void THaploidDistribution::estimate() {
	_pi = TBaseProbabilities::normalize(_piSum);
	_piSum.fill(0.);
}

std::string THaploidDistribution::definition() const noexcept {
	using coretools::str::toString;
	return std::string("AA: ").append(toString(_pi[Base::A]))
		.append(", CC: ").append(toString(_pi[Base::C]))
		.append(", GG: ").append(toString(_pi[Base::G]))
		.append(", TT: ").append(toString(_pi[Base::T]));
}

TGenotypeLikelihoods TDiploidDistribution::getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const {
	return TGenotypeLikelihoods({baseLikelihoods[Base::A], 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::C]),
								 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::G]),
								 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::T]), baseLikelihoods[Base::C],
								 0.5 * (baseLikelihoods[Base::C] + baseLikelihoods[Base::G]),
								 0.5 * (baseLikelihoods[Base::C] + baseLikelihoods[Base::T]), baseLikelihoods[Base::G],
								 0.5 * (baseLikelihoods[Base::G] + baseLikelihoods[Base::T]),
								 baseLikelihoods[Base::T]});
}

coretools::Probability TDiploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// if first == second, then 0.5*first + 0.5*first = first
	return 0.5 * (baseLikelihoods[genometools::first(genotype)] + baseLikelihoods[genometools::second(genotype)]);
}

double TDiploidDistribution::normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base) {
	double sum = 0;
	// all 10
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		likelihoods[g] *= _pi[g];
		sum += likelihoods[g];
	}
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		_piSum[g] += likelihoods[g].scale(sum);
	}
	return sum;
}

void TDiploidDistribution::estimate() {
	// Only update if values make sense, else set to initial value
	const auto sumHomo = _piSum[Genotype::AA] + _piSum[Genotype::CC] + _piSum[Genotype::GG] + _piSum[Genotype::TT];
	const auto sumHet  = _piSum[Genotype::AC] + _piSum[Genotype::AG] + _piSum[Genotype::AT] + _piSum[Genotype::CG] +
						 _piSum[Genotype::CT] + _piSum[Genotype::GT];

	if (sumHomo > sumHet) {
		// estimate looks good enough
		_pi = TGenotypeProbabilities::normalize(_piSum);
	} else {
		// the estimate is definitely wrong. Let's use the default
		_pi = _pi_init;
	}

	_piSum.fill(0.);
}

std::string TDiploidDistribution::definition() const noexcept {
	using coretools::str::toString;
	std::string ret;
	for (auto g = Genotype::min; g < Genotype::max; ++g) {
		ret.append(toString(g)).append(": ").append(toString(_pi[g])).append(", ");
	}
	ret.resize(ret.size() - 2);
	return ret;
}

TGenotypeLikelihoods THKY85::getGenotypeLikelihoods(const TBaseLikelihoods &baseLikelihoods) const {
	return TGenotypeLikelihoods({baseLikelihoods[Base::A], 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::C]),
								 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::G]),
								 0.5 * (baseLikelihoods[Base::A] + baseLikelihoods[Base::T]), baseLikelihoods[Base::C],
								 0.5 * (baseLikelihoods[Base::C] + baseLikelihoods[Base::G]),
								 0.5 * (baseLikelihoods[Base::C] + baseLikelihoods[Base::T]), baseLikelihoods[Base::G],
								 0.5 * (baseLikelihoods[Base::G] + baseLikelihoods[Base::T]),
								 baseLikelihoods[Base::T]});
}

coretools::Probability THKY85::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// if first == second, then 0.5*first + 0.5*first = first
	return 0.5 * (baseLikelihoods[genometools::first(genotype)] + baseLikelihoods[genometools::second(genotype)]);
}

double THKY85::normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base ref) {
	double sum = 0;
	// all 10
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		likelihoods[g] *= _pi[ref][g];
		sum += likelihoods[g];
	}
	for(auto g = Genotype::min; g < Genotype::max; ++g) {
		_likelihoodSum[ref][g] += likelihoods[g].scale(sum);
	}
	return sum;
}


THKY85::THKY85()
	: _pi(impl::piTable(_mu, _theta_r, _theta_g)),
	  _nelderMead([this](auto Vals) { return -impl::Q(std::exp(Vals[0]), std::exp(Vals[1]), std::exp(Vals[2]), _likelihoodSum); }) {
}

void THKY85::estimate() {
	// If likelihoodSum is totally off, it is not worth it
	constexpr size_t AA_CC = 100;
	const bool isWorthIt = _likelihoodSum[Base::A][Genotype::AA] > AA_CC*_likelihoodSum[Base::A][Genotype::CC];

	if (isWorthIt && _nelderMead.minimize({std::log(_mu), std::log(_theta_r), std::log(_theta_g)},
							 10.)) {
		const auto &crds = _nelderMead.coordinates();
		_mu              = std::exp(crds[0]);
		_theta_r         = std::exp(crds[1]);
		_theta_g         = std::exp(crds[2]);
		_pi              = impl::piTable(_mu, _theta_r, _theta_g);
	}
	_likelihoodSum.fill({});
}

std::string THKY85::definition() const noexcept {
	return coretools::str::toString(name, ": mu=", _mu, ", theta_r=", _theta_r, ", theta_g=", _theta_g);
}

}; // namespace GenotypeLikelihoods
