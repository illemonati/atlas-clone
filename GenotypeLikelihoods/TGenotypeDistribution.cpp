/*
 * TGenotypeModel.cpp
 *
 *  Created on: May 12, 2020
 *      Author: phaentu
 */

#include "TGenotypeDistribution.h"
#include "GenotypeData.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Strings/toString.h"
#include "stattools/MLEInference/TNelderMead.h"

#include <armadillo>

namespace GenotypeLikelihoods {
using genometools::Base;
using genometools::Genotype;
using coretools::instances::logfile;

namespace impl {
using coretools::TStrongArray;

TStrongArray<TGenotypeProbabilities, genometools::Base>	piTable(double mu, double theta_r, double theta_g) {
	using coretools::index;

	const auto z                     = (1. - mu) / 2;
	const arma::mat::fixed<4, 4> l   = {{-1., z, mu, z}, {z, -1., z, mu}, {mu, z, -1., z}, {z, mu, z, -1.}};
	const arma::mat::fixed<4, 4> P_g = arma::expmat(theta_g * l);
	const arma::mat::fixed<4, 4> P_r = arma::expmat(theta_r * l);

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
	} catch (...) {
		return std::numeric_limits<double>::lowest();
	}
}

TStrongArray<TBaseProbabilities, genometools::Base>	piTable(double mu, double theta) {
	using coretools::index;

	const auto z                   = (1. - mu) / 2;
	const arma::mat::fixed<4, 4> l = {{-1., z, mu, z}, {z, -1., z, mu}, {mu, z, -1., z}, {z, mu, z, -1.}};
	const arma::mat::fixed<4, 4> P = arma::expmat(theta * l);
	coretools::TStrongArray<TBaseProbabilities, genometools::Base> pi;

	for (auto r = Base::min; r < Base::max; ++r) {
		TBaseData pi_r;
		for (auto b = Base::min; b < Base::max; ++b) {
			pi_r[b] = coretools::Probability(P(index(r), index(b)));
		}
		pi[r] = TBaseProbabilities::normalize(pi_r);
	}
	return pi;
}

double Q(double mu, double theta, const coretools::TStrongArray<TBaseData, genometools::Base> &lkhSum) {
	try {
		const auto pi = impl::piTable(mu, theta);
		double Q      = 0;
		for (auto r = Base::min; r < Base::max; ++r) {
			for (auto b = Base::min; b < Base::max; ++b) { Q += std::log(pi[r][b]) * lkhSum[r][b]; }
		}
		return Q;
	} catch (...) { return std::numeric_limits<double>::lowest(); }
}

double het(double mu, double theta) {
	const auto z                     = (1. - mu) / 2;
	const arma::mat::fixed<4, 4> l   = {{-1., z, mu, z}, {z, -1., z, mu}, {mu, z, -1., z}, {z, mu, z, -1.}};
	const arma::mat::fixed<4, 4> P_g = arma::expmat(theta * l);

	double hom = 0.;
	for (size_t i = 0; i < 4; ++i) {
		hom += P_g(0, i)*P_g(0, i);
	}
	return 1-hom;
}

} // namespace impl

TGenotypeLikelihoods THaploidDistribution::P_dij(const TBaseLikelihoods &baseLikelihoods) const {
	return base2genotype<genometools::Ploidy::haploid>(baseLikelihoods);
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

void THaploidDistribution::log() const {
	logfile().list("AA: ", _pi[Base::A], ", CC: ", _pi[Base::C], ", GG: ", _pi[Base::G], ", TT: ", _pi[Base::T]);
}

void THaploidDistribution::addHeader(std::vector<std::string> &Header, std::string_view Prefix) const {
	using coretools::str::toString;
	Header.push_back(toString(Prefix, "piA"));
	Header.push_back(toString(Prefix, "piC"));
	Header.push_back(toString(Prefix, "piG"));
	Header.push_back(toString(Prefix, "piT"));
}

std::vector<double> THaploidDistribution::pis() const {
	return{_pi[Base::A], _pi[Base::C], _pi[Base::G], _pi[Base::T]};
	
}

void THaploidDistribution::reset() {
	_piSum.fill(0.);
	_pi = TBaseProbabilities::normalize(_piSum);
}

TGenotypeLikelihoods TDiploidDistribution::P_dij(const TBaseLikelihoods &baseLikelihoods) const {
	return base2genotype<genometools::Ploidy::diploid>(baseLikelihoods);
}

coretools::Probability TDiploidDistribution::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// if first == second, then 0.5*first + 0.5*first = first
	return coretools::average(baseLikelihoods[genometools::first(genotype)],
							  baseLikelihoods[genometools::second(genotype)]);
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

void TDiploidDistribution::reset() {
	_piSum.fill(0.);
	_pi = _pi_init;
}

void TDiploidDistribution::log() const {
	using coretools::str::toString;
	std::string ret;
	for (auto g = Genotype::min; g < Genotype::max; ++g) {
		ret.append(toString(g)).append(": ").append(toString(_pi[g])).append(", ");
	}
	ret.resize(ret.size() - 2);
	logfile().list(ret);
}

void TDiploidDistribution::addHeader(std::vector<std::string> &Header, std::string_view Prefix) const {
	using coretools::str::toString;
	Header.push_back(toString(Prefix, "piAA"));
	Header.push_back(toString(Prefix, "piAC"));
	Header.push_back(toString(Prefix, "piAG"));
	Header.push_back(toString(Prefix, "piAT"));
	Header.push_back(toString(Prefix, "piCC"));
	Header.push_back(toString(Prefix, "piCG"));
	Header.push_back(toString(Prefix, "piCT"));
	Header.push_back(toString(Prefix, "piGG"));
	Header.push_back(toString(Prefix, "piGT"));
	Header.push_back(toString(Prefix, "piTT"));
}

std::vector<double> TDiploidDistribution::pis() const {
	return {_pi[Genotype::AA], _pi[Genotype::AC], _pi[Genotype::AG], _pi[Genotype::AT], _pi[Genotype::CC],
			  _pi[Genotype::CG], _pi[Genotype::CT], _pi[Genotype::GG], _pi[Genotype::GT], _pi[Genotype::TT]};
}

TGenotypeLikelihoods THKY85::P_dij(const TBaseLikelihoods &baseLikelihoods) const {
	return base2genotype<genometools::Ploidy::diploid>(baseLikelihoods);
}

coretools::Probability THKY85::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// if first == second, then 0.5*first + 0.5*first = first
	return coretools::average(baseLikelihoods[genometools::first(genotype)],
							  baseLikelihoods[genometools::second(genotype)]);
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
	  _nelderMead([this](auto Vals) { return -impl::Q(coretools::expit(Vals[0]), coretools::expit(Vals[1]), coretools::expit(Vals[2]), _likelihoodSum); }) {
}

void THKY85::estimate() {
	// If likelihoodSum is totally off, it is not worth it
	constexpr size_t AA_CC = 100;
	const bool isWorthIt = _likelihoodSum[Base::A][Genotype::AA] > AA_CC*_likelihoodSum[Base::A][Genotype::CC];

	const auto lMu = coretools::logit(_mu);
	const auto lTheta_r = coretools::logit(_theta_r);
	const auto lTheta_g = coretools::logit(_theta_g);

	const stattools::TSimplex<3> simpl({lMu, lTheta_r, lTheta_g},
									   {std::abs(lMu), std::abs(lTheta_r), std::abs(lTheta_g)});

	if (isWorthIt && _nelderMead.minimize(simpl)) {
		const auto &crds = _nelderMead.coordinates();
		_mu              = coretools::expit(crds[0]);
		_theta_r         = coretools::expit(crds[1]);
		_theta_g         = coretools::expit(crds[2]);
		_pi              = impl::piTable(_mu, _theta_r, _theta_g);
	}
	_likelihoodSum.fill({});
}

void THKY85::log() const {
	logfile().list(name, ": mu=", _mu, ", theta_r=", _theta_r, ", theta_g=", _theta_g);
}

void THKY85::addHeader(std::vector<std::string> &Header, std::string_view Prefix) const {
	using coretools::str::toString;
	Header.push_back(toString(Prefix, "mu"));
	Header.push_back(toString(Prefix, "theta_r"));
	Header.push_back(toString(Prefix, "theta_g"));
	Header.push_back(toString(Prefix, "het"));
}

std::vector<double> THKY85::pis() const {
	return {_mu, _theta_r, _theta_g, impl::het(_mu, _theta_g)};
}

void THKY85::reset() {
	_likelihoodSum.fill({});
	_mu      = _mu_init;
	_theta_r = _theta_r_init;
	_theta_g = _theta_g_init;
	_pi      = impl::piTable(_mu, _theta_r, _theta_g);
}

TGenotypeLikelihoods THKY85fixedMu::P_dij(const TBaseLikelihoods &baseLikelihoods) const {
	return base2genotype<genometools::Ploidy::diploid>(baseLikelihoods);
}

coretools::Probability THKY85fixedMu::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// if first == second, then 0.5*first + 0.5*first = first
	return coretools::average(baseLikelihoods[genometools::first(genotype)],
							  baseLikelihoods[genometools::second(genotype)]);
}

double THKY85fixedMu::normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base ref) {
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


THKY85fixedMu::THKY85fixedMu()
	: _pi(impl::piTable(_mu, _theta_r, _theta_g)),
	  _nelderMead([this](auto Vals) { return -impl::Q(_mu, std::exp(Vals[0]), std::exp(Vals[1]), _likelihoodSum); }) {
}

void THKY85fixedMu::estimate() {
	// If likelihoodSum is totally off, it is not worth it
	constexpr size_t AA_CC = 100;
	const bool isWorthIt = _likelihoodSum[Base::A][Genotype::AA] > AA_CC*_likelihoodSum[Base::A][Genotype::CC];

	if (isWorthIt && _nelderMead.minimize({std::log(_theta_r), std::log(_theta_g)},
							 10.)) {
		const auto &crds = _nelderMead.coordinates();
		_theta_r         = std::exp(crds[0]);
		_theta_g         = std::exp(crds[1]);
		_pi              = impl::piTable(_mu, _theta_r, _theta_g);
	}
	_likelihoodSum.fill({});
}

void THKY85fixedMu::log() const {
	logfile().list(name, ": mu=", _mu, ", theta_r=", _theta_r, ", theta_g=", _theta_g);
}

void THKY85fixedMu::addHeader(std::vector<std::string> &Header, std::string_view Prefix) const {
	using coretools::str::toString;
	Header.push_back(toString(Prefix, "mu"));
	Header.push_back(toString(Prefix, "theta_r"));
	Header.push_back(toString(Prefix, "theta_g"));
	Header.push_back(toString(Prefix, "het"));
}

std::vector<double> THKY85fixedMu::pis() const {
	return {_mu, _theta_r, _theta_g, impl::het(_mu, _theta_g)};
}

void THKY85fixedMu::reset() {
	_likelihoodSum.fill({});
	_theta_r = _theta_r_init;
	_theta_g = _theta_g_init;
	_pi      = impl::piTable(_mu, _theta_r, _theta_g);
}


TGenotypeLikelihoods THKY85_mono::P_dij(const TBaseLikelihoods &baseLikelihoods) const {
	return base2genotype<genometools::Ploidy::haploid>(baseLikelihoods);
}

coretools::Probability THKY85_mono::getGenotypeLikelihood(const TBaseLikelihoods &baseLikelihoods,
																   Genotype genotype) const {
	// if first == second, then 0.5*first + 0.5*first = first
	return coretools::average(baseLikelihoods[genometools::first(genotype)],
							  baseLikelihoods[genometools::second(genotype)]);
}

double THKY85_mono::normalize_add(TGenotypeLikelihoods &likelihoods, genometools::Base ref) {
	double sum = 0;
	// all 10
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g    = genometools::genotype(b, b);
		likelihoods[g] *= _pi[ref][b];
		sum            += likelihoods[g];
	}
	for(auto b = Base::min; b < Base::max; ++b) {
		const auto g = genometools::genotype(b, b);
		_likelihoodSum[ref][b] += likelihoods[g].scale(sum);
	}
	return sum;
}


THKY85_mono::THKY85_mono()
	: _pi(impl::piTable(_mu, _theta)),
	  _nelderMead([this](auto Vals) { return -impl::Q(coretools::expit(Vals[0]), std::exp(Vals[1]), _likelihoodSum); }) {
}

void THKY85_mono::estimate() {
	// If likelihoodSum is totally off, it is not worth it
	constexpr size_t AA_CC = 100;
	const bool isWorthIt = _likelihoodSum[Base::A][Base::A] > AA_CC*_likelihoodSum[Base::A][Base::C];

	if (isWorthIt && _nelderMead.minimize({coretools::logit(_mu), std::log(_theta)},
							 10.)) {
		const auto &crds = _nelderMead.coordinates();
		_mu              = coretools::expit(crds[0]);
		_theta           = std::exp(crds[1]);
		_pi              = impl::piTable(_mu, _theta);
	}
	_likelihoodSum.fill({});
}

void THKY85_mono::log() const {
	logfile().list(name, ": mu=", _mu, ", theta=", _theta);
}

void THKY85_mono::addHeader(std::vector<std::string> &Header, std::string_view Prefix) const {
	using coretools::str::toString;
	Header.push_back(toString(Prefix, "mu"));
	Header.push_back(toString(Prefix, "theta"));
}

std::vector<double> THKY85_mono::pis() const {
	return {_mu, _theta};
}

void THKY85_mono::reset() {
	_likelihoodSum.fill({});
	_mu    = _mu_init;
	_theta = _theta_init;
	_pi    = impl::piTable(_mu, _theta);
}

}; // namespace GenotypeLikelihoods
