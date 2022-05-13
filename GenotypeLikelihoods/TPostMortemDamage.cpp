/*
 * TPostMortemDamage.cpp
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#include "TPostMortemDamage.h"

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>

#include "GenotypeTypes.h"
#include "TFile.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TSequencedBase.h"
#include "probability.h"


namespace GenotypeLikelihoods {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;

using namespace coretools::str;

namespace /* anonymous */ {

std::vector<double> parseParameters(const std::string &string) {
	// expect string of the form NAME[P1,P2,...]
	// extract P1, P2, ... as a vector of doubles
	std::vector<double> ps;
	if (string.find('[') == string.npos) return ps;

	std::string tmp = readAfter(string, '[');
	fillContainerFromString(extractBefore(tmp, ']'), ps, ',', true, true, true);
	return ps;
}

std::unique_ptr<TPMDFunction>initializeFunction(const std::string &pmdString) {
	// parse string to get model. Options are
	//  none
	//  Empiric[0.5,0.3,...]
	//  Exponential[a,b,c]

	const std::string name = readBefore(pmdString, '[');

	if (name == TPMDFunctionNoPMD::name) return std::make_unique<TPMDFunctionNoPMD>(pmdString);
	if (name == TPMDFunctionExponential::name) return std::make_unique<TPMDFunctionExponential>(pmdString);
	if (name == TPMDFunctionEmpiric::name) return std::make_unique<TPMDFunctionEmpiric>(pmdString);

	throw "Cannot initialize PMD function: unknown function '" + name + "'!. Use either " + TPMDFunctionNoPMD::name +
		", " + TPMDFunctionExponential::name + " or " + TPMDFunctionEmpiric::name + ".";
}

std::unique_ptr<TPMDType> createPMDType(const std::string &pmdString) {
	// split by ':'
	std::vector<std::string> details;
	fillContainerFromString(pmdString, details, ":");

	// switch type
	if (details[0] == TPMDTypeNone::name)  return std::make_unique<TPMDTypeNone>();
	if (details[0] == TPMDTypeSingleStrand::name) return std::make_unique<TPMDTypeSingleStrand>(details);
	if (details[0] == TPMDTypeDoubleStrand::name) return std::make_unique<TPMDTypeDoubleStrand>(details);

	throw "Cannot initialize PMD: unknown PMD type '" + details[0] + "'!" + "\nUse " + TPMDTypeNone::name + " or " +
		TPMDTypeSingleStrand::name + " or " + TPMDTypeDoubleStrand::name + ".";
}

} // namespace

//---------------------------------------------------------------
// TPMDFunctionNoPMD
//---------------------------------------------------------------
TPMDFunctionNoPMD::TPMDFunctionNoPMD(const std::string &string) {
	const std::vector<double> params = parseParameters(string);
	if (params.size() != 0) {
		throw "Cannot initialize PMD function '" + name + "': expected 0 but found " +
			toString(params.size()) + " parameters!";
	}
}

//---------------------------------------------------------------
// TPMDFunctionExponential
//--------------------------------------------------------------
TPMDFunctionExponential::TPMDFunctionExponential(const std::string &string) {
	constexpr size_t nParams = 4;
	std::vector<double> params = parseParameters(string);
	if (params.empty()) {
		// parameters missing: set to no PMD
		_lastPosition = 0;
		_a = _b = _c = 0.;
	} else {
		// parameters are provided
		if (params.size() != nParams) {
			throw "Cannot initialize PMD function '" + name + "': expected" +
				toString(nParams) + "(" + example + ") but found " + toString(params.size()) + " parameters!";
		}
		_lastPosition = std::lround(params[0]);
		_a = params[1];
		_b = params[2];
		_c = params[3];

		if (_lastPosition == 0)  throw "Cannot initialize PMD function '" + name + "': last position must be > 0!"; 
		if (_b < 0.0) throw "Cannot initialize PMD function '" + name + "': b must be > 0!";
		if (_c < 0.0) throw "Cannot initialize PMD function '" + name + "': c must be > 0!";
	}
	_fillPMDProbabilities();
}

void TPMDFunctionExponential::_fillPMDProbabilities() {
	_probs.resize(_lastPosition + 1);
	for (size_t p = 0; p < _probs.size(); ++p) {
		_probs[p] = _a * exp(-_b * p) + _c;
	}
}

	void TPMDFunctionExponential::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) {
	if (EstimationParameters.find(epsilon) == EstimationParameters.end()) {
		double eps = parameters().getParameterWithDefault<double>(epsilon, 0.001);
		EstimationParameters.emplace(epsilon, eps);
		logfile().list("Will consider the Newton-Raphson algorithm to have converged if the likelihood difference < " +
			      toString(eps) + ". (parameter '" + epsilon + "')");
	}

	if (EstimationParameters.find(numNR) == EstimationParameters.end()) {
		double numNRIterations = parameters().getParameterWithDefault<int>(numNR, 100);
		EstimationParameters.emplace(numNR, numNRIterations);
		logfile().list("Will run up to " + toString(numNRIterations) + " Newton-Raphson iterations. (parameter '" +
			      numNR + ")");
	}
}

void TPMDFunctionExponential::_initialEstimatesOLS(const countVec &pmdCounts, const countVec &pmdSums,
						   std::array<double, 3> &Parameters) {
	// fill vector y to fit using OLS
	const size_t N = _lastPosition + 1;
	arma::vec y(N);
	double sumYSquared = 0.0;
	for (int p = 0; p <= _lastPosition; ++p) {
		y(p) = (double)pmdCounts[p]/pmdSums[p];
		sumYSquared += y(p) * y(p);
	}

	// some variables
	double gammaStep = 0.01;
	double gammaTmp  = -gammaStep + 0.00000001;
	double SSRold    = N;
	double SSRdiff   = -1.0;
	arma::mat betaHat;
	arma::mat X(N, 2);
	X.ones();

	// do until we get a small alpha
	while (fabs(gammaStep) > 0.00000001) {
		while (SSRdiff < 0.0) {
			// update gamma
			gammaTmp += gammaStep;

			// fill x
			for (size_t p = 0; p < N; ++p) {
				X(p, 1) = exp(-gammaTmp*p);
			}
			betaHat = inv(X.t() * X) * X.t() * y;

			// calc sum of squares
			const arma::mat tmp = betaHat.t() * X.t() * y;
			const double SSRnew = sumYSquared - tmp(0, 0);
			SSRdiff             = SSRnew - SSRold;
			SSRold              = SSRnew;
		}
		// update alpha step
		gammaStep = -0.1 * gammaStep;
		SSRdiff   = -1.0;
	}
	// set parameters
	Parameters = {betaHat(0), betaHat(1), gammaTmp};
}

void TPMDFunctionExponential::_fillFAndJacobian(arma::vec &F, arma::mat &J, const countVec &pmdCounts,
						const countVec &pmdSums, const std::array<double, 3> &Parameters) {
	F.zeros();
	J.zeros();

	for (int p = 0; p <= _lastPosition; ++p) {
		// exp
		const auto expMinusAlphaP  = exp(-Parameters[2] * p);
		const auto dExpMinusAlphaP = Parameters[1] * expMinusAlphaP;

		// first term
		//----------
		const auto tmp1     = Parameters[0] + dExpMinusAlphaP;
		const auto weight1  = pmdCounts[p] / tmp1;
		const auto weightJ1 = weight1 / tmp1;

		// add to F
		F(0) += weight1;
		F(1) += weight1 * expMinusAlphaP;
		F(2) -= weight1 * p * dExpMinusAlphaP;

		// add to J -> only upper triangle, as it is symmetric
		J(0, 0) -= weightJ1;
		J(0, 1) -= weightJ1 * expMinusAlphaP;
		J(0, 2) += weightJ1 * p * dExpMinusAlphaP;
		J(1, 1) -= weightJ1 * expMinusAlphaP * expMinusAlphaP;
		J(1, 2) -= weightJ1 * p * Parameters[0] * expMinusAlphaP;
		J(2, 2) += weightJ1 * p * p * Parameters[0] * dExpMinusAlphaP;

		// second term
		//-----------
		const auto tmp2     = (1.0 - Parameters[0] - dExpMinusAlphaP);
		const auto weight2  = (pmdSums[p] - pmdCounts[p]) / tmp2;
		const auto weightJ2 = weight2 / tmp2;

		// add to F
		F(0) -= weight2;
		F(1) -= weight2 * expMinusAlphaP;
		F(2) += weight2 * p * dExpMinusAlphaP;

		// add to J -> only upper triangle, as it is symmetric
		J(0, 0) -= weightJ2;
		J(0, 1) -= weightJ2 * expMinusAlphaP;
		J(0, 2) += weightJ2 * p * dExpMinusAlphaP;
		J(1, 1) -= weightJ2 * expMinusAlphaP * expMinusAlphaP;
		J(1, 2) += weightJ2 * p * (1.0 - Parameters[0]) * expMinusAlphaP;
		J(2, 2) -= weightJ2 * p * p * (1.0 - Parameters[0]) * dExpMinusAlphaP;
	}

	// now fill in lower triangle of J
	J(1, 0) = J(0, 1);
	J(2, 0) = J(0, 2);
	J(2, 1) = J(1, 2);
}

double TPMDFunctionExponential::_calcLL(const countVec &pmdCounts, const countVec &pmdSums,
					const std::array<double, 3> &Parameters) {
	double LL = 0.0;
	for (int p = 0; p <= _lastPosition; ++p) {
		const auto dExpMinusAlphaP = Parameters[1] * exp(-Parameters[2] * p);
		LL += pmdCounts[p] * log(Parameters[0] + dExpMinusAlphaP) +
		      (pmdSums[p] - pmdCounts[p]) * log(1.0 - Parameters[0] - dExpMinusAlphaP);
	}
	return LL;
}

void TPMDFunctionExponential::_estimateWithNewtonRaphson(const countVec &pmdCounts, const countVec &pmdSums,
							 std::array<double, 3> &Parameters, uint32_t numNRIterations,
							 double epsilon) {
	// Conduct Newton-Raphson to refine
	//----------------------------------
	double oldLL = _calcLL(pmdCounts, pmdSums, Parameters);

	for (size_t _ = 0; _ < numNRIterations; ++_) {
		arma::mat J(3, 3);
		arma::vec F(3);
		arma::mat JxF;
		_fillFAndJacobian(F, J, pmdCounts, pmdSums, Parameters);

		if (!solve(JxF, J, F)) {
			std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << J << std::endl << std::endl;
			throw std::runtime_error("Issue solving JxF in TPMDTable::fitExponentialModel!");
		}

		// estimate new params
		std::array<double, 3> newParams;
		for (size_t i = 0; i < newParams.size(); ++i) newParams[i] = Parameters[i] - JxF(i);

		// calculate LL at new location
		const auto LL = _calcLL(pmdCounts, pmdSums, newParams);

		// check if we accept or backtrack
		if (LL > oldLL) {
			oldLL = LL;
			// store new params
			for (size_t i = 0; i < Parameters.size(); ++i) Parameters[i] = newParams[i];

			// check if we stop NR
			if (LL - oldLL < epsilon) break;
		}
	}
}

void TPMDFunctionExponential::learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
				    const TPMDEstimationParameters &EstimationParameters) {
	// extract counts in PMD direction and the inverse direction
	const countVec &pmdCounts = Table[from][to];
	const countVec &pmdSums   = Table.sums(from);
	const countVec &invCounts = Table[to][from];
	const countVec &invSums   = Table.sums(from);

	// find last entry with counts
	_lastPosition = pmdCounts.size() - 1;
	while (_lastPosition >= 9 && pmdCounts[_lastPosition] < 100) _lastPosition --;

	// Check if we have sufficient data
	if (_lastPosition <= 9)
		throw "Not sufficient data to fit exponential PMD model: less than the ten first positions have > 100 data "
			  "points!\nConsider pooling read groups (parameter poolReadGroups).";

	const auto last = pmdSums.cbegin() + _lastPosition;
	if (std::find(pmdSums.cbegin(), last, 0) != last)
			throw "Not sufficient data to fit exponential PMD model: no observations for some reference "
				  "alleles!<nConsider reducing the relevant length (parameter length).";

	// get initial estimates via OLS
	std::array<double, 3> Parameters;
	_initialEstimatesOLS(pmdCounts, pmdSums, Parameters);


	// run Newton-Raphson
	_estimateWithNewtonRaphson(pmdCounts, pmdSums, Parameters,
				   EstimationParameters.at(numNR),
				   EstimationParameters.at(epsilon));

	if (Parameters[1] < 0) {
		throw "Estimation resulted in a < 0!\nThis is likely due to limited data. Consider pooling read groups "
			  "(parameter poolReadGroups).";
	}

	// transform parameters
	// the exponential PMD model is f(C->T) = mu + (1-mu) *[ a*exp(-b * position) + c ]
	// but we fitted f(C->T) = alpha + beta * exp(-gamma * position).
	// Hence we have:
	//  mu is estimated from T->C transitions
	//  a =  beta / (1 - mu)
	//  b = gamma
	//  c = (alpha - mu) / (1 - mu)
	const double mu = std::accumulate(invCounts.cbegin(), invCounts.cbegin() + _lastPosition + 1, 0.) /
			  std::accumulate(invSums.cbegin(), invSums.cbegin() + _lastPosition + 1, 0);

	// store parameters, including lastPosition
	_a = Parameters[1] / (1.0 - mu);
	_b = Parameters[2];
	_c = (Parameters[0] - mu) / (1.0 - mu);

	_fillPMDProbabilities();


	// check if pattern is negativ
	if (_probs[_lastPosition] < 0) {
		throw "Estimation resulted in negative PMD at high positions!\nThis is likely be due to limited data. Consider "
			  "pooling read groups (parameter poolReadGroups).";
	}
}

double TPMDFunctionExponential::prob(uint16_t pos) const noexcept {
	// Note: distance is zero based!
	// model is fit up to _lastPosition. We assume constant PMD after that
	return pos < _lastPosition ? _probs[pos] : _probs[_lastPosition];
}

//---------------------------------------------------------------
// TPMDFunctionEmpiric
//---------------------------------------------------------------
TPMDFunctionEmpiric::TPMDFunctionEmpiric(const std::string &string) : _parameters(parseParameters(string)) {
	if (_parameters.empty()) {
		// parameters missing: set to no PMD
		_parameters = {0.0};
	} else {
		// parameters are provided
		for (auto &d : _parameters) {
			if (d < 0.0 || d > 1.0) {
				throw "Cannot initialize post mortem damage function '" + name +
					"': some probabilities are outside [0,1]!";
			}
		}
	}
}

void TPMDFunctionEmpiric::learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
				const TPMDEstimationParameters &) {
	// resize parameters
	_parameters.resize(Table.size()); // include extra bin for sites beyond size (available in PMDTables)

	// extract counts in PMD direction and the inverse direction
	const countVec &forwardCounts  = Table[from][to]; //e.g. C -> T
	const countVec &forwardSums    = Table.sums(from);
	const countVec &backwardCounts = Table[to][from]; //e.g. T -> C
	const countVec &backwardSums   = Table.sums(to);

	for (size_t p = 0; p < _parameters.size(); ++p) {
		if (forwardSums[p] == 0 || backwardSums[p] == 0) {
			_parameters[p] = 0.0;
		} else {
			double forward  = (double)forwardCounts[p] / forwardSums[p]; //e.g. C -> T
			double backward = (double)backwardCounts[p] / backwardSums[p]; //e.g. T -> C

			//forward = mu_CT + (1 - mu_CT) * PMD; mu_CT = mu_TC = backward
			_parameters[p] = std::max(0.0, (forward - backward) / (1.0 - backward));
		}
	}
}

double TPMDFunctionEmpiric::prob(uint16_t pos) const noexcept {
	return pos < _parameters.size() ? _parameters[pos] : _parameters.back();
}

//------------------------------------------------------
// TPMDDoubleStrand
//------------------------------------------------------

TPMDTypeDoubleStrand::TPMDTypeDoubleStrand(const std::vector<std::string> &Details) {
	// expect three elements: type, pmdCT, pmdGA
	constexpr size_t nDetails = 3;
	if (Details.size() != nDetails) {
		throw "Cannot initialize PMD type " + name + ": expect " +
			std::to_string(nDetails) + " entries but found " + toString(Details.size()) + "!" + "\nProvided string: '" +
			concatenateString(Details, ':') + "'." + "\nExpect string of the form '" + name +
			"':functionCT:functionGA'.";
	}
	_pmdCT = initializeFunction(Details[1]);
	_pmdGA = initializeFunction(Details[2]);
}

void TPMDTypeDoubleStrand::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters){
	_pmdCT->parseEstimationParameters(EstimationParameters);
	_pmdGA->parseEstimationParameters(EstimationParameters);
}

void TPMDTypeDoubleStrand::estimate(const PMDTable_RG &PMDTable,
				    const TPMDEstimationParameters &EstimationParameters) {
	using genometools::Base;
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: C->T pattern is the same for forward and reverse reads from their respective 5-prime ends.
	TPMDTable from5(PMDTable[forward5]);
	from5.add(PMDTable[reverse5]);
	_pmdCT->learn(from5, Base::C, Base::T, EstimationParameters);

	// Assumption: G->A pattern is the same for forward and reverse reads from their respective 3-prime ends.
	TPMDTable from3(PMDTable[forward3]);
	from3.add(PMDTable[reverse3]);
	_pmdGA->learn(from3, Base::G, Base::A, EstimationParameters);
}

void TPMDTypeDoubleStrand::fillBaseLikelihoods(const BAM::TSequencedBase &base,
					       const TBaseLikelihoods &baseLikelihoodsNoPMD,
					       TBaseLikelihoods &baseLikelihoods) const {
	using genometools::Base;
	// Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	// no PMD for A and C
	baseLikelihoods[Base::A] = baseLikelihoodsNoPMD[Base::A];
	baseLikelihoods[Base::T] = baseLikelihoodsNoPMD[Base::T];

	// get relevant PMD probabilities
	const auto from3  = base.distFrom3Prime < base.distFrom5Prime;
	double pmdProb_CT = 0;
	double pmdProb_GA = 0;
	if (!base.isReverseStrand()) {
		if (from3)
			pmdProb_GA = _pmdGA->prob(base.distFrom3Prime);
		else
			pmdProb_CT = _pmdCT->prob(base.distFrom5Prime);
	} else {
		// ??? Newest insight
		if (from3)
			pmdProb_CT = _pmdGA->prob(base.distFrom3Prime);
		else
			pmdProb_GA = _pmdCT->prob(base.distFrom5Prime);
	}

	// add PMD
	baseLikelihoods[Base::C] =
		(1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[Base::C].get() + pmdProb_CT * baseLikelihoodsNoPMD[Base::T].get();
	baseLikelihoods[Base::G] =
		(1.0 - pmdProb_GA) * baseLikelihoodsNoPMD[Base::G].get() + pmdProb_GA * baseLikelihoodsNoPMD[Base::A].get();
}

void TPMDTypeDoubleStrand::simulate(BAM::TSequencedBase &base) const {
	simulate(base.base, base.distFrom5Prime, base.distFrom3Prime, base.isReverseStrand());
}

void TPMDTypeDoubleStrand::simulate(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
				       const bool &IsReverseStrand) const {
	using genometools::Base;
	// simulate PMD
	if (!IsReverseStrand) {
		// forward strand
		if (base == Base::C && randomGenerator().getRand() < _pmdCT->prob(DistFrom5Prime)) {
			base = Base::T;
		} else if (base == Base::G && randomGenerator().getRand() < _pmdGA->prob(DistFrom3Prime)) {
			base = Base::A;
		}
	} else {
		// reverse strand
		if (base == Base::C && randomGenerator().getRand() < _pmdGA->prob(DistFrom3Prime)) {
			// ??? Newest insight
			base = Base::T;
		} else if (base == Base::G && randomGenerator().getRand() < _pmdCT->prob(DistFrom5Prime)) {
			base = Base::A;
		}
	}
}
//------------------------------------------------------
// TPMDSingleStrand
//------------------------------------------------------

TPMDTypeSingleStrand::TPMDTypeSingleStrand(const std::vector<std::string> &Details) {
	// expect 2 elements: type, pmdCT
	constexpr size_t nDetails = 3;
	if (Details.size() != nDetails) {
		throw "Cannot initialize PMD type " + name + ": expect " +
			std::to_string(nDetails) + " entries but found " + toString(Details.size()) + "!" + "\nProvided string: '" +
			concatenateString(Details, ':') + "'." + "\nExpect string of the form '" + name +
			"':functionCT:functionGA'.";
	}
	_pmdCT3 = initializeFunction(Details[1]);
	_pmdCT5 = initializeFunction(Details[2]);
}

void TPMDTypeSingleStrand::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) {
	_pmdCT3->parseEstimationParameters(EstimationParameters);
	_pmdCT5->parseEstimationParameters(EstimationParameters);
}

void TPMDTypeSingleStrand::estimate(const PMDTable_RG &PMDTable,
				    const TPMDEstimationParameters &EstimationParameters) {
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: 5-prime C->T pattern is the same for forward and reverse reads from their respective
	// 5-prime ends.
	using genometools::Base;
	TPMDTable from5(PMDTable[forward5]);
	from5.add(PMDTable[reverse5]);
	_pmdCT5->learn(from5, Base::C, Base::T, EstimationParameters);

	// Assumption: 3-prime C->T pattern is the same for forward and reverse reads from their
	// respective 3-prime ends.
	TPMDTable from3(PMDTable[forward3]);
	from3.add(PMDTable[reverse3]);
	// ??? G->A  or C->T (reversed gets flipped when read)
	_pmdCT3->learn(from3, Base::C, Base::T, EstimationParameters);
}

void TPMDTypeSingleStrand::fillBaseLikelihoods(const BAM::TSequencedBase &base,
					       const TBaseLikelihoods &baseLikelihoodsNoPMD,
					       TBaseLikelihoods &baseLikelihoods) const {
	using genometools::Base;
	// Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	// no PMD for A, C and G
	baseLikelihoods[Base::A] = baseLikelihoodsNoPMD[Base::A];
	baseLikelihoods[Base::T] = baseLikelihoodsNoPMD[Base::T];
	baseLikelihoods[Base::G] = baseLikelihoodsNoPMD[Base::G];

	// get relevant PMD probabilities
	const double pmdProb_CT = base.distFrom3Prime < base.distFrom5Prime ? _pmdCT3->prob(base.distFrom3Prime)
									    : _pmdCT5->prob(base.distFrom5Prime);

	// add PMD
	baseLikelihoods[Base::C] = (1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[Base::C].get() +
					  pmdProb_CT * baseLikelihoodsNoPMD[Base::T].get();
}

void TPMDTypeSingleStrand::simulate(BAM::TSequencedBase &base) const {
	simulate(base.base, base.distFrom5Prime, base.distFrom3Prime, base.isReverseStrand());
}

void TPMDTypeSingleStrand::simulate(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
				       const bool &) const {
	using genometools::Base;
	if (!(base == Base::C)) return;

	// simulate PMD
	if (DistFrom3Prime < DistFrom5Prime) {
		if (randomGenerator().getRand() < _pmdCT3->prob(DistFrom3Prime)) { base = Base::T; }
	} else {
		if (randomGenerator().getRand() < _pmdCT5->prob(DistFrom5Prime)) { base = Base::T; }
	}
}

//------------------------------------------------------
// TPostMortemDamage
//------------------------------------------------------

void TPostMortemDamage::writeToFile(const BAM::TReadGroups &ReadGroups, const std::string filename) const {
	std::vector<std::string> header = {"ReadGroup", "Type", "Functions"};
	coretools::TOutputFile out(filename, header);

	// write for each read group
	for (auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r)
		out << r->name_ID << _pmdObjects[r->id]->typeString() << _pmdObjects[r->id]->functionString() << std::endl;
}

void TPostMortemDamage::writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
				    const std::string filename) const {
	std::vector<std::string> header = {"ReadGroup", "Type", "Functions"};
	coretools::TOutputFile out(filename, header);

	// write for each read group
	for (auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r)
		out << r->name_ID << _pmdObjects[r->id]->typeString()
			<< _pmdObjects[ReadGroupMap.pooledIndex(r->id)]->functionString() << std::endl;
}

void TPostMortemDamage::_initializeFromString(const std::string &pmdString) {
	// not a file: initialize all read groups have the same pmd
	logfile().startIndent("PMD function used for all read groups:");

	for (auto &p : _pmdObjects) { p = createPMDType(pmdString); }

	// report
	logfile().list(_pmdObjects[0]->functionString());
	logfile().endIndent();
}

void TPostMortemDamage::_initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename,
					    std::vector<uint16_t> &ReadGroupsWithoutPMD) {
	// create an array of TPMD objects for each read group
	// also works if no parameters are provided (e.g. for estimation)
	// read from file for each read group

	logfile().listFlush("Initializing PMD from file '" + filename + "' ...");
	coretools::TInputFile in(filename, {"readGroup", "pmd"}, "/t", "//");

	// parse file that has structure: ReadGroup, Type, Functions
	std::vector<std::string> vec;
	while (in.read(vec)) {
		if (ReadGroups.readGroupExists(vec[0])) { // ignore if it does not exist
			// get read group
			uint16_t readGroupId = ReadGroups.getId(vec[0]);

			// create type
			_pmdObjects[readGroupId] = createPMDType(vec[1]);
		}
	}
	logfile().done();

	// check if we have a function for all read groups
	// create no-PMD types for all remaining ones and return their indexes
	for (uint16_t i = 0; i < ReadGroups.size(); ++i) {
		if (!_pmdObjects[i]) {
			_pmdObjects[i] = createPMDType("non");
			ReadGroupsWithoutPMD.push_back(i);
		}
	}
}

void TPostMortemDamage::_setHasDamage() {
	// check if there is PMD for at least one read group
	_hasPMD = false;
	for (auto &p : _pmdObjects) {
		if (p->hasDamage()) {
			_hasPMD = true;
			break;
		}
	}
}

void TPostMortemDamage::initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups,
				   std::vector<uint16_t> &ReadGroupsWithoutPMD) {
	if (_hasPMD) {
		throw std::runtime_error(
			"void TPostMortemDamage::initialize(const std::string & pmdString, const BAM::TReadGroups & ReadGroups, "
			"TLog* Logfile, std::vector<uint16_t> & ReadGroupsWithoutPMD): Models already initialized!");
	}

	// prepare objects
	ReadGroupsWithoutPMD.clear();
	_pmdObjects.resize(ReadGroups.size());

	// expect either a file name or a type of the form "type:functions"
	// check if it is a file
	if (stringContains(pmdString, ":")) {
		_initializeFromString(pmdString);
	} else {
		// probably a file
		_initializeFromFile(ReadGroups, pmdString, ReadGroupsWithoutPMD);
	}

	// check if there is PMD for at least one read group
	_setHasDamage();
}

void TPostMortemDamage::fillBaseLikelihoods(const BAM::TSequencedBase &base,
                                            const TBaseLikelihoods &baseLikelihoodsNoPMD,
                                            TBaseLikelihoods &baseLikelihoods) const {
	if (_hasPMD) {
		_pmdObjects[base.readGroupID]->fillBaseLikelihoods(base, baseLikelihoodsNoPMD, baseLikelihoods);
	} else {
		// just copy
		baseLikelihoods = baseLikelihoodsNoPMD;
	}
}

}; // namespace GenotypeLikelihoods
