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

namespace GenotypeLikelihoods {

using namespace coretools::str;

namespace /* anonymous */ {

std::vector<double> parseParameters(const std::string &string) {
	// expect string of the form NAME[P1,P2,...]
	// extract P1, P2, ... as a vector of doubles
	std::string tmp = readAfter(string, '[');
	std::vector<double> ps;
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
		_lastPosition = 1;
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

void TPMDFunctionExponential::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters,
							TParameters &Params, TLog *Logfile) {
	if (EstimationParameters.find(epsilon) == EstimationParameters.end()) {
		double eps = Params.getParameterWithDefault<double>(epsilon, 0.001);
		EstimationParameters.emplace(epsilon, eps);
		Logfile->list("Will consider the Newton-Raphson algorithm to have converged if the likelihood difference < " +
			      toString(eps) + ". (parameter '" + epsilon + "')");
	}

	if (EstimationParameters.find(numNR) == EstimationParameters.end()) {
		double numNRIterations = Params.getParameterWithDefault<int>(numNR, 100);
		EstimationParameters.emplace(numNR, numNRIterations);
		Logfile->list("Will run up to " + toString(numNRIterations) + " Newton-Raphson iterations. (parameter '" +
			      numNR + ")");
	}
}

void TPMDFunctionExponential::_initialEstimatesOLS(const countVec &pmdCounts, const countVec &pmdSums,
						   std::vector<double> &Parameters) {
	// fill vector y to fit using OLS
	arma::vec y(_lastPosition + 1);
	double sumYSquared = 0.0;
	for (int p = 0; p <= _lastPosition; ++p) {
		y(p) = (double)pmdCounts[p] / (double)pmdSums[p];
		sumYSquared += y(p) * y(p);
	}

	// some variables
	double gammaStep = 0.01;
	double gammaTmp  = -gammaStep + 0.00000001;
	double SSRold    = pmdCounts.size() + 1;
	double SSRdiff   = -1.0;
	arma::mat betaHat;
	arma::mat X(pmdCounts.size() + 1, 2);
	X.ones();

	// do until we get a small alpha
	while (fabs(gammaStep) > 0.00000001) {
		while (SSRdiff < 0.0) {
			// update gamma
			gammaTmp += gammaStep;

			// fill x
			for (size_t p = 0; p < pmdCounts.size(); ++p) { X(p, 1) = exp(-gammaTmp * (double)p); }

			betaHat = inv(X.t() * X) * X.t() * y;

			// calc sum of squares
			arma::mat tmp = (betaHat.t() * X.t() * y);
			double SSRnew = sumYSquared - tmp(0, 0);
			SSRdiff       = SSRnew - SSRold;
			SSRold        = SSRnew;
		}

		// update alpha step
		gammaStep = -0.1 * gammaStep;
		SSRdiff   = -1.0;
	}

	// set parameters
	Parameters = {betaHat(0), betaHat(1), gammaTmp};
}

void TPMDFunctionExponential::_fillFAndJacobian(arma::vec &F, arma::mat &J, const countVec &pmdCounts,
						const countVec &pmdSums, const std::vector<double> &Parameters) {
	F.zeros();
	J.zeros();

	double weight, weightJ, tmp;
	double expMinusAlphaP;
	double dExpMinusAlphaP;
	for (int p = 0; p <= _lastPosition; ++p) {
		// exp
		expMinusAlphaP  = exp(-Parameters[2] * p);
		dExpMinusAlphaP = Parameters[1] * expMinusAlphaP;

		// first term
		//----------
		tmp     = Parameters[0] + dExpMinusAlphaP;
		weight  = pmdCounts[p] / tmp;
		weightJ = weight / tmp;

		// add to F
		F(0) += weight;
		F(1) += weight * expMinusAlphaP;
		F(2) -= weight * p * dExpMinusAlphaP;

		// add to J -> only upper triangle, as it is symmetric
		J(0, 0) -= weightJ;
		J(0, 1) -= weightJ * expMinusAlphaP;
		J(0, 2) += weightJ * p * dExpMinusAlphaP;
		J(1, 1) -= weightJ * expMinusAlphaP * expMinusAlphaP;
		J(1, 2) -= weightJ * p * Parameters[0] * expMinusAlphaP;
		J(2, 2) += weightJ * p * p * Parameters[0] * dExpMinusAlphaP;

		// second term
		//-----------
		tmp     = (1.0 - Parameters[0] - dExpMinusAlphaP);
		weight  = (pmdSums[p] - pmdCounts[p]) / tmp;
		weightJ = weight / tmp;

		// add to F
		F(0) -= weight;
		F(1) -= weight * expMinusAlphaP;
		F(2) += weight * p * dExpMinusAlphaP;

		// add to J -> only upper triangle, as it is symmetric
		J(0, 0) -= weightJ;
		J(0, 1) -= weightJ * expMinusAlphaP;
		J(0, 2) += weightJ * p * dExpMinusAlphaP;
		J(1, 1) -= weightJ * expMinusAlphaP * expMinusAlphaP;
		J(1, 2) += weightJ * p * (1.0 - Parameters[0]) * expMinusAlphaP;
		J(2, 2) -= weightJ * p * p * (1.0 - Parameters[0]) * dExpMinusAlphaP;
	}

	// now fill in lower triangle of J
	J(1, 0) = J(0, 1);
	J(2, 0) = J(0, 2);
	J(2, 1) = J(1, 2);
}

double TPMDFunctionExponential::_calcLL(const countVec &pmdCounts, const countVec &pmdSums,
					const std::vector<double> &Parameters) {
	double LL = 0.0;
	for (int p = 0; p <= _lastPosition; ++p) {
		double dExpMinusAlphaP = Parameters[1] * exp(-Parameters[2] * p);
		LL += pmdCounts[p] * log(Parameters[0] + dExpMinusAlphaP) +
		      (pmdSums[p] - pmdCounts[p]) * log(1.0 - Parameters[0] - dExpMinusAlphaP);
	}
	return LL;
}

void TPMDFunctionExponential::_estimateWithNewtonRaphson(const countVec &pmdCounts, const countVec &pmdSums,
							 std::vector<double> &Parameters, uint32_t numNRIterations,
							 double epsilon) {
	// variables
	arma::mat J(3, 3);
	arma::vec F(3);
	arma::mat JxF;

	// set starting values
	std::vector<double> newParams(3);

	// Conduct Newton-Raphson to refine
	//----------------------------------
	double oldLL = _calcLL(pmdCounts, pmdSums, Parameters);

	for (size_t i = 0; i < numNRIterations; ++i) {
		_fillFAndJacobian(F, J, pmdCounts, pmdSums, Parameters);
		if (solve(JxF, J, F)) {
			// estimate new params
			for (int x = 0; x < 3; ++x) { newParams[x] = Parameters[x] - JxF(x); }

			// calculate LL at new location
			double LL = _calcLL(pmdCounts, pmdSums, newParams);

			// check if we accept or backtrack
			if (LL > oldLL) {
				// store new params
				for (int x = 0; x < 3; ++x) { Parameters = newParams; }

				// check if we stop NR
				if (LL - oldLL < epsilon) {
					oldLL = LL;
					break;
				}
				oldLL = LL;
			}
		} else {
			std::cout << std::endl << std::endl << "JACOBIAN:" << std::endl << J << std::endl << std::endl;
			throw std::runtime_error("Issue solving JxF in TPMDTable::fitExponentialModel!");
		}
	}
}

void TPMDFunctionExponential::learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
				    const TPMDEstimationParameters &EstimationParameters) {
	// extract counts in PMD direction and the inverse direction
	const countVec &pmdCounts = Table[from][to.get()];
	const countVec &pmdSums   = Table.sums(from);
	const countVec &invCounts = Table[to][from.get()];
	const countVec &invSums   = Table.sums(from);

	// Check if we have sufficient data
	// find last entry with counts
	_lastPosition = -1;
	for (int p = pmdCounts.size() - 1; p >= 0; --p) {
		if (pmdSums[p] > 100) {
			_lastPosition = p;
			break;
		}
	}

	if (_lastPosition < 10)
		throw "Not sufficient data to fit exponential PMD model: less than the ten first positions have > 100 data "
			  "points!\nConsider pooling read groups (parameter poolReadGroups).";
	for (int p = 0; p < _lastPosition; ++p) {
		if (pmdSums[p] == 0) {
			throw "Not sufficient data to fit exponential PMD model: no observations for some reference "
				  "alleles!<nConsider reducing the relevant length (parameter length).";
		}
	}

	// get initial estimates via OLS
	std::vector<double> Parameters;
	_initialEstimatesOLS(pmdCounts, pmdSums, Parameters);

	// run Newton-Raphson
	_estimateWithNewtonRaphson(pmdCounts, pmdSums, Parameters,
				   EstimationParameters.at(epsilon),
				   EstimationParameters.at(numNR));

	// transform parameters
	// the exponential PMD model is f(C->T) = mu + (1-mu) *[ a*exp(-b * position) + c ]
	// but we fitted f(C->T) = alpha + beta * exp(-gamma * position).
	// Hence we have:
	//  mu is estimated from T->C transitions
	//  a =  beta / (1 - mu)
	//  b = gamma
	//  c = (alpha - mu) / (1 - mu)
	double mu  = 0.0;
	double sum = 0.0;
	for (int p = 0; p <= _lastPosition; ++p) {
		mu += invCounts[p];
		sum += invSums[p];
	}
	mu = mu / sum;

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
	if (Parameters[1] < 0) {
		throw "Estimation resulted in a < 0!\nThis is likely due to limited data. Consider pooling read groups "
			  "(parameter poolReadGroups).";
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
	_parameters.resize(Table.size() + 1); // include extra bin for sites beyond size (available in PMDTables)

	// extract counts in PMD direction and the inverse direction
	const countVec &pmdCounts = Table[from][to.get()];
	const countVec &pmdSums   = Table.sums(from);
	const countVec &invCounts = Table[to][from.get()];
	const countVec &invSums   = Table.sums(from);

	for (size_t p = 0; p <= _parameters.size(); ++p) {
		if (pmdSums[p] == 0 || invSums[p] == 0) {
			_parameters[p] = 0.0;
		} else {
			double pmd = (double)pmdCounts[p] / pmdSums[p];
			double inv = (double)invCounts[p] / invSums[p];

			_parameters[p] = std::max(0.0, (pmd - inv) / (1.0 - inv));
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

std::string TPMDTypeDoubleStrand::functionString() const noexcept {
	return name + ":" + _pmdCT->string() + ":" + _pmdGA->string();
}

void TPMDTypeDoubleStrand::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters,
						     TParameters &Params, TLog *Logfile) {
	_pmdCT->parseEstimationParameters(EstimationParameters, Params, Logfile);
	_pmdGA->parseEstimationParameters(EstimationParameters, Params, Logfile);
}

void TPMDTypeDoubleStrand::estimate(const PMDTable_RG &PMDTable,
				    const TPMDEstimationParameters &EstimationParameters) {
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: C->T pattern is the same for forward and reverse reads from their respective 5-prime ends.
	TPMDTable from5(PMDTable[forward5]);
	from5.add(PMDTable[reverse5]);
	_pmdCT->learn(from5, genometools::C, genometools::T, EstimationParameters);

	// Assumption: G->A pattern is the same for forward and reverse reads from their respective 3-prime ends.
	TPMDTable from3(PMDTable[forward3]);
	from3.add(PMDTable[reverse3]);
	_pmdGA->learn(from3, genometools::G, genometools::A, EstimationParameters);
}

void TPMDTypeDoubleStrand::fillBaseLikelihoods(const BAM::TSequencedBase &base,
					       const TBaseProbabilities &baseLikelihoodsNoPMD,
					       TBaseProbabilities &baseLikelihoods) const {
	// Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	// no PMD for A and C
	baseLikelihoods[genometools::A] = baseLikelihoodsNoPMD[genometools::A];
	baseLikelihoods[genometools::T] = baseLikelihoodsNoPMD[genometools::T];

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
	baseLikelihoods[genometools::C] = (1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[genometools::C].get() +
					  pmdProb_CT * baseLikelihoodsNoPMD[genometools::T].get();
	baseLikelihoods[genometools::G] = (1.0 - pmdProb_GA) * baseLikelihoodsNoPMD[genometools::G].get() +
					  pmdProb_GA * baseLikelihoodsNoPMD[genometools::A].get();
}

void TPMDTypeDoubleStrand::simulatePMD(BAM::TSequencedBase &base, TRandomGenerator &RandomGenerator) const {
	simulatePMD(base.base, base.distFrom5Prime, base.distFrom3Prime, base.isReverseStrand(), RandomGenerator);
}

void TPMDTypeDoubleStrand::simulatePMD(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
				       const bool &IsReverseStrand, TRandomGenerator &RandomGenerator) const {
	// simulate PMD
	if (!IsReverseStrand) {
		// forward strand
		if (base == genometools::C) {
			if (RandomGenerator.getRand() < _pmdCT->prob(DistFrom5Prime)) { base = genometools::T; }
		} else if (base == genometools::G) {
			if (RandomGenerator.getRand() < _pmdGA->prob(DistFrom3Prime)) { base = genometools::A; }
		}
	} else {
		// reverse strand
		if (base == genometools::C) {
			// ??? Newest insight
			if (RandomGenerator.getRand() < _pmdGA->prob(DistFrom3Prime)) { base = genometools::T; }
		} else if (base == genometools::G) {
			if (RandomGenerator.getRand() < _pmdCT->prob(DistFrom5Prime)) { base = genometools::A; }
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

std::string TPMDTypeSingleStrand::functionString() const noexcept {
	return name + ":" + _pmdCT3->string() + ":" + _pmdCT5->string();
}

void TPMDTypeSingleStrand::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters,
						     TParameters &Params, TLog *Logfile) {
	_pmdCT3->parseEstimationParameters(EstimationParameters, Params, Logfile);
	_pmdCT5->parseEstimationParameters(EstimationParameters, Params, Logfile);
}

void TPMDTypeSingleStrand::estimate(const PMDTable_RG &PMDTable,
				    const TPMDEstimationParameters &EstimationParameters) {
	// Note: TPMDTables stores bases as during sequencing (not as after mapping)
	// Assumption: 5-prime C->T pattern is the same for forward and reverse reads from their respective
	// 5-prime ends.
	TPMDTable from5(PMDTable[forward5]);
	from5.add(PMDTable[reverse5]);
	_pmdCT5->learn(from5, genometools::C, genometools::T, EstimationParameters);

	// Assumption: 3-prime C->T pattern is the same for forward and reverse reads from their
	// respective 3-prime ends.
	TPMDTable from3(PMDTable[forward3]);
	from3.add(PMDTable[reverse3]);
	// ??? G->A  or C->T (reversed gets flipped when read)
	_pmdCT3->learn(from3, genometools::C, genometools::T, EstimationParameters);
}

void TPMDTypeSingleStrand::fillBaseLikelihoods(const BAM::TSequencedBase &base,
					       const TBaseProbabilities &baseLikelihoodsNoPMD,
					       TBaseProbabilities &baseLikelihoods) const {
	// Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	// no PMD for A, C and G
	baseLikelihoods[genometools::A] = baseLikelihoodsNoPMD[genometools::A];
	baseLikelihoods[genometools::T] = baseLikelihoodsNoPMD[genometools::T];
	baseLikelihoods[genometools::G] = baseLikelihoodsNoPMD[genometools::G];

	// get relevant PMD probabilities
	const double pmdProb_CT = (base.distFrom3Prime < base.distFrom5Prime ? _pmdCT3->prob(base.distFrom3Prime)
									     : _pmdCT5->prob(base.distFrom5Prime));

	// add PMD
	baseLikelihoods[genometools::C] = (1.0 - pmdProb_CT) * baseLikelihoodsNoPMD[genometools::C].get() +
					  pmdProb_CT * baseLikelihoodsNoPMD[genometools::T].get();
}

void TPMDTypeSingleStrand::simulatePMD(BAM::TSequencedBase &base, TRandomGenerator &RandomGenerator) const {
	simulatePMD(base.base, base.distFrom5Prime, base.distFrom3Prime, base.isReverseStrand(), RandomGenerator);
}

void TPMDTypeSingleStrand::simulatePMD(genometools::Base &base, uint16_t DistFrom5Prime, uint16_t DistFrom3Prime,
				       const bool &, TRandomGenerator &RandomGenerator) const {
	if (!(base == genometools::C)) return;

	// simulate PMD
	if (DistFrom3Prime < DistFrom5Prime) {
		if (RandomGenerator.getRand() < _pmdCT3->prob(DistFrom3Prime)) { base = genometools::T; }
	} else {
		if (RandomGenerator.getRand() < _pmdCT5->prob(DistFrom5Prime)) { base = genometools::T; }
	}
}

//------------------------------------------------------
// TPostMortemDamage
//------------------------------------------------------

void TPostMortemDamage::writeToFile(const BAM::TReadGroups &ReadGroups, const std::string filename) const {
	std::vector<std::string> header = {"ReadGroup", "Type", "Functions"};
	coretools::TOutputFile out(filename, header);

	// write for each read group
	for (auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r) {
		out << r->name_ID << _pmdObjects[r->id]->functionString() << std::endl;
	}
}

void TPostMortemDamage::writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
				    const std::string filename) const {
	std::vector<std::string> header = {"ReadGroup", "Type", "Functions"};
	coretools::TOutputFile out(filename, header);

	// write for each read group
	for (auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r) {
		out << r->name_ID << _pmdObjects[ReadGroupMap.pooledIndex(r->id)]->functionString() << std::endl;
	}
}

void TPostMortemDamage::_createPMDType(const std::string &pmdString, std::shared_ptr<TPMDType> &ptr) {
	// split by ':'
	std::vector<std::string> details;
	fillContainerFromString(pmdString, details, ":");

	// switch type
	if (details[0] == TPMDTypeNone::name) {
		ptr = std::make_shared<TPMDTypeNone>();
	} else if (details[0] == TPMDTypeSingleStrand::name) {
		ptr = std::make_shared<TPMDTypeSingleStrand>(details);
	} else if (details[0] == TPMDTypeDoubleStrand::name) {
		ptr = std::make_shared<TPMDTypeDoubleStrand>(details);
	} else {
		throw "Cannot initialize PMD: unknown PMD type '" + details[0] + "'!" + "\nUse " + TPMDTypeNone::name+ " or " +
			TPMDTypeSingleStrand::name + " or " + TPMDTypeDoubleStrand::name+ ".";
	}
}

void TPostMortemDamage::_initializeFromString(const std::string &pmdString, TLog *logfile) {
	// not a file: initialize all read groups have the same pmd
	logfile->startIndent("PMD function used for all read groups:");

	for (auto &p : _pmdObjects) { _createPMDType(pmdString, p); }

	// report
	logfile->list(_pmdObjects[0]->functionString());
	logfile->endIndent();
}

void TPostMortemDamage::_initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename,
					    TLog *logfile, std::vector<uint16_t> &ReadGroupsWithoutPMD) {
	// create an array of TPMD objects for each read group
	// also works if no parameters are provided (e.g. for estimation)
	// read from file for each read group

	logfile->listFlush("Initializing PMD from file '" + filename + "' ...");
	coretools::TInputFile in(filename, {"readGroup", "pmd"}, "/t", "//");

	// parse file that has structure: ReadGroup, Type, Functions
	std::vector<std::string> vec;
	while (in.read(vec)) {
		if (ReadGroups.readGroupExists(vec[0])) { // ignore if it does not exist
			// get read group
			uint16_t readGroupId = ReadGroups.getId(vec[0]);

			// create type
			_createPMDType(vec[1], _pmdObjects[readGroupId]);
		}
	}
	logfile->done();

	// check if we have a function for all read groups
	// create no-PMD types for all remaining ones and return their indexes
	for (uint16_t i = 0; i < ReadGroups.size(); ++i) {
		if (!_pmdObjects[i]) {
			_createPMDType("non", _pmdObjects[i]);
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

void TPostMortemDamage::initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups, TLog *Logfile,
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
		_initializeFromString(pmdString, Logfile);
	} else {
		// probably a file
		_initializeFromFile(ReadGroups, pmdString, Logfile, ReadGroupsWithoutPMD);
	}

	// check if there is PMD for at least one read group
	_setHasDamage();
}

void TPostMortemDamage::fillBaseLikelihoods(const BAM::TSequencedBase &base,
                                            const TBaseProbabilities &baseLikelihoodsNoPMD,
                                            TBaseProbabilities &baseLikelihoods) const {
	if (_hasPMD) {
		_pmdObjects[base.readGroupID]->fillBaseLikelihoods(base, baseLikelihoodsNoPMD, baseLikelihoods);
	} else {
		// just copy
		baseLikelihoods = baseLikelihoodsNoPMD;
	}
}

}; // namespace GenotypeLikelihoods
