/*
 * TEpsilon.cpp
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#include "TEpsilon.h"
#include "armadillo"

namespace GenotypeLikelihoods {
namespace SequencingError {

void TEpsilon::_addToQ(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds,
					   const TGenotypeFloats &P_bbar_I_gds, bool isInvariant) {
	const double eps    = calcErrorRate(data);
	const double eps_c  = 1. - eps;
	const double leps   = std::log(eps);
	const double leps_c = std::log(eps_c);

	for (auto g : Genotypes::homo) {
		const auto P_bbar_I_gd = P_bbar_I_gds[g];
		const auto P_g_I_d     = P_g_I_ds[g];
		_Q += P_g_I_d * (P_bbar_I_gd * leps_c + (1. - P_bbar_I_gd) * leps);
	}
	if (!isInvariant) {
		for (auto g : Genotypes::het) {
			const auto P_bbar_I_gd = P_bbar_I_gds[g];
			const auto P_g_I_d     = P_g_I_ds[g];
			_Q += P_g_I_d * (P_bbar_I_gd * leps_c + (1. - P_bbar_I_gd) * leps);
		}
	}
}

void TEpsilon::_addToQJF(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds,
						 const TGenotypeFloats &P_bbar_I_gds, bool isInvariant) {
	static std::vector<T1stDerivative> der1st;
	static std::vector<T2ndDerivative> der2nd;
	der1st.clear();
	der2nd.clear();
	// get error rate
	const double eps    = _calcErrorRate(data, der1st, der2nd);
	const double eps_c  = 1. - eps;
	const double leps   = std::log(eps);
	const double leps_c = std::log(eps_c);

	double w_ij = 0.;
	for (auto g : Genotypes::homo) {
		const auto P_bbar_I_gd = P_bbar_I_gds[g];
		const auto P_g_I_d     = P_g_I_ds[g];

		_Q += P_g_I_d * (P_bbar_I_gd * leps_c + (1. - P_bbar_I_gd) * leps);
		w_ij += P_g_I_d * (eps_c - P_bbar_I_gd);
	}

	if (!isInvariant) {
		for (auto g : Genotypes::het) {
			const auto P_bbar_I_gd = P_bbar_I_gds[g];
			const auto P_g_I_d     = P_g_I_ds[g];

			_Q += P_g_I_d * (P_bbar_I_gd * leps_c + (1. - P_bbar_I_gd) * leps);
			w_ij += P_g_I_d * (eps_c - P_bbar_I_gd);
		}
	}
	for (auto dm = der1st.begin(); dm != der1st.end(); ++dm) _F(dm->index) += w_ij * dm->derivative;

	// add first derivative products to Jacobian
	const double epsEps_c = eps * eps_c;
	for (auto dm = der1st.begin(); dm != der1st.end(); ++dm) {
		_Jacobian(dm->index, dm->index) -= epsEps_c * dm->derivative * dm->derivative;
		for (auto dn = dm + 1; dn != der1st.end(); ++dn) {
			_Jacobian(dm->index, dn->index) -= epsEps_c * dm->derivative * dn->derivative;
		}
	}
	// add second derivatives to Jacobian
	for (auto &dmn : der2nd) { _Jacobian(dmn.index1, dmn.index2) += w_ij * dmn.derivative; }

	++_numSitesAdded;
}

void TEpsilon::add(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds,
				   const TGenotypeFloats &P_bbar_I_gds, bool updateJF, bool isInvariant) {
	if (_accepted) return;
	if (updateJF)
		_addToQJF(data, P_g_I_ds, P_bbar_I_gds, isInvariant);
	else
		_addToQ(data, P_g_I_ds, P_bbar_I_gds, isInvariant);
}

TEpsilon::TEpsilon(std::string_view Def) : _functions(Def) {
	const size_t numParameters = _functions.numParameters();

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters, numParameters);
	_F.resize(numParameters);
}

TEpsilon::TEpsilon(const BAM::RGInfo::TInfo &Info) : _functions(Info) {
	const size_t numParameters = _functions.numParameters();

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters, numParameters);
	_F.resize(numParameters);
}

void TEpsilon::init(const RecalEstimatorTools::TRecalDataTable &DataTable, size_t MinData) {
	_functions.init(DataTable, MinData);
	const size_t numParameters = _functions.numParameters();

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters, numParameters);
	_F.resize(numParameters);
}

coretools::Probability TEpsilon::calcErrorRate(const BAM::TSequencedData &data) const noexcept {
	return _functions.getEpsilon(data);
}

coretools::Probability TEpsilon::_calcErrorRate(const BAM::TSequencedData &data, std::vector<T1stDerivative> &der1,
											   std::vector<T2ndDerivative> &der2) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	return _functions.getEpsilon(data, der1, der2);
}


void TEpsilon::solveJxF() {
	// scale maxF #sites
	_Jacobian = arma::symmatu(_Jacobian);
	_maxF = std::max(std::abs(_F.max()), std::abs(_F.min())) / _numSitesAdded;
	if (!solve(_JxF, _Jacobian, _F))
		UERROR("Issue solving JxF! This may be due to a lack of data. Consider adding more sites. Jacobian: ",
		       _Jacobian);
	_maxJxF = std::max(std::abs(_JxF.max()), std::abs(_JxF.min()));

	// automatically reset
	_Jacobian.zeros();
	_F.zeros();
	_numSitesAdded = 0;
	_oldQ          = _Q;
}

void TEpsilon::propose(double lambda) {
	if (!_accepted) {
		_functions.propose(lambda, _JxF);
		_Q = 0; // reset to recalculate
	}
}

bool TEpsilon::acceptOrReject() {
	_accepted = _Q > _oldQ;
	if (!_accepted) _functions.reject();
	return _accepted;
}

void TEpsilon::adjust() {
	_Q         = 0.;
	_accepted = false;
	_functions.adjust();
}

void TEpsilon::log() const {
	_functions.log();
}

BAM::RGInfo::TInfo TEpsilon::info() const{
	return _functions.info();
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
