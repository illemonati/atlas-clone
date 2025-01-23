/*
 * TEpsilon.h
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#ifndef TEPSILON_H_
#define TEPSILON_H_

#include <armadillo>

#include "SequencingError/TFunctions.h"
#include "genometools/Genotypes/Containers.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"

namespace RecalEstimatorTools {class TRecalDataTable;}

namespace GenotypeLikelihoods {
namespace SequencingError {

class TEpsilon {
	TFunctions _functions;

	double _Q      = 0.;
	double _oldQ   = 0.;
	double _maxF   = 0.;
	double _maxJxF = 0.;

	arma::mat _Jacobian;
	arma::vec _F;
	arma::mat _JxF;
	size_t _numSitesAdded = 0;
	bool _converged = false;

	coretools::Probability _calcErrorRate(const BAM::TSequencedData &data, std::vector<T1stDerivative> &der1,
												   std::vector<T2ndDerivative> &der2) const noexcept;

	template<bool isInvariant> static constexpr auto _makeGenotype() {
		using genometools::Genotype;
		if constexpr (isInvariant) {
			using genometools::Base;
			std::array<Genotype, 4> ar{};
			for (size_t i = 0; i < ar.size(); ++i) {
				const auto a = Base(i);
				const auto g = genometools::genotype(a, a);
				ar[i] = g;
			}
			return ar;
		} else {
			std::array<Genotype, 10> ar{};
			for (size_t i = 0; i < ar.size(); ++i) {
				const auto g = Genotype(i);
				ar[i] = g;
			}
			return ar;
		}
	}

	template<bool isInvariant>
	void _addToQ(const BAM::TSequencedData &data, const genometools::TGenotypeLikelihoods &P_g_I_ds,
				const genometools::TGenotypeLikelihoods &P_bbar_I_gds) {
		const double eps    = calcErrorRate(data);
		const double eps_c  = 1. - eps;
		const double leps   = std::log(eps);
		const double leps_c = std::log(eps_c);

		for (auto g : _makeGenotype<isInvariant>()) {
			const double P_bbar_I_gd = P_bbar_I_gds[g];
			const double P_g_I_d     = P_g_I_ds[g];
			_Q += P_g_I_d * (P_bbar_I_gd * leps_c + (1. - P_bbar_I_gd) * leps);
		}
	}

	template<bool isInvariant>
	void _addToQJF(const BAM::TSequencedData &data, const genometools::TGenotypeLikelihoods &P_g_I_ds,
				   const genometools::TGenotypeLikelihoods &P_bbar_I_gds) {
		static std::vector<T1stDerivative> der1st;
		static std::vector<T2ndDerivative> der2nd;
		der1st.clear();
		der2nd.clear();
		// get error rate
		const double eps      = _calcErrorRate(data, der1st, der2nd);
		const double eps_c    = 1. - eps;
		const double leps     = std::log(eps);
		const double leps_c   = std::log(eps_c);

		double w_ij = 0.;
		for (auto g : _makeGenotype<isInvariant>()) {
			const double P_bbar_I_gd = P_bbar_I_gds[g];
			const double P_g_I_d     = P_g_I_ds[g];

			_Q   += P_g_I_d * (P_bbar_I_gd * leps_c + (1. - P_bbar_I_gd) * leps);
			w_ij += P_g_I_d * (eps_c - P_bbar_I_gd);

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
		for (auto &dmn : der2nd) {
			_Jacobian(dmn.index1, dmn.index2) += w_ij * dmn.derivative;
		}

		++_numSitesAdded;
	}

public:
	TEpsilon(std::string_view Def);
	TEpsilon(const BAM::RGInfo::TInfo & info);

	void init(const RecalEstimatorTools::TRecalDataTable &DataTable, size_t MinData);

	coretools::Probability calcErrorRate(const BAM::TSequencedData &data) const noexcept; 
	double oldQ() const noexcept {return _oldQ;};
	double Q() const noexcept {return _Q;};
	double maxF() const noexcept {return _maxF;};
	double maxChange() const noexcept {return _maxJxF;}

	template<bool updateJF, bool isInvariant>
	void add(const BAM::TSequencedData &data, const genometools::TGenotypeLikelihoods &P_g_I_ds, const genometools::TGenotypeLikelihoods & P_bbar_I_gds) {
		if (_converged) return;
		if constexpr (updateJF) _addToQJF<isInvariant>(data, P_g_I_ds, P_bbar_I_gds);
		else _addToQ<isInvariant>(data, P_g_I_ds, P_bbar_I_gds);
	}
	void solveJxF();
	void propose(double lambda);
	bool acceptOrReject();
	void adjust();

	void log() const;
	BAM::RGInfo::TInfo info() const;
};
} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
