/*
 * TEpsilon.h
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#ifndef TEPSILON_H_
#define TEPSILON_H_

#include <armadillo>
#include <memory>

#include "TReadGroupInfo.h"
#include "genometools/GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TGenotypeData.h"
#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

class TFunction;

struct T1stDerivative {
	size_t index;
	double derivative;
	T1stDerivative(size_t Index, double Derivative) : index(Index), derivative(Derivative) {}
};

struct T2ndDerivative {
	size_t index1;
	size_t index2;
	double derivative;
	T2ndDerivative(size_t Index1, size_t Index2, double Derivative) : index1(Index1), index2(Index2), derivative(Derivative) {}
};


class TEpsilon {
	std::vector<std::unique_ptr<TFunction>> _functions;
	std::vector<double> _oldBetas;

	double _Q    = 0.;
	double _oldQ = std::numeric_limits<double>::lowest();
	double _maxF = 0.;
	arma::mat _Jacobian;
	arma::vec _F;
	arma::mat _JxF;
	size_t _numSitesAdded = 0;

	coretools::Probability _calcErrorRate(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
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
	void addToQ(const BAM::TSequencedBase &base, const TGenotypeLikelihoods &P_g_I_ds,
				const TGenotypeLikelihoods &P_bbar_I_gds) {
		using genometools::Genotype;
		const double eps    = calcErrorRate(base);
		const double eps_c  = 1. - eps;
		const double leps   = std::log(eps);
		const double leps_c = std::log(eps_c);

		constexpr std::array<Genotype, 4 + (!isInvariant) * 6> gs = _makeGenotype<isInvariant>();
		for (auto g : gs) {
			const double P_bbar_I_gd = P_bbar_I_gds[g];
			const double P_g_I_d     = P_g_I_ds[g];
			_Q += P_g_I_d * (P_bbar_I_gd * leps_c + (1. - P_bbar_I_gd) * leps);
		}
	}

	template<bool isInvariant>
	void _addToQJF(const BAM::TSequencedBase &base, const TGenotypeLikelihoods &P_g_I_ds,
				   const TGenotypeLikelihoods &P_bbar_I_gds) {
		static std::vector<T1stDerivative> der1st;
		static std::vector<T2ndDerivative> der2nd;
		der1st.clear();
		der2nd.clear();
		// get error rate
		const double eps    = _calcErrorRate(base, der1st, der2nd);
		const double eps_c  = 1. - eps;
		const double leps   = std::log(eps);
		const double leps_c = std::log(eps_c);

		for (auto g : _makeGenotype<isInvariant>()) {
			const double P_bbar_I_gd = P_bbar_I_gds[g];
			const double P_g_I_d     = P_g_I_ds[g];

			// add Q
			_Q += P_g_I_d * (P_bbar_I_gd * leps_c + (1 - P_bbar_I_gd) * leps);

			const double w_ij = P_g_I_d * (eps_c - P_bbar_I_gd);

			// add first derivatives
			for (auto dm = der1st.begin(); dm != der1st.end(); ++dm) {
				// add to F
				_F(dm->index) += w_ij * dm->derivative;

				// add to upper half of J
				_Jacobian(dm->index, dm->index) -= eps_c * eps * dm->derivative * dm->derivative;
				for (auto dn = dm + 1; dn != der1st.end(); ++dn) {
					_Jacobian(dm->index, dn->index) -= eps_c * eps * dm->derivative * dn->derivative;
				}
			}

			// add second derivatives to Jacobian (happens to have the same weigth as F!)
			for (auto &dmn : der2nd) _Jacobian(dmn.index1, dmn.index2) += w_ij * dmn.derivative;

			++_numSitesAdded;
		}
	}

public:
	TEpsilon(std::string_view Def);
	~TEpsilon();

	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable);

	coretools::Probability calcErrorRate(const BAM::TSequencedBase &base) const noexcept; 
	double Q() const noexcept {return _Q;};
	double maxF() const noexcept {return _maxF;};

	template<bool updateJF, bool isInvariant>
	void addToEpsilon(const BAM::TSequencedBase &base, const TGenotypeLikelihoods &P_g_I_ds, const TGenotypeLikelihoods & P_bbar_I_gds) {
		if constexpr (updateJF) _addToQJF<isInvariant>(base, P_g_I_ds, P_bbar_I_gds);
		else addToQ<isInvariant>(base, P_g_I_ds, P_bbar_I_gds);
	}
	void solveJxF();
	void propose(double lambda);
	bool acceptOrReject();
	void adjust();

	std::string definition() const noexcept;
	BAM::RGInfo::TInfo info() const;
};
} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
