/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Dec 20, 2022
 *      Author: Andreas
 */

#ifndef GENOTYPELIKELIHOODS_TFUNCTIONS_H_
#define GENOTYPELIKELIHOODS_TFUNCTIONS_H_

#include "TSequencedBase.h"
#include "TFunction.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/probability.h"
#include <memory>
#include <vector>

namespace GenotypeLikelihoods::SequencingError {
namespace impl {
constexpr coretools::Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return coretools::Probability(0.9999999999);
	if (eta < -23.03) return coretools::Probability(0.0000000001);

	return coretools::logistic(eta);
}
}


class TFunctions {
	std::vector<std::unique_ptr<TFunction>> _functions;
public:

	TFunctions(std::string_view Def);
	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable);

	size_t numParameters() const noexcept {
		size_t numParameters = 0;
		for (const auto& f: _functions) {
			numParameters += f->numParameters();
		}
		return numParameters;
	}

	coretools::Probability getEpsilon(const BAM::TSequencedBase &base) const {
		double eta = 0.;
		for (const auto &fn : _functions) eta += fn->getEta(base);
		return impl::calcEpsilon(eta);
	}

	coretools::Probability getEpsilon(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
								  std::vector<T2ndDerivative> &der2) const noexcept {
		double eta = 0.;
		for (const auto &fn : _functions) eta += fn->getEta(base, der1, der2);
		return impl::calcEpsilon(eta);
	}

	auto &front() noexcept {return _functions.front();}
	const auto &front() const noexcept {return _functions.front();}
	auto begin() noexcept { return _functions.begin(); }
	auto begin() const noexcept { return _functions.begin(); }
	auto end() noexcept { return _functions.end(); }
	auto end() const noexcept { return _functions.end(); }
};

inline TFunctions *makeFunctions(std::string_view Def) { return new TFunctions(Def); }
} // namespace GenotypeLikelihoods::SequencingError

#endif
