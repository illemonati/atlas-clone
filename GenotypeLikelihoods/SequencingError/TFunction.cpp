/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#include "SequencingError/TFunction.h"

#include <algorithm>
#include <iterator>

#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/toString.h"
#include "coretools/Math/TProbabilityDistributions.h"

#include "coretools/devtools.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using coretools::str::toString;

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
//--------------------------------------------------------------

void TFunction::_initializeValues(const std::vector<std::string> &betas) {
	if (!betas.empty()) {
		if (betas.size() != numParameters()) {
			throw toString("Failed to initialize recalibration module: wrong number of values (", betas.size(),
						   " instead of ", numParameters(), ")!");
		}
		std::transform(betas.cbegin(), betas.cend(), _begin(), coretools::str::convertStringCheck<double>);
	}
}

void TFunction::proposeNewParameters(const arma::mat &JxF, size_t &index, double lambda) noexcept {
	// update new ones
	std::copy(_cbegin(), _cend(), _obegin());
	for (auto it = _begin(); it != _end(); ++it) {
		*it -= lambda * JxF(index++);
	}
}

void TFunction::rejectProposedParameters() noexcept {
	std::copy(_obegin(), _oend(), _begin());
}

std::string TFunction::modelString() const {
	std::string s = typeString() + "[" + std::accumulate(_cbegin(), _cend(), std::string{}, [](auto tot, auto b){return tot + toString(b) + ",";});
	return s.substr(0, s.size() - 1) + "]";
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
