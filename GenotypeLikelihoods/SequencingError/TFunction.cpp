/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#include "SequencingError/TFunction.h"

#include <algorithm>
#include <iterator>

#include "commonWeakTypes.h"
#include "mathFunctions.h"
#include "probability.h"
#include "stringFunctions.h"
#include "toString.h"
#include "TProbabilityDistributions.h"

#include "devtools.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using coretools::str::toString;

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
//--------------------------------------------------------------
void TFunction::_initializeValues(const BAM::RGInfo::Info &betas) {
	if (!betas.empty()) {
		if (betas.size() != numParameters()) {
			UERROR("Failed to initialize recal function '", typeString(), "': wrong number of values (", betas.size(),
						   " instead of ", numParameters(), ")!");
		}
		if(betas.is_array()){
			for(size_t i = 0; i < numParameters(); ++i){
				_betas[i] = betas[i];
			}
		} else if(numParameters() == 1 && betas.is_number()){
			_betas[0] = betas.get();
		} else {
			UERROR("Failed to initialize recal function '", typeString(), "': expected an array of values!");
		}
	}
}

void TFunction::_initializeValues(const std::vector<std::string> &betas) {
	if (!betas.empty()) {
		if (betas.size() != numParameters()) {
			throw toString("Failed to initialize recal function '", typeString(), "': wrong number of values (", betas.size(),
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

virtual BAM::RGInfo::TInfo TFunction::modelInfo() const {
	return BAM::RGInfo::TInfo{ {typeString(), _betas}  };
}


} // namespace SequencingError
} // namespace GenotypeLikelihoods
