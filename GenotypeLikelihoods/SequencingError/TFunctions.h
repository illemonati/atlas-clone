/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Dec 20, 2022
 *      Author: Andreas
 */

#ifndef GENOTYPELIKELIHOODS_TFUNCTIONS_H_
#define GENOTYPELIKELIHOODS_TFUNCTIONS_H_

#include <armadillo>
#include <memory>
#include <vector>

#include "RecalEstimatorTools.h"
#include "SequencingError/TFunction.h"
#include "TDerivatives.h"
#include "TReadGroupInfo.h"
#include "TSequencedBase.h"

#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods::SequencingError {

class TFunctions {
	TIntercept _intercept;
	coretools::TStrongArray<std::unique_ptr<TFunction>, Covariates> _covariates{
		{nullptr, nullptr, nullptr, nullptr, nullptr}};
	std::vector<double> _oldBetas;

public:
	TFunctions(const BAM::RGInfo::TInfo& info);
	TFunctions(std::string_view Def); 

	void init(const RecalEstimatorTools::TRecalDataTable &DataTable);
	size_t numParameters() const noexcept;

	coretools::Probability getEpsilon(const BAM::TSequencedBase &base) const noexcept;
	coretools::Probability getEpsilon(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
									  std::vector<T2ndDerivative> &der2) const noexcept;
	void reject() noexcept;
	void propose(double lambda, const arma::mat &_JxF) noexcept;
	void adjust() noexcept;

	void log() const;
	BAM::RGInfo::TInfo info() const;
};
} // namespace GenotypeLikelihoods::SequencingError

#endif
