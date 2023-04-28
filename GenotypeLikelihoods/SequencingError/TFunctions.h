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
#include "TDerivatives.h"
#include "TReadGroupInfo.h"
#include "TSequencedBase.h"

#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods::SequencingError {

struct TFunctions {
	virtual ~TFunctions() = default;
	virtual void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable)             = 0;
	virtual size_t numParameters() const noexcept                                               = 0;
	virtual coretools::Probability getEpsilon(const BAM::TSequencedBase &base) const            = 0;
	virtual coretools::Probability getEpsilon(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
											  std::vector<T2ndDerivative> &der2) const noexcept = 0;
	virtual void reject() noexcept                                                              = 0;
	virtual void propose(double lambda, const arma::mat &_JxF) noexcept                         = 0;
	virtual void adjust() noexcept                                                              = 0;
	virtual std::string definition() const noexcept                                             = 0;
	virtual BAM::RGInfo::TInfo info() const                                                     = 0;
};
TFunctions *makeFunctions(std::string_view Def);
TFunctions *makeFunctions(const BAM::RGInfo::TInfo& info);
} // namespace GenotypeLikelihoods::SequencingError

#endif
