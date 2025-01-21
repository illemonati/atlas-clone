/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef SEQUENCINGERROR_TFUNCTION_H_
#define SEQUENCINGERROR_TFUNCTION_H_

#include <string>
#include <vector>

#include "TDerivatives.h"
#include "TReadGroupInfo.h"

namespace GenotypeLikelihoods::RecalEstimatorTools {class TRecalDataTable;}
namespace BAM {struct TSequencedData;}

namespace GenotypeLikelihoods::SequencingError {

//--------------------------------------------------------------
// TCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TFunction {
protected:
	size_t _firstParameterIndex;

public:
	TFunction(size_t FirstParameterIndex) : _firstParameterIndex(FirstParameterIndex) {}
	virtual ~TFunction() = default;

	// non-virtuals
	constexpr size_t firstParameterIndex() const noexcept { return _firstParameterIndex; }

	// virtuals
	virtual double *begin() noexcept              = 0;
	virtual double *end() noexcept                = 0;
	virtual const double *begin() const noexcept  = 0;
	virtual const double *end() const noexcept    = 0;
	virtual size_t numParameters() const noexcept = 0;

	// check value range: to ensure that data can be recalibrated
	virtual void init(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t FirstParameterIndex, size_t MinData) = 0;

	// estimation
	virtual double getEta(const BAM::TSequencedData &data) const noexcept   = 0;
	virtual double getEta(const BAM::TSequencedData &data, std::vector<T1stDerivative> &der1,
						  std::vector<T2ndDerivative> &der2) const noexcept = 0;
	virtual double adjust() noexcept                                        = 0;
	virtual std::string typeString() const noexcept                         = 0;
	virtual void addInfo(BAM::RGInfo::TInfo &info) const                    = 0;
	virtual void log() const; 
};
} // namespace GenotypeLikelihoods::SequencingError

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
