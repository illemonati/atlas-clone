#ifndef SEQUENCINGERROR_TNOFUNCTION_H_
#define SEQUENCINGERROR_TNOFUNCTION_H_

#include "TFunction.h"

namespace GenotypeLikelihoods::SequencingError {

class TNoFunction final : public TFunction{
public:
	TNoFunction() : TFunction(0) {}

	// virtuals
	double *begin() noexcept override {return nullptr;}
	double *end() noexcept override {return nullptr;}
	const double *begin() const noexcept override {return nullptr;}
	const double *end() const noexcept override {return nullptr;}
	virtual size_t numParameters() const noexcept override {return 0;}

	// check value range: to ensure that data can be recalibrated
	void init(const RecalEstimatorTools::TRecalDataTable &, size_t, size_t) override {}

	// estimation
	double getEta(const BAM::TSequencedData &) const noexcept override {return 0.;}
	double getEta(const BAM::TSequencedData &, std::vector<T1stDerivative> &,
						  std::vector<T2ndDerivative> &) const noexcept override {return 0.;}
	double adjust() noexcept override {return 0.;}
	std::string typeString() const noexcept override {return "";}
	void log() const noexcept  override {}; 
	void addInfo(BAM::RGInfo::TInfo &) const override {};
};
}

#endif
