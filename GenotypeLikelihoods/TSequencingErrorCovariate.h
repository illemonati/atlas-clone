/*
 * TRecalibrationEMModules.h
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_

#include "../BAM/TSequencedBase.h"
#include "GenotypeTypes.h"
#include "RecalEstimatorTools.h"
#include "TSequencingErrorCovariateFunction.h"
#include <memory>

namespace GenotypeLikelihoods {
namespace SequencingError {

// define covariate names

//------------------------------------------------------------------------------------
// TCovariate
// This is the base class without any covariate. Not intended to be used in models!
//------------------------------------------------------------------------------------
class TCovariate {
protected:
	std::unique_ptr<TCovariateFunction> _function;


	// extract
	virtual uint16_t _extractCovariate(const BAM::TSequencedBase &) const noexcept = 0;
	TCovariate()                                                    = default;
	virtual void addFunction(const size_t, const std::string &, const RecalEstimatorTools::TRecalDataTable &) = 0;
	virtual void addFunction(const size_t, const std::string &)                                               = 0;

public:
	static inline const std::string name = "none";
	virtual ~TCovariate()                                           = default;

	uint16_t numParameters() const noexcept { return _function ? _function->numParameters() : 0; }
	uint16_t numNonZeroFirstDerivatives() const noexcept {
		return _function ? _function->numNonZeroFirstDerivatives() : 0;
	}
	uint16_t numNonZeroSecondDerivatives() const noexcept {
		return _function ? _function->numNonZeroSecondDerivatives() : 0;
	}
	virtual std::string typeString() const = 0;

	// covariate function
	std::string functionString() const {
		return _function ? _function->getModelString() : TCovariateFunction::name;
	}

	virtual bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable &)  const noexcept = 0;
	virtual bool checkParameterRange(std::vector<uint16_t> &, uint16_t) const noexcept             = 0;
	virtual void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &) = 0;

	TCovariateFunction *getPointerToFunction() const noexcept { return _function.get(); }

	// calculate terms
	double getEtaTerm(const BAM::TSequencedBase &base) const noexcept { return _function->getEtaTerm(_extractCovariate(base)); }

	void fillDerivatives(const BAM::TSequencedBase &base, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) {
		_function->fillDerivatives(_extractCovariate(base), first, second);
	}
	double adjustParametersPostEstimation() { return _function->adjustParametersPostEstimation(); }
};

//-------------------------------------------
// TCovariate_quality
//-------------------------------------------
class TCovariate_quality : public TCovariate {
private:
	TRecalibrationEMQualityTransformationMap _qualityToLogit;

	uint16_t _extractCovariate(const BAM::TSequencedBase &base) const noexcept override {
		return base.originalQuality_phredInt.get();
	}

public:
	static inline const std::string name = "quality";
	TCovariate_quality(size_t FirstParameterIndex, const std::string &functionString,
			   const RecalEstimatorTools::TRecalDataTable &dataTable);
	TCovariate_quality(size_t FirstParameterIndex, const std::string &functionString);

	std::string typeString() const override { return name; }
	void addFunction(size_t FirstParameterIndex, const std::string &functionString,
			 const RecalEstimatorTools::TRecalDataTable &dataTable) override;
	void addFunction(size_t FirstParameterIndex, const std::string &functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override;
	bool checkParameterRange(std::vector<uint16_t> &usedValues, uint16_t maxPos) const noexcept override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) override;
};

//-------------------------------------------
// TCovariate_position
//-------------------------------------------
class TCovariate_position : public TCovariate {
private:
	uint16_t _extractCovariate(const BAM::TSequencedBase &base) const noexcept override { return base.distFrom5Prime; }

public:
	static inline const std::string name = "position";
	TCovariate_position(size_t FirstParameterIndex, const std::string &functionString,
			    const RecalEstimatorTools::TRecalDataTable &dataTable);
	TCovariate_position(size_t FirstParameterIndex, const std::string &functionString);

	std::string typeString() const override { return name; }
	void addFunction(size_t FirstParameterIndex, const std::string &functionString,
			 const RecalEstimatorTools::TRecalDataTable &dataTable) override;
	void addFunction(size_t FirstParameterIndex, const std::string &functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override;
	bool checkParameterRange(std::vector<uint16_t> &usedValues, uint16_t maxPos) const noexcept override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) override;
};

//-------------------------------------------
// TCovariate_context
//-------------------------------------------
class TCovariate_context : public TCovariate {
private:
	static constexpr int numContext = 20;

	uint16_t _extractCovariate(const BAM::TSequencedBase &base) const noexcept override { return genometools::index(base.context); }

public:
	static inline const std::string name = "context";
	TCovariate_context(size_t FirstParameterIndex, const std::string &functionString,
			   const RecalEstimatorTools::TRecalDataTable &dataTable);
	TCovariate_context(size_t FirstParameterIndex, const std::string &functionString);

	std::string typeString() const override { return name; }
	void addFunction(size_t FirstParameterIndex, const std::string &functionString,
			 const RecalEstimatorTools::TRecalDataTable &dataTable) override;
	void addFunction(size_t FirstParameterIndex, const std::string &functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override;
	bool checkParameterRange(std::vector<uint16_t> &usedValues, uint16_t maxPos) const noexcept override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &) override{};
};

//-------------------------------------------
// TCovariate_fragmentLength
//-------------------------------------------

class TCovariate_fragmentLength : public TCovariate {
private:
	uint16_t _extractCovariate(const BAM::TSequencedBase &base) const noexcept override { return base.fragmentLength; }

public:
	static inline std::string name = "fragmentLength";
	TCovariate_fragmentLength(size_t FirstParameterIndex, const std::string &functionString,
				  const RecalEstimatorTools::TRecalDataTable &dataTable);
	TCovariate_fragmentLength(size_t FirstParameterIndex, const std::string &functionString);

	std::string typeString() const override { return name; }
	void addFunction(size_t FirstParameterIndex, const std::string &functionString,
			 const RecalEstimatorTools::TRecalDataTable &dataTable) override;
	void addFunction(size_t FirstParameterIndex, const std::string &functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override;
	bool checkParameterRange(std::vector<uint16_t> &usedValues, uint16_t maxPos) const noexcept override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) override;
};

//-------------------------------------------
// TCovariate_mappingQuality
//-------------------------------------------

class TCovariate_mappingQuality : public TCovariate {
private:
	uint16_t _extractCovariate(const BAM::TSequencedBase &base) const noexcept override { return base.mappingQuality; }

public:
	static inline const std::string name = "mappingQuality";
	TCovariate_mappingQuality(size_t FirstParameterIndex, const std::string &functionString,
				  const RecalEstimatorTools::TRecalDataTable &dataTable);
	TCovariate_mappingQuality(size_t FirstParameterIndex, const std::string &functionString);

	std::string typeString() const override { return name; }
	void addFunction(size_t FirstParameterIndex, const std::string &functionString,
			 const RecalEstimatorTools::TRecalDataTable &dataTable) override;
	void addFunction(size_t FirstParameterIndex, const std::string &functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override;
	bool checkParameterRange(std::vector<uint16_t> &usedValues, uint16_t maxPos) const noexcept override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) override;
};

} // namespace SequencingError
} // end namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
