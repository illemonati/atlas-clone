/*
 * TRecalibrationEMModules.h
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_

#include "RecalEstimatorTools.h"
#include "TSequencedBase.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

// define covariate names

//------------------------------------------------------------------------------------
// TCovariate
// This is the base class without any covariate.
//------------------------------------------------------------------------------------
class TCovariate {
public:
	static inline const std::string name                                                                      = "none";
	virtual ~TCovariate()                                                                                     = default;
	virtual std::string typeString() const                                                                    = 0;
	virtual uint16_t extract(const BAM::TSequencedBase &base) const noexcept                                  = 0;
	virtual std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept = 0;
};

//-------------------------------------------
// TCovariate_quality
//-------------------------------------------
class TCovariate_quality : public TCovariate {
private:
	//TRecalibrationEMQualityTransformationMap _qualityToLogit;
public:
	static inline const std::string name = "quality";

	uint16_t extract(const BAM::TSequencedBase &base) const noexcept override {
		return base.originalQuality_phredInt.get();
	}

	std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override {
		return RecalEstimatorTools::vectorOfUsed(dataTable.qualities());
	}

	std::string typeString() const override { return name; }
};

//-------------------------------------------
// TCovariate_position
//-------------------------------------------
class TCovariate_position : public TCovariate {
public:
	static inline const std::string name = "position";

	uint16_t extract(const BAM::TSequencedBase &base) const noexcept override { return base.distFrom5Prime; }

	std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override {
		const auto N = dataTable.positions().size() - 1;
		std::vector<uint16_t> v(N);
		std::iota(v.begin(), v.end(), uint16_t{});
		return v;
	}

	std::string typeString() const override { return name; }
};

//-------------------------------------------
// TCovariate_context
//-------------------------------------------
class TCovariate_context : public TCovariate {
private:
	static constexpr int numContext = 20;
public:
	static inline const std::string name = "context";

	uint16_t extract(const BAM::TSequencedBase &base) const noexcept override {
		return genometools::index(base.context);
	}

	std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override {
		const auto N = dataTable.positions().size() - 1;
		std::vector<uint16_t> v(N);
		std::iota(v.begin(), v.end(), uint16_t{});
		return v;
	}

	std::string typeString() const override { return name; }
};

//-------------------------------------------
// TCovariate_fragmentLength
//-------------------------------------------

class TCovariate_fragmentLength : public TCovariate {
public:
	static inline std::string name = "fragmentLength";

	uint16_t extract(const BAM::TSequencedBase &base) const noexcept override { return base.fragmentLength; }

	std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override {
		return RecalEstimatorTools::vectorOfUsed(dataTable.fragmentLengths());
	}

	std::string typeString() const override { return name; }
};

//-------------------------------------------
// TCovariate_mappingQuality
//-------------------------------------------

class TCovariate_mappingQuality : public TCovariate {
public:
	static inline const std::string name = "mappingQuality";

	uint16_t extract(const BAM::TSequencedBase &base) const noexcept override { return base.mappingQuality; }

	std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept override {
		return RecalEstimatorTools::vectorOfUsed(dataTable.mappingQualities());
	}

	std::string typeString() const override { return name; }
};

} // namespace SequencingError
} // end namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
