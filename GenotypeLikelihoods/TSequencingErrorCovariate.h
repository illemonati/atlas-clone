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

// define covariate policies

//-------------------------------------------
// TCovariate_quality
//-------------------------------------------
struct TCovariate_quality {
	static inline const std::string name = "quality";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return base.originalQuality_phredInt.get();
	}

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return RecalEstimatorTools::vectorOfUsed(dataTable.qualities());
	}
};

//-------------------------------------------
// TCovariate_position
//-------------------------------------------
class TCovariate_position {
public:
	static inline const std::string name = "position";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.distFrom5Prime; }

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		const auto N = dataTable.positions().size() - 1;
		std::vector<uint16_t> v(N);
		std::iota(v.begin(), v.end(), uint16_t{});
		return v;
	}
};

//-------------------------------------------
// TCovariate_context
//-------------------------------------------
class TCovariate_context {
private:
	static constexpr int numContext = 20;
public:
	static inline const std::string name = "context";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return genometools::index(base.context);
	}

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		const auto N = dataTable.positions().size() - 1;
		std::vector<uint16_t> v(N);
		std::iota(v.begin(), v.end(), uint16_t{});
		return v;
	}
};

//-------------------------------------------
// TCovariate_fragmentLength
//-------------------------------------------

class TCovariate_fragmentLength {
public:
	static inline std::string name = "fragmentLength";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.fragmentLength; }

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return RecalEstimatorTools::vectorOfUsed(dataTable.fragmentLengths());
	}
};

//-------------------------------------------
// TCovariate_mappingQuality
//-------------------------------------------

class TCovariate_mappingQuality {
public:
	static inline const std::string name = "mappingQuality";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.mappingQuality; }

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return RecalEstimatorTools::vectorOfUsed(dataTable.mappingQualities());
	}
};

} // namespace SequencingError
} // end namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
