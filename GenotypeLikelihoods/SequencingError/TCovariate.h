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
	static constexpr std::string_view name = "quality";
	static constexpr bool isIndexed        = true;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return base.originalQuality_phredInt.get();
	}

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return RecalEstimatorTools::indexedRange(dataTable.qualities());
	}
};

//-------------------------------------------
// TCovariate_position
//-------------------------------------------
class TCovariate_position {
public:
	static constexpr std::string_view name = "position";
	static constexpr bool isIndexed        = false;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.distFrom5Prime; }

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return RecalEstimatorTools::fullRange(dataTable.positions());
	}
};

//-------------------------------------------
// TCovariate_context
//-------------------------------------------
class TCovariate_context {
private:
	static constexpr int numContext = 20;
public:
	static constexpr std::string_view name = "context";
	static constexpr bool isIndexed        = false;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return genometools::index(base.context);
	}

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		// There is a context for each position
		return RecalEstimatorTools::fullRange(dataTable.positions());
	}
};

//-------------------------------------------
// TCovariate_fragmentLength
//-------------------------------------------

class TCovariate_fragmentLength {
public:
	static constexpr std::string_view name = "fragmentLength";
	static constexpr bool isIndexed        = true;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.fragmentLength; }

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return RecalEstimatorTools::indexedRange(dataTable.fragmentLengths());
	}
};

//-------------------------------------------
// TCovariate_mappingQuality
//-------------------------------------------

class TCovariate_mappingQuality {
public:
	static constexpr std::string_view name = "mappingQuality";
	static constexpr bool isIndexed        = true;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.mappingQuality; }

	static std::vector<uint16_t> range(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return RecalEstimatorTools::indexedRange(dataTable.mappingQualities());
	}
};

} // namespace SequencingError
} // end namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
