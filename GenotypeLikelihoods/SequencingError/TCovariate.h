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

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return base.originalQuality_phredInt.get();
	}

	static size_t N(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return dataTable.qualities().size();
	}

	static bool isUsed(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t i) noexcept {
		return dataTable.qualities()[i] > 0;
	}
};

//-------------------------------------------
// TCovariate_position
//-------------------------------------------
class TCovariate_position {
public:
	static constexpr std::string_view name = "position";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.distFrom5Prime; }

	static size_t N(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return dataTable.positions().size();
	}

	static bool isUsed(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t i) noexcept {
		return dataTable.positions()[i] > 0;
	}
};

//-------------------------------------------
// TCovariate_context
//-------------------------------------------
class TCovariate_context {
public:
	static constexpr std::string_view name = "context";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return coretools::index(base.previousBase);
	}

	static size_t N(const RecalEstimatorTools::TRecalDataTable &) noexcept {
		return coretools::index(genometools::Base::N) + 1;
	}
	static bool isUsed(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t i) noexcept {
		return i < N(dataTable);
	}
};

//-------------------------------------------
// TCovariate_fragmentLength
//-------------------------------------------

class TCovariate_fragmentLength {
public:
	static constexpr std::string_view name = "fragmentLength";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.fragmentLength; }

	static size_t N(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return dataTable.fragmentLengths().size();
	}

	static bool isUsed(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t i) noexcept {
		return dataTable.fragmentLengths()[i] > 0;
	}
};

//-------------------------------------------
// TCovariate_mappingQuality
//-------------------------------------------

class TCovariate_mappingQuality {
public:
	static constexpr std::string_view name = "mappingQuality";

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.mappingQuality.get(); }

	static size_t N(const RecalEstimatorTools::TRecalDataTable &dataTable) noexcept {
		return dataTable.mappingQualities().size();
	}

	static bool isUsed(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t i) noexcept {
		return dataTable.mappingQualities()[i] > 0;
	}
};

} // namespace SequencingError
} // end namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
