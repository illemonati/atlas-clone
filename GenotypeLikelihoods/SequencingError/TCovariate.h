/*
 * TRecalibrationEMModules.h
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_

#include "TSequencedBase.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

enum class Covariates : size_t { min = 0, Context = min, FragmentLength, MappingQuality, Position, Quality, max };

struct TCovariate_context {
	static constexpr std::string_view name = "context";
	static constexpr Covariates index      = Covariates::Context;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return coretools::index(base.previousBase);
	}

	static size_t N(const std::vector<size_t> &) noexcept {
		return coretools::index(genometools::Base::N) + 1;
	}
	static bool isUsed(const std::vector<size_t> &, size_t) noexcept {
		return true;
	}
};

struct TCovariate_fragmentLength {
	static constexpr std::string_view name = "fragmentLength";
	static constexpr Covariates index      = Covariates::FragmentLength;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.fragmentLength; }

	static size_t N(const std::vector<size_t> &fragmentLengths) noexcept {
		return fragmentLengths.size();
	}

	static bool isUsed(const std::vector<size_t> &fragmentLengths, size_t i) noexcept {
		return fragmentLengths[i];
	}
};

struct TCovariate_mappingQuality {
	static constexpr std::string_view name = "mappingQuality";
	static constexpr Covariates index      = Covariates::MappingQuality;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.mappingQuality.get(); }

	static size_t N(const std::vector<size_t>& mappingQualities) noexcept {
		return mappingQualities.size();
	}

	static bool isUsed(const std::vector<size_t>& mappingQualities, size_t i) noexcept {
		return mappingQualities[i];
	}
};

struct TCovariate_position {
	static constexpr std::string_view name = "position";
	static constexpr Covariates index      = Covariates::Position;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept { return base.distFrom5Prime; }

	static size_t N(const std::vector<size_t>& positions) noexcept {
		return positions.size();
	}

	static bool isUsed(const std::vector<size_t>& positions, size_t i) noexcept {
		return positions[i];
	}
};

struct TCovariate_quality {
	static constexpr std::string_view name = "quality";
	static constexpr Covariates index      = Covariates::Quality;

	static uint16_t extract(const BAM::TSequencedBase &base) noexcept {
		return base.originalQuality_phredInt.get();
	}

	static size_t N(const std::vector<size_t>& qualities) noexcept {
		return qualities.size();
	}

	static bool isUsed(const std::vector<size_t>& qualities, size_t i) noexcept {
		return qualities[i];
	}
};

template<Covariates C>
struct CovariateType {
};
template<> struct CovariateType<Covariates::Context> {
	using type = TCovariate_context;
};
template<> struct CovariateType<Covariates::FragmentLength> {
	using type = TCovariate_fragmentLength;
};
template<> struct CovariateType<Covariates::MappingQuality> {
	using type = TCovariate_mappingQuality;
};
template<> struct CovariateType<Covariates::Position> {
	using type = TCovariate_position;
};
template<> struct CovariateType<Covariates::Quality> {
	using type = TCovariate_quality;
};

template<Covariates C>
using CovariateType_t = typename CovariateType<C>::type;

} // namespace SequencingError
} // end namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
