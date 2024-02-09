
#ifndef GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_
#define GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_

#include "SequencingError/TCovariate.h"
#include "coretools/Containers/TStrongArray.h"
#include <cstdint>
#include <vector>

namespace GenotypeLikelihoods::RecalEstimatorTools {

class TRecalDataTable {
private:
	uint64_t _counts = 0;
	//all vectors are uint16_t, which is used by seq error models for all covariates
// Object to store for which qualities and positions data is available.
	coretools::TStrongArray<std::vector<size_t>, SequencingError::Covariates> _tables;

public:
	void add(const BAM::TSequencedBase & base);

	constexpr size_t size() const noexcept { return _counts; }
	const std::vector<size_t>& operator[](SequencingError::Covariates cov) const {return _tables[cov];}
};
}

#endif
