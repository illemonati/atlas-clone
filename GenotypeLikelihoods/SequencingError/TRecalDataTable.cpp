#include "TRecalDataTable.h"
namespace GenotypeLikelihoods::RecalEstimatorTools {

namespace impl {
void addUsed(std::vector<size_t> &counts, size_t value) {
	if (counts.size() <= value) counts.resize(value + 1, 0);
	++counts[value];
}

} // namespace impl

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
void TRecalDataTable::add(const BAM::TSequencedBase & base){
	++_counts;

	//add quality
	impl::addUsed(_tables[SequencingError::Covariates::Context], coretools::index(base.previousBase));
	impl::addUsed(_tables[SequencingError::Covariates::FragmentLength], base.fragmentLength);
	impl::addUsed(_tables[SequencingError::Covariates::MappingQuality], base.mappingQuality.get());
	impl::addUsed(_tables[SequencingError::Covariates::Position], base.distFrom5Prime);
	impl::addUsed(_tables[SequencingError::Covariates::Quality], base.originalQuality.get());
}
}
