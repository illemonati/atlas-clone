#include "TRecalDataTable.h"
#include "TCovariate.h"
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
	using SequencingError::Covariates;
	++_counts;

	//add quality
	impl::addUsed(_tables[Covariates::Context], coretools::index(base.previousBase));
	impl::addUsed(_tables[Covariates::MappingQuality], base.mappingQuality.get());
	impl::addUsed(_tables[Covariates::Quality], base.originalQuality.get());

	// remain reasonable
	impl::addUsed(_tables[Covariates::FragmentLength], std::min<size_t>(SequencingError::TCovariate_fragmentLength::NMax - 1,base.fragmentLength));
	impl::addUsed(_tables[Covariates::Position], std::min<size_t>(SequencingError::TCovariate_position::NMax - 1, base.distFrom5Prime));
}
}
