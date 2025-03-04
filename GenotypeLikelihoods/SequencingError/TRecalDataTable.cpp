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
void TRecalDataTable::add(const BAM::TSequencedData & data){
	using SequencingError::Covariates;
	using SequencingError::TCovariate_context;
	using SequencingError::TCovariate_fragmentLength;
	using SequencingError::TCovariate_mappingQuality;
	using SequencingError::TCovariate_position;
	using SequencingError::TCovariate_quality;
	++_counts;

	//add quality
	impl::addUsed(_tables[Covariates::Context], TCovariate_context::extract(data));
	impl::addUsed(_tables[Covariates::FragmentLength], TCovariate_fragmentLength::extract(data));
	impl::addUsed(_tables[Covariates::MappingQuality], TCovariate_mappingQuality::extract(data));
	impl::addUsed(_tables[Covariates::Position], TCovariate_position::extract(data));
	impl::addUsed(_tables[Covariates::Quality], TCovariate_quality::extract(data));
}
}
