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
void TRecalDataTable::add(const BAM::TSequencedData & data, size_t ID) {
	using namespace SequencingError;

	if (_lastID != ID) {
		if (_IDcount > 1) ++_counts_g1;
		_IDcount = 0;
		_lastID = ID;
	}
	++_IDcount;
	++_counts;

	//add quality
	impl::addUsed(_tables[Covariates::Context], TCovariate_context::extract(data));
	impl::addUsed(_tables[Covariates::FragmentLength], TCovariate_fragmentLength::extract(data));
	impl::addUsed(_tables[Covariates::MappingQuality], TCovariate_mappingQuality::extract(data));
	impl::addUsed(_tables[Covariates::Position], TCovariate_position::extract(data));
	impl::addUsed(_tables[Covariates::Quality], TCovariate_quality::extract(data));
}
}
