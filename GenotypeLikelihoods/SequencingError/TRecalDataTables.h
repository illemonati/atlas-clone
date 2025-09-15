/*
 * RecalEstimatorTools.h
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#ifndef TRECALDATATABLES_h_
#define TRECALDATATABLES_h_

#include <vector>

#include "TReadGroupMap.h"
#include "TSequencedData.h"
#include "coretools/Containers/TStrongArray.h"
#include "SequencingError/TRecalDataTable.h"

namespace GenotypeLikelihoods { class TSite; }

namespace GenotypeLikelihoods::RecalEstimatorTools {


using TRecalDataTableOneReadGroup = coretools::TStrongArray<TRecalDataTable, BAM::Mate>;

class TRecalDataTables{
private:
	const BAM::TReadGroupMap *_readGroupMap = nullptr;
	std::vector<TRecalDataTableOneReadGroup> _tables; //tables[readGroup][first/second mate]
	size_t _NSites    = 0;
	size_t _NSites_g1 = 0;
	size_t _NBases    = 0;

public:
	TRecalDataTables(const BAM::TReadGroupMap &ReadGroupMapObject)
	    : _readGroupMap(&ReadGroupMapObject), _tables(_readGroupMap->size()){}

	void add(const TSite & site);

	constexpr size_t size() const noexcept {return _NBases;}
	constexpr size_t NSites_g1() const noexcept {return _NSites_g1;}
	const TRecalDataTableOneReadGroup& operator[](size_t readGroupId) const;

	void write(std::string_view Name) const;

	void clear() {
		// empty all memory!
		std::vector<TRecalDataTableOneReadGroup>().swap(_tables);
	};
};

}; // namespace GenotypeLikelihoods::RecalEstimatorTools

#endif /* GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_ */
