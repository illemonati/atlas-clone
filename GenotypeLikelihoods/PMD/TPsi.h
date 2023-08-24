/*
 * PMD/TModels.h
 */

#ifndef GENOTYPELIKELIHOODS_PMD_TPSI_H_
#define GENOTYPELIKELIHOODS_PMD_TPSI_H_

#include "TReadGroupInfo.h"
#include "TSequencedBase.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"
#include <vector>

namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods::PMD {

enum class Type : size_t {min, CT=min, GA, max};
enum class End : size_t {min, from5=min, from3, max};
class TPsi {
	coretools::TStrongArray<coretools::TStrongArray<std::vector<coretools::Probability>, End>, Type> _tables;

	void _fromString(std::string_view Psi);
public:
	TPsi(std::string_view Psi); 
	TPsi(const BAM::RGInfo::TInfo & info);

	BAM::RGInfo::TInfo info() const; 

	template<Type From_To>
	coretools::Probability prob(const BAM::TSequencedBase &data) const noexcept {
		constexpr coretools::TStrongArray<Type, Type> flip{{Type::GA, Type::CT}};

		const auto realType = data.isReverseStrand() ? flip[From_To] : From_To;
		const auto end      = data.distFrom5Prime < data.distFrom3Prime ? End::from5 : End::from3;
		const auto pos      = end == End::from5 ? data.distFrom5Prime : data.distFrom3Prime;

		const auto &table = _tables[realType][end];
		if (pos >= table.size()) return table.back();
		return table[pos];
	}
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
