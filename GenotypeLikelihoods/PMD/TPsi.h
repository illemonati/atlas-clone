/*
 * PMD/TModels.h
 */

#ifndef GENOTYPELIKELIHOODS_PMD_TPSI_H_
#define GENOTYPELIKELIHOODS_PMD_TPSI_H_

#include "TGenotypeData.h"
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
	struct NumDenom {
		double num   = 0.;
		double denom = 0.;
	};
	coretools::TStrongArray<coretools::TStrongArray<std::vector<coretools::Probability>, Type>, End> _tables;
	coretools::TStrongArray<coretools::TStrongArray<std::vector<NumDenom>, Type>, End> _tableSums;

	void _fromString(std::string_view Psi);

	template<Type From_To>
	void _add(const BAM::TSequencedBase &data, const TGenotypeLikelihoods &P_g_I_dij,
			  const TBaseProbabilities &P_bbar_I_From) {
		using genometools::Base;
		constexpr coretools::TStrongArray<Type, Type> Flip{{Type::GA, Type::CT}};
		constexpr coretools::TStrongArray<Base, Type> From{{Base::C, Base::G}};
		constexpr coretools::TStrongArray<Base, Type> To{{Base::T, Base::A}};

		const auto end  = data.distFrom5Prime < data.distFrom3Prime ? End::from5 : End::from3;
		const auto pos  = end == End::from5 ? data.distFrom5Prime : data.distFrom3Prime;
		const auto from = From[From_To];
		const auto to   = To[From_To];

		const auto realType = data.isReverseStrand() ? Flip[From_To] : From_To;
		auto &tSum          = _tableSums[end][realType];
		if (tSum.size() <= pos) tSum.resize(pos + 1);

		for (auto a = Base::min; a < Base::max; ++a) {
			const auto g = genometools::genotype(from, a);
			tSum[pos].num   += P_bbar_I_From[to] * P_g_I_dij[g];
			tSum[pos].denom += (P_bbar_I_From[to] + P_bbar_I_From[from]) * P_g_I_dij[g];
		}
	}

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

		const auto &table = _tables[end][realType];
		if (pos >= table.size()) return table.back();
		return table[pos];
	}

	void add(const BAM::TSequencedBase &data, const TGenotypeLikelihoods &P_g_I_dij,
			 const TBaseProbabilities &P_bbar_I_C, const TBaseProbabilities &P_bbar_I_G) noexcept {
		_add<Type::CT>(data, P_g_I_dij, P_bbar_I_C);
		_add<Type::GA>(data, P_g_I_dij, P_bbar_I_G);
	}

	void estimate() noexcept;

	std::string definition() const noexcept;
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
