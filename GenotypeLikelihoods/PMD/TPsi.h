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
#include "genometools/GenotypeTypes.h"
#include <vector>

namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods::PMD {

enum class Type : size_t {min, CT=min, GA, max};
enum class End : size_t {min, from5=min, from3, max};
using TBaseBaseProbabilities =
	coretools::TStrongArray<coretools::TStrongArray<coretools::Probability, genometools::Base>, genometools::Base>;

class TPsi {
	static constexpr coretools::TStrongArray<Type, Type> _flip{{Type::GA, Type::CT}};
	static constexpr coretools::TStrongArray<genometools::Base, Type> _from{{genometools::Base::C, genometools::Base::G}};
	static constexpr coretools::TStrongArray<genometools::Base, Type> _to{{genometools::Base::T, genometools::Base::A}};

	struct SumType {
		struct NumDenom {
			double num;
			double denom;
		};
		struct FromTo {
			int32_t fromTo;
			int32_t fromSum;
			int32_t toFrom;
			int32_t toSum;
		};
		union {
			FromTo fromTo;
			NumDenom numDenom;
		};
	};

	coretools::TStrongArray<coretools::TStrongArray<std::vector<coretools::Probability>, Type>, End> _tables;
	coretools::TStrongArray<coretools::TStrongArray<std::vector<SumType>, Type>, End> _tableSums;

	void _fromString(std::string_view Psi);

	template<Type From_To>
	void _add(const BAM::TSequencedBase &data, coretools::Probability &P_g_I_di,
			  const TBaseBaseProbabilities &P_b_bbar_I_gdij) {
		using genometools::Base;
		constexpr auto From = _from[From_To];
		constexpr auto To   = _to[From_To];

		const auto end  = data.distFrom5Prime < data.distFrom3Prime ? End::from5 : End::from3;
		const auto pos  = end == End::from5 ? data.distFrom5Prime : data.distFrom3Prime;

		const auto realType = data.isReverseStrand() ? _flip[From_To] : From_To;
		auto &tSum          = _tableSums[end][realType];
		// tSum has already the correct size

		if (pos >= tSum.size()) return; // Not enough data or wrong pattern

		tSum[pos].numDenom.num += P_b_bbar_I_gdij[From][To] * P_g_I_di;
		tSum[pos].numDenom.denom += (P_b_bbar_I_gdij[From][To] + P_b_bbar_I_gdij[From][From]) * P_g_I_di;
	}

	template<Type From_To>
	void _add(const BAM::TSequencedBase &data, genometools::Base ref) {
		using genometools::Base;
		constexpr auto From = _from[From_To];
		constexpr auto To   = _to[From_To];

		const auto end  = data.distFrom5Prime < data.distFrom3Prime ? End::from5 : End::from3;
		const auto pos  = end == End::from5 ? data.distFrom5Prime : data.distFrom3Prime;
		const auto base = data.base;

		const auto realType = data.isReverseStrand() ? _flip[From_To] : From_To;
		auto &tSum          = _tableSums[end][realType];
		if (tSum.size() <= pos) tSum.resize(pos + 1, {{{0,0,0,0}}});

		if (ref == From) {
			if (base == To) tSum[pos].fromTo.fromTo += 1.;
			tSum[pos].fromTo.fromSum += 1.;
		} else if (ref == To) {
			if (base == From) tSum[pos].fromTo.toFrom += 1.;
			tSum[pos].fromTo.toSum += 1.;
		} 
	}

public:

	TPsi(std::string_view Psi); 
	TPsi(const BAM::RGInfo::TInfo & info);

	BAM::RGInfo::TInfo info() const; 

	template<Type From_To>
	coretools::Probability prob(const BAM::TSequencedBase &data) const noexcept {
		const auto realType = data.isReverseStrand() ? _flip[From_To] : From_To;
		const auto end      = data.distFrom5Prime < data.distFrom3Prime ? End::from5 : End::from3;
		const auto pos      = end == End::from5 ? data.distFrom5Prime : data.distFrom3Prime;

		const auto &table = _tables[end][realType];
		if (pos >= table.size()) return table.back();
		return table[pos];
	}

	void addGA(const BAM::TSequencedBase &data, coretools::Probability P_g_I_di,
			 const TBaseBaseProbabilities& bbProbs) noexcept {
		_add<Type::GA>(data, P_g_I_di, bbProbs);
	}

	void addCT(const BAM::TSequencedBase &data, coretools::Probability P_g_I_di,
			 const TBaseBaseProbabilities& bbProbs) noexcept {
		_add<Type::CT>(data, P_g_I_di, bbProbs);
	}

	void add(const BAM::TSequencedBase &data, genometools::Base ref) noexcept {
		_add<Type::CT>(data, ref);
		_add<Type::GA>(data, ref);
	}

	void estimate() noexcept;
	void estimateInit() noexcept;

	void log() const noexcept;
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
