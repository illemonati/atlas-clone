/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRHO_H_
#define TRHO_H_

#include "TSequencedData.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"

#include "genometools/Genotypes/Containers.h"
#include "TReadGroupInfo.h"

namespace GenotypeLikelihoods::SequencingError {

//--------------------------------------------------------------------
// TRho
//--------------------------------------------------------------------
class TRho {
private:
	coretools::TStrongArray<genometools::TBaseProbabilities, genometools::Base> _rho{
		{genometools::TBaseProbabilities::normalize({0., 1. / 3, 1. / 3, 1. / 3}),
		 genometools::TBaseProbabilities::normalize({1. / 3, 0., 1. / 3, 1. / 3}),
		 genometools::TBaseProbabilities::normalize({1. / 3, 1. / 3, 0., 1. / 3}),
		 genometools::TBaseProbabilities::normalize({1. / 3, 1. / 3, 1. / 3, 0.})}}; //[from k][to l]

	coretools::TStrongArray<coretools::TStrongArray<double, genometools::Base>, genometools::Base> _rhoSum{};
public:
	TRho(std::string_view Def);
	TRho(const BAM::RGInfo::TInfo & info);

	const genometools::TBaseProbabilities& operator[](genometools::Base from) const noexcept {
		return _rho[from];
	}

	const coretools::Probability prob(const BAM::TSequencedData& from, genometools::Base to) const noexcept {
		if (from.get<BAM::Flags::ReversedStrand>()) {
			using genometools::flipped;
			return _rho[flipped(from.base)][flipped(to)];
		} else {
			return _rho[from.base][to];
		}
	}

	const coretools::Probability prob(genometools::Base from, const BAM::TSequencedData& to) const noexcept {
		if (to.get<BAM::Flags::ReversedStrand>()) {
			using genometools::flipped;
			return _rho[flipped(from)][flipped(to.base)];
		} else {
			return _rho[from][to.base];
		}
	}

	// functions used to estimate
	void add(const BAM::TSequencedData& data, coretools::Probability P_g_I_d, const genometools::TBaseProbabilities &P_bbar_I_d) noexcept;
	void estimate() noexcept;

	void log() const;
	BAM::RGInfo::TInfo info() const;
};

	
}

#endif
