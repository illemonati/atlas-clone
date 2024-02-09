/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include <array>
#include <string>
#include <vector>

#include "coretools/Math/TSubsamplePicker.h"

#include "genometools/GenotypeTypes.h"

#include "GenotypeData.h"
#include "TSequencedBase.h"

namespace GenotypeLikelihoods {

// TODO: write as templated classes
//----------------------------------------------------------------------------------------------------------------------------------
//  TSite
//  Class that stores bases.
//----------------------------------------------------------------------------------------------------------------------------------
class TSite {
private:
	std::vector<BAM::TSequencedBase> _bases;
public:
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;
	genometools::Base refBase      = genometools::Base::N;
	genometools::Genotype genotype = genometools::Genotype::NN;

	void clear() noexcept;

	// access
	BAM::TSequencedBase &operator[](size_t i) noexcept { return _bases[i]; };
	const BAM::TSequencedBase &operator[](size_t i) const noexcept{ return _bases[i]; };

	// add
	void add(const BAM::TSequencedBase &base);
	TBaseData baseFrequencies() const noexcept;
	void downsample(size_t maxDepth, const coretools::TSubsamplePicker &picker);

	// getters
	bool empty() const noexcept { return _bases.empty(); };
	size_t depth() const noexcept { return _bases.size(); };
	size_t refDepth() const;
	std::string getBases() const;
	std::vector<genometools::Base> sampleBases() const;
	std::string getQualities() const;

	TBaseCounts countAlleles() const;
	coretools::TStrongArray<size_t, BAM::Mate> countMates() const;
	std::array<int, 2> countFwdRev() const;

	// loop
	using iterator       = std::vector<BAM::TSequencedBase>::iterator;
	using const_iterator = std::vector<BAM::TSequencedBase>::const_iterator;

	iterator begin() noexcept { return _bases.begin(); };
	iterator end() noexcept{ return _bases.end(); };
	const_iterator begin() const noexcept { return _bases.cbegin(); };
	const_iterator end() const noexcept { return _bases.cend(); };
	const_iterator cbegin() const noexcept { return _bases.cbegin(); };
	const_iterator cend() const noexcept { return _bases.cend(); };
};

}; // namespace GenotypeLikelihoods

#endif /* TSITE_H_ */
