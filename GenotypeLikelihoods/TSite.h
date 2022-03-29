/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include <vector>

#include "TSequencedBase.h"
#include "GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TSubsamplePicker.h"

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
	genometools::Base refBase = genometools::Base::N;
	genometools::Genotype genotype  = genometools::Genotype::NN;

	void clear() noexcept;

	// access
	BAM::TSequencedBase &operator[](size_t i) noexcept { return _bases[i]; };
	const BAM::TSequencedBase &operator[](size_t i) const noexcept{ return _bases[i]; };

	// add
	void add(const BAM::TSequencedBase &base);
	void add(BAM::TSequencedBase &&base);
	void addToBaseFrequencies(TBaseData &frequencies) const noexcept;
	void downsample(uint32_t maxDepth, const coretools::TSubsamplePicker &picker);

	// getters
	bool empty() const noexcept { return _bases.empty(); };
	uint32_t depth() const noexcept { return _bases.size(); };
	uint32_t refDepth() const;
	std::string getBases() const;
	std::string getQualities() const;

	void countAlleles(TBaseCounts &alleleCounts) const;
	void countMates(std::array<int, 2>& mateCounts) const;
	void countFwdRev(std::array<int, 2> &) const;

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
