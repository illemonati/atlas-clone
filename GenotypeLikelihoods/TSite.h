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

#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"

#include "genometools/Genotypes/Containers.h"
#include "TSequencedData.h"

namespace GenotypeLikelihoods {

// TODO: write as templated classes
//----------------------------------------------------------------------------------------------------------------------------------
//  TSite
//  Class that stores bases.
//----------------------------------------------------------------------------------------------------------------------------------
class TSite {
public:
	std::vector<BAM::TSequencedData> _data;
public:
	genometools::TGenotypeLikelihoods genotypeLikelihoods;
	genometools::Base refBase      = genometools::Base::N;
	genometools::Genotype genotype = genometools::Genotype::NN;

	void clear() noexcept;

	// access
	BAM::TSequencedData &operator[](size_t i) noexcept { return _data[i]; };
	const BAM::TSequencedData &operator[](size_t i) const noexcept{ return _data[i]; };

	// add
	void add(const BAM::TSequencedData &data);
	genometools::TBaseData baseFrequencies() const noexcept;
	void shuffle();

	void downsample(coretools::Probability p);
	void downsample(size_t UpToDepth);

	// getters
	bool empty() const noexcept { return _data.empty(); };
	size_t depth() const noexcept { return _data.size(); };
	size_t refDepth() const;
	std::string getBases() const;
	std::vector<genometools::Base> sampleBases() const;
	std::string getQualities() const;

	genometools::TBaseCounts countAlleles() const;
	coretools::TStrongArray<size_t, BAM::Mate> countMates() const;
	std::array<int, 2> countFwdRev() const;

	// loop
	using iterator       = std::vector<BAM::TSequencedData>::iterator;
	using const_iterator = std::vector<BAM::TSequencedData>::const_iterator;

	iterator begin() noexcept { return _data.begin(); };
	iterator end() noexcept{ return _data.end(); };
	const_iterator begin() const noexcept { return _data.cbegin(); };
	const_iterator end() const noexcept { return _data.cend(); };
	const_iterator cbegin() const noexcept { return _data.cbegin(); };
	const_iterator cend() const noexcept { return _data.cend(); };
};

}; // namespace GenotypeLikelihoods

#endif /* TSITE_H_ */
