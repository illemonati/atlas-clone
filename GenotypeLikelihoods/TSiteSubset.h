/*
 * TSiteSubset.h
 *
 *  Created on: Nov 16, 2015
 *      Author: wegmannd
 */

#ifndef TSITESUBSET_H_
#define TSITESUBSET_H_

#include <stddef.h>
#include <stdint.h>
#include <functional>
#include <set>
#include <string>
#include <vector>

#include "GenotypeTypes.h"
#include "TGenomePosition.h"

namespace BAM { class TFastaBuffer; }
namespace coretools { class TOutputFile; }
namespace genometools { class TChromosomes; }

namespace GenotypeLikelihoods{

//-----------------------------------------
// A class to store site per window used for masking and filtering
// Positions are 0-based
//-----------------------------------------

//-----------------------------------------------
// TSiteSubsetSite
//-----------------------------------------------
class TSiteSubsetSite:public genometools::TGenomePosition{
private:
	genometools::Base _ref, _alt;
public:
	constexpr TSiteSubsetSite(uint32_t refID, uint32_t position, genometools::Base Ref, genometools::Base Alt)
		: TGenomePosition(refID, position), _ref(Ref), _alt(Alt){};
	constexpr TSiteSubsetSite(const genometools::TGenomePosition &Position, genometools::Base Ref,
							  genometools::Base Alt)
		: TGenomePosition(Position), _ref(Ref), _alt(Alt){};
	void write(coretools::TOutputFile & out) const;

	constexpr genometools::Base ref() const noexcept { return _ref; };
	constexpr genometools::Base alt() const noexcept { return _alt; };
};

//-----------------------------------------------
// TSiteSubset
//-----------------------------------------------
class TSiteSubset{
private:
	std::set<uint32_t> _refIDUsed;
	std::set<TSiteSubsetSite, std::less<>> _sites;

	bool _storesInvariantSites;

	void _readFile(const std::string &Filename, const genometools::TChromosomes & Chromosomes);
	void _readFile(const std::string &Filename, const genometools::TChromosomes & Chromosomes, BAM::TFastaBuffer & Reference);

public:
	TSiteSubset(const std::string &Filename, const genometools::TChromosomes & Chromosomes, bool InvariantSites);
	TSiteSubset(const std::string &Filename, const genometools::TChromosomes & Chromosomes, bool InvariantSites, BAM::TFastaBuffer & Reference);
	void write(const std::string &Filename) const;
	bool hasPositionsInWindow(const genometools::TGenomeWindow & Window) const;
	std::set<TSiteSubsetSite> getPositionInWindow(const genometools::TGenomeWindow & Window) const;
	size_t size() const noexcept { return _sites.size(); }
};

}; //end namespace

#endif /* TSITESUBSET_H_ */
