/*
 * TSiteSubset.h
 *
 *  Created on: Nov 16, 2015
 *      Author: wegmannd
 */

#ifndef TSITESUBSET_H_
#define TSITESUBSET_H_

#include <set>
#include <string>
#include <vector>

#include "coretools/Containers/TView.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/algorithms.h"

#include "genometools/GenomePositions/TChromosomes.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods {

//-----------------------------------------
// A class to store site per window used for masking and filtering
// Positions are 0-based
//-----------------------------------------

namespace SiteSubset {
class TSitePolymorphic : public genometools::TGenomePosition {
private:
	genometools::TGenomePosition _pos;
	genometools::Base _ref, _alt;

public:
	TSitePolymorphic(uint32_t refID, uint32_t position,
                   char ref, char alt,
                   const genometools::TChromosomes &Chromosomes);

	constexpr genometools::Base ref() const noexcept { return _ref; };
	constexpr genometools::Base alt() const noexcept { return _alt; };
	constexpr const genometools::TGenomePosition& pos() const noexcept { return *this; };

	friend bool operator<(const TSitePolymorphic& lhs, const TSitePolymorphic& rhs) noexcept {return lhs.pos() < rhs.pos();}
};

class TSiteMonomorphic : public genometools::TGenomePosition {
private:
  genometools::Base _ref;

public:
  TSiteMonomorphic(uint32_t refID, uint32_t position, char ref,
                   const genometools::TChromosomes &Chromosomes);

  constexpr genometools::Base ref() const noexcept { return _ref; };
};

//-----------------------------------------------
// TSiteSubset
//-----------------------------------------------
template <typename SiteType> class TSiteSubset {
private:
  std::vector<SiteType> _sites;
  std::string _filename;

  void _readFile(const std::string &Filename,
                 const genometools::TChromosomes &Chromosomes) {
    coretools::instances::logfile().listFlushTime(
        "Reading sites to be used from '", Filename, "' ...");

    // open file
    _filename = Filename;
    coretools::TInputFile in(Filename, coretools::FileType::Header);

    // read file and add sites
    std::set<uint32_t> refIDUsed;
    std::vector<size_t> idx;
    if constexpr((std::is_same<SiteType, TSitePolymorphic>::value)){
      idx = in.indices<std::vector<size_t>, std::vector<std::string>>({"Chr", "Pos", "Ref", "Alt"});
    } else {
      idx = in.indices<std::vector<size_t>, std::vector<std::string>>({"Chr", "Pos", "Ref"});
    }

    for (; !in.empty(); in.popFront()) {
      // get chromosome and position: throws error if chromosome does not exist
      const genometools::TChromosome &chr =
          Chromosomes.getChromosome(in.get(idx[0]));
      refIDUsed.emplace(chr.refID());
      uint32_t pos = coretools::str::fromString<uint32_t, true>(in.get(idx[1])) -
                     1; // make 0-based

      // add site
      if constexpr((std::is_same<SiteType, TSitePolymorphic>::value)){
        _sites.emplace_back(chr.refID(), pos, in.get(idx[2])[0], in.get(idx[3])[0], Chromosomes);
      } else {
        _sites.emplace_back(chr.refID(), pos, in.get(idx[2])[0], Chromosomes);
      }
    }

    // sort sites
    std::sort(_sites.begin(), _sites.end());

    // check for duplicates
    auto res = coretools::findDuplicate(_sites);
    if (res.first) {
      UERROR("Duplicates in sites file '", Filename,
             "': ", _sites[res.second].asFormattedString(Chromosomes),
             " is present multiple times!");
    }

    // report
    coretools::instances::logfile().doneTime();
    coretools::instances::logfile().conclude("Parsed ", size(), " sites on ",
                                             refIDUsed.size(), " chromosomes.");
  };

public:
  TSiteSubset(const std::string &Filename,
              const genometools::TChromosomes &Chromosomes) {
    _readFile(Filename, Chromosomes);
  };

  void write(const std::string &Filename,
             const genometools::TChromosomes &Chromosomes) const {
    coretools::TOutputFile out(Filename, _sites[0].getHeader());
    for (auto &s : _sites) {
      s.write(out, Chromosomes);
    }
  };
  bool hasPositionsInWindow(const genometools::TGenomeWindow &Window) const {
    auto it = std::lower_bound(_sites.begin(), _sites.end(), Window);
    return !(it == _sites.end() || *it < Window);
  };

  coretools::TConstView<SiteType>
  getPositionInWindow(const genometools::TGenomeWindow &Window) const {
    coretools::TConstView<SiteType> view(_sites);
    const auto start = std::lower_bound(_sites.begin(), _sites.end(), Window);
    auto end = start;
    while (end != _sites.end() && Window.within(*end)) {
      ++end;
    }
    return view.subview(std::distance(_sites.begin(), start), end - start);
  }

  [[nodiscard]] size_t size() const noexcept { return _sites.size(); }
  [[nodiscard]] std::string filename() const { return _filename; }
};

} // namespace SiteSubset

using TSiteSubsetPolymorphic =
    SiteSubset::TSiteSubset<SiteSubset::TSitePolymorphic>;
using TSiteSubsetMonomorphic =
    SiteSubset::TSiteSubset<SiteSubset::TSiteMonomorphic>;

}; // namespace GenotypeLikelihoods

#endif /* TSITESUBSET_H_ */
