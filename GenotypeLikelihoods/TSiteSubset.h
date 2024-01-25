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

#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/algorithms.h"
#include "coretools/Containers/TView.h"

#include "genometools/GenomePositions/TChromosomes.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods{

//-----------------------------------------
// A class to store site per window used for masking and filtering
// Positions are 0-based
//-----------------------------------------

//-----------------------------------------------
// TSitePolymorphic / TSiteMonomorphic
//-----------------------------------------------
namespace SiteSubset {
	class TSitePolymorphic:public genometools::TGenomePosition{
	private:
		genometools::Base _ref, _alt;
	public:
		TSitePolymorphic(uint32_t refID, uint32_t position, const std::vector<std::string_view> & Line, const genometools::TChromosomes & Chromosomes);		

		static std::vector<std::string> getHeader() {
			return {"Chr", "Pos", "Allele1", "Allele2"};
		};
		void write(coretools::TOutputFile & out, const genometools::TChromosomes & Chromosomes) const;		
		
		constexpr genometools::Base ref() const noexcept { return _ref; };
		constexpr genometools::Base alt() const noexcept { return _alt; };		
	};

	class TSiteMonomorphic:public genometools::TGenomePosition{
	private:
		genometools::Base _ref;
	public:
		TSiteMonomorphic(uint32_t refID, uint32_t position, const std::vector<std::string_view> & Line, const genometools::TChromosomes & Chromosomes);
			
		static std::vector<std::string> getHeader() {
			return {"Chr", "Pos", "Allele1"};
		};
		void write(coretools::TOutputFile & out, const genometools::TChromosomes & Chromosomes) const;

		constexpr genometools::Base ref() const noexcept { return _ref; };		
	};

//-----------------------------------------------
// TSiteSubset
//-----------------------------------------------
template <typename SiteType>
class TSiteSubset{
private:
	std::vector<SiteType> _sites;
	std::string _filename;

	void _readFile(const std::string &Filename, const genometools::TChromosomes & Chromosomes){
		coretools::instances::logfile().listFlushTime("Reading sites to be used from '" + Filename + "' ...");

		// open file
		_filename = Filename;
		coretools::TInputFile in(Filename, coretools::FileType::Header);
		const auto indices = in.indices(SiteType::getHeader());

		// read file and add sites
		std::set<uint32_t> refIDUsed;		
		for (;!in.empty(); in.popFront()) {
			// get chromosome and position: throws error if chromosome does not exist
			const genometools::TChromosome &chr = Chromosomes.getChromosome(in.get(indices[0]));
			refIDUsed.emplace(chr.refID());
			uint32_t pos = coretools::str::fromString<uint32_t, true>(in.get(indices[1])) - 1; //make 0-based

			// add site			
			_sites.emplace_back(chr.refID(), pos, in.front(), Chromosomes);
		}

		//sort sites
		std::sort(_sites.begin(), _sites.end());
		
		//check for duplicates
		auto res = coretools::findDuplicate(_sites);
		if(res.first){
			UERROR("Duplicates in sites file '", Filename, "': ", _sites[res.second].asFormattedString(Chromosomes), " is present multiple times!");
		}

		// report
		coretools::instances::logfile().doneTime();
		coretools::instances::logfile().conclude("Parsed ", size(), " sites on ", refIDUsed.size(), " chromosomes.");
	};

public:
	TSiteSubset(const std::string &Filename, const genometools::TChromosomes & Chromosomes){
		_readFile(Filename, Chromosomes);
	};
	
	void write(const std::string &Filename, const genometools::TChromosomes & Chromosomes) const {		
		coretools::TOutputFile out(Filename, _sites[0].getHeader());
		for (auto &s : _sites) { s.write(out, Chromosomes); }
	};
	bool hasPositionsInWindow(const genometools::TGenomeWindow & Window) const {
		auto it = std::lower_bound(_sites.begin(), _sites.end(), Window);
		return !(it == _sites.end() || *it < Window);	
	};

	coretools::TConstView<SiteType> getPositionInWindow(const genometools::TGenomeWindow & Window) const {
		coretools::TConstView<SiteType> view(_sites);
		const auto start = std::lower_bound(_sites.begin(), _sites.end(), Window);
		auto end         = start;
		while (end != _sites.end() && Window.within(*end)){
 			++end;
		}
		return view.subview(std::distance(_sites.begin(), start), end - start);
	}
	
	[[nodiscard]] size_t size() const noexcept { return _sites.size(); }
	[[nodiscard]] std::string filename() const { return _filename; }
}; 

} //end namespace sitesubset

using TSiteSubsetPolymorphic = SiteSubset::TSiteSubset<SiteSubset::TSitePolymorphic>;
using TSiteSubsetMonomorphic = SiteSubset::TSiteSubset<SiteSubset::TSiteMonomorphic>;

}; //end namespace

#endif /* TSITESUBSET_H_ */
