/*
 * TSiteSubset.cpp
 *
 *  Created on: Nov 25, 2019
 *      Author: vivian
 */

#include "TSiteSubset.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <map>

#include "coretools/Containers/TView.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringFunctions.h"
#include "genometools/TFastaReader.h"

namespace GenotypeLikelihoods{

using coretools::str::toString;
using coretools::instances::logfile;

namespace impl {

void checkAlleles(const std::string &chr, uint32_t pos, const genometools::Base &ref, const genometools::Base &alt,
				  const std::string &refAllele, const std::string &altAllele, bool storesInvariantSites) {
	if (ref == genometools::Base::N) {
		throw "Unknown reference allele '" + refAllele + "' on chr " + chr + " at " + toString(pos + 1) + "!";
	}
	if (alt == genometools::Base::N) {
		throw "Unknown alternative allele '" + altAllele + "' on chr " + chr + " at " + toString(pos + 1) + "!";
	}

	if (ref == alt && !storesInvariantSites) {
		throw "Site on chr " + chr + " at " + toString(pos + 1) + " is invariant!";
	}
	if (ref != alt && storesInvariantSites) {
		throw "Site on chr " + chr + " at " + toString(pos + 1) + " is polymorphic!";
	}
};
} // namespace impl

//-------------------------------------------------
// TSiteSubsetSite
//-------------------------------------------------

void TSiteSubsetSite::write(coretools::TOutputFile & out) const{
	out.writeln(_refID, _position, _ref, _alt);
};

//-------------------------------------------------
// TSiteSubset
//-------------------------------------------------

void TSiteSubset::_readFile(const std::string &Filename, const genometools::TChromosomes & Chromosomes){
	logfile().listFlushTime("Reading sites to be used from '" + Filename + "' ...");

	//open file
	coretools::TInputFile in(Filename, 4);

	//read file and add sites
	std::vector<std::string> line;
	std::set<uint32_t> refIDUsed;
	std::set<TSiteSubsetSite, std::less<>> setSites;
	while(in.read(line)){
		//get chromosome: throws error if chromosome does not exist
		const genometools::TChromosome& chr = Chromosomes.getChromosome(line[0]);
		refIDUsed.emplace(chr.refID());

		//extract positions
		uint32_t pos = coretools::str::fromString<uint32_t, true>(line[1]) - 1; //make 0-based
		const genometools::Base ref = genometools::char2base(line[2][0]);
		const genometools::Base alt = genometools::char2base(line[3][0]);

		//check alleles
		impl::checkAlleles(chr.name, pos, ref, alt, line[2], line[3], _storesInvariantSites);

		//add site
		setSites.emplace(chr.refID(), pos, ref, alt);
	}

	_sites.clear();
	std::copy(setSites.begin(), setSites.end(), std::back_inserter(_sites));

	//report
	logfile().doneTime();
	logfile().conclude("Parsed " + toString(size()) + " sites on " + toString(refIDUsed.size()) + " chromosomes.");
};

void TSiteSubset::_readFile(const std::string &Filename, const genometools::TChromosomes & Chromosomes, const genometools::TFastaReader & Reference) {
	logfile().listFlushTime("Reading sites to be used from '" + Filename + "' ...");

	//open file
	coretools::TInputFile in(Filename, 4);

	//conflicts with ref
	bool conflictsFound = false;
	coretools::TOutputFile conflicts;

	//read file and add sites
	std::vector<std::string> line;
	std::set<uint32_t> refIDUsed;
	std::set<TSiteSubsetSite, std::less<>> setSites;
	while(in.read(line)){
		//get chromosome: throws error if chromosome does not exist
		const genometools::TChromosome& chr = Chromosomes.getChromosome(line[0]);
		refIDUsed.emplace(chr.refID());

		//extract positions
		uint32_t pos = coretools::str::fromString<uint32_t, true>(line[1]) - 1; //make 0-based
		const genometools::Base ref = genometools::char2base(line[2][0]);
		const genometools::Base alt = genometools::char2base(line[3][0]);

		//check alleles
		impl::checkAlleles(chr.name, pos, ref, alt, line[2], line[3], _storesInvariantSites);

		//check with reference
		genometools::Base trueRef = Reference(chr.refID(), pos);
		if(trueRef != ref && trueRef != alt){
			//conflict with fasta
			if(!conflictsFound){
				conflicts.open(Filename + ".conflicts.gz");
				conflicts.writeHeader({"chr", "position", "reference", "allele1", "allele2"});
				conflictsFound = true;
			}
			conflicts.writeln(chr.name, pos+1, trueRef, ref, alt);
		}

		//add site
		setSites.emplace(chr.refID(), pos, ref, alt);
	}

	_sites.clear();
	std::copy(setSites.begin(), setSites.end(), std::back_inserter(_sites));

	//report
	logfile().doneTime();
	logfile().conclude("Parsed " + toString(size()) + " sites on " + toString(refIDUsed.size()) + " chromosomes.");

	//write conflicts, if any
	if(conflictsFound){
		logfile().conclude("Reference conflicted with provided alleles at " + toString(conflicts.curLine() - 1) + " positions!");
		logfile().conclude("These positions were written to " + conflicts.name());
	}
};

TSiteSubset::TSiteSubset(const std::string &Filename, const genometools::TChromosomes & Chromosomes, bool InvariantSites){
	_storesInvariantSites = InvariantSites;
	_readFile(Filename, Chromosomes);
};

TSiteSubset::TSiteSubset(const std::string &Filename, const genometools::TChromosomes & Chromosomes, bool InvariantSites, const genometools::TFastaReader & Reference) {
	_storesInvariantSites = InvariantSites;
	_readFile(Filename, Chromosomes, Reference);
};

void TSiteSubset::write(const std::string &Filename) const{
	coretools::TOutputFile out(Filename, 4);
	for(auto& s : _sites){
		s.write(out);
	}
};

bool TSiteSubset::hasPositionsInWindow(const genometools::TGenomeWindow & Window) const{
	auto it = std::lower_bound(_sites.begin(), _sites.end(), Window);
	return !(it == _sites.end() || *it < Window);
};

coretools::TConstView<TSiteSubsetSite> TSiteSubset::getPositionInWindow(const genometools::TGenomeWindow & Window) const {
	coretools::TConstView<TSiteSubsetSite> view(_sites);
	const auto start = std::lower_bound(_sites.begin(), _sites.end(), Window);
	auto end = start;
	while(end != _sites.end() && *end < Window){
		++end;
	}

	return view.subview(std::distance(_sites.begin(), start), end-start);
};

}; //end namespace
