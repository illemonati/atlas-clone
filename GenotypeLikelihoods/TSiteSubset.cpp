/*
 * TSiteSubset.cpp
 *
 *  Created on: Nov 25, 2019
 *      Author: vivian
 */

#include "TSiteSubset.h"

#include <cstdint>
#include <fstream>
#include <map>

#include "genometools/GenotypeTypes.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "TFastaBuffer.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringFunctions.h"

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
	while(in.read(line)){
		//get chromosome: throws error if chromosome does not exist
		const genometools::TChromosome& chr = Chromosomes.getChromosome(line[0]);
		_refIDUsed.emplace(chr.refID());

		//extract positions
		uint32_t pos = coretools::str::convertStringCheck<uint32_t>(line[1]) - 1; //make 0-based
		const genometools::Base ref = genometools::char2base(line[2][0]);
		const genometools::Base alt = genometools::char2base(line[3][0]);

		//check alleles
		impl::checkAlleles(chr.name, pos, ref, alt, line[2], line[3], _storesInvariantSites);

		//add site
		_sites.emplace(chr.refID(), pos, ref, alt);
	}

	//report
	logfile().doneTime();
	logfile().conclude("Parsed " + toString(_sites.size()) + " sites on " + toString(_refIDUsed.size()) + " chromosomes.");
};

void TSiteSubset::_readFile(const std::string &Filename, const genometools::TChromosomes & Chromosomes, BAM::TFastaBuffer & Reference){ //version that checks witth fasta reference
	logfile().listFlushTime("Reading sites to be used from '" + Filename + "' ...");

	//open file
	coretools::TInputFile in(Filename, 4);

	//conflicts with ref
	bool conflictsFound = false;
	coretools::TOutputFile conflicts;

	//read file and add sites
	std::vector<std::string> line;
	while(in.read(line)){
		//get chromosome: throws error if chromosome does not exist
		const genometools::TChromosome& chr = Chromosomes.getChromosome(line[0]);
		_refIDUsed.emplace(chr.refID());

		//extract positions
		uint32_t pos = coretools::str::convertStringCheck<uint32_t>(line[1]) - 1; //make 0-based
		const genometools::Base ref = genometools::char2base(line[2][0]);
		const genometools::Base alt = genometools::char2base(line[3][0]);

		//check alleles
		impl::checkAlleles(chr.name, pos, ref, alt, line[2], line[3], _storesInvariantSites);

		//check with reference
        genometools::TGenomePosition genoPos(chr.refID(), pos);
		genometools::Base trueRef = Reference.refAt(genoPos);
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
		_sites.emplace(chr.refID(), pos, ref, alt);
	}

	//report
	logfile().doneTime();
	logfile().conclude("Parsed " + toString(_sites.size()) + " sites on " + toString(_refIDUsed.size()) + " chromosomes.");

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

TSiteSubset::TSiteSubset(const std::string &Filename, const genometools::TChromosomes & Chromosomes, bool InvariantSites, BAM::TFastaBuffer & Reference){
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
	auto it = _sites.lower_bound(Window);
	return !(it == _sites.end() || *it < Window);
};

std::set<TSiteSubsetSite> TSiteSubset::getPositionInWindow(const genometools::TGenomeWindow & Window) const{
	std::set<TSiteSubsetSite> set;
	auto it = _sites.lower_bound(Window);
	while(it != _sites.end() && *it < Window){
		set.emplace(*it);
		++it;
	}

	return set;
};

}; //end namespace
