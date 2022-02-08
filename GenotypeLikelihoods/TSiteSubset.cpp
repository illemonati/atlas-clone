/*
 * TSiteSubset.cpp
 *
 *  Created on: Nov 25, 2019
 *      Author: vivian
 */

#include "TSiteSubset.h"
#include "GenotypeTypes.h"

namespace GenotypeLikelihoods{

using coretools::str::toString;

//-------------------------------------------------
// TSiteSubsetSite
//-------------------------------------------------
TSiteSubsetSite::TSiteSubsetSite(uint32_t refID, uint32_t position, const genometools::Base & Ref, const genometools::Base & Alt):TGenomePosition(refID, position){
	_ref = Ref;
	_alt = Alt;
};

TSiteSubsetSite::TSiteSubsetSite(const BAM::TGenomePosition & Position, const genometools::Base & Ref, const genometools::Base & Alt):TGenomePosition(Position){
	_ref = Ref;
	_alt = Alt;
};

void TSiteSubsetSite::write(coretools::TOutputFile & out) const{
	out << _refID << _position << _ref << _alt << std::endl;
};

//-------------------------------------------------
// TSiteSubset
//-------------------------------------------------
void TSiteSubset::_checkAlleles(const std::string & chr, uint32_t pos, const genometools::Base & ref, const genometools::Base & alt, const std::string & refAllele, const std::string & altAllele){
	if(ref == genometools::N){
		throw "Unknown reference allele '" + refAllele + "' on chr " + chr + " at " + toString(pos+1) + "!";
	}
	if(alt == genometools::N){
		throw "Unknown alternative allele '" + altAllele + "' on chr " + chr + " at " + toString(pos+1) + "!";
	}

	if(ref == alt && !_storesInvariantSites){
		throw "Site on chr " + chr + " at " + toString(pos+1) + " is invariant!";
	}
	if(ref != alt && _storesInvariantSites){
		throw "Site on chr " + chr + " at " + toString(pos+1) + " is polymorphic!";
	}
};

void TSiteSubset::_readFile(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile){
	Logfile->listFlushTime("Reading sites to be used from '" + Filename + "' ...");

	//open file
	coretools::TInputFile in(Filename, 4);

	//read file and add sites
	std::vector<std::string> line;
	while(in.read(line)){
		//get chromosome: throws error if chromosome does not exist
		const BAM::TChromosome& chr = Chromosomes.getChromosome(line[0]);
		_refIDUsed.emplace(chr.refID());

		//extract positions
		uint32_t pos = coretools::str::convertStringCheck<uint32_t>(line[1]) - 1; //make 0-based
		const genometools::Base ref = genometools::fromChar(line[2][0]);
		const genometools::Base alt = genometools::fromChar(line[3][0]);

		//check alleles
		_checkAlleles(chr.name, pos, ref, alt, line[2], line[3]);

		//add site
		_sites.emplace(chr.refID(), pos, ref, alt);
	}

	//report
	Logfile->doneTime();
	Logfile->conclude("Parsed " + toString(_sites.size()) + " sites on " + toString(_refIDUsed.size()) + " chromosomes.");
};

void TSiteSubset::_readFile(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile, BAM::TFastaBuffer & Reference){ //version that checks witth fasta reference
	Logfile->listFlushTime("Reading sites to be used from '" + Filename + "' ...");

	//open file
	coretools::TInputFile in(Filename, 4);

	//conflicts with ref
	bool conflictsFound = false;
	coretools::TOutputFile conflicts;

	//read file and add sites
	std::vector<std::string> line;
	while(in.read(line)){
		//get chromosome: throws error if chromosome does not exist
		const BAM::TChromosome& chr = Chromosomes.getChromosome(line[0]);
		_refIDUsed.emplace(chr.refID());

		//extract positions
		uint32_t pos = coretools::str::convertStringCheck<uint32_t>(line[1]) - 1; //make 0-based
		const genometools::Base ref = genometools::fromChar(line[2][0]);
		const genometools::Base alt = genometools::fromChar(line[3][0]);

		//check alleles
		_checkAlleles(chr.name, pos, ref, alt, line[2], line[3]);

		//check with reference
		BAM::TGenomePosition genoPos(chr.refID(), pos);
		genometools::Base trueRef = Reference.refAt(genoPos);
		if(trueRef != ref && trueRef != alt){
			//conflict with fasta
			if(!conflictsFound){
				conflicts.open(Filename + ".conflicts.gz");
				conflicts.writeHeader({"chr", "position", "reference", "allele1", "allele2"});
			}
			conflicts << chr.name << pos+1 << trueRef << ref << alt << std::endl;
		}

		//add site
		_sites.emplace(chr.refID(), pos, ref, alt);
	}

	//report
	Logfile->doneTime();
	Logfile->conclude("Parsed " + toString(_sites.size()) + " sites on " + toString(_refIDUsed.size()) + " chromosomes.");

	//write conflicts, if any
	if(conflictsFound){
		Logfile->conclude("Reference conflicted with provided alleles at " + toString(conflicts.lineNumber() - 1) + " positions!");
		Logfile->conclude("These positions were written to " + conflicts.name());
	}
};

TSiteSubset::TSiteSubset(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile, bool InvariantSites){
	_storesInvariantSites = InvariantSites;
	_readFile(Filename, Chromosomes, Logfile);
};

TSiteSubset::TSiteSubset(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile, bool InvariantSites, BAM::TFastaBuffer & Reference){
	_storesInvariantSites = InvariantSites;
	_readFile(Filename, Chromosomes, Logfile, Reference);
};

void TSiteSubset::write(const std::string Filename) const{
	coretools::TOutputFile out(Filename, 4);
	for(auto& s : _sites){
		s.write(out);
	}
};

bool TSiteSubset::hasPositionsInWindow(const BAM::TGenomeWindow & Window) const{
	auto it = _sites.lower_bound(Window);
	if(it == _sites.end() || *it < Window){
		return false;
	} else {
		return true;
	}
};

std::set<TSiteSubsetSite> TSiteSubset::getPositionInWindow(const BAM::TGenomeWindow & Window) const{
	std::set<TSiteSubsetSite> set;
	auto it = _sites.lower_bound(Window);
	while(it != _sites.end() && *it < Window){
		set.emplace(*it);
		++it;
	}

	return set;
};

size_t TSiteSubset::size(){
	return _sites.size();
};

}; //end namespace
