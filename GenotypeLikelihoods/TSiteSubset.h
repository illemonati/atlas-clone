/*
 * TSiteSubset.h
 *
 *  Created on: Nov 16, 2015
 *      Author: wegmannd
 */

#ifndef TSITESUBSET_H_
#define TSITESUBSET_H_

#include <fstream>
#include <set>
#include <map>
#include "TLog.h"
#include "stringFunctions.h"

#include "TFastaBuffer.h"
#include "TChromosomes.h"
#include "TGenotypeMap.h"
#include "TFile.h"
#include "TGenomePosition.h"

namespace GenotypeLikelihoods{

//-----------------------------------------
// A class to store site per window used for masking and filtering
// Positions are 0-based
//-----------------------------------------

//-----------------------------------------------
// TSiteSubsetSite
//-----------------------------------------------
class TSiteSubsetSite:public BAM::TGenomePosition{
private:
	genometools::Base _ref, _alt;

public:
	TSiteSubsetSite(const uint32_t & refID, const uint32_t & position, const genometools::Base & Ref, const genometools::Base & Alt);
	TSiteSubsetSite(const BAM::TGenomePosition & Position, const genometools::Base & Ref, const genometools::Base & Alt);
	TSiteSubsetSite(const TSiteSubsetSite & other) = default;
	void write(coretools::TOutputFile & out) const;

	genometools::Base ref() const{ return _ref; };
	genometools::Base alt() const{ return _alt; };
};

//-----------------------------------------------
// TSiteSubset
//-----------------------------------------------
class TSiteSubset{
private:
	std::set<uint32_t> _refIDUsed;
	std::set<TSiteSubsetSite, std::less<>> _sites;

	std::vector<TSiteSubsetSite> empty; //an empty vector to be returned in case there are no positions
	bool _storesInvariantSites;

	void _checkAlleles(const std::string & chr, const uint32_t & pos, const genometools::Base & ref, const genometools::Base & alt, const std::string & refAllele, const std::string & altAllele);
	void _readFile(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile);
	void _readFile(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile, BAM::TFastaBuffer & Reference);

public:
	TSiteSubset(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile, bool InvariantSites);
	TSiteSubset(const std::string Filename, const BAM::TChromosomes & Chromosomes, coretools::TLog* Logfile, bool InvariantSites, BAM::TFastaBuffer & Reference);
	void write(const std::string Filename) const;
	bool hasPositionsInWindow(const BAM::TGenomeWindow & Window) const;
	std::set<TSiteSubsetSite> getPositionInWindow(const BAM::TGenomeWindow & Window) const;
	size_t size();
};

}; //end namespace

#endif /* TSITESUBSET_H_ */
