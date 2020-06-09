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
#include "commonutilities/TLog.h"
#include "commonutilities/stringFunctions.h"

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
class TSiteSubsetSite{
private:
	BAM::TGenomePosition _genomicPosition;
	Base _ref, _alt;

public:
	TSiteSubsetSite(const BAM::TGenomePosition Position, const Base Ref, const Base Alt);
	TSiteSubsetSite(const TSiteSubsetSite & other) = default;
	void write(TOutputFile & out) const;

	uint32_t refID() const{ return _genomicPosition.refId(); };
	uint32_t position() const{ return _genomicPosition.position(); };
	Base ref() const{ return _ref; };
	Base alt() const{ return _alt; };

	//operators: needed for sorting and finding
	bool operator ==(const TSiteSubsetSite & other) const{
		return _genomicPosition == other._genomicPosition;
	};
	bool operator<(const TSiteSubsetSite & other) const{
		return _genomicPosition < other._genomicPosition;
	};
	bool operator<(const BAM::TGenomePosition pos) const{
		return _genomicPosition < pos;
	};
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

	void _checkAlleles(const std::string & chr, const uint32_t & pos, const Base & ref, const Base & alt, const std::string & refAllele, const std::string & altAllele);
	void _readFile(const std::string Filename, const BAM::TChromosomes & Chromosomes, const TGenotypeMap & GenoMap, TLog* Logfile);
	void _readFile(const std::string Filename, const BAM::TChromosomes & Chromosomes, const TGenotypeMap & GenoMap, TLog* Logfile, BAM::TFastaBuffer & Reference);

public:
	TSiteSubset(const std::string Filename, const BAM::TChromosomes & Chromosomes, const TGenotypeMap & GenoMap, TLog* Logfile, bool InvariantSites);
	TSiteSubset(const std::string Filename, const BAM::TChromosomes & Chromosomes, const TGenotypeMap & GenoMap, TLog* Logfile, bool InvariantSites, BAM::TFastaBuffer & Reference);
	void setChr(const std::string chr);
	void write(const std::string Filename) const;
	bool hasPositionsInWindow(const BAM::TGenomePosition & Start, const BAM::TGenomePosition & End) const;
	std::set<TSiteSubsetSite> getPositionInWindow(const BAM::TGenomePosition & Start, const BAM::TGenomePosition & End) const;
	size_t size();
};

}; //end namespace

#endif /* TSITESUBSET_H_ */
