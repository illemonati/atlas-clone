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


namespace GenotypeLikelihoods{

//-----------------------------------------
// A class to store site per window used for masking and filtering
// Positions are 0-based
//-----------------------------------------

//-----------------------------------------------
// TSiteSubsetSite
//-----------------------------------------------
class TSiteSubsetSite{
public:
	uint32_t position;
	Base ref;
	Base alt;

	TSiteSubsetSite(const uint32_t Position, const Base Ref, const Base Alt);
	void print() const;

	//operators: needed for sorting and finding
	bool operator ==(const TSiteSubsetSite & other) const{
		return position == other.position;
	};
	bool operator <(const TSiteSubsetSite & other) const{
		return position < other.position;
	};
	bool operator <(const uint32_t & pos) const{
		return position < pos;
	};
};

//-----------------------------------------------
// TSiteSubsetWindow
//-----------------------------------------------
class TSiteSubsetWindow{
private:
	uint32_t _start, _end;
	std::set<TSiteSubsetSite> _positions; //stores reference and alternative allele

public:
	TSiteSubsetWindow(uint32_t Start, uint32_t End);
	~TSiteSubsetWindow();

	void addPosition(const uint32_t Position, const Base Ref, const Base Alt);
	void print();
	size_t size();
};

//-----------------------------------------------
// TSiteSubsetChr
//-----------------------------------------------
class TSiteSubsetChr{
private:
	std::map<uint32_t, TSiteSubsetWindow>::iterator _findWindow(const unsigned int & pos);
	std::map<uint32_t, TSiteSubsetWindow>::iterator _findOrCreateWindow(const unsigned int & pos);

public:
	std::map<uint32_t, TSiteSubsetWindow> windows;
	uint32_t windowSize;

	TSiteSubsetChr(unsigned int & WindowSize);
	~TSiteSubsetChr();

	void addPosition(const uint32_t Position, const Base Ref, const Base Alt);

	void addPosition(std::vector<std::string> & tmp, const std::string & chr, bool invariantSites);
	bool addPosition(std::vector<std::string> & tmp, const std::string & chr, const char refBase, std::string & error, const bool invariantSites);

	void print(const std::string chrName);
	bool hasPositionsInWindow(const unsigned int & windowStart);
	std::set<TSiteSubsetSite>& getPositionInWindow(const unsigned int & windowStart);
	size_t size();
};

//-----------------------------------------------
// TSiteSubset
//-----------------------------------------------
class TSiteSubset{
private:
	std::string _filename;
	std::map<std::string, TSiteSubsetChr> _chromosomes;
	std::map<std::string, TSiteSubsetChr>::iterator curChrIt;
	std::vector<TSiteSubsetSite> empty; //an empty vector to be returned in case there are no positions
	uint32_t _windowSize;
	bool _storesInvariantSites;

	std::map<std::string, TSiteSubsetChr>::iterator _findOrCreateChromosome(const std::string chrName);
	void _checkAlleles(const std::string & chr, const uint32_t & pos, const Base & ref, const Base & alt, const std::string & refAllele, const std::string & altAllele);
	void _readFile(TLog* logfile);
	void _readFile(BAM::TFastaBuffer & reference, BAM::TChromosomes & Chromosomes, TLog* logfile);

public:
	TSiteSubset(std::string Filename, uint32_t & WindowSize, TLog* logfile, bool InvariantSites);
	TSiteSubset(std::string Filename, uint32_t WindowSize, TLog* logfile, bool InvariantSites, BAM::TFastaBuffer & reference, BAM::TChromosomes & chromosomes);
	void setChr(const std::string chr);
	void print();
	bool hasPositionsInWindow(const unsigned int & windowStart);
	std::set<TSiteSubsetSite>& getPositionInWindow(const unsigned int & windowStart);
	size_t size();
};

}; //end namespace

#endif /* TSITESUBSET_H_ */
