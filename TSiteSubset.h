/*
 * TSiteSubset.h
 *
 *  Created on: Nov 16, 2015
 *      Author: wegmannd
 */

#ifndef TSITESUBSET_H_
#define TSITESUBSET_H_

#include <fstream>
#include <vector>
#include <map>
#include "bamtools/utils/bamtools_fasta.h"
#include "TLog.h"
#include "bamtools/api/BamReader.h"
#include "stringFunctions.h"



//store sites 0-based!

class TSiteSubsetWindow{
public:
	bool hasData;
	unsigned int start, end;
	std::map< unsigned int, std::pair<char,char> > positions; //stores reference and alternative allele
	std::map< int, char > positions2; //stores reference and alternative allele

	TSiteSubsetWindow(unsigned int Start, unsigned int End);
	~TSiteSubsetWindow();

	void addPosition(unsigned int pos, char & ref, char & alt);
	void print();
	size_t size();
};

class TSiteSubsetChr{
public:
	std::string name;
	std::map<unsigned int, TSiteSubsetWindow*> windows;
	std::map<unsigned int, TSiteSubsetWindow*>::iterator windowIt;
	unsigned int windowSize;
	int chrNumberInFasta;

	TSiteSubsetChr(std::string & Name, unsigned int & WindowSize);
	TSiteSubsetChr(std::string & Name, unsigned int & WindowSize, BamTools::SamHeader bamHeader);
	~TSiteSubsetChr();

	void findWindow(const unsigned int & pos);
	void findOrCreateWindow(const unsigned int & pos);
	void addPosition(std::vector<std::string> & tmp, const std::string & chr, bool invariantSites);
	bool addPosition(std::vector<std::string> & tmp, const std::string & chr, BamTools::Fasta & reference, std::string & error, bool invariantSites);
	void print();
	bool hasPositionsInWindow(const unsigned int & windowStart);
	std::map<unsigned int,std::pair<char,char> >& getPositionInWindow(const unsigned int & windowStart);
	size_t size();
};

class TSiteSubset{
private:
	std::map<std::string, TSiteSubsetChr*> chromosomes;
	std::map<std::string, TSiteSubsetChr*>::iterator chrIt;
	unsigned int windowSize;
	std::string curChr;
	bool invariantSites;

	void readFile(TLog* logfile);

	void readFile(BamTools::Fasta & reference, BamTools::SamHeader bamHeader, TLog* logfile);

public:
	std::string filename;
	TSiteSubset(std::string Filename, int & WindowSize, TLog* logfile, bool InvariantSites);
	TSiteSubset(std::string Filename, BamTools::Fasta & reference, BamTools::SamHeader bamHeader, const unsigned int WindowSize, TLog* logfile, bool InvariantSites);
	~TSiteSubset();
	void setChr(const std::string & chr);
	void print();
	bool hasPositionsInWindow(const unsigned int & windowStart);
	std::map<unsigned int,std::pair<char,char> >& getPositionInWindow(const unsigned int & windowStart);
	size_t size();
};


#endif /* TSITESUBSET_H_ */
