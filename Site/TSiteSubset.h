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
#include <iostream>
#include "../IOTools/IOAbstractClasses/Fasta.h"
#include "../IOTools/IOAbstractClasses/SamHeader.h"
#include "TLog.h"
#include "stringFunctions.h"

//store sites 0-based!

class TSiteSubsetWindow{
public:
	bool hasData;
	long start, end;
	std::map< long, std::pair<char,char> > positions; //stores reference and alternative allele
	std::map< int, char > positions2; //stores reference and alternative allele

    TSiteSubsetWindow(long Start, long End);
    ~TSiteSubsetWindow();

    void addPosition(long pos, char & ref, char & alt);
    void print();
    long size();
};

class TSiteSubsetChr{
public:
	std::string name;
	std::map<int, TSiteSubsetWindow*> windows;
	std::map<int, TSiteSubsetWindow*>::iterator windowIt;
	int windowSize;
	int chrNumberInFasta;


    TSiteSubsetChr(std::string & Name, int & WindowSize);
    TSiteSubsetChr(std::string & Name, int & WindowSize, SamHeader & bamHeader);
    ~TSiteSubsetChr();

    void findWindow(const long & pos);
    void findOrCreateWindow(const long & pos);
    void addPosition(std::vector<std::string> & tmp, const std::string & chr, bool invariantSites);
    bool addPosition(std::vector<std::string> & tmp, const std::string & chr, Fasta & reference, std::string & error, bool invariantSites);
    void print();
    bool hasPositionsInWindow(const long & windowStart);
    std::map<long,std::pair<char,char> >& getPositionInWindow(const long & windowStart);
    long size();
};

class TSiteSubset{
private:
	std::map<std::string, TSiteSubsetChr*> chromosomes;
	std::map<std::string, TSiteSubsetChr*>::iterator chrIt;
	int windowSize;
	std::string curChr;
	bool invariantSites;

    void readFile(TLog* logfile);
    void readFile(Fasta & reference, SamHeader & bamHeader, TLog* logfile);

public:
	std::string filename;

    TSiteSubset(std::string Filename, int & WindowSize, TLog* logfile, bool InvariantSites);
    TSiteSubset(std::string Filename, Fasta & reference, SamHeader & bamHeader, int & WindowSize, TLog* logfile, bool InvariantSites);
    ~TSiteSubset();

    void setChr(const std::string & chr);
    void print();
    bool hasPositionsInWindow(const long & windowStart);
    std::map<long,std::pair<char,char> >& getPositionInWindow(long & windowStart);
    long size();
};


#endif /* TSITESUBSET_H_ */
