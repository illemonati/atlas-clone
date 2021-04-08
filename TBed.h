/*
 * TBed.h
 *
 *  Created on: Mar 16, 2017
 *      Author: phaentu
 */

#ifndef TBED_H_
#define TBED_H_

#include <map>
#include "stringFunctions.h"
#include "gzstream.h"
#include "TFile.h"

typedef std::map<uint64_t, uint64_t> TBedWindowMap;

class TBedChromosome{
private:
	TBedWindowMap _windows;
	TBedWindowMap::iterator _windowIt;
	bool _parsed;

	void _fuseNeighboringWindowsIfNeeded(TBedWindowMap::iterator & it);

public:
	TBedChromosome(){
		rewind();
	};

	~TBedChromosome(){
		//delete all windows
		_windows.clear();
	};

	void addWindow(uint64_t start, uint64_t end);
	void addOrExtendWindow(const uint64_t start, const uint64_t end);
	void print(const std::string & chrName);
	size_t size(){ return _windows.size(); };
	uint64_t length();
	void setAsParsed(){ _parsed = true; };
	void rewind(){ _parsed = false; };
	bool parsed(){ return _parsed; };
	bool begin();
	bool next();
	bool reachedEnd();
	uint64_t curStart(){ return _windowIt->first; };
	uint64_t curEnd(){ return _windowIt->second; };
	uint64_t curLength(){ return _windowIt->second - _windowIt->first; };
	void write(TOutputFile & out, const std::string & chrName);
};

typedef std::map<std::string, TBedChromosome> TBedChrMap;

class TBed{
private:
	TBedChrMap _chromosomes;
	TBedChrMap::iterator _chrIt;
	std::string _curChr;
	bool _allChrParsed;

	void readFile(const std::string filename);
	void _setCurChrAsParsed();

public:
	TBed(){
		_curChr = "";
		rewind();
	};

	TBed(const std::string filename){
		readFile(filename);
		_curChr = "";
		rewind();
	};

	~TBed(){
		clear();
	};

	void addWindowCurChr(const uint64_t start, const uint64_t end);
	void addSite(const uint64_t pos);
	void rewind();
	void clear();
	bool setChr(const std::string & chr);
	void setChrCreateIfNew(const std::string & chr);
	std::string curChr();
	bool nextWindow();
	bool reachedEnd(){ return _allChrParsed; };
	bool reachedEndOfChr();
	uint64_t curWindowStart();
	uint64_t curWindowEnd();
	uint64_t curWindowSize();
	void print();
	size_t size();
	uint64_t length();
	int getNumChromosomes(){ return _chromosomes.size(); };
	int getNumWindowsOnCurChr();
	bool test();
	void write(const std::string filename);
};



#endif /* TBED_H_ */
