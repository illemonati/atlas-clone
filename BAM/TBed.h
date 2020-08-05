/*
 * TBed.h
 *
 *  Created on: Mar 16, 2017
 *      Author: phaentu
 */

#ifndef TBED_H_
#define TBED_H_

#include <set>
#include "TFile.h"
#include "TGenomePosition.h"
#include "TChromosomes.h"

namespace BAM{


//-------------------------------------------------------------
// TBedChromosome
// a class to store the chromosomes parsed
// also used to write
//-------------------------------------------------------------
class TBedChromosome{
public:
	uint32_t refId;
	std::string name;

	TBedChromosome(const uint32_t RefId, const std::string Name){
		refId = RefId;
		name = Name;
	};
	bool operator<(const TBedChromosome & Other) const{
		return name < Other.name;
	};
	bool operator<(const std::string & Name) const{
		return name < Name;
	};
};

bool operator<(const std::string & Name, const TBedChromosome & Chr);

//-------------------------------------------------------------
// TBed_base
// base class to store a list of windows
//-------------------------------------------------------------
class TBed_base{
protected:
	std::set<TBedChromosome, std::less<>> _chromosomes; //ordered by name
	std::map<uint32_t, std::string> _chromosomeNames; //ordered by ID

	template <typename T> uint32_t _numChromosomesWithWindows(const T List) const{
		if(List.empty()){
			return 0;
		}

		uint32_t c = 1;
		auto it = List.begin();
		uint32_t curChr = it->refID();

		for(; it!=List.end(); ++it){
			if(it->refID() > curChr){
				++c;
				curChr = it->refID();
			}
		}
		return c;
	};

public:

	std::set<TBedChromosome, std::less<>>::iterator addChromosome(const std::string & Chr);
	void addChromosomes(const TChromosomes & Chromosomes);
	virtual uint32_t numChromosomesWithWindows() const;
	bool hasWindowsOnChr(const std::string& Chr);
	bool hasWindowsOnChr(const uint32_t refId);
	uint32_t getRefID(const std::string& Chr);
	std::string getChromosomeName(const uint32_t refId);
};


//-------------------------------------------------------------
// TBed
// a list of non-overlapping windows
//-------------------------------------------------------------
class TBed:public TBed_base{
private:
	std::set<TGenomeWindow, std::less<>> _bed;

public:
	TBed()= default;
	TBed(const std::string Filename){ add(Filename); };
	TBed(const std::string Filename, const TChromosomes & Chromosomes){ add(Filename, Chromosomes); };

	void add(TGenomeWindow Window);
	void add(const TGenomePosition & Position);
	void add(const uint32_t Chr, const uint32_t Pos);
	void add(const std::string& Chr, const uint32_t Pos);

	void add(const std::string& Filename);
	void add(const std::string& Filename, const TChromosomes & Chromosomes);
	void write(const std::string& Filename) const;

	uint64_t size() const;
	uint64_t length() const;
	uint32_t numChromosomesWithWindows() const;

	bool exists(const TGenomeWindow& Window) const;

	//loop
	std::set<TGenomeWindow>::iterator lower_bound(const TGenomeWindow & Window){ return _bed.lower_bound(Window); };
	std::set<TGenomeWindow>::iterator begin(){ return _bed.begin(); };
	std::set<TGenomeWindow>::iterator end(){ return _bed.end(); };

	std::set<TGenomeWindow>::iterator begin(const uint32_t refId){ return _bed.lower_bound(TGenomePosition(refId, 0)); };
	std::set<TGenomeWindow>::iterator end(const uint32_t refId){ return _bed.lower_bound(TGenomePosition(refId+1, 0)); };
};


//-----------------------------------------------------
// TGenomeWindowList
// a list of potentially overlapping windows
//-----------------------------------------------------
class TGenomeWindowList:public TBed_base{
private:
	std::multiset<TGenomeWindow> _list;

public:
	TGenomeWindowList()= default;

	void add(const TGenomeWindow& Window);
	void add(const std::string Filename,  const TChromosomes & Chromosomes);
	void write(const std::string Filename, const TChromosomes & Chromosomes) const;

	bool empty() const{ return _list.empty(); };
	uint64_t size() const;
	uint64_t length() const;
	uint32_t numChromosomesWithWindows() const;

	bool exists(const TGenomeWindow& Window) const;
	bool hasWindowsOnChr(uint32_t refId) const;
	uint32_t numWindowsOnChr(const uint32_t refId) const;

	//loop
	std::multiset<TGenomeWindow>::iterator begin(){ return _list.begin(); };
	std::multiset<TGenomeWindow>::iterator end(){ return _list.end(); };
};



/*


//------------------ OLD ----------------------------
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
*/

}; //end namespace

#endif /* TBED_H_ */
