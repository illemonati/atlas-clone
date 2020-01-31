/*
 * TBedReader.h
 *
 *  Created on: Oct 6, 2015
 *      Author: wegmannd
 */

#ifndef TBEDREADER_H_
#define TBEDREADER_H_

#include <fstream>
#include <vector>
#include <map>
#include "bamtools/api/SamSequenceDictionary.h"
#include "stringFunctions.h"
#include "TLog.h"
#include "gzstream.h"


//read sorted bed files window by window
//store all data in chr / window combinations using vectors
//Store all positions 0-based, as in TWindow

class TBedReaderWindow{
public:
	bool hasData;
	unsigned int start, end;
	std::vector<unsigned int> positions;

	TBedReaderWindow(unsigned int Start, unsigned int End);
	~TBedReaderWindow();
	void addPosition(unsigned int & pos);
	void print();
	unsigned int size();
};

class TBedReaderChromosome{
public:
	std::string name;
	std::map<unsigned int, TBedReaderWindow*> windows;
	std::map<unsigned int, TBedReaderWindow*>::iterator windowIt;
	int windowSize;

	TBedReaderChromosome(std::string & Name, unsigned int & WindowSize);
	~TBedReaderChromosome();
	void findWindow(const unsigned int & pos);
	void findOrCreateWindow(const unsigned int & pos);
	void addPosition(std::vector<std::string> & tmp, unsigned int & numPositionsAdded);
	void print();
	bool hasPositionsInWindow(const unsigned int & windowStart);
	std::vector<unsigned int>& getPositionInWindow(unsigned int windowStart);
	unsigned int size();
};

class TBedReader{
private:
	std::map<std::string, TBedReaderChromosome*> chromosomes;
	std::map<std::string, TBedReaderChromosome*>::iterator chrIt;
	unsigned int windowSize;
	std::string curChr;
	unsigned int numPositionsAdded;

	void readFile(BamTools::SamSequenceDictionary & Sequences, unsigned int siteLimit, TLog* logfile);

public:
	std::string filename;
	TBedReader(std::string Filename, const unsigned int & WindowSize, BamTools::SamSequenceDictionary & Sequences, unsigned int siteLimit, TLog* logfile);
	~TBedReader();
	void setChr(const std::string & chr);
	void print();
	bool hasPositionsInWindow(const unsigned int & windowStart);
	std::vector<unsigned int>& getPositionInWindow(unsigned int & windowStart);
	unsigned int size();
	unsigned int getNumChromosomes();

};


#endif /* TBEDREADER_H_ */
