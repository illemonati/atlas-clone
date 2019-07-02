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
#include <iostream>
#include "stringFunctions.h"
#include "gzstream.h"
#include "TLog.h"
#include "IOTools/IOAbstractClasses/SamHeader.h"

//read sorted bed files window by window
//store all data in chr / window combinations using vectors
//Store all positions 0-based, as in TWindow

class TBedReaderWindow{
public:
    bool hasData;
    long start, end;
    std::vector<long> positions;

    TBedReaderWindow(long Start, long End);
    ~TBedReaderWindow();

    void addPosition(long & pos);
    void print();
    long size();
};

class TBedReaderChromosome{
public:
    std::string name;
    std::map<int, TBedReaderWindow*> windows;
    std::map<int, TBedReaderWindow*>::iterator windowIt;
    int windowSize;


    TBedReaderChromosome(std::string & Name, int & WindowSize);
    ~TBedReaderChromosome();

    void findWindow(const long & pos);
    void findOrCreateWindow(const long & pos);
    void addPosition(std::vector<std::string> & tmp);
    void print();
    bool hasPositionsInWindow(const long & windowStart);
    std::vector<long>& getPositionInWindow(long windowStart);
    long size();
};

class TBedReader{
private:
    std::map<std::string, TBedReaderChromosome*> chromosomes;
    std::map<std::string, TBedReaderChromosome*>::iterator chrIt;
    int windowSize;
    std::string curChr;

    void readFile(SamSequenceDictionary & Sequences, TLog* logfile);

public:
    std::string filename;

    TBedReader(std::string Filename, int & WindowSize, SamSequenceDictionary & Sequences, TLog* logfile);
    ~TBedReader();

    void setChr(const std::string & chr);
    void print();
    bool hasPositionsInWindow(const long & windowStart);
    std::vector<long>& getPositionInWindow(long & windowStart);
    long size();
    int getNumChromosomes();
};


#endif /* TBEDREADER_H_ */
