/*
 * TBed.h
 *
 *  Created on: Mar 16, 2017
 *      Author: phaentu
 */

#ifndef TBED_H_
#define TBED_H_

#include <string>
#include <map>
#include <stringFunctions.h>
#include "gzstream.h"

class TBedChromosome{
public:
    std::string name;
    std::map<long, long> windows;
    std::map<long, long>::iterator windowIt;

    TBedChromosome(std::string & Name);
    ~TBedChromosome();

    void addWindow(long start, long end);
    void print();
    long size();
    bool begin();
    bool next();
    long curStart();
    long curEnd();
};

class TBed{
private:
    std::map<std::string, TBedChromosome*> chromosomes;
    std::map<std::string, TBedChromosome*>::iterator chrIt;
    std::string curChr;

    void readFile();

public:
    std::string filename;

    TBed(std::string Filename);
    ~TBed();

    void setChr(const std::string & chr);
    bool nextWindow();
    long curWindowStart();
    long curWindowEnd();
    void print();
    long size();
    int getNumChromosomes();
    int getNumWindowsOnCurChr();
};



#endif /* TBED_H_ */
