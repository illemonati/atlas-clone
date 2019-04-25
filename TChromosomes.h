/*
 * TChromosomes.h
 *
 *  Created on: Apr 25, 2019
 *      Author: linkv
 */

#ifndef TCHROMOSOMES_H_
#define TCHROMOSOMES_H_

#include <string>
#include <vector>
#include "bamtools/api/BamReader.h"
#include "stringFunctions.h"
#include "TLog.h"


class TChromosomes{
private:
	BamTools::SamHeader* bamHeader;
	int numChromosomes;
	std::vector<std::string> names;
	std::vector<long> lengths;
	std::vector<bool> diploid;
	std::vector<bool> inUse;
	std::map<std::string, int> nameMap;

	int curChrNumber;
	BamTools::SamSequenceIterator curChrIterator;

public:
	TChromosomes();
	TChromosomes(BamTools::SamHeader* BamHeader);
	void limitChr(std::string & limitName);
	void useSpecifiedChr(std::vector<std::string> & chrNames, TLog* logfile);

	//move
	void begin();
	void next();
	bool end();
	void jumpToEnd();
	void jumpToBeginningOfLastChr();

	//getters
	int getNumberOfChr();
	int getCurIndex();
	long calcReferenceLength();
	long getCurLength();
	std::string getCurName();
	std::string getNameAtIndex(int & index);
	int getIndexFromName(const std::string & chrName);
	bool chrInUse(int & index);
	bool curChrInUse();
};



#endif /* TCHROMOSOMES_H_ */
