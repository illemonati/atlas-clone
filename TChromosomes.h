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
	std::vector<int> ploidies;
	std::vector<bool> isInUse;
	std::map<std::string, int> nameMap;

	int curChrNumber;
	BamTools::SamSequenceIterator curChrIterator;

public:
	TChromosomes();
	TChromosomes(BamTools::SamHeader* BamHeader);
	void limitChr(std::string & limitName);
	void useSpecifiedChr(std::vector<std::string> & chrNames, TLog* logfile);
	void specifyPloidy(std::ifstream & ploidyFile, TLog* logfile);
	void setToHaploid(std::vector<std::string> chrNames, TLog* logfile);

	//move
	void begin();
	void next();
	bool end();
	void jumpToEnd();
	void jumpToBeginningOfLastChr();

	//getters
	int size();
	int curIndex();
	long referenceLength();
	int getIndexFromName(const std::string chrName);

	long length(const int index);
	long curLength();
	std::string name(const int index);
	std::string curName();
	bool inUse(const int index);
	bool curInUse();
	int ploidy(const int index);
	int curPloidy();
};



#endif /* TCHROMOSOMES_H_ */
