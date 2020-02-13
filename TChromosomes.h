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
	std::vector<uint32_t> lengths;
	std::vector<uint8_t> ploidies;
	std::vector<bool> isInUse;
	std::map<std::string, uint16_t> nameMap;

	uint16_t curChrNumber;
	BamTools::SamSequenceIterator curChrIterator;

public:
	TChromosomes();
	TChromosomes(BamTools::SamHeader* BamHeader);
	int limitChr(std::string & limitName);
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
	uint16_t size();
	uint16_t curIndex();
	uint32_t referenceLength();
	uint16_t getIndexFromName(const std::string chrName);

	uint32_t length(const uint16_t index);
	uint32_t curLength();
	std::string name(const uint16_t index);
	std::string curName();
	bool inUse(const uint16_t index);
	bool curInUse();
	uint8_t ploidy(const uint16_t index);
	uint8_t curPloidy();
};



#endif /* TCHROMOSOMES_H_ */
