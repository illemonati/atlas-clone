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


class TChromosome{
public:
	std::string name;
	uint32_t length;
	uint8_t ploidy;
	bool inUse;

	TChromosome(const std::string Name, const uint32_t Length){
		name = Name;
		length = Length;
		ploidy = 2; //default: diploid
		inUse = true;
	};
};

class TChromosomes{
private:
	std::vector<TChromosome> chromosomes;

	std::vector<TChromosome>::iterator curChr;
	//uint16_t curChrNumber;
	//uint16_t

	TChromosome& _find(const std::string chrName);

public:
	TChromosomes();
	TChromosomes(BamTools::SamHeader* BamHeader);
	void readChromosomes(BamTools::SamHeader* BamHeader);
	void limitChr(const std::string limitName);
	void useSpecifiedChr(const std::vector<std::string> & chrNames);
	void specifyPloidy(std::ifstream & ploidyFile, TLog* logfile);
	void setToHaploid(std::vector<std::string> chrNames, TLog* logfile);
	void writeUsedChromosomes(TLog* logfile);

	//move
	void begin();
	bool next();
	bool end() const;
	void jumpToEnd();
	void jumpToBeginningOfLastChr();

	//getters
	uint16_t size() const;
	uint16_t curIndex() const;
	uint32_t referenceLength() const;
	uint16_t getIndexFromName(const std::string chrName) const;

	uint32_t length(const uint16_t index) const;
	uint32_t curLength() const;
	std::string name(const uint16_t index) const;
	std::string curName() const;
	bool exists(const std::string name) const;
	bool inUse(const uint16_t index) const;
	bool curInUse() const;
	uint8_t ploidy(const uint16_t index) const;
	uint8_t curPloidy() const;
};



#endif /* TCHROMOSOMES_H_ */
