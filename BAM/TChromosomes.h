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

namespace BAM{

//---------------------------------------------------------
// TChromosome
//---------------------------------------------------------
class TChromosome{
public:
	uint32_t refID;
	std::string name;
	uint32_t length;
	uint8_t ploidy;
	bool inUse;

	TChromosome(const uint32_t RefID, const std::string Name, const uint32_t Length){
		refID = RefID;
		name = Name;
		length = Length;
		ploidy = 2; //default: diploid
		inUse = true;
	};
};

//---------------------------------------------------------
// TChromosomes
//---------------------------------------------------------

class TChromosomes{
private:
	std::vector<TChromosome> _chromosomes;
	std::vector<TChromosome>::iterator _curChr;

	TChromosome& _find(const std::string chrName) const;

public:
	TChromosomes(){};
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
	uint32_t size() const;
	uint32_t referenceLength() const;
	bool exists(const std::string name) const;

	TChromosome& getChromosome(const std::string chrName);
	TChromosome& curChromosome();

	uint32_t refID(const std::string chrName) const;
	uint32_t curRefID() const;

	uint32_t length(const uint16_t index) const;
	uint32_t curLength() const;

	std::string name(const uint16_t index) const;
	std::string curName() const;

	bool inUse(const uint16_t index) const;
	bool curInUse() const;

	uint8_t ploidy(const uint16_t index) const;
	uint8_t curPloidy() const;
};

}; //end namespace

#endif /* TCHROMOSOMES_H_ */
