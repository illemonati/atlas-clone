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
#include "stringFunctions.h"
#include "TParameters.h"
#include "TLog.h"
#include "TFile.h"
#include "TGenomePosition.h"

namespace BAM{


//---------------------------------------------------------
// TChromosome
//---------------------------------------------------------
class TChromosome{
private:
	void _initialize(const uint32_t & RefID, const std::string & Name, const uint32_t & Length, const uint8_t & Ploidy);

public:
	std::string name; //SN field
	uint32_t length; //LS field
	uint8_t ploidy;
	TGenomePosition chrStart, chrEnd;
	bool inUse;

	//the following fields are not parsed.
	std::string alternateLocus; //AH field
	std::string alternativeNames; //AN field
	std::string assemblyIdentifier; //AS field
	std::string description; //DS field
	std::string MD5; //M5 field
	std::string species; //SP field
	std::string topology; //TP field;
	std::string uri; //UR field;

	TChromosome(const uint32_t & RefID, const std::string & Name, const uint32_t & Length);
	TChromosome(const uint32_t & RefID, const std::string & Name, const uint32_t & Length, const uint8_t & Ploidy);


	uint32_t refID() const { return chrStart.refID(); };

	bool operator<(const TChromosome & other){
		return chrStart < other.chrStart;
	};

	std::string compileSamHeader() const;
};

inline bool operator<(const TGenomePosition & pos, const TChromosome & chr){
	return pos < chr.chrStart;
};
inline bool operator<(const TChromosome & chr, const TGenomePosition & pos){
	return chr.chrEnd < pos;
};
inline bool operator<(const TChromosome & chr, const TGenomeWindow & win){
	return chr.chrEnd < win;
};
inline bool operator<(const TGenomeWindow & win, const TChromosome & chr){
	return win < chr.chrStart;
};

//---------------------------------------------------------
// TChromosomes
//---------------------------------------------------------

class TChromosomes{
private:
	std::vector<TChromosome> _chromosomes;
	mutable std::vector<TChromosome>::iterator _curChr;

	TChromosome& _find(const std::string chrName);
    const TChromosome& _find(const std::string chrName) const;
    void _setAllNotInUse();
	void _useSpecifiedChr(const std::vector<std::string> & chrNames);
	void _useUpToAndIncluding(const std::string & name);
	void _specifyPloidy(const std::string & ploidyFileName, TLog* logfile);
	void _setToHaploid(const std::vector<std::string> & chrNames, TLog* logfile);

public:
	TChromosomes(){};

	void clear();
	void appendChromosome(const std::string & name, const uint32_t & length);
	void appendChromosome(const std::string & name, const uint32_t & length, const uint8_t & ploidy);

	void limitAndSetPloidy(TParameters & params, TLog* logfile);
	void limitChr(TParameters & params, TLog* logfile);
	void setPloidy(TParameters & params, TLog* logfile);

	void limitChr(const std::string limitName);
	void useSpecifiedChr(const std::vector<std::string> & chrNames);
	void writeUsedChromosomes(TLog* logfile);

	//iterate
	std::vector<TChromosome>::iterator begin(){ return _chromosomes.begin(); };
	std::vector<TChromosome>::iterator begin(const uint32_t & RefID){ return _chromosomes.begin() + RefID; };
	std::vector<TChromosome>::iterator end(){ return _chromosomes.end(); };
	std::vector<TChromosome>::const_iterator cbegin() const{ return _chromosomes.cbegin(); };
	std::vector<TChromosome>::const_iterator cbegin(const uint32_t & RefID) const{ return _chromosomes.begin() + RefID; };
	std::vector<TChromosome>::const_iterator cend() const{ return _chromosomes.cend(); };

	//getters
	uint32_t size() const { return _chromosomes.size(); };
	uint32_t referenceLength() const;
	uint32_t minLength() const;
	uint32_t maxLength() const;

	//getters by name
	bool exists(const std::string name) const;
	const TChromosome& getChromosome(const std::string chrName) const;
	uint32_t refID(const std::string chrName) const;

	//getters by refID
	bool exists(const uint32_t & RefID) const{ return RefID < size(); };
	TChromosome& operator[](const uint32_t & RefID){ return _chromosomes[RefID]; };
	const TChromosome& at(const uint32_t & RefID){ return _chromosomes[RefID]; };
	uint32_t length(const uint32_t & RefID) const;
	std::string name(const uint32_t & RefID) const;
	bool inUse(const uint32_t & RefID) const;
	uint8_t ploidy(const uint32_t & RefID) const;
	const TGenomePosition& chrStart(const uint32_t & RefID) const;
	const TGenomePosition& chrEnd(const uint32_t & RefID) const;

	std::string compileSamHeader() const;
};

}; //end namespace

#endif /* TCHROMOSOMES_H_ */
