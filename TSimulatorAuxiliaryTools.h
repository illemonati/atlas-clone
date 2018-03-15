/*
 * TSimulatorAuxiliaryTools.h
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATORAUXILIARYTOOLS_H_
#define TSIMULATORAUXILIARYTOOLS_H_

#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamHeader.h"
#include "bamtools/api/BamAlignment.h"
#include "TLog.h"
#include "TRandomGenerator.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#include "gzstream.h"

//---------------------------------------------------------
//TSimulatorChromosome
//---------------------------------------------------------
class TSimulatorChromosome{
public:
	std::string name;
	long length;
	bool haploid;
	int refID;

	TSimulatorChromosome(std::string Name, int RefID, long Length, bool Haploid){
		name = Name;
		refID = RefID;
		length = Length;
		haploid = Haploid;
	};
};

//---------------------------------------------------------
//TSimulatorReference
//---------------------------------------------------------
class TSimulatorReference{
private:
	TLog* logfile;

	//fasta file
	std::ofstream fasta;
	std::ofstream fastaIndex;
	long oldOffset;
	bool fastaOpen;
	std::string filename;
	TGenotypeMap genoMap;

	//reference storage
	Base* ref;
	bool storageInitialized;
	long storageLength;
	long chrLength;
	std::string chrName;

	void allocateStorage(long length);
	void freeStorage();

	void openFastaFile();
	void closeFastaFile();
	void writeRefToFasta();

public:
	TSimulatorReference();
	TSimulatorReference(std::string Filename, TLog* Logfile);
	~TSimulatorReference(){
		close();
	};

	void initialize(std::string Filename, TLog* Logfile);
	void close(){
		if(chrName != "")
			writeRefToFasta();
		closeFastaFile();
		freeStorage();
	};

	void setChr(std::string ChrName, long ChrLength);
	void simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float* cumulBaseFreq);
	Base* getPointerToRef(){ return ref; };
};


//---------------------------------------------------------
//TSimulatorBamFile
//---------------------------------------------------------
class TSimulatorBamFile{
private:
	bool isOpen;
	std::string filename;
	BamTools::RefVector references;
	BamTools::BamWriter bamWriter;
	TLog* logfile;

	void init(TLog* Logfile){
		logfile = Logfile;
		filename = "";
		isOpen = false;
	};

public:
	TSimulatorBamFile(std::string Filename, std::vector<std::string> & readGroupNames, std::vector<TSimulatorChromosome> & chromosomes, TLog* Logfile){
		init(Logfile);
		open(Filename, readGroupNames, chromosomes);
	};
	~TSimulatorBamFile(){
		close();
	};

	void open(std::string Filename, std::vector<std::string> & readGroupNames, std::vector<TSimulatorChromosome> & chromosomes);
	bool saveAlignment(const BamTools::BamAlignment & bamAlignment){
		return bamWriter.SaveAlignment(bamAlignment);
	};
	void close();
	void indexBamFile();
};

//---------------------------------------------------------
//TSimulatorHaplotypes
//---------------------------------------------------------
class TSimulatorHaplotypes{
private:
	int numInd;
	long curLength;
	long storageLength;
	bool initialized;
	int ind;
	Base*** haplotypes;

	void allocateStorage(long length);
	void freeStorage();

public:
	TSimulatorHaplotypes(int NumIndividuals);
	~TSimulatorHaplotypes(){
		freeStorage();
	};

	void setLength(long length);
	Base** getHaplotypesOfIndividual(int i);
	Base** getHaplotypesFirstIndividual(){
		return haplotypes[0];
	};
	void writeGenotypes(gz::ogzstream & out, std::string & chrName, TGenotypeMap & genoMap);
	int size(){ return numInd; };
	Base& operator()(int ind, int hap, long site){
		return haplotypes[ind][hap][site];
	};
};







#endif /* TSIMULATORAUXILIARYTOOLS_H_ */
