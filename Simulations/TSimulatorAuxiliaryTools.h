/*
 * TSimulatorAuxiliaryTools.h
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATORAUXILIARYTOOLS_H_
#define TSIMULATORAUXILIARYTOOLS_H_

#include "TLog.h"
#include "TRandomGenerator.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#include "gzstream.h"
#include "TBamFile.h"
#include "globalConstants.h"

namespace Simulations{

//---------------------------------------------------------
//TSimulatorChromosome
//---------------------------------------------------------
class TSimulatorChromosome{
public:
	std::string name;
	uint64_t length;
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
	bool needsWriting;

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
	void close();

	void setChr(std::string ChrName, long ChrLength);
//	void simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float* cumulBaseFreq);
	Base* getPointerToRef(){ return ref; };
};

//---------------------------------------------------------
//TSimulatorBamFile
//---------------------------------------------------------
class TSimulatorBamFile{
private:
	BAM::TOutputBamFile _outBam;

public:
	TSimulatorBamFile(){};
	TSimulatorBamFile(const std::string Filename, const std::string SampleName, const std::vector<std::string> & ReadGroupNames, const std::vector<TSimulatorChromosome> & Chromosomes, TLog* Logfile, TGenotypeMap & GenoMap, TQualityMap & QualMap){
		open(Filename, SampleName, ReadGroupNames, Chromosomes, Logfile, GenoMap, QualMap);
	};
	~TSimulatorBamFile();

	void open(const std::string Filename, const std::string SampleName, const std::vector<std::string> & ReadGroupNames, const std::vector<TSimulatorChromosome> & Chromosomes, TLog* Logfile, TGenotypeMap & GenoMap, TQualityMap & QualMap);

	void saveAlignment(const BAM::TAlignment & Alignment){
		_outBam.writeAlignment(Alignment);
	};
	void close(TLog* Logfile);
	void indexBamFile();
};

class TSimulatorBamFiles{
private:
	std::vector<TSimulatorBamFile> _files;
	TLog* _logfile;

public:
	TSimulatorBamFiles(uint32_t NumFiles, const std::string Outname, const std::vector<std::string> & ReadGroupNames, const std::vector<TSimulatorChromosome> & Chromosomes, TLog* Logfile, TGenotypeMap & GenoMap, TQualityMap & QualMap);
	~TSimulatorBamFiles();

	void close();
	TSimulatorBamFile& operator[](size_t i);
};

//---------------------------------------------------------
//TSimulatorAlleleIndex
//---------------------------------------------------------
class TSimulatorAlleleIndex{
private:
	int nextIndex;
	Base refBase;
	Base indexToBase[4];

public:
	int index[4];
	bool used[4];

	TSimulatorAlleleIndex(){
		nextIndex = 0;
		refBase = N;
	};

	void clear(const Base & ref){
		used[0] = false; used[1] = false; used[2] = false; used[3] = false;
		used[ref] = true;
		index[ref] = 0;
		nextIndex = 1;
		refBase = ref;
	};

	void add(const Base & base){
		if(!used[base]){
			used[base] = true;
			index[base] = nextIndex;
			indexToBase[nextIndex] = base;
			++nextIndex;
		}
	};

	void writeRefAltToVCF(gz::ogzstream & VCF, TGenotypeMap & genoMap){
		VCF << genoMap.baseToChar[refBase] << '\t';
		if(nextIndex == 1) //no alt
			VCF << ".";
		else {
			VCF << genoMap.baseToChar[indexToBase[1]];
			for(int i=2; i<nextIndex; ++i)
				VCF << ',' << genoMap.baseToChar[indexToBase[i]];
		}
	};
};

//---------------------------------------------------------
//TSimulatorHaplotypes
//---------------------------------------------------------
class TSimulatorHaplotypes{
private:
	int numInd;
	uint64_t curLength;
	uint64_t storageLength;
	bool initialized;
	Base*** haplotypes;

	//write true genotypes to VCF
	gz::ogzstream trueGenoVCF;
	bool trueGenoVCFOpend;


	void allocateStorage(uint64_t length);
	void freeStorage();

public:
	TSimulatorHaplotypes(int NumIndividuals);
	~TSimulatorHaplotypes(){
		freeStorage();
	};

	void setLength(uint64_t length);
	void openTrueGenotypeVCF(std::string filename);
	void closeTrueGenotypeVCF();
	Base** getHaplotypesOfIndividual(int i);
	Base** getHaplotypesFirstIndividual(){
		return haplotypes[0];
	};
	void writeTrueGenotypes(TSimulatorChromosome & chromosome, Base* ref, TGenotypeMap & genoMap);
	int size(){ return numInd; };
	Base& operator()(int ind, int hap, uint64_t site){
		return haplotypes[ind][hap][site];
	};
	bool isPolymoprhic(uint64_t pos);
};

//---------------------------------------------------------
//TSimulatorMutationtable
//---------------------------------------------------------
class TSimulatorMutationtable{
private:
	float** mutTable;
	bool tableAllocated;

	void allocateTable();
	void deleteTable();
	void normalizeAndMakeCumulative();

public:
	TSimulatorMutationtable();
	TSimulatorMutationtable(float* baseFreq);
	TSimulatorMutationtable(float* baseFreq, const double theta);
	~TSimulatorMutationtable(){ deleteTable(); };
	void fill(float* baseFreq);
	void fill(float* baseFreq, double theta);
	float* operator[](int i){ return mutTable[i]; };
	void print();
};

//---------------------------------------------------------
//TSimulatorVariantInvariantBedFiles
//---------------------------------------------------------
class TSimulatorVariantInvariantBedFiles{
private:
	gz::ogzstream variantSitesFile;
	gz::ogzstream invariantSitesFile;
	bool filesOpend;

	void openFile(gz::ogzstream & file, const std::string filename);
	void close();

public:
	TSimulatorVariantInvariantBedFiles();
	TSimulatorVariantInvariantBedFiles(std::string outname);
	~TSimulatorVariantInvariantBedFiles();

	void open(std::string outname);

	void write(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome);
};

}; //end namespace

#endif /* TSIMULATORAUXILIARYTOOLS_H_ */
