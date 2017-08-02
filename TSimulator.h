/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include "TLog.h"
#include "TRandomGenerator.h"
#include "SFS.h"
#include "TPostMortemDamage.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamHeader.h"
#include "bamtools/api/BamAlignment.h"
#include <math.h>

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
	char* toBase;

	//reference storage
	short* ref;
	bool storageInitialized;
	long storageLength;
	long chrLength;
	std::string chrName;

	void allocateStorage(long length){
		freeStorage();

		//allocate storage
		ref = new short[length];
		storageInitialized = true;
		storageLength = length;
	};

	void freeStorage(){
		if(storageInitialized){
			delete[] ref;
		}
	};

	void openFastaFile();
	void closeFastaFile();
	void writeRefToFasta();

public:
	TSimulatorReference(std::string Filename, char* ToBase, TLog* Logfile){
		filename = Filename;
		logfile = Logfile;
		toBase = ToBase;
		chrName = "";
		ref = NULL;
		storageInitialized = false;
		storageLength = 0;

		openFastaFile();
	};
	~TSimulatorReference(){
		if(chrName != "")
			writeRefToFasta();
		closeFastaFile();
		freeStorage();
	};

	void setChr(std::string ChrName, long ChrLength){
		//write if not yet written
		if(chrName != "")
			writeRefToFasta();

		//move to new chr
		chrName = ChrName;
		if(ChrLength > storageLength)
			allocateStorage(ChrLength);
		chrLength = ChrLength;
	};

	void simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float* cumulBaseFreq);
	short* getPointerToRef(){ return ref; };
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
	TSimulatorBamFile(std::string Filename, const std::string & readGroupName, std::vector<TSimulatorChromosome> & chromosomes, TLog* Logfile){
		init(Logfile);
		open(Filename, readGroupName, chromosomes);
	};

	void open(std::string Filename, const std::string & readGroupName, std::vector<TSimulatorChromosome> & chromosomes);
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
	short*** haplotypes;

	void allocateStorage(long length){
		freeStorage();
		//allocate storage



		haplotypes = new short**[numInd];
		for(ind=0; ind<numInd; ++ind){
			haplotypes[ind] = new short*[2];
			haplotypes[ind][0] = new short[length];
			haplotypes[ind][1] = new short[length];
		}
		initialized = true;
		storageLength = length;
	};

	void freeStorage(){
		if(initialized){
			for(ind=0; ind<numInd; ++ind){
				delete[] haplotypes[ind][0];
				delete[] haplotypes[ind][1];
				delete[] haplotypes[ind];
			}
			delete[] haplotypes;
		}
	};

public:
	TSimulatorHaplotypes(int NumIndividuals){
		numInd = NumIndividuals;
		ind = 0;
		haplotypes = NULL;
		initialized = false;
		curLength = 0;
		storageLength = 0;
	};

	~TSimulatorHaplotypes(){
		freeStorage();
	};

	void setLength(long length){
		if(length > storageLength){
			allocateStorage(length);
		}
		curLength = length;
	};

	short** getHaplotypesOfIndividual(int i){
		if(i >= numInd)
			throw "Haplotypes of individual " + toString(i+1) + " requested, but defined for only " + toString(numInd) + " individuals!";
		return haplotypes[i];
	};

	short** getHaplotypesFirstIndividual(){
		return haplotypes[0];
	};

	void writeGenotypes(std::ofstream & out, std::string & chrName, char* toBase){
		for(int l=0; l<curLength; ++l){
			out << chrName << "\t" << l+1;
			for(ind=0; ind < numInd; ++ind)
				out << "\t" << toBase[haplotypes[ind][0][l]] << "/" << toBase[haplotypes[ind][1][l]];
			out << "\n";
		}
	};
};

//---------------------------------------------------------
//TSimulatorGenotypecombination
//---------------------------------------------------------
class TSimulatorGenotypeCombination{
private:
	double cumulGenoCaseFrequencies[9];
	int numGenotypeCombinations[9];
	double** cumulGenoCombinationFreq;
	short*** genoTrans;
	short** orderLookup;
	bool tablesInitialized;;

	void fillTables(std::vector<double> & phis, std::vector<float> & baseFreq);
	void deleteTables();

public:
	TSimulatorGenotypeCombination(std::vector<double> & phis, std::vector<float> & baseFreq){
		tablesInitialized = false;
		fillTables(phis, baseFreq);
	};
	~TSimulatorGenotypeCombination(){
		deleteTables();
	}

	void simulateHaplotypes(short** haplotypesInd0, short** haplotypesInd1, short* ref, float referenceDivergence, long length, TRandomGenerator* randomGenerator);

};

//---------------------------------------------------------
//TSimulator
//---------------------------------------------------------
class TSimulator{
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	bool bamFileOpen;

	//general simulation parameters
	double meanQual, sdQual;
	int maxQual;
	float seqDepth;
	int readLength;
	std::vector<TSimulatorChromosome> chromosomes;
	std::vector<TSimulatorChromosome>::iterator chrIt;
	std::string readGroupName;

	//Qual to error table
	double* qualToErroTable;
	bool qualToErroTableInitialized;

	//PMD
	TPMD* pmdObject;
	bool pmdInitialized;

	//Quality transformation
	TGenotypeMap genoMap;
	double* beta;
	double* qualTermForTransformation;
	double* posTermForTransformation;
	bool qualTransformationInitialized;

	//helper tools
	BamTools::BamAlignment bamAlignment;
	char toBase[4];
	std::vector<float> baseFreq;
	float cumulBaseFreq[4];
	bool refInitialized;

	int sampleQuality();
	double dePhred(double x);
	void initializeQualToErrorTable();
	int transformQuality(int & qual, int pos, int context);
	void fillMutationTable(float** & mutTable, double theta);
	void simulateDiploidHaplotypesCurChromosome(short** haplotypes, float** & mutTable, short* ref, const double & referenceDivergence);
	void writeInvariantSites(short** haplotypes, std::ofstream & genoFile);
	//void simulateReads(int & numReads, long & pos, float* & altFreq);
	void simulateReadsFromHaplotypes(std::vector<TSimulatorChromosome>::iterator & thisChr, short** haplotypes, TSimulatorBamFile & bamFile);
	void writeRead(long & pos, short* haplotype, TSimulatorBamFile & bamFile);

public:
	TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator, std::vector<float> & baseFreq);
	~TSimulator(){
		if(qualToErroTableInitialized)
			delete[] qualToErroTable;
	}

	//functions to set general parameters
	void setQualityDistribution(double mean, double sd);
	void setMaxQual(int maxQual);
	void setReadLength(int length);
	void setDepth(float depth);
	void setBaseFreq();
	void setReadGroupName(std::string name);
	void setPMD(TPMD* PmdObject);
	void setQualityTransformation(std::vector<double> & Betas);
	void initializeChromosomes(int numChr, long chrLength, bool haploid);
	void initializeChromosomes(std::vector<long> & chrLength, std::vector<bool> haploid);

	void simulatePooledData(int sampleSize, SFS & sfs, std::string outname);
	void simulateSingleIndividual(double theta, double referenceDivergence, std::string outname);

	void simulateIndividualPair(std::vector<double> & phis, double referenceDivergence, std::string outname);

};



#endif /* TSIMULATOR_H_ */
