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
	~TSimulatorBamFile(){
		close();
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
			initialized = false;
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

	void writeGenotypes(gz::ogzstream & out, std::string & chrName, char* toBase){
		for(int l=0; l<curLength; ++l){
			out << chrName << "\t" << l+1;
			for(ind=0; ind < numInd; ++ind){


				//std::cout << ind << " @ " << l << ": " << std::flush;
				//std::cout << "\t" << toBase[haplotypes[ind][0][l]] << "/" << toBase[haplotypes[ind][1][l]] << std::endl;

				out << "\t" << toBase[haplotypes[ind][0][l]] << "/" << toBase[haplotypes[ind][1][l]];
			}
			out << "\n";
		}
	};

	int size(){
		return numInd;
	};

	short& operator()(int ind, int hap, long site){
		return haplotypes[ind][hap][site];
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

	void fillTables(std::vector<double> & phis, float* baseFreq);
	void deleteTables();

public:
	TSimulatorGenotypeCombination(std::vector<double> & phis, float* baseFreq){
		tablesInitialized = false;
		fillTables(phis, baseFreq);
	};
	~TSimulatorGenotypeCombination(){
		deleteTables();
	}

	void simulateHaplotypes(short** haplotypesInd0, short** haplotypesInd1, short* ref, float referenceDivergence, long length, TRandomGenerator* randomGenerator);

};

//---------------------------------------------------------
//TSimulatorReadLength
//---------------------------------------------------------
struct readLengthContainer{
	int fragmentLength;
	int readLength;
};

class TSimulatorReadLength{
protected:
	TRandomGenerator* randomGenerator;
	double meanLength;
	double cumulAtMin;

public:
	TSimulatorReadLength(TRandomGenerator* RandomGenerator, std::string & s){
		randomGenerator = RandomGenerator;

		//is a fixed length
		meanLength = stringToDouble(s);
		if(meanLength < 5 || meanLength > 10000)
			throw "Read length must be between 5 and 10,000!";
		cumulAtMin = 0.0;
	};
	TSimulatorReadLength(TRandomGenerator* RandomGenerator){
		randomGenerator = RandomGenerator;
		meanLength = -1;
		cumulAtMin = 0.0;
	};
	virtual ~TSimulatorReadLength(){};

	virtual void sample(readLengthContainer & rl){
		rl.fragmentLength = meanLength;
		rl.readLength = meanLength;
	};
	virtual int max(){return meanLength;};
	virtual double mean(){return meanLength;};
	virtual double probAcceptance(){return 1.0 - cumulAtMin;};
	virtual std::string getFunctionString(){ return "Will simulate reads of fixed length " + toString(meanLength) + ".";};
};

class TSimulatorReadLengthGamma:public TSimulatorReadLength{
protected:
	double alpha, beta;
	int _min, _max;

	void parseFunctionString(std::string & s, double & param1, double & param2);
	void calculateAverageLength();

public:
	TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator, std::string & s);
	TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator);
	void sample(readLengthContainer & rl);
	virtual int max(){return _max;};
	virtual std::string getFunctionString(){ return "Will simulate reads of gamma distributed length with alpha=" + toString(alpha) + " and beta=" + toString(beta) + ".";};
};

class TSimulatorReadLengthGammaMode:public TSimulatorReadLengthGamma{
protected:
	double mode, var;

public:
	TSimulatorReadLengthGammaMode(TRandomGenerator* RandomGenerator, std::string & s);
	std::string getFunctionString(){ return "Will simulate reads of gamma distributed length with mode=" + toString(mode) + " and variance=" + toString(var) + ".";};
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
	double referenceDivergence;
	double meanQual, sdQual;
	int maxQual;
	float seqDepth;
	TSimulatorReadLength* readLengthDist;
	bool readLengthDistInitialized;
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
	float baseFreq[4];
	float cumulBaseFreq[4];
	bool refInitialized;

	void initializeQualityTransform(TParameters & params);
	int sampleQuality();
	double dePhred(double x);
	void initializeQualToErrorTable();
	int transformQuality(int & qual, int pos, int context);
	void fillMutationTable(float** & mutTable, double theta);
	void simulateDiploidHaplotypesCurChromosome(short** haplotypes, float** & mutTable, short* ref);
	void writeBEDFiles(short** haplotypes, gz::ogzstream & invariantSitesFile, gz::ogzstream & variantSitesFile);
	//void simulateReads(int & numReads, long & pos, float* & altFreq);
	void simulateReadsFromHaplotypes(std::vector<TSimulatorChromosome>::iterator & thisChr, short** haplotypes, TSimulatorBamFile & bamFile, std::string extraProgressText);
	void writeRead(long & pos, short* haplotype, TSimulatorBamFile & bamFile);

	//from SFS
	void fillMutationTable(float** & mutTable);
	void simulateHaplotypes(TSimulatorHaplotypes & haplotypes, SFS* sfs, float** & mutTable, short* ref);

public:
	TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator, TParameters & params);
	~TSimulator(){
		if(qualToErroTableInitialized)
			delete[] qualToErroTable;
		if(readLengthDistInitialized)
			delete readLengthDist;
	}

	//functions to set general parameters
	void setQualityDistribution(double mean, double sd, int maxQual);
	void setReadLength(std::string s);
	void setDepth(float depth);
	void setBaseFreq(std::vector<float> & freq);
	void setReadGroupName(std::string name);
	void setPMD(TPMD* PmdObject);
	void setQualityTransformation(std::vector<double> & Betas);
	void initializeChromosomes(TParameters & params, TLog* logfile);
	void initializeChromosomes(int numChr, long chrLength, bool haploid);
	void initializeChromosomes(std::vector<long> & chrLength, std::vector<bool> haploid);

	void simulatePooledData(int sampleSize, SFS & sfs, std::string outname);
	void simulateSingleIndividual(double theta, std::string & outname);
	void simulateSingleIndividual(std::vector<double> theta, std::string & outname);
	void simulateIndividualPair(std::vector<double> & phis, std::string & outname);
	void simulatePopulationFromSFS(double theta, int numIndividuals, std::string & outname);
	void simulatePopulationFromSFS(std::vector<double> & thetas, int numIndividuals, std::string & outname);
	void simulatePopulationFromSFS(std::vector<std::string> & sfsFileNames, bool folded, int numIndividuals, std::string & outname);
	void simulatePopulationFromSFS(std::vector<SFS*> sfs, int numIndividuals, std::string & outname);
};



#endif /* TSIMULATOR_H_ */
