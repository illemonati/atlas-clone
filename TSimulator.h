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

class TSimulatorChromosome{
public:
	std::string name;
	long length;
	bool haploid;

	TSimulatorChromosome(std::string Name, long Length, bool Haploid){
		name = Name;
		length = Length;
		haploid = Haploid;
	};
};

class TSimulator{
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	BamTools::BamWriter bamWriter;
	BamTools::RefVector references;
	std::string bamFileName;
	bool bamFileOpen;
	std::ofstream fasta;
	std::ofstream fastaIndex;
	long oldOffset;
	bool fastaOpen;

	//general simulation parameters
	double meanQual, sdQual;
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
	long refLength;
	short* ref;
	float baseFreq[4];
	float cumulBaseFreq[4];
	bool refInitialized;

	void openBamFile(std::string Filename);
	void closeBamFile();
	void indexBamFile(std::string & filename);
	void openFastaFile(std::string filename);
	void closeFastaFile();
	void writeRefToFasta();
	void simulateReferenceSequenceCurChromosome();
	void initializeRefStorage();
	void clearRefStorage();
	int sampleQuality();
	double dePhred(double x);
	void initializeQualToErrorTable();
	int transformQuality(int & qual, int pos, int context);
	void fillMutationTable(float** & mutTable, double theta);
	void simulateDiploidHaplotypesCurChromosome(short** haplotypes, float** & mutTable, const double & referenceDivergence);
	void writeTrueGenotypes(short** haplotypes, std::ofstream & genoFile);
	void simulateReads(int & numReads, long & pos, float* & altFreq);
	void writeRead(long & pos, short* haplotype);

public:
	TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TSimulator(){
		closeBamFile();
		closeFastaFile();
		clearRefStorage();
		if(qualToErroTableInitialized)
			delete[] qualToErroTable;
	}

	//functions to set general parameters
	void setQualityDistribution(double mean, double sd);
	void setReadLength(int length);
	void setDepth(float depth);
	void setBaseFreq(float* freq);
	void setReadGroupName(std::string name);
	void setPMD(TPMD* PmdObject);
	void setQualityTransformation(std::vector<double> & Betas);
	void initializeChromosomes(int numChr, long chrLength, bool haploid);
	void initializeChromosomes(std::vector<long> & chrLength, std::vector<bool> haploid);

	void simulatePooledData(int sampleSize, SFS & sfs, std::string outname);
	void simulateSingleIndividual(double theta, double referenceDivergence, std::string outname);

};



#endif /* TSIMULATOR_H_ */
