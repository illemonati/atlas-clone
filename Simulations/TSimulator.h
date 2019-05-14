/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include "SFS.h"
#include "../stringFunctions.h"
#include <math.h>
#include <numeric>
#include <algorithm>
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorQualityTransformation.h"
#include "TSimulatorRead.h"

//---------------------------------------------------------
//TSimulator
//---------------------------------------------------------
class TSimulator{
protected:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	std::string outname;

	//general simulation parameters
	int sampleSize;
	double referenceDivergence;
	float cumulRef[4];
	double seqDepth;
	double averageReadLength;
	double maxReadLength;

	//chromosomes
	std::vector<TSimulatorChromosome> chromosomes;
	bool writeTrueGenotypes;
	bool writeVariantInvariantBedFiles;
	TSimulatorReference referenceObj;

	//read simulator
	std::vector<TSimulatorSingleEndRead*> readSimulators;
	std::vector<std::string> readGroupNames;
	std::vector<double> simGroupFrequencies;
	std::vector<double> cumulSimGroupFrequenies;

	//Quality transformation
	TGenotypeMap genoMap;

	//helper tools
	char toBase[4] = {'A', 'C', 'G', 'T'};
	float baseFreq[4];
	float cumulBaseFreq[4];
	bool refInitialized;

	//function to initialize read groups
	void initializeCommonSettings(TParameters & params);
	void saveToMap(std::string & name, std::string args, std::map<std::string, std::string> & map, std::string & filename);
	void initializeReadLengthDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & readLengthMap);
	void initializeQualityDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & qualityDistMap);
	void initializeQualityTransformations(TParameters & params, bool & perReadGroup, std::map<std::string, TSimulatorQualityTransformParameters > & qualTransformMap);
	void initializePMD(TParameters & params, bool & perReadGroup, std::map<std::string, std::pair<std::string, std::string> > & pmdMap);
	void initializeContamination(TParameters & params, bool & perReadGroup, std::map<std::string, double> & contaminationMap);
	void addToReadGroupVector(std::vector<std::string> & vec, const std::string & rg);
	void initializeReadGroup(const std::string & readLengthString, std::string & readGroupName, int rgNumber, int maxPrintQual);
	void initializeReadSimulator(TParameters & params);
	void initializeReadGroupFrequencies(TParameters & params);

	//functions to simulate
	virtual void simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){ throw "simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref) not implemented for base class TSimulator!"; };
	virtual void simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref){ throw "simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref) not implemented for base class TSimulator!"; };
	void simulateReadsFromHaplotypes(std::vector<TSimulatorChromosome>::iterator & thisChr, Base** haplotypes, TSimulatorBamFile & bamFile, std::string extraProgressText);
	//void writeRead(const long & pos, short* haplotype, TSimulatorBamFile & bamFile);

	//from SFS
	//void fillMutationTable(float** & mutTable);
	//void simulateHaplotypes(TSimulatorHaplotypes & haplotypes, SFS* sfs, float** & mutTable, Base* ref);

public:
	TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TSimulator(){
		for(TSimulatorSingleEndRead* readSimIt: readSimulators)
			delete readSimIt;
	};

	//functions to set general parameters
	void setQualityDistribution(double mean, double sd, int maxQual);
	void setReadLength(std::string s);
	void setDepth(float depth);
	void setBaseFreq(std::vector<float> & freq);
	void setQualityTransformation(std::vector<double> & Betas);
	void initializeChromosomes(TParameters & params, TLog* logfile);
	void initializeChromosomes(int numChr, long chrLength, bool haploid);
	void initializeChromosomes(std::vector<long> & chrLength, std::vector<bool> haploid);

	void runSimulations();
};

//---------------------------------------------------------
//TSimulatorOneIndividual
//---------------------------------------------------------
class TSimulatorOneIndividual:public TSimulator{
private:
	std::vector<double> thetas;
	TSimulatorMutationtable mutTable;

	void simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);
	void simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);


public:
	TSimulatorOneIndividual(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorOneIndividual();

	void runSimulations();
};


//---------------------------------------------------------
//TSimulatorPairOfIndividuals
//---------------------------------------------------------
class TSimulatorPairOfIndividuals:public TSimulator{
private:
	std::vector<double> phis;
	double cumulGenoCaseFrequencies[9];
	int numGenotypeCombinations[9];
	double** cumulGenoCombinationFreq;
	Base*** genoTrans;
	short** orderLookup;
	bool tablesInitialized;

	void fillTables();
	void deleteTables();

	void simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);
	void simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);

public:

	TSimulatorPairOfIndividuals(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorPairOfIndividuals(){ deleteTables(); };

	void runSimulations();
};

//---------------------------------------------------------
//TSimulatorSFS
//---------------------------------------------------------
class TSimulatorSFS:public TSimulator{
private:
	std::vector<SFS*> sfs;
	TSimulatorMutationtable mutTable;

	void initializeSFS(std::vector<double> & thetas);
	void initializeSFS(std::vector<std::string> & sfsFileNames, bool folded);
	void simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);
	void simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);

public:
	TSimulatorSFS(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorSFS();

	void runSimulations();
};

//---------------------------------------------------------
//TSimulatorHardyWeinberg
//---------------------------------------------------------
class TSimulatorHardyWeinberg:public TSimulator{
private:
	double fracPoly, alpha, beta, F;
	double cumulGenoProb[3];
	TSimulatorMutationtable mutTable;
	bool writeTrueAlleleFreq;
	std::string alleleFreqFile;
	std::string alleleFreqFileMAF;

	void fillCumulGenoProb(const double & f);
	void fillhaplotypesMonomoprhic(TSimulatorHaplotypes & haplotypes, int & locus, Base* ref);
	void simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);
	void simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref);

public:
	TSimulatorHardyWeinberg(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorHardyWeinberg(){};

	void runSimulations();
};


#endif /* TSIMULATOR_H_ */
