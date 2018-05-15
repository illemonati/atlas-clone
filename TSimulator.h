/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include "SFS.h"
#include "stringFunctions.h"
#include <math.h>
#include <numeric>
#include <algorithm>
#include "TSimulatorRead.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorQualityTransformation.h"

//---------------------------------------------------------
//TSimulatorGenotypecombination
//---------------------------------------------------------
class TSimulatorGenotypeCombination{
private:
	double cumulGenoCaseFrequencies[9];
	int numGenotypeCombinations[9];
	double** cumulGenoCombinationFreq;
	Base*** genoTrans;
	short** orderLookup;
	bool tablesInitialized;

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

	void simulateHaplotypes(Base** haplotypesInd0, Base** haplotypesInd1, Base* ref, float referenceDivergence, long length, TRandomGenerator* randomGenerator);
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
	double seqDepth;
	double averageReadLength;
	double maxReadLength;

	//chromosomes
	std::vector<TSimulatorChromosome> chromosomes;
	std::vector<TSimulatorChromosome>::iterator chrIt;
	bool writeTrueGenotypes;
	TSimulatorReference referenceObj;

	//read simulator
	std::vector<TSimulatorRead*> readSimulators;
	std::vector<TSimulatorRead*>::iterator readSimsIt;
	std::vector<std::string> readGroupNames;
	std::vector<double> simGroupFrequencies;
	std::vector<double> cumulSimGroupFrequenies;

	//Quality transformation
	TGenotypeMap genoMap;
	double* beta;
	double* qualTermForTransformation;
	double* posTermForTransformation;

	//helper tools
	char toBase[4] = {'A', 'C', 'G', 'T'};
	float baseFreq[4];
	float cumulBaseFreq[4];
	bool refInitialized;

	//function to initialize read groups
	void saveToMap(std::string & name, std::string args, std::map<std::string, std::string> & map, std::string & filename);
	void initializeReadLengthDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & readLengthMap);
	void initializeQualityDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & qualityDistMap);
	void initializeQualityTransformations(TParameters & params, bool & perReadGroup, std::map<std::string, std::pair<std::string, std::string> > & qualTransformMap);
	void initializePMD(TParameters & params, bool & perReadGroup, std::map<std::string, std::pair<std::string, std::string> > & pmdMap);
	void initializeContamination(TParameters & params, bool & perReadGroup, std::map<std::string, double> & contaminationMap);
	void addToReadGroupVector(std::vector<std::string> & vec, const std::string & rg);
	void initializeReadSimulator(TParameters & params);
	void initializeReadGroupFrequencies(TParameters & params);

	//function sto simulate
	void fillMutationTable(float** & mutTable, double theta);
	void simulateDiploidHaplotypesCurChromosome(Base** haplotypes, float** & mutTable, Base* ref);
	void writeBEDFiles(Base** haplotypes, Base* ref, gz::ogzstream & invariantSitesFile, gz::ogzstream & variantSitesFile);
	void writeVCFwithVariantPositions(Base** haplotypes, Base* ref, const std::string & chr, gz::ogzstream & vcf);
	void simulateReadsFromHaplotypes(std::vector<TSimulatorChromosome>::iterator & thisChr, Base** haplotypes, TSimulatorBamFile & bamFile, std::string extraProgressText);
	void writeRead(const long & pos, short* haplotype, TSimulatorBamFile & bamFile);

	//from SFS
	void fillMutationTable(float** & mutTable);
	void simulateHaplotypes(TSimulatorHaplotypes & haplotypes, SFS* sfs, float** & mutTable, Base* ref);

public:
	TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator, TParameters & params);
	~TSimulator(){
		for(readSimsIt=readSimulators.begin(); readSimsIt!=readSimulators.end(); ++readSimsIt)
			delete *readSimsIt;

	}

	//functions to set general parameters
	void setQualityDistribution(double mean, double sd, int maxQual);
	void setReadLength(std::string s);
	void setDepth(float depth);
	void setBaseFreq(std::vector<float> & freq);
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
