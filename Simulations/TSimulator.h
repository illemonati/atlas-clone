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
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorQualityTransformation.h"
#include "TSimulatorRead.h"
#include "TFile.h"
#include "TTask.h"

namespace Simulations{

//---------------------------------------------------------
//TSimulator
//---------------------------------------------------------
class TSimulator{
protected:
	TLog* _logfile;
	TRandomGenerator* _randomGenerator;
	std::string _outname;

	//general simulation parameters
	int _sampleSize;
	double _referenceDivergence;
	float _cumulRef[4];
	double _seqDepth;
	double _averageReadLength;
	double _maxReadLength;

	//chromosomes
	BAM::TChromosomes _chromosomes;
	bool _writeTrueGenotypes;
	bool _writeVariantInvariantBedFiles;
	TSimulatorReference _referenceObj;

	//read simulator
	std::vector<TSimulatorSingleEndRead*> _readSimulators;
	std::vector<std::string> _readGroupNames;
	std::vector<double> _simGroupFrequencies;
	std::vector<double> _cumulSimGroupFrequenies;

	//Quality transformation
	GenotypeLikelihoods::TGenotypeMap _genoMap;
	BAM::TQualityMap _qualMap;

	//helper tools
	char _toBase[4] = {'A', 'C', 'G', 'T'};
	float _baseFreq[4];
	float _cumulBaseFreq[4];
	bool _refInitialized;

	//function to initialize read groups
	void _initializeCommonSettings(TParameters & params);
	void _saveToMap(std::string & name, std::string args, std::map<std::string, std::string> & map, std::string & filename);
	void _initializeReadLengthDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & readLengthMap);
	void _initializeQualityDistribution(TParameters & params, bool & perReadGroup, std::map<std::string, std::string> & qualityDistMap);
	void _initializeQualityTransformations(TParameters & params, bool & perReadGroup, std::map<std::string, TSimulatorQualityTransformParameters > & qualTransformMap);
	void _initializePMD(TParameters & params, bool & perReadGroup, std::map<std::string, std::pair<std::string, std::string> > & pmdMap);
	void _initializeContamination(TParameters & params, bool & perReadGroup, std::map<std::string, double> & contaminationMap);
	void _addToReadGroupVector(std::vector<std::string> & vec, const std::string & rg);
	void _initializeReadGroup(const std::string & readLengthString, std::string & readGroupName, int rgNumber, int maxPrintQual);
	void _initializeReadSimulator(TParameters & params);
	void _initializeReadGroupFrequencies(TParameters & params);

	//functions to simulate
	virtual void _simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){ throw "_simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref) not implemented for base class TSimulator!"; };
	virtual void _simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref){ throw "_simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome, Base* ref) not implemented for base class TSimulator!"; };
	void _simulateReadsFromHaplotypes(const BAM::TChromosome & thisChr, Base** haplotypes, TSimulatorBamFile & bamFile, std::string extraProgressText);

public:
	TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TSimulator(){
		for(TSimulatorSingleEndRead* readSimIt: _readSimulators)
			delete readSimIt;
	};

	//functions to set general parameters
	void setQualityDistribution(double mean, double sd, int maxQual);
	void setReadLength(std::string s);
	void setDepth(float depth);
	void setBaseFreq(std::vector<float> & freq);
	void setQualityTransformation(std::vector<double> & Betas);
	void initializeChromosomes(TParameters & params, TLog* logfile);
	void initializeChromosomes(const uint32_t & numChr, const uint32_t & chrLength, const uint8_t & ploidy);
	void initializeChromosomes(std::vector<uint32_t> & chrLength, std::vector<uint8_t> haploid);

	void runSimulations();
};

//---------------------------------------------------------
//TSimulatorOneIndividual
//---------------------------------------------------------
class TSimulatorOneIndividual:public TSimulator{
private:
	std::vector<double> thetas;
	TSimulatorMutationtable mutTable;

	void _simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);


public:
	TSimulatorOneIndividual(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorOneIndividual();
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

	void _simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);

public:

	TSimulatorPairOfIndividuals(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorPairOfIndividuals(){ deleteTables(); };
};

//---------------------------------------------------------
//TSimulatorSFS
//---------------------------------------------------------
class TSimulatorSFS:public TSimulator{
private:
	std::vector<SFS*> sfs;
	TSimulatorMutationtable mutTable;

	void _initializeSFS(std::vector<double> & thetas);
	void _initializeSFS(std::vector<std::string> & sfsFileNames, bool folded);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);

public:
	TSimulatorSFS(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorSFS();
};

//---------------------------------------------------------
//TSimulatorHardyWeinberg
//---------------------------------------------------------
struct TSimulatorHardyWeinbergSite{
	bool isPolymorphic;
	Base reference;
	Base alternative;
	double f;
};

class TSimulatorHardyWeinberg:public TSimulator{
private:
	double fracPoly, alpha, beta, F;
	double cumulGenoProb[3];
	TSimulatorMutationtable mutTable;
	bool writeTrueAlleleFreq;
	TOutputFile trueFreqFile;

	void _fillCumulGenoProb(const double & f);
	void _simulateSite(TSimulatorHardyWeinbergSite & site, const std::string & chr, const uint64_t & pos, TSimulatorReference & ref);
	void _fillhaplotypesMonomoprhic(TSimulatorHaplotypes & haplotypes, const uint64_t & locus, TSimulatorHardyWeinbergSite & site);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes & haplotypes, BAM::TChromosome & chromosome, TSimulatorReference & ref);

public:
	TSimulatorHardyWeinberg(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TSimulatorHardyWeinberg(){};
};


//--------------------------------------
// Tasks
//--------------------------------------
class TTask_simulate:public TTask{
public:
	TTask_simulate(){ _explanation = "Generating simulations"; };

	void run(TParameters & Parameters, TLog* Logfile){
		//initialize simulator
		TSimulator* simulator;
		std::string method = Parameters.getParameterStringWithDefault("type", "one");
		if(method == "one")
			simulator = new TSimulatorOneIndividual(Logfile, Parameters, _randomGenerator);
		else if(method == "pair")
			simulator = new TSimulatorPairOfIndividuals(Logfile, Parameters, _randomGenerator);
		else if(method == "SFS")
			simulator = new TSimulatorSFS(Logfile, Parameters, _randomGenerator);
		else if(method == "HW")
			simulator = new TSimulatorHardyWeinberg(Logfile, Parameters, _randomGenerator);
		else throw "Unknown simulation method '" + method + "'!";

		//now run simulations
		simulator->runSimulations();

		//clean up
		delete simulator;
	};
};

}; //end namespace

#endif /* TSIMULATOR_H_ */
