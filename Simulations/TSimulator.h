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
#include "algorithmsAndVectors.h"
#include "progressTools.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorRead.h"
#include "TFile.h"
#include "TTask.h"
#include <filesystem>
#include <functional>

namespace Simulations{

//TODO: add cross-contamination between samples or RGs? That would be easier to model contamination that the way it is done now as it would allow for contaminated reads to have different characteristsics.

using coretools::TParameters;
using coretools::TLog;
using coretools::TRandomGenerator;
using genometools::Base;

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
	std::array<double, 4> _cumulRef;
	double _seqDepth;
	double _averageReadLength;
	double _maxReadLength;

	//chromosomes
	BAM::TChromosomes _chromosomes;
	bool _writeTrueGenotypes;
	bool _writeVariantInvariantBedFiles;
	TSimulatorReference _referenceObj;

	//simulation tools
	BAM::TReadGroups _readGroups;
	GenotypeLikelihoods::TPostMortemDamage _PMD;
	//BAM::TReadGroupMap _readGroupMap; //needed by recal REALLYY??????
	GenotypeLikelihoods::TSequencingErrorModels _recal;

	//read simulator
	std::vector<TSimulatorSingleEndRead*> _readSimulators;
	std::vector<double> _simGroupFrequencies;
	std::vector<double> _cumulSimGroupFrequenies;

	//helper tools
	GenotypeLikelihoods::TBaseProbabilities _baseFreq;
	std::array<double, 4> _cumulBaseFreq;
	bool _refInitialized;

	//function to initialize read groups
	void _initializeCommonSettings(TParameters & params);
	std::vector<std::string>  _readSimInfoPerReadGroup(const std::string & Filename, const std::string & Column, const std::string & Name);
	void _initializeReadGroup(const std::string & readLengthString, const BAM::TReadGroup & ReadGroup);
	void _initializeReadGroupsFromReadLengthDistribution(TParameters & params, const std::string & ParameterName, const std::string & DefaultValue, const std::string & Name);
	void _initializeDistribution(TParameters & params, const std::string & ParameterName, const std::string & DefaultValue, const std::string & Name, std::function<void(TSimulatorSingleEndRead&, std::string)> function);
	void _initializePMD(TParameters & params, const std::string & ParameterName, const std::string & Name);
	void _initializeQualityTransformations(TParameters & params, const std::string & ParameterName, const std::string & Name);

	void _initializeContamination(TParameters & params, bool & perReadGroup, std::map<std::string, double> & contaminationMap);
	void _addToReadGroupVector(std::vector<std::string> & vec, const std::string & rg);
	void _addReadGroupsIfFile(const std::string & ParameterName, TParameters & Parameters, BAM::TReadGroups & ReadGroups);
	void _initializeReadSimulator(TParameters & params);
	void _initializeReadGroupFrequencies(TParameters & params);

	//functions to simulate
	Base _sampleBase(const std::array<double, 4> & cumulProbs);
	Base _mutateBase(const Base & base, const std::array<double, 4> & cumulProbs);
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
	std::array< std::array<uint8_t, 4>, 4> orderLookup;
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
	genometools::Base reference;
	genometools::Base alternative;
	double f;
};

class TSimulatorHardyWeinberg:public TSimulator{
private:
	double fracPoly, alpha, beta, F;
	double cumulGenoProb[3];
	TSimulatorMutationtable mutTable;
	bool writeTrueAlleleFreq;
	coretools::TOutputFile trueFreqFile;

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
class TTask_simulate:public coretools::TTask{
public:
	TTask_simulate(){ _explanation = "Generating simulations"; };

	void run(TParameters & Parameters, TLog* Logfile){
		//initialize simulator
		TSimulator* simulator;
		std::string method = Parameters.getParameterWithDefault<std::string>("type", "one");
		if(method == "one"){
			Logfile->startIndent("Simulating a single individual (parameter type=one):");
			simulator = new TSimulatorOneIndividual(Logfile, Parameters, _randomGenerator);
		} else if(method == "pair"){
			Logfile->startIndent("Simulating a pair of individual (parameter type=pair):");
			simulator = new TSimulatorPairOfIndividuals(Logfile, Parameters, _randomGenerator);
		} else if(method == "SFS"){
			Logfile->startIndent("Simulating individuals from an SFS (parameter type=SFS):");
			simulator = new TSimulatorSFS(Logfile, Parameters, _randomGenerator);
		} else if(method == "HW"){
			Logfile->startIndent("Simulating a individuals under Hardy-Weinberg (parameter type=HW):");
			simulator = new TSimulatorHardyWeinberg(Logfile, Parameters, _randomGenerator);
		}else throw "Unknown simulation method '" + method + "'!";

		//now run simulations
		simulator->runSimulations();

		//clean up
		Logfile->endIndent();
		delete simulator;
	};
};

}; //end namespace

#endif /* TSIMULATOR_H_ */
