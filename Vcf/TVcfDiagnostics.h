/*
 * TAnnotator.h
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */

#ifndef TVCFDIAGNOSTICS_H_
#define TVCFDIAGNOSTICS_H_

#include "../stringFunctions.h"
#include <vector>
#include "../TParameters.h"
#include <iostream>
#include <fstream>
#include "TVcfFile.h"
#include "../TRandomGenerator.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "../TQualityMap.h"
#include "../TFile.h"
#include "../PopulationTools/TGLF.h"
#include "../PopulationTools/TGenotypeFrequencies.h"

class TGenotype{
	public:
		char first;
		char second;
		TGenotype(char & f, char & s){
			first=f;
			second=s;
		};
		TGenotype(char f){
			first=f;
			second=f;
		};
};

class TCountTable{
private:
	int nrows;
	int ncols;
	int initializationValue;
	std::ofstream out;
	std::string outname;

public:
	int** table;

	TCountTable(int Nrows, int Ncols, std::string Outname, TLog* Logfile){
		nrows = Nrows;
		ncols = Ncols;
		initializationValue = 0;
		outname = Outname;
		initialize();
		openOut(Outname, Logfile);
	}

	void initialize(){
		table = new int*[nrows];
		for(int i=0; i<nrows; ++i){
			table[i] = new int[nrows];
		}

		for(int i=0; i<ncols; ++i){
			for(int j=0; j<ncols; ++j)
				table[i][j] = initializationValue;
		}
	}

	void openOut(std::string & outname, TLog* logfile){
		logfile->list("Writing count table to '" + outname + "'.");
		out.open(outname.c_str());
		if(!out)
			throw "Failed to open file '" + outname + " for writing!";
	}

	void writeTable(std::string & description, std::string & rowPrefix, std::string & colPrefix){
		//header
		out << description; //this goes in top left corner
		for(int i=0; i<ncols; ++i){
			out << "\t" << colPrefix << i;
		}
		out << "\n";

		//write rows
		for(int i=0; i<nrows; ++i){
			out << rowPrefix << i;
			for(int j=0; j<ncols; ++j)
				out << "\t" << table[i][j];
			out << "\n";
		}
	}

	~TCountTable(){
		//delete table
		for(int i = 0; i < nrows; ++i)
		    delete[] table[i];
		delete[] table;

		//close file
		out.close();
	}
};

class TVcfFilters{
private:
	void specifyChromosomesToKeep(TLog* logfile, TParameters & Parameters);
	void initialize(TParameters & Parameters, TLog* logfile);

public:
	//settings
	long limitLines;
	int minDepth;
	int minNumSamplesWithData;
	//double freqFilter;
	//double epsilonF; //F for EM algorithm to estimate allele frequencies
	int minVariantQuality;
	//bool estimateGenotypeFrequencies;
	long progressFrequency;
	std::string trueAlleleFreqFile;
	std::vector<std::string> chromosomesToKeep;

	TVcfFilters(TParameters & Parameters, TLog* logfile);
};

class VcfDiagnostics{
private:
	TLog* logfile;
	int chr;
	std::ifstream vcfFileStream;
	std::ofstream vcfOutFilestream;
	bool verbose;
	TRandomGenerator* randomGenerator;
	bool randomGeneratorInitialized;
	TVcfFileSingleLine vcfFile;
	std::string outname;
	bool isZipped;

	void openVCF(std::string filename, TVcfFile_base & vcfFile);
	void initializeRandomGenerator();
	std::pair<char, char> getGenotypeFromIndex(int index);

	// beagle
	void writeBeagleHeader(TOutputFileZipped & beagleFile);
	void printProgressFrequencyFiltering(const struct timeval & startTime, long lineCounter, long numAcceptedLoci);
	void concludeFilters(TVcfFilters & vcfFilters, const struct timeval & startTime,
			long lineCounter, long notBialleleicCounter,
			long missingSNPCounter, long lowFreqSNPCounter,
			long lowVariantQualityCounter,
			long noPLCounter, long notOnChrCounter, long numAcceptedLoci);
	int filterOnDepth(int numSamples, double * data, TVcfFilters & vcfFilters);
	bool readVcfAndWriteBeagle(int numSamples, TVcfFilters & vcfFilters,
			TOutputFileZipped & beagleFile,
			long & lineCounter, long & notBialleleicCounter,
			long & missingSNPCounter, long & lowFreqSNPCounter,
			long & lowVariantQualityCounter,
			long & noPLCounter, long & notOnChrCounter, long & numAcceptedLoci, struct timeval & startTime);
	void writeLocusToBeagleFile(TOutputFileZipped & beagleFile,
			double * data, int numSamples);
	inline double phredToError(double phred);
public:
	VcfDiagnostics(TParameters & Params, TLog* Logfile);
	~VcfDiagnostics(){if(randomGeneratorInitialized) delete randomGenerator;};
	int baseToNumber(char base, std::string & marker);
	void vcfToBeagle();
	int findLastPassedFilterIndex(int obsValue, std::vector<int> & filtersAscendingOrder);
	void assessAllelicImbalance(TParameters & Params);
	//void filterAllelicImbalance();
	void vcfToInvariantBed();

	// beagle
	void vcfToBeagle(TParameters & Params);

	void fixIntAsFloat();
};




#endif /* TVCFDIAGNOSTICS_H_ */

