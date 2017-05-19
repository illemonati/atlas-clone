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
#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamHeader.h"
#include "bamtools/api/BamAlignment.h"
#include <math.h>

class TSimulator{
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	BamTools::BamWriter bamWriter;
	BamTools::RefVector references;
	std::string bamFileName;
	bool bamFileOpen;
	std::ofstream fasta;
	bool fastaOpen;

	//general simulation parameters
	double meanQual, sdQual;
	float seqDepth;
	int readLength;
	std::map<std::string, long> chromosomes;
	std::map<std::string, long>::iterator chrIt;
	std::string readGroupName;

	//Qual to error table
	double* qualToErroTable;
	bool qualToErroTableInitialized;

	//PMD
	TPMD* pmdObject;
	bool pmdInitialized;

	//Quality transformation
	double* beta;
	double* qualTermForTransformation;
	double* posTermForTransformation;
	bool qualTransformationInitialized;

	//helper tools
	BamTools::BamAlignment bamAlignment;
	char toBase[4];
	short* ref;
	short* alt;
	bool refInitialized;

	void openBamFile(std::string Filename);
	void closeBamFile();
	void indexBamFile(std::string & filename);
	void openFastaFile(std::string filename);
	void closeFastaFile();
	void simulateReferenceAndAlternativeSequenceCurChromosome();
	void clearRefStorage();
	int sampleQuality();
	double dePhred(double x);
	void initializeQualToErrorTable();
	void simulateReads(int & numReads, long & pos, float* & altFreq);

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
	void setReadGroupName(std::string name);
	void setPMD(TPMD* PmdObject);
	void setQualityTransformation(std::vector<double> & Betas);
	void initializeChromosomes(int numChr, long chrLength);
	void initializeChromosomes(std::map<std::string, long> & chr);


	void simulatePooledData(int sampleSize, SFS & sfs, std::string outname);
	void simulateSingleIndividual(double theta, std::string outname);

};



#endif /* TSIMULATOR_H_ */
