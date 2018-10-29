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

class VcfDiagnostics{
private:
	TParameters* params;
	TLog* logfile;
	int chr;
	std::ifstream vcfFileStream;
	std::ofstream vcfOutFilestream;
	bool verbose;
	TRandomGenerator* randomGenerator;
	bool randomGeneratorInitialized;
	TVcfFileSingleLine vcfFile;

	void openVCF(TVcfFile_base & vcfFile);
	void initializeRandomGenerator();
	void initializeCountsTable(int** table, int nrows, int ncols);
	void deleteCountsTable(int** table, int nrows);
	std::pair<char, char> getGenotypeFromIndex(int index);

public:
	VcfDiagnostics(TParameters* Params, TLog* Logfile);
	~VcfDiagnostics(){if(randomGeneratorInitialized) delete randomGenerator;};
	int baseToNumber(char base, std::string & marker);
	void vcfToBeagle();
	void assessAllelicImbalance();
	void filterAllelicImbalance();
	void vcfToInvariantBed();
};




#endif /* TVCFDIAGNOSTICS_H_ */

