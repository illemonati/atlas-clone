/*
 * TAnnotator.h
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */

#ifndef TVCFDIAGNOSTICS_H_
#define TVCFDIAGNOSTICS_H_

#include <string>
#include <vector>

#include "coretools/Main/TParameters.h"
#include "coretools/Main/TLog.h"
#include "genometools/VCF/TVcfFile.h"

namespace VCF{

class TCountTable{
private:
	int nrows;
	int ncols;
	int initializationValue;
	std::ofstream out;
	std::string outname;

public:
	int** table;

	TCountTable(int Nrows, int Ncols, std::string_view Outname){
		nrows = Nrows;
		ncols = Ncols;
		initializationValue = 0;
		outname = Outname;
		initialize();
		openOut(Outname);
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

	void openOut(std::string_view outname){
		coretools::instances::logfile().list("Writing count table to '", outname, "'.");
		out.open(std::string(outname).c_str());
		coretools::user_assert(out.is_open(), "Failed to open file '", outname, " for writing!");
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

class TVcfDiagnostics{
private:
	int chr = -1;
	std::ifstream _vcfFileStream;
	std::ofstream _vcfOutFilestream;
	genometools::TVcfFileSingleLine _vcfFile;
	std::string _outName;

	void _openVCF(const std::string &VCFName, bool isZipped);

public:
	TVcfDiagnostics();

	int findLastPassedFilterIndex(int obsValue, std::vector<int> & filtersAscendingOrder);
	void assessAllelicImbalance();
	void vcfToInvariantBed();
	void fixIntAsFloat();
	void run() {
		using coretools::instances::parameters; 
		if (parameters().exists("fixInt")) {
			fixIntAsFloat();
		} else if (parameters().exists("writeBED")) {
			vcfToInvariantBed();
		} else {
			assessAllelicImbalance();
		}
	}
};

}; //end namespace


#endif /* TVCFDIAGNOSTICS_H_ */

