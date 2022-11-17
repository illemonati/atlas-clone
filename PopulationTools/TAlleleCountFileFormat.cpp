/*
 * TAlleleCountFileFormat.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: linkv
 */

#include "TAlleleCountFileFormat.h"

#include <stddef.h>
#include <ostream>

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/VCF/TPopulation.h"
#include "coretools/Strings/stringFunctions.h"

namespace PopulationTools{
using coretools::instances::logfile;
using coretools::instances::parameters;

//------------------------------
// TAlleleCountFile
//------------------------------

TAlleleCountFile::TAlleleCountFile(std::string Filename){
	filename = Filename;
	sep = '\t';
};

void TAlleleCountFile::openFileToWrite(std::string filename){
	outFile.open(filename.c_str());;
	if(!outFile)
		throw "Failed to open file '" + filename + "' for writing!";
};

void TAlleleCountFile::writeHeader(genometools::TPopulationSamples & samples){
	bool useLocusName = parameters().parameterExists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile << "Locus";
		sep = '_';
	}
	outFile << "chr" << sep << "pos";
	for(size_t p=0; p<samples.numPopulations(); p++)
		outFile << "\t" << samples.getPopulationName(p);
	outFile << "\n";
};

void TAlleleCountFile::writeHeader(std::vector<std::string> populationNames){
	bool useLocusName = parameters().parameterExists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile << "Locus";
		sep = '_';
	}
	outFile << "chr" << sep << "pos";
	for(size_t p=0; p<populationNames.size(); p++)
		outFile << "\t" << populationNames[p];
	outFile << "\n";
};

void TAlleleCountFile::writePosition(std::string chr, long pos){
	outFile << chr << sep << pos;
};

void TAlleleCountFile::writePosition(std::string chr, std::string pos){
	outFile << chr << sep << pos;
};

void TAlleleCountFile::writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader){
	outFile << reader.chr() << sep << reader.position();
}

void TAlleleCountFile::writeCounts(int count, int numAlleles, int){
	outFile << "\t" << count << "/" << numAlleles;
};

void TAlleleCountFile::writeCounts(std::string count, std::string numAlleles, int){
	outFile << "\t" << count << "/" << numAlleles;
};

void TAlleleCountFile::endl(){
	outFile << "\n";
};

//------------------------------
// TAlleleCountFileWithAlleles
//------------------------------

void TAlleleCountFileWithAlleles::writeHeader(genometools::TPopulationSamples & samples){
	bool useLocusName = parameters().parameterExists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile << "Locus";
		sep = '_';
	}
	outFile << "chr" << sep << "pos" << sep <<"ref" << sep << "alt";
	for(size_t p=0; p<samples.numPopulations(); p++)
		outFile << "\t" << samples.getPopulationName(p);
	outFile << "\n";
};

void TAlleleCountFileWithAlleles::writeHeader(std::vector<std::string> populationNames){
	bool useLocusName = parameters().parameterExists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile << "Locus";
		sep = '_';
	}
	outFile << "chr" << sep << "pos" << sep <<"ref" << sep << "alt";
	for(size_t p=0; p<populationNames.size(); p++)
		outFile << "\t" << populationNames[p];
	outFile << "\n";
};

void TAlleleCountFileWithAlleles::writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader){
	outFile << reader.chr() << sep << reader.position() << sep << reader.refAllele() << sep << reader.altAllele();
}

//------------------
// TTreeMixFile
//------------------

void TTreeMixFile::writeHeader(genometools::TPopulationSamples & samples){
	outFile << samples.getPopulationName(0);
	for(size_t p=1; p<samples.numPopulations(); p++)
		outFile << " " << samples.getPopulationName(p);
	outFile << "\n";
};

void TTreeMixFile::writeHeader(std::vector<std::string> populationNames){
	outFile << populationNames[0];
		for(size_t p=1; p<populationNames.size(); p++)
			outFile << " " << populationNames[p];
	outFile << "\n";
}

void TTreeMixFile::writeCounts(int count, int numAlleles, int populationNum){
	if(populationNum == 0)
		outFile << count << "," << numAlleles - count;
	else
		outFile << " " << count << "," << numAlleles - count;
};

void TTreeMixFile::writeCounts(std::string count, std::string numAlleles, int populationNum){
	if(populationNum == 0)
		outFile << count << "," << coretools::str::fromString<int>(numAlleles) - coretools::str::fromString<int>(count);
	else
		outFile << " " << count << "," << coretools::str::fromString<int>(numAlleles) - coretools::str::fromString<int>(count);
};

//------------------
// TFlinkFile
//------------------

void TFlinkFile::writeHeader(genometools::TPopulationSamples & samples){
	outFile << "-\t-";
	for(size_t g=0; g<samples.numPopulations(); ++g){
		outFile << "\tGroup_A";
	}
	outFile << "\n";
	outFile << "-\t-";
	for(size_t p=0; p<samples.numPopulations(); p++)
		outFile << "\t" << samples.getPopulationName(p);
	outFile << "\n";

};

void TFlinkFile::writeHeader(std::vector<std::string> populationNames){
	outFile << "-\t-";
		for(unsigned int g=0; g<populationNames.size(); ++g){
			outFile << "\tGroup_A";
		}
		outFile << "\n";
		outFile << "-\t-";
		for(unsigned int p=0; p<populationNames.size(); p++)
			outFile << "\t" << populationNames[p];
		outFile << "\n";
}

void TFlinkFile::writePosition(std::string chr, long pos){
	outFile << chr << sep << pos;
};

void TFlinkFile::writePosition(std::string chr, std::string pos){
	outFile << chr << sep << pos;
};

void TFlinkFile::writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader){
	outFile << reader.chr() << sep << reader.position();
}

void TFlinkFile::writeCounts(int count, int numAlleles, int){
		outFile << "\t" << count << "/" << numAlleles;
};

void TFlinkFile::writeCounts(std::string count, std::string numAlleles, int){
		outFile << "\t" << count << "/" << numAlleles;
};

}; //end namespace
