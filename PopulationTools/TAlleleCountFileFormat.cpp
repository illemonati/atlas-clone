/*
 * TAlleleCountFileFormat.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: linkv
 */


#include "TAlleleCountFileFormat.h"

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

void TAlleleCountFile::writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile){
	bool useLocusName = params.parameterExists("useLocusName");
	if(useLocusName){
		logfile->list("Will print locus names (rather than chromosome and position).");
		outFile << "Locus";
		sep = '_';
	}
	outFile << "chr" << sep << "pos";
	for(int p=0; p<samples.numPopulations(); p++)
		outFile << "\t" << samples.getPopulationName(p);
	outFile << "\n";
};

void TAlleleCountFile::writeHeader(std::vector<std::string> populationNames, TParameters & params, TLog* logfile){
	bool useLocusName = params.parameterExists("useLocusName");
	if(useLocusName){
		logfile->list("Will print locus names (rather than chromosome and position).");
		outFile << "Locus";
		sep = '_';
	}
	outFile << "chr" << sep << "pos";
	for(int p=0; p<populationNames.size(); p++)
		outFile << "\t" << populationNames[p];
	outFile << "\n";
};

void TAlleleCountFile::writePosition(std::string chr, long pos){
	outFile << chr << sep << pos;
};

void TAlleleCountFile::writePosition(std::string chr, std::string pos){
	outFile << chr << sep << pos;
};

void TAlleleCountFile::writeCounts(int count, int numAlleles, int populationNum){
	outFile << "\t" << count << "/" << numAlleles;
};

void TAlleleCountFile::writeCounts(std::string count, std::string numAlleles, int populationNum){
	outFile << "\t" << count << "/" << numAlleles;
};

void TAlleleCountFile::endl(){
	outFile << "\n";
};

//------------------
// TTreeMixFile
//------------------

TTreeMixFile::TTreeMixFile(std::string Filename):TAlleleCountFile(Filename){
	filename = Filename;
};

void TTreeMixFile::writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile){
	outFile << samples.getPopulationName(0);
	for(int p=1; p<samples.numPopulations(); p++)
		outFile << " " << samples.getPopulationName(p);
	outFile << "\n";
};

void TTreeMixFile::writeHeader(std::vector<std::string> populationNames, TParameters & params, TLog* logfile){
	outFile << populationNames[0];
		for(int p=1; p<populationNames.size(); p++)
			outFile << " " << populationNames[p];
	outFile << "\n";
}


void TTreeMixFile::writePosition(std::string chr, long pos){
	//do nothing, treemix does not need position
};

void TTreeMixFile::writePosition(std::string chr, std::string pos){
	//do nothing, treemix does not need position
};

void TTreeMixFile::writeCounts(int count, int numAlleles, int populationNum){
	if(populationNum == 0)
		outFile << count << "," << numAlleles - count;
	else
		outFile << " " << count << "," << numAlleles - count;
};

void TTreeMixFile::writeCounts(std::string count, std::string numAlleles, int populationNum){
	if(populationNum == 0)
		outFile << count << "," << stringToInt(numAlleles) - stringToInt(count);
	else
		outFile << " " << count << "," << stringToInt(numAlleles) - stringToInt(count);
};


//------------------
// TFlinkFile
//------------------

TFlinkFile::TFlinkFile(std::string Filename):TAlleleCountFile(Filename){
	filename = Filename;
	sep = "\t";
};

void TFlinkFile::writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile){
	outFile << "-\t-";
	for(int g=0; g<samples.numPopulations(); ++g){
		outFile << "\tGroup_A";
	}
	outFile << "\n";
	outFile << "-\t-";
	for(int p=0; p<samples.numPopulations(); p++)
		outFile << "\t" << samples.getPopulationName(p);
	outFile << "\n";

};

void TFlinkFile::writeHeader(std::vector<std::string> populationNames, TParameters & params, TLog* logfile){
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

void TFlinkFile::writeCounts(int count, int numAlleles, int populationNum){
		outFile << "\t" << count << "/" << numAlleles;
};

void TFlinkFile::writeCounts(std::string count, std::string numAlleles, int populationNum){
		outFile << "\t" << count << "/" << numAlleles;
};
