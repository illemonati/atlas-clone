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

	//open file
	outFile.open(filename.c_str());;
	if(!outFile)
		throw "Failed to open file '" + filename + "' for writing!";

}

void TAlleleCountFile::writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile){
	bool useLocusName = params.parameterExists("useLocusName");
	if(useLocusName){
		logfile->list("Will print locus names (rather than chromosome and position).");
		outFile << "Locus";
		sep = '_';
	}
	outFile << "chr" << sep << "pos\tnumAlleles";
	for(int p=0; p<samples.numPopulations(); p++)
		outFile << "\t" << samples.getPopulationName(p);
	outFile << "\n";

}

void TAlleleCountFile::writePosition(std::string chr, long pos){
	outFile << chr << sep << pos;
}

void TAlleleCountFile::writeCounts(int count, int numAlleles){
	outFile << "\t" << count << "/" << numAlleles;
}

void TAlleleCountFile::endl(){
	outFile << std::endl;
}
