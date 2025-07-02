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

namespace PopulationTools{
using coretools::instances::logfile;
using coretools::instances::parameters;

//------------------------------
// TAlleleCountFile
//------------------------------

TAlleleCountFile::TAlleleCountFile(std::string Filename){
	filename = Filename;
	sep      = '\t';
}

void TAlleleCountFile::openFileToWrite(std::string filename) { outFile.open(filename); }

void TAlleleCountFile::writeHeader(genometools::TPopulationSamples & samples){
	bool useLocusName = parameters().exists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile.write("Locus");
		sep = '_';
	}
	outFile.write("chr",sep,"pos");
	for(size_t p=0; p<samples.numPopulations(); p++)
		outFile.write("\t",samples.getPopulationName(p));
	outFile.endln();
}

void TAlleleCountFile::writeHeader(std::vector<std::string> populationNames){
	bool useLocusName = parameters().exists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile.write("Locus");
		sep = '_';
	}
	outFile.write("chr",sep,"pos");
	for(size_t p=0; p<populationNames.size(); p++)
		outFile.write("\t",populationNames[p]);
	outFile.endln();
}

void TAlleleCountFile::writePosition(std::string chr, long pos){
	outFile.write(chr,sep,pos);
}

void TAlleleCountFile::writePosition(std::string chr, std::string pos){
	outFile.write(chr,sep,pos);
}

void TAlleleCountFile::writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader){
	outFile.write(reader.chr(),sep,reader.position());
}

void TAlleleCountFile::writeCounts(int count, int numAlleles, int){
	outFile.write("\t",count,"/",numAlleles);
}

void TAlleleCountFile::writeCounts(std::string count, std::string numAlleles, int){
	outFile.write("\t",count,"/",numAlleles);
}

void TAlleleCountFile::endl(){
	outFile.endln();
}

//------------------------------
// TAlleleCountFileWithAlleles
//------------------------------

void TAlleleCountFileWithAlleles::writeHeader(genometools::TPopulationSamples & samples){
	bool useLocusName = parameters().exists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile.write("Locus");
		sep = '_';
	}
	outFile.write("chr",sep,"pos",sep,"ref",sep,"alt");
	for(size_t p=0; p<samples.numPopulations(); p++)
		outFile.write("\t",samples.getPopulationName(p));
	outFile.endln();
}

void TAlleleCountFileWithAlleles::writeHeader(std::vector<std::string> populationNames){
	bool useLocusName = parameters().exists("useLocusName");
	if(useLocusName){
		logfile().list("Will print locus names (rather than chromosome and position).");
		outFile.write("Locus");
		sep = '_';
	}
	outFile.write("chr",sep,"pos",sep, "ref",sep,"alt");
	for(size_t p=0; p<populationNames.size(); p++)
		outFile.write("\t",populationNames[p]);
	outFile.endln();
}

void TAlleleCountFileWithAlleles::writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader){
	outFile.write(reader.chr(),sep,reader.position(),sep,reader.refAllele(),sep,reader.altAllele());
}

//------------------
// TTreeMixFile
//------------------

void TTreeMixFile::writeHeader(genometools::TPopulationSamples & samples){
	outFile.write(samples.getPopulationName(0));
	for(size_t p=1; p<samples.numPopulations(); p++)
		outFile.write(" ",samples.getPopulationName(p));
	outFile.endln();
}

void TTreeMixFile::writeHeader(std::vector<std::string> populationNames){
	outFile.write(populationNames[0]);
		for(size_t p=1; p<populationNames.size(); p++)
			outFile.write(" ",populationNames[p]);
		outFile.endln();
}

void TTreeMixFile::writeCounts(int count, int numAlleles, int populationNum){
	if(populationNum == 0)
		outFile.write(count,",",numAlleles - count);
	else
		outFile.write(" ",count,",",numAlleles - count);
}

void TTreeMixFile::writeCounts(std::string count, std::string numAlleles, int populationNum){
	if(populationNum == 0)
		outFile.write(count,",",coretools::str::fromString<int>(numAlleles) - coretools::str::fromString<int>(count));
	else
		outFile.write(" ",count,",",coretools::str::fromString<int>(numAlleles) - coretools::str::fromString<int>(count));
}

//------------------
// TFlinkFile
//------------------

void TFlinkFile::writeHeader(genometools::TPopulationSamples & samples){
	outFile.write("-\t-");
	for(size_t g=0; g<samples.numPopulations(); ++g){
		outFile.write("\tGroup_A");
	}
	outFile.endln();
	outFile.write("-\t-");
	for(size_t p=0; p<samples.numPopulations(); p++)
		outFile.write("\t",samples.getPopulationName(p));
	outFile.endln();
}

void TFlinkFile::writeHeader(std::vector<std::string> populationNames){
	outFile.write("-\t-");
		for(unsigned int g=0; g<populationNames.size(); ++g){
			outFile.write("\tGroup_A");
		}
		outFile.endln();
		outFile.write("-\t-");
		for(unsigned int p=0; p<populationNames.size(); p++)
			outFile.write("\t",populationNames[p]);
		outFile.endln();
}

void TFlinkFile::writePosition(std::string chr, long pos) { outFile.write(chr, sep, pos); }

void TFlinkFile::writePosition(std::string chr, std::string pos){
	outFile.write(chr,sep,pos);
}

void TFlinkFile::writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader){
	outFile.write(reader.chr(),sep,reader.position());
}

void TFlinkFile::writeCounts(int count, int numAlleles, int){
	outFile.write("\t",count,"/",numAlleles);
}

void TFlinkFile::writeCounts(std::string count, std::string numAlleles, int){
	outFile.write("\t",count,"/",numAlleles);
}

} //end namespace
