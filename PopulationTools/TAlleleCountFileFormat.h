/*
 * TAlleleCountFileFormat.h
 *
 *  Created on: Oct 23, 2019
 *      Author: vivian
 */

#ifndef POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_
#define POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_

#include <string>
#include <vector>

#include "coretools/Files/gzstream.h"
#include "genometools/VCF/TPopulationLikelihoods.h"

namespace genometools { class TPopulationSamples; }

namespace PopulationTools{

class TAlleleCountFile{
protected:
	std::string filename;
	gz::ogzstream outFile;
	std::string sep;


public:
	void openFileToWrite(std::string filename);
	virtual void writeHeader(genometools::TPopulationSamples & samples);
	virtual void writeHeader(std::vector<std::string> populationNames);
	virtual void writePosition(std::string chr, long pos);
	virtual void writePosition(std::string chr, std::string pos);
	virtual void writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader);
	virtual void writeCounts(int count, int numAlleles, int populationNum);
	virtual void writeCounts(std::string count, std::string numAlleles, int populationNum);
	void endl();

	TAlleleCountFile(std::string Filename);
	virtual ~TAlleleCountFile(){};
};

class TAlleleCountFileWithAlleles:public TAlleleCountFile{
public:
	virtual void writeHeader(genometools::TPopulationSamples & samples);
	virtual void writeHeader(std::vector<std::string> populationNames);
	virtual void writePosition(std::string, long){ DEVERROR("Need to provide alleles for this format!"); }
	virtual void writePosition(std::string, std::string){ DEVERROR("Need to provide alleles for this format!"); }
	virtual void writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader);

	TAlleleCountFileWithAlleles(std::string Filename):TAlleleCountFile(Filename){};
	~TAlleleCountFileWithAlleles(){};
};

class TTreeMixFile:public TAlleleCountFile{
public:
	void writeHeader(genometools::TPopulationSamples & samples);
	void writeHeader(std::vector<std::string> populationNames);

	// do nothing, treemix does not need position
	void writePosition(std::string, long) {}
	void writePosition(std::string, std::string) {}
	void writePosition(const genometools::TPopulationLikelihoodReaderLocus &){};

	void writeCounts(int count, int numAlleles, int populationNum);
	void writeCounts(std::string count, std::string numAlleles, int populationNum);

	TTreeMixFile(std::string Filename):TAlleleCountFile(Filename){};
	~TTreeMixFile(){};

};

class TFlinkFile:public TAlleleCountFile{
public:
	void writeHeader(genometools::TPopulationSamples & samples);
	void writeHeader(std::vector<std::string> populationNames);
	void writePosition(std::string chr, long pos);
	void writePosition(std::string chr, std::string pos);
	void writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader);
	void writeCounts(int count, int numAlleles, int populationNum);
	void writeCounts(std::string count, std::string numAlleles, int populationNum);


	TFlinkFile(std::string Filename):TAlleleCountFile(Filename){};
	~TFlinkFile(){};

};

}; //end namespace

#endif /* POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_ */
