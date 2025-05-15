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

#include "coretools/Files/TLineWriter.h"
#include "genometools/VCF/TPopulationLikelihoods.h"

namespace genometools { class TPopulationSamples; }

namespace PopulationTools{

class TAlleleCountFile{
protected:
	std::string filename;
	coretools::TLineWriter outFile;
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
	virtual ~TAlleleCountFile() = default;
};

class TAlleleCountFileWithAlleles:public TAlleleCountFile{
public:
	void writeHeader(genometools::TPopulationSamples & samples) override;
	void writeHeader(std::vector<std::string> populationNames) override;
	void writePosition(std::string, long) override{ DEVERROR("Need to provide alleles for this format!"); }
	void writePosition(std::string, std::string) override{ DEVERROR("Need to provide alleles for this format!"); }
	void writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader) override;

	TAlleleCountFileWithAlleles(std::string Filename):TAlleleCountFile(Filename){};
};

class TTreeMixFile:public TAlleleCountFile{
public:
	void writeHeader(genometools::TPopulationSamples & samples) override;
	void writeHeader(std::vector<std::string> populationNames) override;

	// do nothing, treemix does not need position
	void writePosition(std::string, long) override {}
	void writePosition(std::string, std::string) override {}
	void writePosition(const genometools::TPopulationLikelihoodReaderLocus &) override{};

	void writeCounts(int count, int numAlleles, int populationNum) override;
	void writeCounts(std::string count, std::string numAlleles, int populationNum) override;

	TTreeMixFile(std::string Filename):TAlleleCountFile(Filename){};

};

class TFlinkFile:public TAlleleCountFile{
public:
	void writeHeader(genometools::TPopulationSamples & samples) override;
	void writeHeader(std::vector<std::string> populationNames) override;
	void writePosition(std::string chr, long pos) override;
	void writePosition(std::string chr, std::string pos) override;
	void writePosition(const genometools::TPopulationLikelihoodReaderLocus & reader) override;
	void writeCounts(int count, int numAlleles, int populationNum) override;
	void writeCounts(std::string count, std::string numAlleles, int populationNum) override;


	TFlinkFile(std::string Filename):TAlleleCountFile(Filename){};

};

}; //end namespace

#endif /* POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_ */
