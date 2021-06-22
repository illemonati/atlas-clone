/*
 * TAlleleCountFileFormat.h
 *
 *  Created on: Oct 23, 2019
 *      Author: vivian
 */

#ifndef POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_
#define POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_

#include "mathFunctions.h"
#include "TPopulationLikelihoods.h"
#include "TFile.h"
#include <string>

namespace PopulationTools{

class TAlleleCountFile{
protected:
	std::string filename;
	gz::ogzstream outFile;
	std::string sep;


public:
	void openFileToWrite(std::string filename);
	virtual void writeHeader(TPopulationSamples & samples, coretools::TParameters & params, coretools::TLog* logfile);
	virtual void writeHeader(std::vector<std::string> populationNames, coretools::TParameters & params, coretools::TLog* logfile);
	virtual void writePosition(std::string chr, long pos);
	virtual void writePosition(std::string chr, std::string pos);
	virtual void writeCounts(int count, int numAlleles, int populationNum);
	virtual void writeCounts(std::string count, std::string numAlleles, int populationNum);
	void endl();

	TAlleleCountFile(std::string Filename);
	virtual ~TAlleleCountFile(){};
};

class TTreeMixFile:public TAlleleCountFile{
public:
	void writeHeader(TPopulationSamples & samples, coretools::TParameters & params, coretools::TLog* logfile);
	void writeHeader(std::vector<std::string> populationNames, coretools::TParameters & params, coretools::TLog* logfile);
	void writePosition(std::string chr, long pos);
	void writePosition(std::string chr, std::string pos);
	void writeCounts(int count, int numAlleles, int populationNum);
	void writeCounts(std::string count, std::string numAlleles, int populationNum);

	TTreeMixFile(std::string Filename);
	~TTreeMixFile(){};

};

class TFlinkFile:public TAlleleCountFile{
public:
	void writeHeader(TPopulationSamples & samples, coretools::TParameters & params, coretools::TLog* logfile);
	void writeHeader(std::vector<std::string> populationNames, coretools::TParameters & params, coretools::TLog* logfile);
	void writePosition(std::string chr, long pos);
	void writePosition(std::string chr, std::string pos);
	void writeCounts(int count, int numAlleles, int populationNum);
	void writeCounts(std::string count, std::string numAlleles, int populationNum);


	TFlinkFile(std::string Filename);
	~TFlinkFile(){};

};

}; //end namespace

#endif /* POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_ */
