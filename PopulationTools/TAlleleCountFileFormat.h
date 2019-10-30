/*
 * TAlleleCountFileFormat.h
 *
 *  Created on: Oct 23, 2019
 *      Author: vivian
 */

#ifndef POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_
#define POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_

#include "../mathFunctions.h"
#include "TPopulationLikelihoods.h"
#include "../TFile.h"

class TAlleleCountFile{
protected:
	std::string filename;
	gz::ogzstream outFile;
	std::string sep;

	void openFileToWrite(std::string filename);

public:
	virtual void writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile);
	virtual void writePosition(std::string chr, long pos);
	virtual void writeCounts(int count, int numAlleles, int populationNum);
	void endl();

	TAlleleCountFile(std::string Filename);
	virtual ~TAlleleCountFile(){};
};

class TTreeMixFile:TAlleleCountFile{
public:
	void writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile);
	void writePosition(std::string chr, long pos);
	void writeCounts(int count, int numAlleles, int populationNum);

	TTreeMixFile(std::string Filename);
	~TTreeMixFile(){};

};

class TFlinkFile:TAlleleCountFile{
public:
	void writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile);
	void writePosition(std::string chr, long pos);
	void writeCounts(int count, int numAlleles, int populationNum);

	TFlinkFile(std::string Filename);
	~TFlinkFile(){};

};


#endif /* POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_ */
