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
private:
	std::string filename;
	gz::ogzstream outFile;
	char sep;

public:
	TAlleleCountFile(std::string Filename);
	void writeHeader(TPopulationSamples & samples, TParameters & params, TLog* logfile);
	void writePosition(std::string chr, long pos);
	void writeCounts(int count, int numAlleles);
	void endl();

};



#endif /* POPULATIONTOOLS_TALLELECOUNTFILEFORMAT_H_ */
