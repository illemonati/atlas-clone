/*
 * TAlleleCountEstimator.h
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */

#ifndef TALLELECOUNTESTIMATOR_H_
#define TALLELECOUNTESTIMATOR_H_

#include "TAlleleCountFileFormat.h"

//-------------------------------------------------
// TSAFChooseStorage
//-------------------------------------------------
class TSAFChooseStorage{
private:
	int _k;
	double* log_choose_k_j;

public:
	TSAFChooseStorage(int k);
	~TSAFChooseStorage();

	double logChoose(int j);
	double operator[](int j);
};

//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
class TSiteAlleleFrequencyLikelihoods{
private:
	double logOf2;
	int numAlleleCounts; //from 0 to 2k
	int storageSize;
	double* log_alleleFrequencyLikelihoods_h;
	std::map<int, TSAFChooseStorage*> log_choose;

	void updateStorage(int numRequiredAlleleCounts);
	double protectedSumInLog(double a, double b);
	double protectedSumInLog(double a, double b, double c);
	void normalize();

	TSAFChooseStorage* getLogChoose(int counts);
	void fillLog(const TSampleLikelihoods* data, int numSamples, TGlfConverter & glfConverter);
	void fillNatural(const TSampleLikelihoods* data, int numSamples, TGlfConverter & glfConverter);


public:
	TSiteAlleleFrequencyLikelihoods(int numIndividuals);
	~TSiteAlleleFrequencyLikelihoods();
	void fill(const TSampleLikelihoods* data, int numSamples, TGlfConverter & glfConverter);
	void print();
	void write(gz::ogzstream & file);
	int getMLAlleleCount(TRandomGenerator & randomGenerator);
	int getNumAlleles(){ return numAlleleCounts; };
};

//-------------------------------------------------
// TAlleleCountEstimator
//-------------------------------------------------
class TAlleleCountEstimator{
private:
	TLog* logfile;
	//TRandomGenerator* randomGenerator;
	TGlfConverter glfConverter;

public:
	TAlleleCountEstimator(TParameters & params, TLog* Logfile);
	~TAlleleCountEstimator();

	void estimateAlleleCounts(TParameters & params, TRandomGenerator* randomGenerator);
	void writeAlleleFrequencyLikelihoods(TParameters & params);
};


#endif /* TALLELECOUNTESTIMATOR_H_ */

