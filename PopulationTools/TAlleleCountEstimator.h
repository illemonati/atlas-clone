/*
 * TAlleleCountEstimator.h
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */

#ifndef TALLELECOUNTESTIMATOR_H_
#define TALLELECOUNTESTIMATOR_H_

#include "TAlleleCountFileFormat.h"
#include "probability.h"
#include <string>
#include <map>

namespace PopulationTools{

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
	uint32_t numAlleleCounts; //from 0 to 2k
	std::vector<coretools::LogProbability> log_alleleFrequencyLikelihoods_h;
	std::map<int, TSAFChooseStorage*> log_choose;

	coretools::LogProbability _protectedSumInLog(const coretools::LogProbability a, const coretools::LogProbability c);
	coretools::LogProbability _protectedSumInLog(const coretools::LogProbability a, const coretools::LogProbability b, const coretools::LogProbability c);
	void normalize();

	TSAFChooseStorage* _getLogChoose(int counts);

protected:
    void _fillLog(TSampleLikelihoods* data, uint32_t numSamples);
    void _fillNatural(TSampleLikelihoods* data, uint32_t numSamples);

public:
	TSiteAlleleFrequencyLikelihoods(int numIndividuals);
	~TSiteAlleleFrequencyLikelihoods();
	void fill(TSampleLikelihoods* data, uint32_t numSamples);
	void print();
	void write(gz::ogzstream & file);
	int getMLAlleleCount(coretools::TRandomGenerator & randomGenerator);
	int getNumAlleles(){ return numAlleleCounts; };
	const std::vector<coretools::LogProbability> & getLogAlleleFrequencyLikelihoods() const;
};

//-------------------------------------------------
// TAlleleCountEstimator
//-------------------------------------------------
class TAlleleCountEstimator{
private:
	coretools::TLog* logfile;

public:
	TAlleleCountEstimator(coretools::TParameters & params, coretools::TLog* Logfile);
	~TAlleleCountEstimator();

	void estimateAlleleCounts(coretools::TParameters & params, coretools::TRandomGenerator* randomGenerator);
	void writeAlleleFrequencyLikelihoods(coretools::TParameters & params);
	TAlleleCountFile* prepareOutputFile(std::string type, std::string filePrefix, coretools::TParameters& params);
	void transformFormat(coretools::TParameters & params);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimateAlleleCounts:public coretools::TTask{
public:
	TTask_estimateAlleleCounts(){ _explanation = "Estimating population allele counts"; };

	void run(){
		using namespace coretools::instances;
		TAlleleCountEstimator alleleCountEst(parameters(), &logfile());
		alleleCountEst.estimateAlleleCounts(parameters(), &randomGenerator());
	};
};

class TTask_writeAlleleCountLikelihoods:public coretools::TTask{
public:
	TTask_writeAlleleCountLikelihoods(){ _explanation = "Writing sample allele count likelihoods"; };

	void run(){
		using namespace coretools::instances;
		TAlleleCountEstimator alleleCountEst(parameters(), &logfile());
		alleleCountEst.writeAlleleFrequencyLikelihoods(parameters());
	};
};

class TTask_transformAlleleCountFormat:public coretools::TTask{
public:
	TTask_transformAlleleCountFormat(){ _explanation = "Transforming allele counts file format"; };

	void run(){
		using namespace coretools::instances;
		TAlleleCountEstimator alleleCountEst(parameters(), &logfile());
		alleleCountEst.transformFormat(parameters());
	};
};

}; //end namespace

#endif /* TALLELECOUNTESTIMATOR_H_ */
