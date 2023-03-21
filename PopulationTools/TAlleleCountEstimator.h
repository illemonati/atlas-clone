/*
 * TAlleleCountEstimator.h
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */

#ifndef TALLELECOUNTESTIMATOR_H_
#define TALLELECOUNTESTIMATOR_H_

#include <stdint.h>
#include <map>
#include <string>
#include <vector>

#include "coretools/Main/TParameters.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/TSampleLikelihoods.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"

namespace PopulationTools { class TAlleleCountFile; }
namespace gz { class ogzstream; }

namespace PopulationTools{

using TSampleLikelihoods = genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability>;


//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
class TSiteAlleleFrequencyLikelihoods{
private:
	static constexpr double logOf2 = 0.6931471805599453;
	uint32_t numAlleleCounts; //from 0 to 2k
	std::vector<coretools::LogProbability> log_alleleFrequencyLikelihoods_h;
	std::map<int, std::vector<double>> log_choose;

	void normalize();

	const std::vector<double>& _getLogChoose(int counts);

protected:
    void _fillLog(TSampleLikelihoods* data, uint32_t numSamples);
    void _fillNatural(TSampleLikelihoods* data, uint32_t numSamples);

public:
	TSiteAlleleFrequencyLikelihoods(int numIndividuals);
	void fill(TSampleLikelihoods* data, uint32_t numSamples);
	void print();
	void write(gz::ogzstream & file);
	int getMLAlleleCount();
	int getNumAlleles(){ return numAlleleCounts; };
	const std::vector<coretools::LogProbability> & getLogAlleleFrequencyLikelihoods() const;
};

struct TAlleleCounter {
	void run();
};

}; //end namespace

#endif /* TALLELECOUNTESTIMATOR_H_ */
