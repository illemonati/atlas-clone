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
	std::vector<coretools::LogProbability> log_alleleFrequencyLikelihoods_h;
	std::vector<std::vector<double>> log_choose;

	const std::vector<double>& _getLogChoose(size_t counts);

protected:
    void _fillLog(TSampleLikelihoods* data, uint32_t numSamples);
    void _fillNatural(TSampleLikelihoods* data, uint32_t numSamples);

public:
	void fill(TSampleLikelihoods* data, uint32_t numSamples);
	void write(gz::ogzstream & file);
	size_t MLAlleleCount();
	size_t Nalleles(){ return log_alleleFrequencyLikelihoods_h.size() - 1; };
	const std::vector<coretools::LogProbability> & getLogAlleleFrequencyLikelihoods() const;
};

struct TAlleleCounter {
	void run();
};

}; //end namespace

#endif /* TALLELECOUNTESTIMATOR_H_ */
