#ifndef TSAFESTIMATOR_H_
#define TSAFESTIMATOR_H_

#include "genometools/GLF/TGLFMultiReader.h"

namespace PopulationTools {

class TSafEstimator {
	std::vector<double> _logProbs;
	genometools::TGLFMultiReader _glfReader;
	size_t _lower = 0;

	void _iterate(const std::vector<genometools::TGLFEntry>& data, genometools::Base major);
public:
	TSafEstimator();
	void run();
};

} // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
