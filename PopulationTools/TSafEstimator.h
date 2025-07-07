#ifndef TSAFESTIMATOR_H_
#define TSAFESTIMATOR_H_

#include "coretools/Containers/TView.h"
#include "genometools/GLF/TMultiGLFTraverser.h"
#include "genometools/TFastaReader.h"

namespace PopulationTools {

class TSafEstimator {
	std::vector<double> _logProbs;
	genometools::TMultiGLFTraverser _glfTraverser;
	genometools::TFastaReader _fasta;

	size_t _lower = 0;

	void _iterate(coretools::TConstView<genometools::TGLFEntry> data, genometools::Base major);
public:
	TSafEstimator();
	void run();
};

} // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
