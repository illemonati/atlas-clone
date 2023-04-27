/*
 * TAlleleCountEstimator.h
 *
 *  Created on: Jan 25, 2023
 *      Author: raphael
 */

#ifndef TANCESTRALALLELETESTIMATOR_H_
#define TANCESTRALALLELETESTIMATOR_H_

#include <stdint.h>
#include <map>
#include <string>
#include <vector>

#include "coretools/Main/TTask.h"
#include "TAlleleCountReader.h"
#include "genometools/TFastaIndex.h"

namespace PopulationTools{

//-------------------------------------------------
// TAncestralAllelEestimator
//-------------------------------------------------    
class TAncestralAlleleEstimator{
private:
    TAlleleCountReader _alleleCounts;
    genometools::TFastaIndex _fastaIndex;	
	std::string filename;

	size_t _minorCountMax;
	size_t _totalCountMin;
public:
	TAncestralAlleleEstimator();
    void run();
};

} // namespace PopulationTools

#endif /* TANCESTRALALLELETESTIMATOR_H_ */
