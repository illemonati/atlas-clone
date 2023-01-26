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

	uint32_t _minorCountMin;
	uint32_t _totalCountMin;
public:
	TAncestralAlleleEstimator();
    void printFasta();
};


class TTask_estimateAncestralAlleles:public coretools::TTask{
public:
	TTask_estimateAncestralAlleles(){ _explanation = "Writing FASTA-file with ancestral alleles"; };

	void run(){
		TAncestralAlleleEstimator ancestralAlleleEst;
		ancestralAlleleEst.printFasta();
	};
};

}; //end namespace

#endif /* TANCESTRALALLELETESTIMATOR_H_ */