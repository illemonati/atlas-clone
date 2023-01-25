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
#include "genometools/TFastaWriter.h"

namespace PopulationTools{

//-------------------------------------------------
// TAncestralAllelEstimator
//-------------------------------------------------    
class TAncestralAllelEstimator{
private:
    TAlleleCountReader alleleCounts;
    TFastaIndex fastaIndex;
public:
	void readAlleleCounts();
    void readFastaIndex();
    void printFasta();
};


class TTask_estimateAncestralAlleles:public coretools::TTask{
public:
	TTask_estimateAncestralAlleles(){ _explanation = "Writing FASTA-file with ancestral alleles"; };

	void run(){
		TAncestralAlleleEstimator ancestralAlleleEst;
		//ancestralAlleleEst.function();
	};
};

}; //end namespace

#endif /* TANCESTRALALLELETESTIMATOR_H_ */