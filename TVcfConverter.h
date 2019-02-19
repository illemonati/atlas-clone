/*
 * TVcfConverter.h
 *
 *  Created on: Feb 19, 2019
 *      Author: caduffm
 */

#ifndef TVCFCONVERTER_H_
#define TVCFCONVERTER_H_

#include "TPopulationLikelihoods.h"

class TVcfConverter {
private:
	TLog* logfile;
	int findMaxGenotype(uint8_t* phred);
    void convertToLfmm(TParameters & parameters);

public:
	TVcfConverter(TParameters & Parameters, TLog* Logfile);
    void convertVcf(TParameters & parameters);
};



#endif /* TVCFCONVERTER_H_ */
