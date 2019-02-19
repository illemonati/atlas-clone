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

    void convertToLfmm(TParameters & parameters, TPopulationSamples & samples, TPopulationLikelihoodReader & reader);
    void prepareReadingVcf(TParameters & parameters, TPopulationSamples & samples, TPopulationLikelihoodReader & reader);
    std::ofstream openOutputFile(TParameters & parameters, std::string fileExtension);


public:
	TVcfConverter(TParameters & Parameters, TLog* Logfile);
    void convertVcf(TParameters & parameters);
};



#endif /* TVCFCONVERTER_H_ */
