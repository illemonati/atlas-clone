/*
 * TSimulatedOutputFiles.h
 *
 *  Created on: Marc 22, 2023
 *      Author: Michael Jopiti
 */

/**
 * 
 * Create structure for working with BAM or FASTQ files:
 *  - create vector populated from one of those two file types
 *  - create empty constructor
 *  - create two different open() methods (openFASTQ() and openBAM())
 *  - still have operators for the twos
 * 
*/

#ifndef TSimulatedOutputFiles_H_
#define TSimulatedOutputFiles_H_

#include <vector>
#include "TSimulatedOutputFile.h"
//#include "TAlignment.h"


namespace Simulations{

class TSimulatedOutputFiles{

    private:
        //std::vector<TSimulatedBamOutputFile> _BamFileVector;          //Don't need it since BAM is a single file of multiple alignment  + already used
                                                                        //in TBAMSimulator!!          
        std::vector<FASTQ::TFastqFile> _FastqFileVector;        //need to be a pointer to single files
    
    public:
        TSimulatedOutputFiles();
        void openBAMFile();
        void openFASTQFile();

        //void operator();

};

};      //namespace Simulations

#endif

