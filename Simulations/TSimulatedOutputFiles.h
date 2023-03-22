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

#include <vector>
#include "TSimulatedOutputFile.h"


namespace Simulations{

class TSimulatedOutputFiles{

    private:
        //std::vector<TSimulatedOutputFile> _fileVector;        //need to be a pointer to single files
    
    public:
        TSimulatedOutputFiles();
        void openBAMFile();
        void openFASTQFile();

        //void operator();

};

 }