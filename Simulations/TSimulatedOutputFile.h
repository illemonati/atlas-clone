/*
 * TSimulatedOutputFile.h
 *
 *  Created on: Marc 22, 2023
 *      Author: Michael Jopiti
 */

#ifndef TSIMULATEDOUTPUTFILE_H
#define TSIMULATEDOUTPUTFILE_H

#include <string>
#include "../coretools/core/coretools/Files/TOutputFile.h"
#include "TAlignment.h"

namespace Simulations {

class TSimulatedOutputFile{
    public:
        TSimulatedOutputFile();   //since it is pure abstract a constructor does not need to exist

        coretools::TOutputFile simulatedOutputFile;

        virtual void open(std::string_view filename) = 0;       //following TFile instructions
        virtual void close() = 0;
        virtual void writeAlignment(const BAM::TAlignment &alignment) = 0;
        
        //virtual void writeAlignmentLater(const BAM::TAlignment & alignment) = 0; 
            //not pure virtual because of some random compile time error
};

};      //namespace Simulations

#endif