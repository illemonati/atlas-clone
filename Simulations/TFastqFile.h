/*
 * TFastqFile.h
 *
 *  Created on: Marc 21, 2023
 *      Author: Michael Jopiti
 */

/*
 * 
 *  Create basic File set up for FASTQ file
 * 
*/

#ifndef TFASTQFILE_H_
#define TFASTQFILE_H_

#include <string>
#include <vector>
#include <queue>
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Math/counters.h"
#include "TAlignment.h"
#include "TReadSimulator.h"
#include "TSimulatedOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Files/TFile.h"
#include "coretools/Files/TOutputFile.h"

namespace Simulations { class TSimulatedOutputFile; }

namespace FASTQ{

    class TFastqFile : public Simulations::TSimulatedOutputFile {

    private:
        std::string _fileName;
        coretools::TOutputFile _file;

        bool _open = false;

    public:
        //constructor
        TFastqFile(std::string_view fileName);

        //methods
        void open(std::string_view filename);
        void close();
        void writeAlignment(const BAM::TAlignment &alignment);

        std::string_view getName();
    };

};      //namespace FASTQ

#endif