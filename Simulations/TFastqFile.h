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
        std::string _outputName;
        std::vector<coretools::TOutputFile*> _files;
        //coretools::TOutputFile _file;

        uint16_t _readGroupId;
        int _exists;
        std::vector<uint16_t> _readGroups;

        bool _open = false;

        void _writeAlignment(coretools::TOutputFile& file, const BAM::TAlignment &alignment);
        std::string _newName(int readGroupId);
        bool _idExists(uint16_t iD);

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