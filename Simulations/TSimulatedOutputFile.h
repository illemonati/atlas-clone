/*
 * TSimulatedOutputFile.h
 *
 *  Created on: Marc 22, 2023
 *      Author: Michael Jopiti
 */

#ifndef TSIMULATEDOUTPUTFILE_H
#define TSIMULATEDOUTPUTFILE_H

#include <string>
#include "TFastqFile.h"
#include "../coretools/core/coretools/Files/TOutputFile.h"

#include "TBamFile.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <exception>
#include <filesystem>
#include "coretools/Main/TLog.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/Main/TParameters.h"
#include "TSamFlags.h"
#include "api/BamIndex.h"
#include "api/SamProgram.h"
#include "api/SamProgramChain.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"
#include "coretools/Main/globalConstants.h"
#include "coretools/Types/strongTypes.h"

#include "TAlignment.h"

namespace Simulations{

class TSimulatedOutputFile{
    private:
        coretools::TOutputFile simulatedOutputFile;

    public:
        virtual void open(std::string_view filename) = 0;       //following TFile instructions
        virtual void close() = 0;
        virtual void writeAlignment(BAM::TAlignment &alignment) = 0;
        virtual void writeAlignmentLater(const TAlignment & alignment) = 0;

};
}

#endif