#include <string>
#include <vector>
#include "TFastqFile.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Math/counters.h"
#include "TAlignment.h"
#include "TReadSimulator.h"


namespace FASTQ{

std::string TFastqFile::_generateIlluminaHeader(std::string machineID, short runID, short FlowCellID, short lane, short tile, unsigned short xCoordinate, unsigned short yCoordinate,
								bool readDirection, bool passFilter, short controlBits, std::string barCodeSequence){
    _cPassFilter = (passFilter) ? _cPassFilter = 'Y': _cPassFilter = 'N'; 

	return "@" + machineID + ":" + toString(runID) + ":" + toString(FlowCellID) + ":" + toString(lane) + ":" + toString(tile) + ":" + toString(xCoordinate) + ":" + toString(yCoordinate) + " " +
			_cPassFilter + ":" + toString(controlBits) + ":" + toString(barCodeSequence);
}

void TOutputFASTQFile::_writeSequence(){
}



}