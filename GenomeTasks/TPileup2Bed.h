#ifndef GENOMETASKS_TPILEUP2BED_H_
#define GENOMETASKS_TPILEUP2BED_H_

#include "coretools/Files/TInputFile.h"
#include "coretools/Math/TNumericRange.h"
#include "genometools/TBed.h"

namespace GenomeTasks {
	
class TPileup2Bed {
	coretools::TInputFile _pileup;
	std::string _bedName;
	coretools::TNumericRange<size_t> _depthRange;

public:
	TPileup2Bed();
	void run();
		
};

}

#endif
