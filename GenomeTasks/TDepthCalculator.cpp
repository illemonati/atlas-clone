#include "TDepthCalculator.h"
#include "TCigar.h"
#include "api/BamReader.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {
using coretools::instances::logfile;

void TDepthCalculator::run() {
	const auto depth = _genome.bamFile().averageDepth();
	logfile().list("Average Depth = ", depth);
}


}
