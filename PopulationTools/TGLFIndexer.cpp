#include "TGLFIndexer.h"

#include "genometools/GLF/TGLFReader.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace PopulationTools {
using coretools::instances::parameters;
using coretools::instances::logfile;

void TGLFIndexer::run() {
	// open GLF for reading, do not open index
	const auto glf = parameters().get<std::string>("glf");
	if (parameters().exists("check")) {
		logfile().list("Checking index-file of ", glf);

		for (genometools::TGLFReader reader(glf, true); !reader.empty(); reader.popFront()) {}
		logfile().list("GLF-file is OK!");
	} else {
		logfile().list("Indexing GLF-file: ", glf);

		genometools::TGLFReader reader(glf, false);
		reader.writeIndex();
	}
}
} // namespace PopulationTools
