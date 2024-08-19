#include "TGLFPrinter.h"
#include "genometools/GLF/TGLFReader.h"
#include "coretools/Main/TParameters.h"

namespace PopulationTools {
void TGLFPrinter::run() {
	genometools::TGLFReader reader(coretools::instances::parameters().get<std::string>("glf"));
	reader.printToEnd();
}
} // namespace PopulationTools
