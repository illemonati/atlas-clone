#include "TSafEstimator.h"
#include "TGlfMultiReader.h"
#include "TSafFile.h"
#include "coretools/Main/TParameters.h"

namespace PopulationTools {
void TSafEstimator::run() {
	GLF::TGlfMultiReader glfReader;
	glfReader.openGLFs();
	glfReader.setAllActive();
	const auto oName = coretools::instances::parameters().get("out", "saf");
	TSafFile safFile(oName, glfReader.numActiveSamples());
	for (auto ids = glfReader.readWindow(); !ids.empty(); ids = glfReader.readWindow()) {
	}
}
}
