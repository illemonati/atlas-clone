#ifndef TGLFPRINTER_H_
#define TGLFPRINTER_H_

#include "coretools/Main/TParameters.h"
#include "genometools/GLF/TSingleGLFTraverser.h"
namespace PopulationTools {
class TGLFPrinter {
	genometools::TSingleGLFTraverser _traverser{};
public:
	void run();
};
} // namespace PopulationTools

#endif
