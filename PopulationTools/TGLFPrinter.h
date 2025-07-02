#ifndef TGLFPRINTER_H_
#define TGLFPRINTER_H_

#include "coretools/Main/TParameters.h"
namespace PopulationTools {
class TGLFPrinter {
	std::string _glfName;
	std::string _bedName;
public:
	TGLFPrinter()
		: _glfName(coretools::instances::parameters().get("glf")),
		  _bedName(coretools::instances::parameters().get("bed", "")) {}
	void run();
};
} // namespace PopulationTools

#endif
