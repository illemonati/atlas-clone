/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"
#include "TGLFReader.h"
#include "coretools/Main/TParameters.h"

namespace GLF {

void TGLFPrinter::run() {
	TGLFReader reader(coretools::instances::parameters().get<std::string>("glf"));
	reader.printToEnd();
}

void TGLFIndexer::run(){
	// open GLF for reading, do not open index
	TGLFReader reader(coretools::instances::parameters().get<std::string>("glf"), false);
	reader.writeIndex();
}

}; // namespace GLF
