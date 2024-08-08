/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"
#include "TGLFReader.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace GLF {
using coretools::instances::parameters;
using coretools::instances::logfile;

void TGLFPrinter::run() {
	TGLFReader reader(parameters().get<std::string>("glf"));
	reader.printToEnd();
}

void TGLFIndexer::run() {
	// open GLF for reading, do not open index
	const auto glf = parameters().get<std::string>("glf");
	if (parameters().exists("check")) {
		logfile().list("Checking index-file of ", glf);
		TGLFReader reader(glf, true);
		while (reader.readNext()) {}
		logfile().list("GLF-file is OK!");
	} else {
		logfile().list("Indexing GLF-file: ", glf);
		TGLFReader reader(glf, false);
		reader.writeIndex();
	}
}

}; // namespace GLF
