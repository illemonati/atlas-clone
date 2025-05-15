#ifndef TSIMULATORVARIANTINVARIANTBEDFILES_H_
#define TSIMULATORVARIANTINVARIANTBEDFILES_H_

#include "TSimulatorHaplotypes.h"
#include "coretools/Files/TLineWriter.h"

namespace Simulations {

class TSimulatorVariantInvariantBedFiles {
private:
	coretools::TLineWriter variantSitesFile;
	coretools::TLineWriter invariantSitesFile;

	void openFile(coretools::TLineWriter &file, const std::string filename);
	void close();

public:
	TSimulatorVariantInvariantBedFiles(){}
	TSimulatorVariantInvariantBedFiles(std::string outname) { open(outname); }
	~TSimulatorVariantInvariantBedFiles() { close(); }

	void open(std::string outname);
	void write(TSimulatorHaplotypes &haplotypes, const std::string &chrName);
};
}

#endif
