#ifndef TSIMULATORVARIANTINVARIANTBEDFILES_H_
#define TSIMULATORVARIANTINVARIANTBEDFILES_H_

#include "TSimulatorHaplotypes.h"
#include "coretools/Files/gzstream.h"

namespace Simulations {

class TSimulatorVariantInvariantBedFiles {
private:
	gz::ogzstream variantSitesFile;
	gz::ogzstream invariantSitesFile;

	void openFile(gz::ogzstream &file, const std::string filename);
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
