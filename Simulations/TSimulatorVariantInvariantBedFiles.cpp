#include "TSimulatorVariantInvariantBedFiles.h"
#include "coretools/Main/TError.h"

namespace Simulations {
	
void TSimulatorVariantInvariantBedFiles::openFile(coretools::TLineWriter &file, const std::string filename) {
	file.open(filename);
}

void TSimulatorVariantInvariantBedFiles::open(std::string outname) {
	// make sure files are closed
	close();

	openFile(variantSitesFile, outname + "_variantSites.bed.gz");
	openFile(invariantSitesFile, outname + "_invariantSites.bed.gz");
}

void TSimulatorVariantInvariantBedFiles::close() {
	if (variantSitesFile.isOpen()) variantSitesFile.close();
	if (invariantSitesFile.isOpen()) invariantSitesFile.close();
}

void TSimulatorVariantInvariantBedFiles::write(TSimulatorHaplotypes &haplotypes, const std::string &chrName) {
	// 0-based
	size_t invariantStart = 0;
	for (size_t l = 0; l < haplotypes.length(); ++l) {
		if (haplotypes.isPolymoprhic(l)) {
			// write invariant
			if (invariantStart < l) invariantSitesFile.writeln(chrName,"\t",invariantStart,"\t",l);
			invariantStart = l + 1;

			// write variant
			variantSitesFile.writeln(chrName,"\t",l,"\t",l + 1);
		}
	}

	// write last invariant interval
	if (invariantStart < haplotypes.length())
		invariantSitesFile.writeln(chrName,"\t",invariantStart,"\t",haplotypes.length());
}
}
