#include "TSimulatorVariantInvariantBedFiles.h"

namespace Simulations {
	
void TSimulatorVariantInvariantBedFiles::openFile(gz::ogzstream &file, const std::string filename) {
	file.open(filename.c_str());
	if (!file) UERROR("Failed to open file '", filename, "' for writing!");
}

void TSimulatorVariantInvariantBedFiles::open(std::string outname) {
	// make sure files are closed
	close();

	openFile(variantSitesFile, outname + "_variantSites.bed.gz");
	openFile(invariantSitesFile, outname + "_invariantSites.bed.gz");
}

void TSimulatorVariantInvariantBedFiles::close() {
	if (variantSitesFile) variantSitesFile.close();
	if (invariantSitesFile) invariantSitesFile.close();
}

void TSimulatorVariantInvariantBedFiles::write(TSimulatorHaplotypes &haplotypes, const std::string &chrName) {
	// 0-based
	size_t invariantStart = 0;
	for (size_t l = 0; l < haplotypes.length(); ++l) {
		if (haplotypes.isPolymoprhic(l)) {
			// write invariant
			if (invariantStart < l) invariantSitesFile << chrName << "\t" << invariantStart << "\t" << l << "\n";
			invariantStart = l + 1;

			// write variant
			variantSitesFile << chrName << "\t" << l << "\t" << l + 1 << "\n";
		}
	}

	// write last invariant interval
	if (invariantStart < haplotypes.length())
		invariantSitesFile << chrName << "\t" << invariantStart << "\t" << haplotypes.length() << "\n";
}
}
