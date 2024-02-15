#include "TSimulatorHaplotypes.h"

#include "TSimulatorAlleleIndex.h"
#include "TSimulatorReference.h"
#include "coretools/Strings/toString.h"
#include "genometools/GenotypeTypes.h"

namespace Simulations {
using genometools::Base;
using coretools::str::toString;

void TSimulatorHaplotypes::allocateStorage() {
	// allocate storage
	haplotypes.resize(numInd);
	for (size_t ind = 0; ind < numInd; ++ind) {
		haplotypes[ind][0].resize(_length);
		haplotypes[ind][1].resize(_length);
	}
}

void TSimulatorHaplotypes::setLength(size_t length) noexcept {
	if (length > _length) {
		_length = length;
		allocateStorage();
	}
}

void TSimulatorHaplotypes::openTrueGenotypeVCF(std::string filename) {
	// open file
	trueGenoVCF.open(filename.c_str());
	if (!trueGenoVCF) UERROR("Failed to open VCF file '", filename, "' for writing!");

	// write header
	trueGenoVCF << "##fileformat=VCFv4.3\n";
	trueGenoVCF << "##source=ATLAS_Simulator\n";
	trueGenoVCF << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	trueGenoVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for (size_t ind = 0; ind < numInd; ++ind) trueGenoVCF << "\tInd" << ind + 1;
	trueGenoVCF << '\n';
}

const std::array<std::vector<Base>,2>& TSimulatorHaplotypes::get(size_t i) {
	if (i >= numInd)
		UERROR("Haplotypes of individual ", i + 1, " requested, but defined for only ", numInd, " individuals!");
	return haplotypes[i];
}

void TSimulatorHaplotypes::writeTrueGenotypes(const std::string &chrName, const TSimulatorReference &ref) {
	// prepare allele storage
	TSimulatorAlleleIndex index;
	std::string genoString;

	for (size_t l = 0; l < _length; ++l) {
		// chromosome name, position and ID
		trueGenoVCF << chrName << '\t' << l + 1 << "\t.\t";

		// assemble alleles and genotypes
		genoString.clear();
		index.clear(ref[l]);

		// loop over all individuals to figure out which alleles are used
		for (size_t ind = 0; ind < numInd; ++ind) {
			// homozygous or heterozygous?
			if (haplotypes[ind][0][l] == haplotypes[ind][1][l]) {
				// make sure allele exists
				index.add(haplotypes[ind][0][l]);

				// add genotype
				genoString += '\t' + toString(index.index[haplotypes[ind][0][l]]) + '/' +
					      toString(index.index[haplotypes[ind][0][l]]);
			} else {
				// make sure allele exists
				index.add(haplotypes[ind][0][l]);
				index.add(haplotypes[ind][1][l]);

				if (index.index[haplotypes[ind][0][l]] < index.index[haplotypes[ind][1][l]])
					genoString += '\t' + toString(index.index[haplotypes[ind][0][l]]) + '/' +
						      toString(index.index[haplotypes[ind][1][l]]);
				else
					genoString += '\t' + toString(index.index[haplotypes[ind][1][l]]) + '/' +
						      toString(index.index[haplotypes[ind][0][l]]);
			}
		}

		// write ref allele
		index.writeRefAltToVCF(trueGenoVCF);

		// write (no) quality of variant, (no) filter, (no) info and format
		trueGenoVCF << "\t.\t.\t.\tGT";

		// now write genotypes
		trueGenoVCF << genoString << '\n';
	}
}

bool TSimulatorHaplotypes::isPolymoprhic(size_t pos) const noexcept {
	// count how many allele match that of first individual
	const Base testBase = haplotypes[0][0][pos];
	size_t counts    = 0;
	for (auto& h: haplotypes) {
		if (h[0][pos] == testBase) ++counts;
		if (h[1][pos] == testBase) ++counts;
	}
	return counts != 2*numInd;
}
}
