#include "TSimulatorHaplotypes.h"

#include "TSimulatorAlleleIndex.h"
#include "coretools/Strings/toString.h"

namespace Simulations {
using genometools::Base;
using coretools::str::toString;

void TSimulatorHaplotypes::openTrueGenotypeVCF(std::string_view Filename) {
	// open file
	_trueGenoVCF.open(Filename);

	// write header
	_trueGenoVCF.writeln("##fileformat=VCFv4.3");
	_trueGenoVCF.writeln("##source=ATLAS_Simulator");
	_trueGenoVCF.writeln("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	_trueGenoVCF.write("#CHROM","POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT");
	for (size_t ind = 0; ind < size(); ++ind) _trueGenoVCF.writeNoDelim("Ind", ind + 1);

	_trueGenoVCF.endln();
}

void TSimulatorHaplotypes::writeTrueGenotypes(std::string_view ChrName, coretools::TView<genometools::Base> Ref) {
	// prepare allele storage
	TSimulatorAlleleIndex index;
	std::string genoString;

	for (size_t l = 0; l < length(); ++l) {
		if (Ref[l] == Base::N) continue; // skip
		// chromosome name, position and ID
		_trueGenoVCF.write(ChrName, l + 1,".");

		// assemble alleles and genotypes
		genoString.clear();
		index.clear(Ref[l]);

		// loop over all individuals to figure out which alleles are used
		for (size_t ind = 0; ind < size(); ++ind) {
			// homozygous or heterozygous?
			const auto haplo1 = first(_haplotypes[ind][l]);
			const auto haplo2 = second(_haplotypes[ind][l]);
			if (haplo1 == haplo2) {
				// make sure allele exists
				index.add(haplo1);

				// add genotype
				genoString += toString(index.index[haplo1], "/", index.index[haplo1]);
			} else {
				// make sure allele exists
				index.add(haplo1);
				index.add(haplo2);

				if (index.index[haplo1] < index.index[haplo2])
					genoString += toString(index.index[haplo1], "/", index.index[haplo2]);
				else
					genoString += toString(index.index[haplo2], "/", index.index[haplo1]);
			}
		}

		// write ref allele
		index.writeRefAltToVCF(_trueGenoVCF);

		// write (no) quality of variant, (no) filter, (no) info and format, genoString
		_trueGenoVCF.writeln(".", ".",".", "GT", genoString);
	}
}

bool TSimulatorHaplotypes::isPolymoprhic(size_t pos) const noexcept {
	// count how many allele match that of first individual
	const Base testBase = first(_haplotypes[0][pos]);
	for (auto& h: _haplotypes) {
		if (first(h[pos]) != testBase) return true;
		if (second(h[pos]) != testBase) return true;
	}
	return false;
}
}
