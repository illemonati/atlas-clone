/*
 * TSimulatorAuxiliaryTools.cpp
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */

#include "TSimulatorAuxiliaryTools.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>

#include "TChromosomes.h"
#include "TLog.h"
#include "TReadGroups.h"
#include "TSamHeader.h"
#include "globalConstants.h"
#include "probability.h"
#include "stringFunctions.h"
#include "weakTypes.h"

namespace Simulations {

using coretools::str::toString;
using coretools::instances::logfile;
using genometools::Base;

//---------------------------------------------------
// TSimulatorReference
//---------------------------------------------------

TSimulatorReference::TSimulatorReference(std::string Filename){
	open(Filename);
}

void TSimulatorReference::open(std::string Filename){
	_filename = Filename;
	// open FASTA file for reference sequences
	logfile().list("Will write reference sequence to '" + _filename + "'.");
	_fasta.open(_filename.c_str());
	if (!_fasta) throw "Failed to open file '" + _filename + "' for writing!";
	_filename += ".fai";
	_fastaIndex.open(_filename.c_str());
	if (!_fastaIndex) throw "Failed to open file '" + _filename + "' for writing!";
	_oldOffset = 0;
}

TSimulatorReference::~TSimulatorReference() {
	if (_chrName != "" && _needsWriting) _writeRefToFasta();
	_fasta.close();
	_fastaIndex.close();
}

void TSimulatorReference::_writeRefToFasta() {
	// write to fasta
	_fasta << ">" << _chrName;
	for (int l = 0; l < _chrLength; ++l) {
		if (l % 70 == 0) _fasta << "\n";
		_fasta << _ref[l];
	}
	_fasta << "\n";

	// add to index
	std::string tmp = _chrName;
	_oldOffset += _chrName.size() + 2;
	_fastaIndex << coretools::str::extractBeforeWhiteSpace(tmp) << "\t" << _chrLength << "\t" << _oldOffset
				<< "\t70\t71\n";
	_oldOffset += _chrLength + (int)(_chrLength / 70);
	if (_chrLength % 70 != 0) _oldOffset += 1;

	_needsWriting = false;
}

void TSimulatorReference::_allocateStorage(long length) {
	// allocate storage
	_ref.resize(length);
}

void TSimulatorReference::setChr(std::string ChrName, long ChrLength) {
	if(!_fasta){
		DEVERROR("Fasta file not opened yet!");
	}
	// write if not yet written
	if (_chrName != "" && _needsWriting) _writeRefToFasta();

	// move to new chr
	_chrName   = ChrName;
	_chrLength = ChrLength;
	_ref.resize(_chrLength);
	_needsWriting = true;
}

// void TSimulatorReference::simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float*
// cumulBaseFreq){ 	logfile->listFlush("Simulating reference alleles ...");
//
//	if(!storageInitialized)
//		throw "Can not simulate reference sequence, no chromosome set!";
//
//	//simulate reference sequence
//	for(int l=0; l<chrLength; ++l){
//		ref[l] = static_cast<Base>(randomGenerator->pickOne(4, cumulBaseFreq));
//	}
//
//	if(fastaOpen){
//		writeRefToFasta();
//	};
// };

//---------------------------------------------------
// TSimulatorBamFile
//---------------------------------------------------
void TSimulatorBamFile::open(const std::string &Filename, const std::string &SampleName, BAM::TReadGroups ReadGroups,
			     const genometools::TChromosomes &Chromosomes) {
	logfile().listFlush("Opening BAM file '" + Filename + "' ...");

	if (_outBam.isOpen()) throw "A BAM file is already open for writing!";

	if (Chromosomes.size() < 1) throw "Can not open a BAM file without specified chromosomes!";

	// create header, read group and chromosome objects
	for (auto &rg : ReadGroups) {
		rg.sequencingCenter_CN =
			coretools::__GLOBAL_APPLICATION_NAME__ + " " + coretools::__GLOBAL_APPLICATION_VERSION__;
		rg.description_DS          = "Simulated with commit " + coretools::__GLOBAL_APPLICATION_COMMIT__;
		rg.sample_SM               = SampleName;
		rg.sequencingTechnology_PL = "ILLUMINA";
	}
	const BAM::TSamHeader header("1.6", "coordinate", "none");
	_outBam.open(Filename, header, Chromosomes, ReadGroups);
	logfile().done();
}
TSimulatorBamFile::~TSimulatorBamFile() { _outBam.closeNoIndex(); }

void TSimulatorBamFile::close() {
	_outBam.close(&logfile());
}

TSimulatorBamFiles::TSimulatorBamFiles(uint32_t NumFiles, const std::string & Outname, std::vector<TReadSimulators> & ReadSimulators,
				       const genometools::TChromosomes &Chromosomes) {
	if (NumFiles < 1) DEVERROR("Can not open less than one BAM file!");
	_files.resize(NumFiles);

	if(ReadSimulators.size() > 1 && ReadSimulators.size() != NumFiles){
		DEVERROR("Numbe rof read simulators does not match number of files!");
	}

	// open BAM files
	if (_files.size() == 1) {
		_files[0].open(Outname + ".bam", "Ind1", ReadSimulators.front().readGroups(), Chromosomes);
	} else {
		logfile().startIndent("Opening ", _files.size(), " BAM files:");
		for (size_t i = 0; i < _files.size(); ++i) {
			std::string filename = Outname + "_ind" + toString(i + 1) + ".bam";
			_files[i].open(filename, "Ind" + toString(i + 1), ReadSimulators[i].readGroups(), Chromosomes);
		}
		logfile().endIndent();
	}
}

TSimulatorBamFiles::~TSimulatorBamFiles() {
	logfile().startIndent("Indexing BAM files:");
	for (auto &f : _files) { f.close(); }
	logfile().endIndent();
}

TSimulatorBamFile &TSimulatorBamFiles::operator[](size_t i) {
	if (i >= _files.size()) throw "BAM file " + toString(i) + " does not exist!";
	return _files[i];
}

//---------------------------------------------------------
// TSimulatorHaplotypes
//---------------------------------------------------------
void TSimulatorHaplotypes::allocateStorage() {
	// allocate storage
	haplotypes.resize(numInd);
	for (size_t ind = 0; ind < numInd; ++ind) {
		haplotypes[ind][0].resize(_length);
		haplotypes[ind][1].resize(_length);
	}
}

void TSimulatorHaplotypes::setLength(uint32_t length) noexcept {
	if (length > _length) {
		_length = length;
		allocateStorage();
	}
}

void TSimulatorHaplotypes::openTrueGenotypeVCF(std::string filename) {
	// open file
	trueGenoVCF.open(filename.c_str());
	if (!trueGenoVCF) throw "Failed to open VCF file '" + filename + "' for writing!";

	// write header
	trueGenoVCF << "##fileformat=VCFv4.3\n";
	trueGenoVCF << "##source=ATLAS_Simulator\n";
	trueGenoVCF << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	trueGenoVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for (size_t ind = 0; ind < numInd; ++ind) trueGenoVCF << "\tInd" << ind + 1;
	trueGenoVCF << '\n';
}

std::array<std::vector<Base>,2> TSimulatorHaplotypes::getHaplotypesOfIndividual(size_t i) {
	if (i >= numInd)
		throw "Haplotypes of individual " + toString(i + 1) + " requested, but defined for only " + toString(numInd) +
			" individuals!";
	return haplotypes[i];
}

void TSimulatorHaplotypes::writeTrueGenotypes(const std::string &chrName, const TSimulatorReference &ref) {
	// prepare allele storage
	TSimulatorAlleleIndex index;
	std::string genoString;

	for (uint64_t l = 0; l < _length; ++l) {
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

bool TSimulatorHaplotypes::isPolymoprhic(uint64_t pos) const noexcept {
	// count how many allele match that of first individual
	const Base testBase = haplotypes[0][0][pos];
	size_t counts    = 0;
	for (auto& h: haplotypes) {
		if (h[0][pos] == testBase) ++counts;
		if (h[1][pos] == testBase) ++counts;
	}
	return counts != 2*numInd;
}

//---------------------------------------------------------
// TSimulatorMutationtable
//---------------------------------------------------------
TSimulatorMutationtable::TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq) {
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			_mutTable[a][b] = baseFreq[a] * baseFreq[b];
		}
		_mutTable[a][a] = 0.0;
	}
	_normalizeAndMakeCumulative();
}

TSimulatorMutationtable::TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq,
						 double theta) {
	const double exp_t   = exp(-theta);
	const double o_exp_t = 1 - exp_t;
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			_mutTable[a][b] = baseFreq[a] * baseFreq[b] * o_exp_t;
		}
		_mutTable[a][a] += baseFreq[a]*exp_t;
	}
	_normalizeAndMakeCumulative();
}

void TSimulatorMutationtable::_normalizeAndMakeCumulative() {
	// normalize within row
	for (auto & mu : _mutTable) {
		const double sum = std::accumulate(mu.begin(), mu.end(), 0.);
		mu[Base::A] /= sum;
		mu[Base::C] = mu[Base::C] / sum + mu[Base::A];
		mu[Base::G] = mu[Base::G] / sum + mu[Base::C];
		mu[Base::T] = 1.0;
	}
}

//---------------------------------------------------------
// TSimulatorVariantInvariantBedFiles
//---------------------------------------------------------
void TSimulatorVariantInvariantBedFiles::openFile(gz::ogzstream &file, const std::string filename) {
	file.open(filename.c_str());
	if (!file) throw "Failed to open file '" + filename + "' for writing!";
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
	uint64_t invariantStart = 0;
	for (uint64_t l = 0; l < haplotypes.length(); ++l) {
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

} // namespace Simulations
