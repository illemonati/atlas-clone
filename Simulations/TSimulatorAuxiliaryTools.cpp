/*
 * TSimulatorAuxiliaryTools.cpp
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */

#include "TSimulatorAuxiliaryTools.h"

namespace Simulations {

using coretools::str::toString;

//---------------------------------------------------
// TSimulatorReference
//---------------------------------------------------

TSimulatorReference::TSimulatorReference(std::string Filename, TLog *Logfile) : _logfile(Logfile), _filename(Filename) {
	// open FASTA file for reference sequences
	_logfile->list("Will write reference sequence to '" + _filename + "'.");
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
void TSimulatorBamFile::open(const std::string Filename, const std::string SampleName, BAM::TReadGroups ReadGroups,
			     const BAM::TChromosomes &Chromosomes, TLog *Logfile) {
	Logfile->listFlush("Opening BAM file '" + Filename + "' ...");

	if (_outBam.isOpen()) throw "A BAM file is already open for writing!";

	if (Chromosomes.size() < 1) throw "Can not open a BAM file without specified chromosomes!";

	// create header, read group and chromosome objects
	_header.set("1.6", "coordinate", "none");
	for (auto &rg : ReadGroups) {
		rg.sequencingCenter_CN =
			coretools::__GLOBAL_APPLICATION_NAME__ + " " + coretools::__GLOBAL_APPLICATION_VERSION__;
		rg.description_DS          = "Simulated with commit " + coretools::__GLOBAL_APPLICATION_COMMIT__;
		rg.sample_SM               = SampleName;
		rg.sequencingTechnology_PL = "ILLUMINA";
	}

	_outBam.open(Filename, _header, Chromosomes, _readGroups);

	Logfile->done();
};

TSimulatorBamFile::~TSimulatorBamFile() { _outBam.closeNoIndex(); };

void TSimulatorBamFile::close(TLog *Logfile) { _outBam.close(Logfile); };

TSimulatorBamFiles::TSimulatorBamFiles(uint32_t NumFiles, const std::string Outname, const BAM::TReadGroups ReadGroups,
				       const BAM::TChromosomes &Chromosomes, TLog *Logfile) {
	if (NumFiles < 1) throw "Can not open less than one BAM file!";
	_logfile = Logfile;
	_files.resize(NumFiles);

	// open BAM files
	if (_files.size() == 1) {
		_files[0].open(Outname + ".bam", "Ind1", ReadGroups, Chromosomes, Logfile);
	} else {
		Logfile->startIndent("Opening ", _files.size(), " BAM files:");
		for (size_t i = 0; i < _files.size(); ++i) {
			_files[i].open(Outname + "_ind" + toString(i + 1) + ".bam", "Ind" + toString(i + 1), ReadGroups,
				       Chromosomes, Logfile);
		}
		Logfile->endIndent();
	}
};

void TSimulatorBamFiles::close() {
	_logfile->startIndent("Indexing BAM files:");
	for (auto &f : _files) { f.close(_logfile); }
	_logfile->endIndent();
};

TSimulatorBamFile &TSimulatorBamFiles::operator[](size_t i) {
	if (i >= _files.size()) throw "BAM file " + toString(i) + " does not exist!";
	return _files[i];
};

//---------------------------------------------------------
// TSimulatorHaplotypes
//---------------------------------------------------------
void TSimulatorHaplotypes::allocateStorage(uint64_t length) {
	// allocate storage
	haplotypes.resize(numInd);
	for (int ind = 0; ind < numInd; ++ind) {
		haplotypes[ind][0].resize(length);
		haplotypes[ind][1].resize(length);
	}
	initialized   = true;
	storageLength = length;
};

void TSimulatorHaplotypes::setLength(uint32_t length) {
	if (length > storageLength) { allocateStorage(length); }
	_length = length;
};

void TSimulatorHaplotypes::openTrueGenotypeVCF(std::string filename) {
	// open file
	trueGenoVCF.open(filename.c_str());
	if (!trueGenoVCF) throw "Failed to open VCF file '" + filename + "' for writing!";

	trueGenoVCFOpend = true;

	// write header
	trueGenoVCF << "##fileformat=VCFv4.3\n";
	trueGenoVCF << "##source=ATLAS_Simulator\n";
	trueGenoVCF << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	trueGenoVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for (int ind = 0; ind < numInd; ++ind) trueGenoVCF << "\tInd" << ind + 1;
	trueGenoVCF << '\n';
};

void TSimulatorHaplotypes::closeTrueGenotypeVCF() {
	if (trueGenoVCFOpend) {
		trueGenoVCF.close();
		trueGenoVCFOpend = false;
	}
};

std::array<std::vector<Base>,2> TSimulatorHaplotypes::getHaplotypesOfIndividual(int i) {
	if (i >= numInd)
		throw "Haplotypes of individual " + toString(i + 1) + " requested, but defined for only " + toString(numInd) +
			" individuals!";
	return haplotypes[i];
};

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
		for (int ind = 0; ind < numInd; ++ind) {
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
};

bool TSimulatorHaplotypes::isPolymoprhic(uint64_t pos) {
	// count how many allele match that of first individual
	Base testBase = haplotypes[0][0][pos];
	int counts    = 0;
	for (int ind = 0; ind < numInd; ++ind) {
		if (haplotypes[ind][0][pos] == testBase) ++counts;
		if (haplotypes[ind][1][pos] == testBase) ++counts;
	}

	if (counts == 2 * numInd)
		return false;
	else
		return true;
};

//---------------------------------------------------------
// TSimulatorMutationtable
//---------------------------------------------------------
TSimulatorMutationtable::TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq) {
	fill(baseFreq);
};

TSimulatorMutationtable::TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq,
						 double theta) {
	fill(baseFreq, theta);
};

void TSimulatorMutationtable::_normalizeAndMakeCumulative() {
	// normalize within row
	for (int i = 0; i < 4; ++i) {
		double sum = 0.0;
		for (int j = 0; j < 4; ++j) { sum += _mutTable[i][j]; }
		for (int j = 0; j < 4; ++j) { _mutTable[i][j] /= sum; }

		// make cumulative
		_mutTable[i][1] += _mutTable[i][0];
		_mutTable[i][2] += _mutTable[i][1];
		_mutTable[i][3] = 1.0;
	}
};

void TSimulatorMutationtable::fill(const GenotypeLikelihoods::TBaseProbabilities &baseFreq) {
	for (uint8_t i = 0; i < 4; ++i) {
		for (uint8_t j = 0; j < 4; ++j) {
			_mutTable[i][j] = baseFreq[static_cast<Base>(i)] * baseFreq[static_cast<Base>(j)];
		}
		_mutTable[i][i] = 0.0;
	}

	_normalizeAndMakeCumulative();
};

void TSimulatorMutationtable::fill(const GenotypeLikelihoods::TBaseProbabilities &baseFreq, double theta) {
	double expMinusTheta = exp(-theta);
	for (uint8_t i = 0; i < 4; ++i) {
		for (uint8_t j = 0; j < 4; ++j) {
			_mutTable[i][j] = baseFreq[static_cast<Base>(i)] * baseFreq[static_cast<Base>(j)] * (1.0 - expMinusTheta);
		}
		_mutTable[i][i] += baseFreq[static_cast<Base>(i)] * expMinusTheta;
	}

	_normalizeAndMakeCumulative();
};

void TSimulatorMutationtable::print() {
	for (int i = 0; i < 4; ++i) {
		std::cout << "Mutation table " << i << ":";
		for (int j = 0; j < 4; ++j) { std::cout << " " << _mutTable[i][j]; }
		std::cout << std::endl;
	}
};

//---------------------------------------------------------
// TSimulatorVariantInvariantBedFiles
//---------------------------------------------------------
TSimulatorVariantInvariantBedFiles::TSimulatorVariantInvariantBedFiles() { filesOpend = false; };

TSimulatorVariantInvariantBedFiles::TSimulatorVariantInvariantBedFiles(std::string outname) {
	filesOpend = false;
	open(outname);
};

TSimulatorVariantInvariantBedFiles::~TSimulatorVariantInvariantBedFiles() { close(); };

void TSimulatorVariantInvariantBedFiles::openFile(gz::ogzstream &file, const std::string filename) {
	file.open(filename.c_str());
	if (!file) throw "Failed to open file '" + filename + "' for writing!";
}

void TSimulatorVariantInvariantBedFiles::open(std::string outname) {
	// make sure files are closed
	close();

	// now open files
	openFile(variantSitesFile, outname + "_variantSites.bed.gz");
	openFile(invariantSitesFile, outname + "_invariantSites.bed.gz");
	filesOpend = true;
};

void TSimulatorVariantInvariantBedFiles::close() {
	if (filesOpend) {
		variantSitesFile.close();
		invariantSitesFile.close();
		filesOpend = false;
	}
};

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
};

}; // namespace Simulations
