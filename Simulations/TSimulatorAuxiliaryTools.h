/*
 * TSimulatorAuxiliaryTools.h
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATORAUXILIARYTOOLS_H_
#define TSIMULATORAUXILIARYTOOLS_H_

#include "TBamFile.h"
#include "globalConstants.h"
#include "gzstream.h"
#include "stringFunctions.h"

namespace Simulations {


//---------------------------------------------------------
// TSimulatorReference
//---------------------------------------------------------
class TSimulatorReference {
private:

	// fasta file
	std::ofstream _fasta;
	std::ofstream _fastaIndex;
	long _oldOffset = 0;
	std::string _filename;

	// reference storage
	std::vector<genometools::Base> _ref;
	long _chrLength      = 0;
	std::string _chrName = "";
	bool _needsWriting   = false;

	void _allocateStorage(long length);
	void _closeFastaFile();
	void _writeRefToFasta();

public:
	TSimulatorReference(std::string Filename);
	~TSimulatorReference();
	void setChr(std::string ChrName, long ChrLength);
	//	void simulateReferenceSequenceCurChromosome(float* cumulBaseFreq);

	genometools::Base &operator[](size_t index) { return _ref[index]; }
	const genometools::Base &operator[](size_t index) const { return _ref[index]; }
	const std::vector<genometools::Base> & reference() const {return _ref;}
};

//---------------------------------------------------------
// TSimulatorBamFile
//---------------------------------------------------------
class TSimulatorBamFile {
private:
	BAM::TOutputBamFile _outBam;

public:
	TSimulatorBamFile(){}
	TSimulatorBamFile(const std::string &Filename, const std::string &SampleName, BAM::TReadGroups &ReadGroups,
			  const genometools::TChromosomes &Chromosomes) {
		open(Filename, SampleName, ReadGroups, Chromosomes);
	}
	~TSimulatorBamFile();

	void open(const std::string &Filename, const std::string &SampleName, BAM::TReadGroups &ReadGroups,
		  const genometools::TChromosomes &Chromosomes);

	void saveAlignment(const BAM::TAlignment &Alignment) {_outBam.writeAlignment(Alignment); }
	void close();
	void indexBamFile();
};

class TSimulatorBamFiles {
private:
	std::vector<TSimulatorBamFile> _files;

public:
	TSimulatorBamFiles(uint32_t NumFiles, const std::string & Outname, BAM::TReadGroups & ReadGroups,
			   const genometools::TChromosomes &Chromosomes);
	~TSimulatorBamFiles();
	TSimulatorBamFile &operator[](size_t i);
};

//---------------------------------------------------------
// TSimulatorAlleleIndex
//---------------------------------------------------------
class TSimulatorAlleleIndex {
private:
	int nextIndex = 0;
	genometools::Base refBase = genometools::Base::N;
	genometools::Base indexToBase[4];

public:
	GenotypeLikelihoods::TBaseData_base<int> index;
	GenotypeLikelihoods::TBaseData_base<bool> used;

	void clear(const genometools::Base &ref) noexcept {
		used.set(false);
		used[ref]  = true;
		index[ref] = 0;
		nextIndex  = 1;
		refBase    = ref;
	}

	void add(const genometools::Base &base) noexcept {
		if (!used[base]) {
			used[base]             = true;
			index[base]            = nextIndex;
			indexToBase[nextIndex] = base;
			++nextIndex;
		}
	}

	void writeRefAltToVCF(gz::ogzstream &VCF) {
		VCF << refBase << '\t';
		if (nextIndex == 1) // no alt
			VCF << ".";
		else {
			VCF << indexToBase[1];
			for (int i = 2; i < nextIndex; ++i) VCF << ',' << indexToBase[i];
		}
	}
};

//---------------------------------------------------------
// TSimulatorHaplotypes
//---------------------------------------------------------
class TSimulatorHaplotypes {
private:
	int numInd;
	uint32_t _length = 0;
	std::vector<std::array<std::vector<genometools::Base>,2>> haplotypes;

	// write true genotypes to VCF
	gz::ogzstream trueGenoVCF;

	void allocateStorage();
public:
	TSimulatorHaplotypes(int NumIndividuals): numInd(NumIndividuals){};
	~TSimulatorHaplotypes() {
		if (trueGenoVCF) trueGenoVCF.close();
	}
	void setLength(uint32_t length) noexcept;
	uint32_t length() const { return _length; };
	void openTrueGenotypeVCF(std::string filename);
	std::array<std::vector<genometools::Base>,2> getHaplotypesOfIndividual(int i);
	std::array<std::vector<genometools::Base>,2> getHaplotypesFirstIndividual() { return haplotypes[0]; };
	void writeTrueGenotypes(const std::string &chrName, const TSimulatorReference &ref);
	int size() const noexcept { return numInd; };
	genometools::Base &operator()(int ind, int hap, uint64_t site) noexcept { return haplotypes[ind][hap][site]; };
	const genometools::Base &operator()(int ind, int hap, uint64_t site) const noexcept { return haplotypes[ind][hap][site]; };
	bool isPolymoprhic(uint64_t pos) const noexcept;
};

//---------------------------------------------------------
// TSimulatorMutationtable
//---------------------------------------------------------
class TSimulatorMutationtable {
private:
	std::array<std::array<double, 4>, 4> _mutTable;
	void _normalizeAndMakeCumulative();
public:
	TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq);
	TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq, double theta);
	~TSimulatorMutationtable() = default;
	const std::array<double, 4> operator[](genometools::Base base) const noexcept { return _mutTable[genometools::index(base)]; }
	void print();
};

//---------------------------------------------------------
// TSimulatorVariantInvariantBedFiles
//---------------------------------------------------------
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

} // namespace Simulations

#endif /* TSIMULATORAUXILIARYTOOLS_H_ */
