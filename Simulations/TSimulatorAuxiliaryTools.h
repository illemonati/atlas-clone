/*
 * TSimulatorAuxiliaryTools.h
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATORAUXILIARYTOOLS_H_
#define TSIMULATORAUXILIARYTOOLS_H_

#include "TLog.h"
#include "TRandomGenerator.h"
#include "stringFunctions.h"
#include "gzstream.h"
#include "TBamFile.h"
#include "globalConstants.h"

namespace Simulations{

using coretools::TLog;
using genometools::Base;

//---------------------------------------------------------
//TSimulatorReference
//---------------------------------------------------------
class TSimulatorReference{
private:
	TLog* _logfile;

	//fasta file
	std::ofstream _fasta;
	std::ofstream _fastaIndex;
	long _oldOffset;
	bool _fastaOpen;
	std::string _filename;

	//reference storage
	std::vector<Base> _ref;
	long _chrLength;
	std::string _chrName;
	bool _needsWriting;

	void _allocateStorage(long length);
	void _openFastaFile();
	void _closeFastaFile();
	void _writeRefToFasta();

public:
	TSimulatorReference();
	TSimulatorReference(std::string Filename, TLog* Logfile);
	~TSimulatorReference(){
		close();
	};

	void initialize(std::string Filename, TLog* Logfile);
	void close();

	void setChr(std::string ChrName, long ChrLength);
//	void simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float* cumulBaseFreq);
	std::vector<Base>& reference(){ return _ref; };

	Base& operator[](const uint32_t & index){ return _ref[index]; };
	const Base& operator[](const uint32_t & index) const { return _ref[index]; };
};

//---------------------------------------------------------
//TSimulatorBamFile
//---------------------------------------------------------
class TSimulatorBamFile{
private:
	BAM::TSamHeader _header;
	BAM::TReadGroups _readGroups;
	BAM::TOutputBamFile _outBam;

public:
	TSimulatorBamFile(){};
	TSimulatorBamFile(const std::string Filename, const std::string SampleName, BAM::TReadGroups ReadGroups, const BAM::TChromosomes & Chromosomes, TLog* Logfile){
		open(Filename, SampleName, ReadGroups, Chromosomes, Logfile);
	};
	~TSimulatorBamFile();

	void open(const std::string Filename, const std::string SampleName, BAM::TReadGroups ReadGroups, const BAM::TChromosomes & Chromosomes, TLog* Logfile);

	void saveAlignment(const BAM::TAlignment & Alignment){
		_outBam.writeAlignment(Alignment);
	};
	void close(TLog* Logfile);
	void indexBamFile();
};

class TSimulatorBamFiles{
private:
	std::vector<TSimulatorBamFile> _files;
	TLog* _logfile;

public:
	TSimulatorBamFiles(uint32_t NumFiles, const std::string Outname, BAM::TReadGroups ReadGroups, const BAM::TChromosomes & Chromosomes, TLog* Logfile);

	void close();
	TSimulatorBamFile& operator[](size_t i);
};

//---------------------------------------------------------
//TSimulatorAlleleIndex
//---------------------------------------------------------
class TSimulatorAlleleIndex{
private:
	int nextIndex;
	Base refBase;
	Base indexToBase[4];

public:
	GenotypeLikelihoods::TBaseData_base<int> index;
	GenotypeLikelihoods::TBaseData_base<bool> used;

	TSimulatorAlleleIndex(){
		nextIndex = 0;
		refBase = genometools::N;
	};

	void clear(const Base & ref){
		used.set(false);
		used[ref] = true;
		index[ref] = 0;
		nextIndex = 1;
		refBase = ref;
	};

	void add(const Base & base){
		if(!used[base]){
			used[base] = true;
			index[base] = nextIndex;
			indexToBase[nextIndex] = base;
			++nextIndex;
		}
	};

	void writeRefAltToVCF(gz::ogzstream & VCF){
		VCF << refBase << '\t';
		if(nextIndex == 1) //no alt
			VCF << ".";
		else {
			VCF << indexToBase[1];
			for(int i=2; i<nextIndex; ++i)
				VCF << ',' << indexToBase[i];
		}
	};
};

//---------------------------------------------------------
//TSimulatorHaplotypes
//---------------------------------------------------------
class TSimulatorHaplotypes{
private:
	int numInd;
	uint32_t _length;
	uint32_t storageLength;
	bool initialized;
	Base*** haplotypes;

	//write true genotypes to VCF
	gz::ogzstream trueGenoVCF;
	bool trueGenoVCFOpend;

	void allocateStorage(uint64_t length);
	void freeStorage();

public:
	TSimulatorHaplotypes(int NumIndividuals);
	~TSimulatorHaplotypes(){
		freeStorage();
	};

	void setLength(uint32_t length);
	uint32_t length() const{ return _length; };
	void openTrueGenotypeVCF(std::string filename);
	void closeTrueGenotypeVCF();
	Base** getHaplotypesOfIndividual(int i);
	Base** getHaplotypesFirstIndividual(){
		return haplotypes[0];
	};
	void writeTrueGenotypes(const std::string & chrName, const TSimulatorReference & ref);
	int size(){ return numInd; };
	Base& operator()(int ind, int hap, uint64_t site){
		return haplotypes[ind][hap][site];
	};
	bool isPolymoprhic(uint64_t pos);
};

//---------------------------------------------------------
//TSimulatorMutationtable
//---------------------------------------------------------
class TSimulatorMutationtable{
private:
	std::array< std::array<double, 4>, 4> _mutTable;

	void _normalizeAndMakeCumulative();

public:
	TSimulatorMutationtable() = default;
	TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities & baseFreq);
	TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities & baseFreq, const double & theta);
	~TSimulatorMutationtable() = default;
	void fill(const GenotypeLikelihoods::TBaseProbabilities & baseFreq);
	void fill(const GenotypeLikelihoods::TBaseProbabilities & baseFreq, const double & theta);
	const std::array<double, 4> operator[](const genometools::Base & base){ return _mutTable[base.get()]; };
	void print();
};

//---------------------------------------------------------
//TSimulatorVariantInvariantBedFiles
//---------------------------------------------------------
class TSimulatorVariantInvariantBedFiles{
private:
	gz::ogzstream variantSitesFile;
	gz::ogzstream invariantSitesFile;
	bool filesOpend;

	void openFile(gz::ogzstream & file, const std::string filename);
	void close();

public:
	TSimulatorVariantInvariantBedFiles();
	TSimulatorVariantInvariantBedFiles(std::string outname);
	~TSimulatorVariantInvariantBedFiles();

	void open(std::string outname);

	void write(TSimulatorHaplotypes & haplotypes, const std::string & chrName);
};

}; //end namespace

#endif /* TSIMULATORAUXILIARYTOOLS_H_ */
