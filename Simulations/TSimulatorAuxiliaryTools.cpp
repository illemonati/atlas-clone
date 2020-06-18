/*
 * TSimulatorAuxiliaryTools.cpp
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */



#include "TSimulatorAuxiliaryTools.h"

namespace Simulations{

//---------------------------------------------------
//TSimulatorReference
//---------------------------------------------------
TSimulatorReference::TSimulatorReference(){
	logfile = NULL;
	chrName = "";
	ref = NULL;
	storageInitialized = false;
	storageLength = 0;
	fastaOpen = false;
	oldOffset = 0;
	chrLength = 0;
	needsWriting = false;
};

TSimulatorReference::TSimulatorReference(std::string Filename, TLog* Logfile){
	chrName = "";
	ref = NULL;
	storageInitialized = false;
	storageLength = 0;
	fastaOpen = false;
	chrLength = 0;
	needsWriting = false;

	initialize(Filename, Logfile);
};

void TSimulatorReference::initialize(std::string Filename, TLog* Logfile){
	filename = Filename;
	logfile = Logfile;

	openFastaFile();
};

void TSimulatorReference::close(){
	if(chrName != "" && needsWriting)
		writeRefToFasta();
	closeFastaFile();
	freeStorage();
};

void TSimulatorReference::openFastaFile(){
	//open FASTA file for reference sequences
	logfile->list("Will write reference sequence to '" + filename + "'.");
	fasta.open(filename.c_str());
	if(!fasta)
		throw "Failed to open file '" + filename + "' for writing!";
	filename += ".fai";
	fastaIndex.open(filename.c_str());
	if(!fastaIndex)
		throw "Failed to open file '" + filename + "' for writing!";
	oldOffset = 0;
	fastaOpen = true;
};

void TSimulatorReference::closeFastaFile(){
	if(fastaOpen){
		fasta.close();
		fastaIndex.close();
	}
	fastaOpen = false;
};

void TSimulatorReference::writeRefToFasta(){
	if(fastaOpen){
		//write to fasta
		fasta << ">" << chrName;
		for(int l=0; l<chrLength; ++l){
			if(l % 70 == 0)
				fasta << "\n";
			fasta << genoMap.baseToChar[ref[l]];

			//std::cout << "Writing base " << ref[l] << " -> " << toBase[ref[l]] << " to fasta..." << std::endl;

		}
		fasta << "\n";

		//add to index
		std::string tmp = chrName;
		oldOffset += chrName.size() + 2;
		fastaIndex << extractBeforeWhiteSpace(tmp) << "\t" << chrLength << "\t" << oldOffset << "\t70\t71\n";
		oldOffset += chrLength + (int) (chrLength / 70);
		if(chrLength % 70 != 0) oldOffset += 1;

		needsWriting = false;
	} else
		throw "Can not write to FASTA file: file was never opened!";
};

void TSimulatorReference::allocateStorage(long length){
	freeStorage();

	//allocate storage
	ref = new Base[length];
	storageInitialized = true;
	storageLength = length;
};

void TSimulatorReference::freeStorage(){
	if(storageInitialized){
		delete[] ref;
		storageInitialized = false;
	}
};

void TSimulatorReference::setChr(std::string ChrName, long ChrLength){
	//write if not yet written
	if(chrName != "" && needsWriting)
		writeRefToFasta();

	//move to new chr
	chrName = ChrName;
	if(ChrLength > storageLength)
		allocateStorage(ChrLength);
	chrLength = ChrLength;
	needsWriting = true;
};

//void TSimulatorReference::simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float* cumulBaseFreq){
//	logfile->listFlush("Simulating reference alleles ...");
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
//};


//---------------------------------------------------
//TSimulatorBamFile
//---------------------------------------------------
void TSimulatorBamFile::open(const std::string Filename,
		                     const std::string SampleName,
							 const std::vector<std::string> & ReadGroupNames,
							 const std::vector<TSimulatorChromosome> & Chromosomes,
							 TLog* Logfile,
							 TGenotypeMap & GenoMap,
							 TQualityMap & QualMap){
	Logfile->listFlush("Opening BAM file '" + Filename + "' ...");

	if(_outBam.isOpen())
		throw "A BAM file is already open for writing!";

	if(Chromosomes.size() < 1)
		throw "Can not open a BAM file without specified chromosomes!";

	//create header, read group and chromosome objects
	BAM::TSamHeader header("1.6", "coordinate", "none");
	BAM::TChromosomes chr;
	for(auto it = Chromosomes.cbegin(); it != Chromosomes.cend(); ++it){
		chr.appendChromosome(it->name, it->length);
	}
	BAM::TReadGroups readGroups;
	for(auto it = ReadGroupNames.cbegin(); it != ReadGroupNames.cend(); ++it){
		BAM::TReadGroup& rg = readGroups.add(*it);
		rg.sequencingCenter_CN = __GLOBAL_APPLICATION_NAME__ + " " + __GLOBAL_APPLICATION_VERSION__;
		rg.description_DS = "Simulated with commit " + __GLOBAL_APPLICATION_COMMIT__;
		rg.sample_SM = SampleName;
		rg.sequencingTechnology_PL = "ILLUMINA";
	}

	_outBam.open(Filename, header, chr, readGroups, &GenoMap, &QualMap);

	Logfile->done();
};

TSimulatorBamFile::~TSimulatorBamFile(){
	_outBam.closeNoIndex();
};

void TSimulatorBamFile::close(TLog* Logfile){
	_outBam.close(Logfile);
};

TSimulatorBamFiles::TSimulatorBamFiles(uint32_t NumFiles, const std::string Outname, const std::vector<std::string> & ReadGroupNames, const std::vector<TSimulatorChromosome> & Chromosomes, TLog* Logfile, TGenotypeMap & GenoMap, TQualityMap & QualMap){
	if(NumFiles < 1) throw "Can not open less than one BAM file!";
	_logfile = Logfile;
	_files.resize(NumFiles);

	//open BAM files
	if(_files.size() == 1){
		_files[0].open(Outname + ".bam", "Ind1", ReadGroupNames, Chromosomes, Logfile, GenoMap, QualMap);
	} else {
		Logfile->startIndent("Opening " + toString(_files.size()) + " BAM files:");
		for(int i=0; i<_files.size(); ++i){
			_files[i].open(Outname + "_ind" + toString(i+1) + ".bam", "Ind" + toString(i+1), ReadGroupNames, Chromosomes, Logfile, GenoMap, QualMap);
		}
		Logfile->endIndent();
	}
};

TSimulatorBamFiles::~TSimulatorBamFiles(){
	delete[] _files;
};

void TSimulatorBamFiles::close(){
	_logfile->startIndent("Indexing BAM files:");
	for(auto& f : _files){
		f.close(_logfile);
	}
	_logfile->endIndent();
};

TSimulatorBamFile& TSimulatorBamFiles::operator[](size_t i){
	if(i >= _files.size()) throw "BAM file " + toString(i) + " does not exist!";
	return _files[i];
};

//---------------------------------------------------------
//TSimulatorHaplotypes
//---------------------------------------------------------
TSimulatorHaplotypes::TSimulatorHaplotypes(int NumIndividuals){
	numInd = NumIndividuals;
	haplotypes = NULL;
	initialized = false;
	curLength = 0;
	storageLength = 0;
	trueGenoVCFOpend = false;
};

void TSimulatorHaplotypes::allocateStorage(uint64_t length){
	freeStorage();
	//allocate storage
	haplotypes = new Base**[numInd];
	for(int ind=0; ind<numInd; ++ind){
		haplotypes[ind] = new Base*[2];
		haplotypes[ind][0] = new Base[length];
		haplotypes[ind][1] = new Base[length];
	}
	initialized = true;
	storageLength = length;
};

void TSimulatorHaplotypes::freeStorage(){
	if(initialized){
		for(int ind=0; ind<numInd; ++ind){
			delete[] haplotypes[ind][0];
			delete[] haplotypes[ind][1];
			delete[] haplotypes[ind];
		}
		delete[] haplotypes;
		initialized = false;
	}
};

void TSimulatorHaplotypes::setLength(uint64_t length){
	if(length > storageLength){
		allocateStorage(length);
	}
	curLength = length;
};

void TSimulatorHaplotypes::openTrueGenotypeVCF(std::string filename){
	//open file
	trueGenoVCF.open(filename.c_str());
	if(!trueGenoVCF)
		throw "Failed to open VCF file '" + filename + "' for writing!";

	trueGenoVCFOpend = true;

	//write header
	trueGenoVCF << "##fileformat=VCFv4.3\n";
	trueGenoVCF << "##source=ATLAS_Simulator\n";
	trueGenoVCF << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	trueGenoVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for(int ind=0; ind<numInd; ++ind)
		trueGenoVCF << "\tInd" << ind + 1;
	trueGenoVCF << '\n';
};

void TSimulatorHaplotypes::closeTrueGenotypeVCF(){
	if(trueGenoVCFOpend){
		trueGenoVCF.close();
		trueGenoVCFOpend = false;
	}
};

Base** TSimulatorHaplotypes::getHaplotypesOfIndividual(int i){
	if(i >= numInd)
		throw "Haplotypes of individual " + toString(i+1) + " requested, but defined for only " + toString(numInd) + " individuals!";
	return haplotypes[i];
};

void TSimulatorHaplotypes::writeTrueGenotypes(TSimulatorChromosome & chromosome, Base* ref, TGenotypeMap & genoMap){
	//prepare allele storage
	TSimulatorAlleleIndex index;
	std::string genoString;

	for(uint64_t l=0; l<curLength; ++l){
		//chromosome name, position and ID
		trueGenoVCF << chromosome.name << '\t' << l+1 << "\t.\t";

		//assemble alleles and genotypes
		genoString.clear();
		index.clear(ref[l]);

		//loop over all individuals to figure out which alleles are used
		for(int ind=0; ind < numInd; ++ind){
			//homozygous or heterozygous?
			if(haplotypes[ind][0][l] == haplotypes[ind][1][l]){
				//make sure allele exists
				index.add(haplotypes[ind][0][l]);

				//add genotype
				genoString += '\t' + toString(index.index[haplotypes[ind][0][l]]) + '/' + toString(index.index[haplotypes[ind][0][l]]);
			} else {
				//make sure allele exists
				index.add(haplotypes[ind][0][l]);
				index.add(haplotypes[ind][1][l]);

				if(index.index[haplotypes[ind][0][l]] < index.index[haplotypes[ind][1][l]])
					genoString += '\t' + toString(index.index[haplotypes[ind][0][l]]) + '/' + toString(index.index[haplotypes[ind][1][l]]);
				else
					genoString += '\t' + toString(index.index[haplotypes[ind][1][l]]) + '/' + toString(index.index[haplotypes[ind][0][l]]);
			}
		}

		//write ref allele
		index.writeRefAltToVCF(trueGenoVCF, genoMap);

		//write (no) quality of variant, (no) filter, (no) info and format
		trueGenoVCF << "\t.\t.\t.\tGT";

		//now write genotypes
		trueGenoVCF <<  genoString << '\n';
	}
};

bool TSimulatorHaplotypes::isPolymoprhic(uint64_t pos){
	//count how many allele match that of first individual
	Base testBase = haplotypes[0][0][pos];
	int counts = 0;
	for(int ind=0; ind < numInd; ++ind){
		if(haplotypes[ind][0][pos] == testBase) ++counts;
		if(haplotypes[ind][1][pos] == testBase) ++counts;
	}

	if(counts == 2*numInd) return false;
	else return true;
};

//---------------------------------------------------------
//TSimulatorMutationtable
//---------------------------------------------------------
TSimulatorMutationtable::TSimulatorMutationtable(){
	mutTable = NULL;
	tableAllocated = false;
};

TSimulatorMutationtable::TSimulatorMutationtable(float* baseFreq){
	tableAllocated = false;
	fill(baseFreq);
};

TSimulatorMutationtable::TSimulatorMutationtable(float* baseFreq, double theta){
	tableAllocated = false;
	fill(baseFreq, theta);
};

void TSimulatorMutationtable::allocateTable(){
	if(!tableAllocated){
		mutTable = new float*[4];
		for(int i=0; i<4; ++i)
			mutTable[i] = new float[4];
		tableAllocated = true;
	}
};

void TSimulatorMutationtable::deleteTable(){
	if(tableAllocated){
		for(int i=0; i<4; ++i)
			delete[] mutTable[i];
		delete[] mutTable;
		tableAllocated = false;
	}
};

void TSimulatorMutationtable::normalizeAndMakeCumulative(){
	//normalize within row
	for(int i=0; i<4; ++i){
		double sum = 0.0;
		for(int j=0; j<4; ++j){
			sum += mutTable[i][j];
		}
		for(int j=0; j<4; ++j){
			mutTable[i][j] /= sum;
		}

		//make cumulative
		mutTable[i][1] += mutTable[i][0];
		mutTable[i][2] += mutTable[i][1];
		mutTable[i][3] = 1.0;
	}
};

void TSimulatorMutationtable::fill(float* baseFreq){
	//create storage
	allocateTable();

	//fill table
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			mutTable[i][j] = baseFreq[i] * baseFreq[j];
		}
		mutTable[i][i] = 0.0;
	}

	normalizeAndMakeCumulative();
};

void TSimulatorMutationtable::fill(float* baseFreq, const double theta){
	//create storage
	allocateTable();

	//fill table
	double expMinusTheta = exp(-theta);
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			mutTable[i][j] = baseFreq[i] * baseFreq[j] * (1.0 - expMinusTheta);
		}
		mutTable[i][i] += baseFreq[i] * expMinusTheta;
	}

	normalizeAndMakeCumulative();
};

void TSimulatorMutationtable::print(){
	for(int i=0; i<4; ++i){
		std::cout << "Mutation table " << i << ":";
		for(int j=0; j<4; ++j){
			std::cout << " " << mutTable[i][j];
		}
		std::cout << std::endl;
	}
};

//---------------------------------------------------------
//TSimulatorVariantInvariantBedFiles
//---------------------------------------------------------
TSimulatorVariantInvariantBedFiles::TSimulatorVariantInvariantBedFiles(){
	filesOpend = false;
};

TSimulatorVariantInvariantBedFiles::TSimulatorVariantInvariantBedFiles(std::string outname){
	filesOpend = false;
	open(outname);
};

TSimulatorVariantInvariantBedFiles::~TSimulatorVariantInvariantBedFiles(){
	close();
};

void TSimulatorVariantInvariantBedFiles::openFile(gz::ogzstream & file, const std::string filename){
	file.open(filename.c_str());
	if(!file)
		throw "Failed to open file '" + filename + "' for writing!";
}

void TSimulatorVariantInvariantBedFiles::open(std::string outname){
	//make sure files are closed
	close();

	//now open files
	openFile(variantSitesFile, outname + "_variantSites.bed.gz");
	openFile(invariantSitesFile, outname + "_invariantSites.bed.gz");
	filesOpend = true;
};

void TSimulatorVariantInvariantBedFiles::close(){
	if(filesOpend){
		variantSitesFile.close();
		invariantSitesFile.close();
		filesOpend = false;
	}
};

void TSimulatorVariantInvariantBedFiles::write(TSimulatorHaplotypes & haplotypes, TSimulatorChromosome & chromosome){
	//0-based
	uint64_t invariantStart = 0;
	for(uint64_t l=0; l<chromosome.length; ++l){
		if(haplotypes.isPolymoprhic(l)){
			//write invariant
			if(invariantStart < l)
				invariantSitesFile << chromosome.name << "\t" << invariantStart << "\t" << l << "\n";
			invariantStart = l+1;

			//write variant
			variantSitesFile << chromosome.name << "\t" << l << "\t" << l+1 << "\n";
		}
	}

	//write last invariant interval
	if(invariantStart < chromosome.length)
		invariantSitesFile << chromosome.name << "\t" << invariantStart << "\t" << chromosome.length << "\n";
};

}; //end namespace
