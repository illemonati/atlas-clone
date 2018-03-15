/*
 * TSimulatorAuxiliaryTools.cpp
 *
 *  Created on: Nov 27, 2017
 *      Author: phaentu
 */



#include "TSimulatorAuxiliaryTools.h"


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
	if(chrName != "")
		writeRefToFasta();

	//move to new chr
	chrName = ChrName;
	if(ChrLength > storageLength)
		allocateStorage(ChrLength);
	chrLength = ChrLength;
	needsWriting = true;
};

void TSimulatorReference::simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float* cumulBaseFreq){
	logfile->listFlush("Simulating reference alleles ...");

	if(!storageInitialized)
		throw "Can not simulate reference sequence, no chromosome set!";

	//simulate reference sequence
	for(int l=0; l<chrLength; ++l){
		ref[l] = static_cast<Base>(randomGenerator->pickOne(4, cumulBaseFreq));
	}

	if(fastaOpen){
		writeRefToFasta();
	};
};


//---------------------------------------------------
//TSimulatorBamFile
//---------------------------------------------------
void TSimulatorBamFile::open(std::string Filename, std::vector<std::string> & readGroupNames, std::vector<TSimulatorChromosome> & chromosomes){
	logfile->listFlush("Opening BAM file '" + Filename + "' ...");

	if(isOpen)
		throw "A BAM file is already open for writing!";

	filename = Filename;

	if(chromosomes.size() < 1)
		throw "Can not open a BAM file without specified chromosomes!";

	BamTools::SamHeader header("");
	header.Version = "1.4";
	header.GroupOrder = "none";
	header.SortOrder = "coordinate";

	//add read group names
	for(std::vector<std::string>::iterator it=readGroupNames.begin(); it!=readGroupNames.end(); ++it)
		header.ReadGroups.Add(*it + "\tPU:UNKNOWN\tLB:UNKNOWN\tSM:Sim1\tCN:UNKNOWN\tPL:ILLUMINA");

	for(std::vector<TSimulatorChromosome>::iterator chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
		references.push_back(BamTools::RefData(chrIt->name, chrIt->length));
		header.Sequences.Add(BamTools::SamSequence(chrIt->name, chrIt->length));
	}

	if (!bamWriter.Open(filename, header, references))
		throw "Failed to open BAM file '" + filename + "'!";
	isOpen = true;
	logfile->done();
}

void TSimulatorBamFile::close(){
	if(isOpen){
		bamWriter.Close();

		//now generate bam index
		indexBamFile();
	}
	references.clear();
	isOpen = false;
}

void TSimulatorBamFile::indexBamFile(){
	logfile->listFlush("Creating index of BAM file '" + filename + "' ...");
	BamTools::BamReader reader;
	if(!reader.Open(filename))
		throw "Failed to open BAM file '" + filename + "' for indexing!";

	// create index for BAM file
	reader.CreateIndex(BamTools::BamIndex::STANDARD);

	//close BAM file
	reader.Close();
	logfile->done();
};

//---------------------------------------------------------
//TSimulatorHaplotypes
//---------------------------------------------------------
void TSimulatorHaplotypes::allocateStorage(long length){
	freeStorage();
	//allocate storage
	haplotypes = new Base**[numInd];
	for(ind=0; ind<numInd; ++ind){
		haplotypes[ind] = new Base*[2];
		haplotypes[ind][0] = new Base[length];
		haplotypes[ind][1] = new Base[length];
	}
	initialized = true;
	storageLength = length;
};

void TSimulatorHaplotypes::freeStorage(){
	if(initialized){
		for(ind=0; ind<numInd; ++ind){
			delete[] haplotypes[ind][0];
			delete[] haplotypes[ind][1];
			delete[] haplotypes[ind];
		}
		delete[] haplotypes;
		initialized = false;
	}
};

TSimulatorHaplotypes::TSimulatorHaplotypes(int NumIndividuals){
	numInd = NumIndividuals;
	ind = 0;
	haplotypes = NULL;
	initialized = false;
	curLength = 0;
	storageLength = 0;
};

void TSimulatorHaplotypes::setLength(long length){
	if(length > storageLength){
		allocateStorage(length);
	}
	curLength = length;
};

Base** TSimulatorHaplotypes::getHaplotypesOfIndividual(int i){
	if(i >= numInd)
		throw "Haplotypes of individual " + toString(i+1) + " requested, but defined for only " + toString(numInd) + " individuals!";
	return haplotypes[i];
};

void TSimulatorHaplotypes::writeGenotypes(gz::ogzstream & out, std::string & chrName, TGenotypeMap & genoMap){
	for(int l=0; l<curLength; ++l){
		out << chrName << "\t" << l+1;
		for(ind=0; ind < numInd; ++ind){


			//std::cout << ind << " @ " << l << ": " << std::flush;
			//std::cout << "\t" << toBase[haplotypes[ind][0][l]] << "/" << toBase[haplotypes[ind][1][l]] << std::endl;

			out << "\t" << genoMap.baseToChar[haplotypes[ind][0][l]] << "/" << genoMap.baseToChar[haplotypes[ind][1][l]];
		}
		out << "\n";
	}
};




