/*
 * TSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#include "TSimulator.h"

TSimulator::TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	randomGenerator = RandomGenerator;
	bamFileOpen = false;
	fastaOpen = false;

	//set default parameters
	setReadLength(100);
	setQualityDistribution(30.0, 10.0);
	setDepth(10.0);
	bamAlignment.Name = "*";
	bamAlignment.MapQuality = 50;
	setReadGroupName("SimReadGroup");

	//helper tools
	toBase[0] = 'A'; toBase[1] = 'C'; toBase[2] = 'G'; toBase[3] = 'T';
	toBase[0] = 'A'; toBase[1] = 'C'; toBase[2] = 'G'; toBase[3] = 'T';
	ref = NULL;
	alt = NULL;
	refInitialized = false;
	qualToErroTableInitialized = false;
};

void TSimulator::setQualityDistribution(double mean, double sd){
	meanQual = mean + 33.0; //add 33 to mean quality to get in in char
	sdQual = sd;
}

void TSimulator::initializeChromosomes(int numChr, long chrLength){
	chromosomes.clear();
	for(int i=0; i<numChr; ++i){
		chromosomes.insert(std::pair<std::string, long>("chr" + toString(i+1), chrLength));
	}
}

void TSimulator::initializeChromosomes(std::map<std::string, long> & chr){
	chromosomes.clear();
	chromosomes = chr;
}

void TSimulator::setReadLength(int length){
	readLength = length;
	bamAlignment.Length = readLength;
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', readLength));
}

void TSimulator::setDepth(float depth){
	seqDepth = depth;
}

void TSimulator::setReadGroupName(std::string name){
	if(bamFileOpen)
		throw "Can not change read group name after opening BAM file!";
	readGroupName = name;
	bamAlignment.AddTag("RG", "Z", readGroupName);
}

void TSimulator::openBamFile(std::string filename){
	logfile->listFlush("Opening BAM file '" + filename + "' ...");

	if(bamFileOpen)
		throw "A BAM file is already open for writing!";

	bamFileName = filename;

	if(chromosomes.size() < 1)
		throw "Can not open a BAM file without specififed chromosomes!";

	BamTools::SamHeader header("");
	header.Version = "1.4";
	header.GroupOrder = "none";
	header.SortOrder = "coordinate";
	header.ReadGroups.Add(readGroupName);
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
		references.push_back(BamTools::RefData(chrIt->first, chrIt->second));
		header.Sequences.Add(BamTools::SamSequence(chrIt->first, chrIt->second));
	}

	if (!bamWriter.Open(bamFileName, header, references))
		throw "Failed to open BAM file '" + bamFileName + "'!";
	bamFileOpen = true;
	logfile->write(" done!");
}

void TSimulator::closeBamFile(){
	if(bamFileOpen){
		bamWriter.Close();

		//now generate bam index
		indexBamFile(bamFileName);
	}
	references.clear();
	bamFileOpen = false;
}

void TSimulator::indexBamFile(std::string & filename){
	logfile->listFlush("Creating index of BAM file '" + filename + "' ...");
	BamTools::BamReader reader;
	if(!reader.Open(filename))
		throw "Failed to open BAM file '" + filename + "' for indexing!";

	// create index for BAM file
	reader.CreateIndex(BamTools::BamIndex::STANDARD);

	//close BAM file
	reader.Close();
	logfile->write(" done!");
}

void TSimulator::openFastaFile(std::string filename){
	//open FASTA file for reference sequences
	logfile->list("Will write reference sequence to '" + filename + "'.");
	fasta.open(filename.c_str());
	if(!fasta)
		throw "Failed to open file '" + filename + "' for writing!";
	fastaOpen = true;
}

void TSimulator::closeFastaFile(){
	fasta.close();
	fastaOpen = false;
}

void TSimulator::simulateReferenceAndAlternativeSequenceCurChromosome(){
	logfile->listFlush("Simulating reference and alternative alleles ...");

	//initialize storage
	clearRefStorage();
	ref = new short[chrIt->second];
	alt = new short[chrIt->second];
	refInitialized = true;

	fasta << ">" << chrIt->first;
	for(int l=0; l<chrIt->second; ++l){
		ref[l] = randomGenerator->pickOne(4);
		alt[l] = (ref[l] + randomGenerator->pickOne(3)) % 4;
		//add to fasta
		if(l % 70 == 0)
			fasta << std::endl;
		fasta << ref[l];
	}
	fasta << std::endl;
	logfile->write(" done!");
}

void TSimulator::clearRefStorage(){
	if(refInitialized){
		delete[] ref;
		delete[] alt;
	}
}

int TSimulator::sampleQuality(){
	int qual = round(randomGenerator->getNormalRandom(meanQual, sdQual));
	if(qual > 126) qual = 126;
	if(qual < 33) qual = 33;
	return qual;
};

double TSimulator::dePhred(double x){
	return pow(10, -(x-33.0) / 10.0);
}

void TSimulator::initializeQualToErrorTable(){
	if(!qualToErroTableInitialized){
		qualToErroTable = new double[127];
		for(int i=0; i<33; ++i)
			qualToErroTable[i] = 1.0;
		for(int i=33; i<127; ++i)
			qualToErroTable[i] = dePhred(i);
	}
	qualToErroTableInitialized = true;
};

void TSimulator::simulateReads(int & numReads, long & pos, float* & altFreq){
	//TODO: Add PMD as general feature to be set prior to simulations.
	static short base;
	static int qual;
	//Whole sim process, numReads pooled and individual can be seen as frequency, no?
	for(int r=0; r<numReads; ++r){
		//simulate a read starting here
		bamAlignment.Position = pos;
		bamAlignment.QueryBases = "";
		bamAlignment.Qualities = "";
		for(long p=pos; p<pos+readLength; ++p){
			//sample base
			if(randomGenerator->getRand() < altFreq[p])
				base = alt[p];
			else
				base = ref[p];

			//sample quality and add error
			qual = sampleQuality();
			if(randomGenerator->getRand() < qualToErroTable[qual])
				base = (base + randomGenerator->pickOne(3) + 1) % 4;

			//add to bam alignment
			bamAlignment.Qualities += (char) qual;
			bamAlignment.QueryBases += toBase[base];

		}
		bamWriter.SaveAlignment(bamAlignment);
	}
}

void TSimulator::simulatePooledData(int sampleSize, SFS & sfs, std::string outname){
	//open BAM file
	openBamFile(outname + ".bam");

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	openFastaFile(filename);

	//prepare variables
	float* altFreq = NULL;
	long numReads;
	long chrLengthForStart;
	double probReadPerSite;
	int numReadsHere;
	long numReadsSimulated;
	initializeQualToErrorTable();

	//open frequency file
	filename = outname + "_frequencies.txt";
	std::ofstream freqFile(filename.c_str());

	//simulate sequences
	int refId = 0;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId){
		logfile->startIndent("Simulating chromosome " + chrIt->first + ":");

		//simulate reference and alternative sequence
		simulateReferenceAndAlternativeSequenceCurChromosome();

		//simulate alternative frequencies (and write to file)
		logfile->listFlush("Simulating alternative allele frequencies ...");
		delete[] altFreq;
		altFreq = new float[chrIt->second];
		for(int l=0; l<chrIt->second; ++l){
			altFreq[l] = sfs.getRandomFrequency(randomGenerator);
			freqFile << chrIt->first << "\t" << altFreq[l] << "\n";
		}
		logfile->write(" done!");

		//simulating reads
		numReads = chrIt->second * seqDepth / readLength;
		chrLengthForStart = chrIt->second - readLength;
		probReadPerSite = 1.0 / (double) chrLengthForStart;
		numReadsSimulated = 0;
		bamAlignment.RefID = refId;
		int prog;
		int oldProg = 0;
		std::string progressString = "Simulating about " + toString(numReads) + " reads ...";
		logfile->listFlush(progressString);
		for(long l=0; l<chrLengthForStart; ++l){
			//draw random number to get number of reads starting at this position
			numReadsHere = randomGenerator->getBiomialRand(probReadPerSite, numReads);

			//now simulate
			if(numReadsHere > 0){
				simulateReads(numReadsHere, l, altFreq);
				numReadsSimulated += numReadsHere;

				//report progress
				prog = 100.0 * (float) numReadsSimulated / (float) numReads;
				if(prog > oldProg){
					oldProg = prog;
					logfile->listOverFlush(progressString + "(" + toString(prog) + "%)");
				}
			}
		}
		logfile->overList(progressString + " done!  ");
		logfile->conclude("Simulated a total of " + toString(numReadsSimulated) + " reads.");
		logfile->endIndent();
	}

	//close stuff
	closeBamFile();
	closeFastaFile();
	freqFile.close();

	//clear memory
	delete[] altFreq;
}

void TSimulator::simulateSingleIndividual(double theta, std::string outname){
	//open BAM file
	openBamFile(outname + ".bam");

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	openFastaFile(filename);

	//prepare cumulative frequencies for genotypes
	float cumulFreq[3];
	cumulFreq[1] = theta;
	cumulFreq[2] = theta / 2.0;
	cumulFreq[0] = 1.0 - cumulFreq[1] - cumulFreq[2];
	if(cumulFreq[0] < 0.0)
		throw "Error when simulating data: current theta value too big, leads to too many mutations!";

	//prepare variables
	float* altFreq = NULL;
	long numReads;
	long chrLengthForStart;
	double probReadPerSite;
	int numReadsHere;
	long numReadsSimulated;
	initializeQualToErrorTable();

	//open frequency file
	filename = outname + "_frequencies.txt";
	std::ofstream freqFile(filename.c_str());

	//simulate sequences
	int refId = 0;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId){
		logfile->startIndent("Simulating chromosome " + chrIt->first + ":");

		//simulate reference and alternative sequence
		simulateReferenceAndAlternativeSequenceCurChromosome();

		//simulate genotypes according to theta
		logfile->listFlush("Simulating genotypes ...");
		delete[] altFreq;
		altFreq = new float[chrIt->second];
		for(int l=0; l<chrIt->second; ++l){
			altFreq[l] = randomGenerator->pickOne(3, cumulFreq) / 2.0;
		}
		logfile->write(" done!");

		//simulating reads
		numReads = chrIt->second * seqDepth / readLength;
		chrLengthForStart = chrIt->second - readLength;
		probReadPerSite = 1.0 / (double) chrLengthForStart;
		numReadsSimulated = 0;
		bamAlignment.RefID = refId;
		int prog;
		int oldProg = 0;
		std::string progressString = "Simulating about " + toString(numReads) + " reads ...";
		logfile->listFlush(progressString);
		for(long l=0; l<chrLengthForStart; ++l){
			//draw random number to get number of reads starting at this position
			numReadsHere = randomGenerator->getBiomialRand(probReadPerSite, numReads);

			//now simulate
			if(numReadsHere > 0){
				simulateReads(numReadsHere, l, altFreq);

				//report progress
				prog = 100.0 * (float) numReadsSimulated / (float) numReads;
				if(prog > oldProg){
					oldProg = prog;
					logfile->listOverFlush(progressString + "(" + toString(prog) + "%)");
				}
			}
		}
		logfile->overList(progressString + " done!  ");
		logfile->conclude("Simulated a total of " + toString(numReadsSimulated) + " reads.");
		logfile->endIndent();
	}

	//close stuff
	closeBamFile();
	closeFastaFile();
	freqFile.close();

	//clear memory
	delete[] altFreq;
}
