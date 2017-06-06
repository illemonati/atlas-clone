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
	oldOffset = 0;
	refInitialized = false;
	qualToErroTableInitialized = false;
	pmdInitialized = false;
	qualTransformationInitialized = false;

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
	float tmp[4];
	tmp[0] = 0.25; tmp[1] = 0.25; tmp[2] = 0.25; tmp[3] = 0.25;
	setBaseFreq(tmp);
	refLength = 0;
	ref = NULL;
};

void TSimulator::setQualityDistribution(double mean, double sd){
	meanQual = mean + 33.0; //add 33 to mean quality to get in in char
	sdQual = sd;
}

void TSimulator::initializeChromosomes(int numChr, long chrLength, bool haploid){
	chromosomes.clear();
	for(int i=0; i<numChr; ++i){
		chromosomes.push_back(TSimulatorChromosome("chr" + toString(i+1), chrLength, haploid));
	}
}

void TSimulator::initializeChromosomes(std::vector<long> & chrLength, std::vector<bool> haploid){
	chromosomes.clear();
	for(unsigned int i=0; i<chrLength.size(); ++i){
		chromosomes.push_back(TSimulatorChromosome("chr" + toString(i+1), chrLength[i], haploid[i]));
	}
}

void TSimulator::setReadLength(int length){
	if(qualTransformationInitialized && length != readLength)
		throw "TSimulator: Can not change read length after quality transformation was initialized!";
	readLength = length;
	bamAlignment.Length = readLength;
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', readLength));
}

void TSimulator::setDepth(float depth){
	seqDepth = depth;
}

void TSimulator::setBaseFreq(float* freq){
	float sum = 0.0;
	for(int i=0; i<4; ++i){
		baseFreq[i] = freq[i];
		sum += freq[i];
	}
	for(int i=0; i<4; ++i){
		baseFreq[i] /= sum;
	}
	cumulBaseFreq[0] = baseFreq[0];
	cumulBaseFreq[1] = cumulBaseFreq[0] + baseFreq[1];
	cumulBaseFreq[2] = cumulBaseFreq[1] + baseFreq[2];
	cumulBaseFreq[3] = 1.0;
}

void TSimulator::setReadGroupName(std::string name){
	if(bamFileOpen)
		throw "Can not change read group name after opening BAM file!";
	readGroupName = name;
	bamAlignment.AddTag("RG", "Z", readGroupName);
}

void TSimulator::setPMD(TPMD* PmdObject){
	pmdObject = PmdObject;
	pmdInitialized = true;
}


void TSimulator::setQualityTransformation(std::vector<double> & Betas){
	if(Betas.size() != 24)
		throw "Wrong size of beta vector when initializing quality transformation: need 24 values (quality, quality^2, pos, pos^2 and 20 contexts).";

	//copy betas
	beta = new double[24];
	for(int i=0; i<24; ++i)
		beta[i] = Betas[i];

	//precalculate stuff
	qualTermForTransformation = new double[127];
	for(int i=0; i<33; ++i)
		qualTermForTransformation[i] = 1.0;
	double tmp;
	for(int i=33; i<127; ++i){
		tmp = pow(10.0, -(double) (i - 33) / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}

	posTermForTransformation = new double[readLength];
	for(int i=0; i<readLength; ++i){
		posTermForTransformation[i] = beta[2] * i + beta[3] * i*i;
	}

	qualTransformationInitialized = true;
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
		references.push_back(BamTools::RefData(chrIt->name, chrIt->length));
		header.Sequences.Add(BamTools::SamSequence(chrIt->name, chrIt->length));
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
	filename += ".fai";
	fastaIndex.open(filename.c_str());
	if(!fastaIndex)
		throw "Failed to open file '" + filename + "' for writing!";
	oldOffset = 0;
	fastaOpen = true;
}

void TSimulator::closeFastaFile(){
	fasta.close();
	fastaIndex.close();
	fastaOpen = false;
}

void TSimulator::writeRefToFasta(){
	if(!fastaOpen)
		throw "Can not write reference to fasta: fasta file was never opend!";

	//write to fasta
	fasta << ">" << chrIt->name;
	for(int l=0; l<chrIt->length; ++l){
		if(l % 70 == 0)
			fasta << "\n";
		fasta << toBase[ref[l]];
	}
	fasta << "\n";

	//add to index
	std::string tmp = chrIt->name;
	oldOffset += chrIt->name.size() + 2;
	fastaIndex << extractBeforeWhiteSpace(tmp) << "\t" << chrIt->length << "\t" << oldOffset << "\t70\t71\n";
	oldOffset += chrIt->length + (int) (chrIt->length / 70);
	if(chrIt->length % 70 != 0) oldOffset += 1;
}

void TSimulator::simulateReferenceSequenceCurChromosome(){
	logfile->listFlush("Simulating reference alleles ...");

	//initialize storage
	initializeRefStorage();

	//simulate reference sequence
	for(int l=0; l<chrIt->length; ++l){
		ref[l] = randomGenerator->pickOne(4, cumulBaseFreq);
	}

	if(fastaOpen){
		writeRefToFasta();
	}
	logfile->write(" done!");
}

void TSimulator::initializeRefStorage(){
	if(refInitialized && refLength != chrIt->length)
		clearRefStorage();
	if(!refInitialized){
		ref = new short[chrIt->length];
		refInitialized = true;
		refLength = chrIt->length;
	}
}

void TSimulator::clearRefStorage(){
	if(refInitialized){
		delete[] ref;
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

int TSimulator::transformQuality(int & qual, int pos, int context){
	static double constant;
	static double tmp;
	static double q;

	constant = posTermForTransformation[pos] + beta[context+4] - qualTermForTransformation[qual];
	if(beta[1] == 0.0){
		q = -constant / beta[0];
	} else {
		tmp = sqrt(beta[0] * beta[0] - 4.0 * beta[1] * constant);
		q = (-tmp - beta[0]) / 2.0 / beta[1];
		//if(q < 0) q = (-tmp - beta[0]) / 2.0 / beta[1];
	}
	tmp = exp(q);
	q = -10.0 * log10(tmp / (1.0 - tmp));
	return (int) round(q) + 33;
}

/*
void TSimulator::simulateReads(int & numReads, long & pos, float* & altFreq){
	//TODO: switch pooled data to haplotype model as well!
	static short base;
	static int qual;	static int r;
	static long p;
	static int previousBase;

	for(r=0; r<numReads; ++r){
		//simulate a read starting here
		bamAlignment.Position = pos;
		bamAlignment.QueryBases = "";
		bamAlignment.Qualities = "";

		previousBase = 4; //means N
		for(p=0; p<readLength; ++p){
			//sample base
			if(randomGenerator->getRand() < altFreq[p+pos])
				base = alt[p+pos];
			else
				base = ref[p+pos];

			//apply PMD
			if(pmdInitialized){
				if(base == 1 ){ //means is C
					if(randomGenerator->getRand() < pmdObject->getProbCT(p))
						base = 3; //means T
				} else if(base == 2){ //means is G
					if(randomGenerator->getRand() < pmdObject->getProbGA(readLength - p - 1))
						base = 0; //means A
				}
			}

			//sample quality and add error
			qual = sampleQuality();
			if(randomGenerator->getRand() < qualToErroTable[qual])
				base = (base + randomGenerator->pickOne(3) + 1) % 4;

			//add to bam alignment
			if(qualTransformationInitialized){
				bamAlignment.Qualities += (char) transformQuality(qual, p, genoMap.getContext(previousBase, base));
				previousBase = base;
			} else
				bamAlignment.Qualities += (char) qual;
			bamAlignment.QueryBases += toBase[base];
		}
		bamWriter.SaveAlignment(bamAlignment);
	}
}
*/

void TSimulator::writeRead(long & pos, short* haplotype){
	static short base;
	static int qual;
	static long p;
	static int previousBase;

	//simulate a read starting here
	bamAlignment.Position = pos;
	bamAlignment.QueryBases = "";
	bamAlignment.Qualities = "";

	//choose haployt
	previousBase = 4; //means N
	for(p=0; p<readLength; ++p){
		//get true nucleotide
		base = haplotype[p + pos];

		//apply PMD
		if(pmdInitialized){
			if(base == 1 ){ //means is C
				if(randomGenerator->getRand() < pmdObject->getProbCT(p))
					base = 3; //means T
			} else if(base == 2){ //means is G
				if(randomGenerator->getRand() < pmdObject->getProbGA(readLength - p - 1))
					base = 0; //means A
			}
		}

		//sample quality and add error
		qual = sampleQuality();
		if(randomGenerator->getRand() < qualToErroTable[qual])
			base = (base + randomGenerator->pickOne(3) + 1) % 4;

		//add to bam alignment
		if(qualTransformationInitialized){
			bamAlignment.Qualities += (char) transformQuality(qual, p, genoMap.getContext(previousBase, base));
			previousBase = base;
		} else
			bamAlignment.Qualities += (char) qual;
		bamAlignment.QueryBases += toBase[base];
	}
	bamWriter.SaveAlignment(bamAlignment);
}

/*
//TODO: Need to switch to haplotype model
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
		logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

		//simulate reference and alternative sequence
		simulateReferenceAndAlternativeSequenceCurChromosome();

		//simulate alternative frequencies (and write to file)
		logfile->listFlush("Simulating alternative allele frequencies ...");
		delete[] altFreq;
		altFreq = new float[chrIt->length];
		for(int l=0; l<chrIt->length; ++l){
			altFreq[l] = sfs.getRandomFrequency(randomGenerator);
			freqFile << chrIt->name << "\t" << l+1 << altFreq[l] << "\n";
		}
		logfile->write(" done!");

		//simulating reads
		numReads = chrIt->length * seqDepth / readLength;
		chrLengthForStart = chrIt->length - readLength;
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
*/

void TSimulator::fillMutationTable(float** & mutTable, double theta){
	double expMinusTheta = exp(-theta);
	double sum;
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			mutTable[i][j] = baseFreq[i] * baseFreq[j] * (1.0 - expMinusTheta);
			sum += mutTable[i][j];
		}
		mutTable[i][i] += baseFreq[i] * expMinusTheta;

		//normalize within row
		sum = 0.0;
		for(int j=0; j<4; ++j){
			sum += mutTable[i][j];
		}
		for(int j=0; j<4; ++j){
			mutTable[i][j] /= sum;
		}

		//make cumulative
		mutTable[i][2] += mutTable[i][1];
		mutTable[i][3] += mutTable[i][2];
		mutTable[i][3] = 1.0;
	}
}

void TSimulator::simulateDiploidHaplotypesCurChromosome(short** haplotypes, float** & mutTable, const double & referenceDivergence){
	//initialize storage
	initializeRefStorage();

	//prepare cumulative probability distribution for reference
	float cumulRef[4];
	cumulRef[0] = 1.0 - referenceDivergence;
	cumulRef[1] = cumulRef[0] + referenceDivergence / 3.0;
	cumulRef[2] = cumulRef[1] + referenceDivergence / 3.0;
	cumulRef[3] = 1.0;

	//now simulate haplotypes
	if(chrIt->haploid){
		for(int l=0; l<chrIt->length; ++l){
			haplotypes[0][l] = randomGenerator->pickOne(4, cumulBaseFreq);
			haplotypes[1][l] = haplotypes[0][l];

			//decide on ref
			ref[l] = (haplotypes[0][l] + randomGenerator->pickOne(4, cumulRef)) % 4;
		}
	} else {
		for(int l=0; l<chrIt->length; ++l){
			haplotypes[0][l] = randomGenerator->pickOne(4, cumulBaseFreq);
			haplotypes[1][l] = randomGenerator->pickOne(4, mutTable[haplotypes[0][l]]);

			//decide on reference sequence
			if(haplotypes[0][l] == haplotypes[1][l])
				ref[l] = (haplotypes[0][l] + randomGenerator->pickOne(4, cumulRef)) % 4;
			else
				ref[l] = haplotypes[randomGenerator->pickOne(2)][l];
		}
	}

	if(fastaOpen)
		writeRefToFasta();
}

void TSimulator::writeTrueGenotypes(short** haplotypes, std::ofstream & genoFile){
	for(int l=0; l<chrIt->length; ++l)
		genoFile << chrIt->name << "\t" << l+1 << "\t" << toBase[haplotypes[0][l]] << "/" << toBase[haplotypes[1][l]] << "\n";
}

void TSimulator::writeInvariantSites(short** haplotypes, std::ofstream & genoFile){
	//0-based
	for(int l=0; l<chrIt->length; ++l){
		if(haplotypes[0][l] == haplotypes[1][l]){
			genoFile << chrIt->name << "\t" << l << "\t" << l+1 << "\t" << toBase[haplotypes[0][l]] << "\t" << toBase[haplotypes[1][l]] << "\n";
		}
	}
}



void TSimulator::simulateSingleIndividual(double theta, double referenceDivergence, std::string outname){
	//open BAM file
	openBamFile(outname + ".bam");

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	openFastaFile(filename);

	//prepare cumulative frequencies for genotypes
	float cumulFreq[3];
	cumulFreq[0] = 1.0 - theta;
	cumulFreq[1] = 1.0;
	if(cumulFreq[0] < 0.0)
		throw "Error when simulating data: current theta value too big, leads to too many mutations!";

	//prepare variables
	short** haplotypes = new short*[2];
	haplotypes[0] = NULL; haplotypes[1] = NULL;
	long numReads;
	long chrLengthForStart;
	double probReadPerSite;
	int numReadsHere;
	long numReadsSimulated;
	initializeQualToErrorTable();
	long oldSize = 0;
	int r;

	//open file for true genotypes
	filename = outname + "_trueGenotypes.txt";
	std::ofstream genoFile(filename.c_str());

	//open file for invariant positions
	filename = outname + "_invariantSites.txt";
	std::ofstream invariantSitesFile(filename.c_str());


	//prepare mutation table
	float** mutTable;
	mutTable = new float*[4];
	for(int i=0; i<4; ++i)
		mutTable[i] = new float[4];
	fillMutationTable(mutTable, theta);

	//simulate sequences
	int refId = 0;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId){
		logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

		//simulate genotypes according to theta
		logfile->listFlush("Simulating genotypes ...");

		if(chrIt->length != oldSize){
			delete[] haplotypes[0];
			delete[] haplotypes[1];
			haplotypes[0] = new short[chrIt->length];
			haplotypes[1] = new short[chrIt->length];
		}

		simulateDiploidHaplotypesCurChromosome(haplotypes, mutTable, referenceDivergence);
		writeTrueGenotypes(haplotypes, genoFile);
		writeInvariantSites(haplotypes, invariantSitesFile);

		logfile->write(" done!");

		//simulating reads
		numReads = chrIt->length * seqDepth / readLength;
		chrLengthForStart = chrIt->length - readLength;
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
				numReadsSimulated += numReadsHere;
				for(r=0; r<numReadsHere; ++r)
					writeRead(l, haplotypes[randomGenerator->pickOne(2)]);

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

		//set old size
		oldSize = chrIt->length;
	}

	//close stuff
	closeBamFile();
	closeFastaFile();
	genoFile.close();

	//clear memory
	for(int i=0; i<4; ++i)
		delete[] mutTable[i];
	delete[] mutTable;
	delete[] haplotypes[0];
	delete[] haplotypes[1];
	delete[] haplotypes;
}
