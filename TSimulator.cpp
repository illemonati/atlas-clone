/*
 * TSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#include "TSimulator.h"

//---------------------------------------------------
//TSimulatorBamFile
//---------------------------------------------------
void TSimulatorBamFile::open(std::string Filename, const std::string & readGroupName, std::vector<TSimulatorChromosome> & chromosomes){
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
	header.ReadGroups.Add(readGroupName);
	for(std::vector<TSimulatorChromosome>::iterator chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
		references.push_back(BamTools::RefData(chrIt->name, chrIt->length));
		header.Sequences.Add(BamTools::SamSequence(chrIt->name, chrIt->length));
	}

	if (!bamWriter.Open(filename, header, references))
		throw "Failed to open BAM file '" + filename + "'!";
	isOpen = true;
	logfile->write(" done!");
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
	logfile->write(" done!");
}

//---------------------------------------------------
//TSimulatorReference
//---------------------------------------------------
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
}

void TSimulatorReference::closeFastaFile(){
	if(fastaOpen){
		fasta.close();
		fastaIndex.close();
	}
	fastaOpen = false;
}

void TSimulatorReference::writeRefToFasta(){
	if(fastaOpen){
		//write to fasta
		fasta << ">" << chrName;
		for(int l=0; l<chrLength; ++l){
			if(l % 70 == 0)
				fasta << "\n";
			fasta << toBase[ref[l]];
		}
		fasta << "\n";

		//add to index
		std::string tmp = chrName;
		oldOffset += chrName.size() + 2;
		fastaIndex << extractBeforeWhiteSpace(tmp) << "\t" << chrLength << "\t" << oldOffset << "\t70\t71\n";
		oldOffset += chrLength + (int) (chrLength / 70);
		if(chrLength % 70 != 0) oldOffset += 1;
	}
}


void TSimulatorReference::simulateReferenceSequenceCurChromosome(TRandomGenerator * randomGenerator, float* cumulBaseFreq){
	logfile->listFlush("Simulating reference alleles ...");

	if(!storageInitialized)
		throw "Can not simulate reference sequence, no chromosome set!";

	//simulate reference sequence
	for(int l=0; l<chrLength; ++l){
		ref[l] = randomGenerator->pickOne(4, cumulBaseFreq);
	}

	if(fastaOpen){
		writeRefToFasta();
	}
	logfile->write(" done!");
}


//---------------------------------------------------
//TSimulator
//---------------------------------------------------

TSimulator::TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator, std::vector<float> & Freq){
	logfile = Logfile;
	randomGenerator = RandomGenerator;
	bamFileOpen = false;
	baseFreq = Freq;
	qualToErroTableInitialized = false;
	pmdInitialized = false;
	qualTransformationInitialized = false;
	initializeQualToErrorTable();

	//set default parameters
	setReadLength(100);
	setQualityDistribution(30.0, 10.0);
	setDepth(10.0);
	bamAlignment.Name = "*";
	bamAlignment.MapQuality = 50;
	setReadGroupName("SimReadGroup");

	//helper tools
	toBase[0] = 'A'; toBase[1] = 'C'; toBase[2] = 'G'; toBase[3] = 'T';
	setBaseFreq();
};



void TSimulator::setQualityDistribution(double mean, double sd){
	meanQual = mean + 33.0; //add 33 to mean quality to get in in char
	sdQual = sd;
}

void TSimulator::setMaxQual(int MaxQual){
	maxQual = MaxQual;
}

void TSimulator::initializeChromosomes(int numChr, long chrLength, bool haploid){
	chromosomes.clear();
	for(int i=0; i<numChr; ++i){
		chromosomes.push_back(TSimulatorChromosome("chr" + toString(i+1), i, chrLength, haploid));
	}
}

void TSimulator::initializeChromosomes(std::vector<long> & chrLength, std::vector<bool> haploid){
	chromosomes.clear();
	for(unsigned int i=0; i<chrLength.size(); ++i){
		chromosomes.push_back(TSimulatorChromosome("chr" + toString(i+1), i, chrLength[i], haploid[i]));
	}
}

void TSimulator::setReadLength(int length){
	if(qualTransformationInitialized && length != readLength)
		throw "TSimulator: Can not change read length after quality transformation was initialized!";
	readLength = length;
	bamAlignment.Length = readLength;
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', readLength));
}

void TSimulator::setDepth(float depth){
	seqDepth = depth;
}

void TSimulator::setBaseFreq(){
	float sum = 0.0;
	for(int i=0; i<4; ++i){
		sum += baseFreq[i];
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
	for(int i=0; i<34; ++i)
		qualTermForTransformation[i] = 100.0;
	double tmp;
	for(int i=34; i<127; ++i){
		tmp = pow(10.0, -(double) (i - 33) / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}

	posTermForTransformation = new double[readLength];
	for(int i=0; i<readLength; ++i){
		posTermForTransformation[i] = beta[2] * i + beta[3] * i*i;
	}

	qualTransformationInitialized = true;
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
	return -10.0 * log10(tmp / (1.0 + tmp)) + 33.0;
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


//simulating reads
void TSimulator::simulateReadsFromHaplotypes(std::vector<TSimulatorChromosome>::iterator & thisChr, short** haplotypes, TSimulatorBamFile & bamFile){
	//Initialize probabilities to simulate reads
	long numReads = thisChr->length * seqDepth / readLength;
	long chrLengthForStart = thisChr->length - readLength;
	double probReadPerSite = 1.0 / (double) chrLengthForStart;
	long numReadsSimulated = 0;
	int numReadsHere;
	int r;

	//prepare bam alignment
	bamAlignment.RefID = thisChr->refID;

	//initialize progress reporting
	int prog;
	int oldProg = 0;
	std::string progressString = "Simulating about " + toString(numReads) + " reads ...";
	logfile->listFlush(progressString);

	//now simulate
	for(long l=0; l<chrLengthForStart; ++l){
		//draw random number to get number of reads starting at this position
		numReadsHere = randomGenerator->getBiomialRand(probReadPerSite, numReads);

		//now simulate
		if(numReadsHere > 0){
			numReadsSimulated += numReadsHere;
			for(r=0; r<numReadsHere; ++r)
				writeRead(l, haplotypes[randomGenerator->pickOne(2)], bamFile);

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
}

void TSimulator::writeRead(long & pos, short* haplotype, TSimulatorBamFile & bamFile){
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
			int transQual = transformQuality(qual, p, genoMap.getContext(previousBase, base));
			//cap at highest possible illumina score
			if(transQual > maxQual)	bamAlignment.Qualities += (char) maxQual;
			else bamAlignment.Qualities += (char) transQual;
			previousBase = base;
		} else {
			if(qual > maxQual) bamAlignment.Qualities += (char) maxQual;
			else bamAlignment.Qualities += (char) qual;
		}
		bamAlignment.QueryBases += toBase[base];
	}
	bamFile.saveAlignment(bamAlignment);
}


//-------------------------------------------------------
//Functions to simulate a single individual
//-------------------------------------------------------
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

void TSimulator::simulateDiploidHaplotypesCurChromosome(short** haplotypes, float** & mutTable, short* ref, const double & referenceDivergence){
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
	TSimulatorBamFile bamFile(outname + ".bam", readGroupName, chromosomes, logfile);
	bamFileOpen = true;

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	TSimulatorReference referenceObj(filename, toBase, logfile);

	//prepare cumulative frequencies for genotypes
	float cumulFreq[3];
	cumulFreq[0] = 1.0 - theta;
	cumulFreq[1] = 1.0;
	if(cumulFreq[0] < 0.0)
		throw "Error when simulating data: current theta value too big, leads to too many mutations!";

	//prepare haplotypes and
	TSimulatorHaplotypes haplotypes(1);

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
		//update reference storage and update haplotype lengths
		referenceObj.setChr(chrIt->name, chrIt->length);
		haplotypes.setLength(chrIt->length);

		//simulate genotypes
		logfile->listFlush("Simulating genotypes ...");
		simulateDiploidHaplotypesCurChromosome(haplotypes.getHaplotypesFirstIndividual(), mutTable, referenceObj.getPointerToRef(), referenceDivergence);
		haplotypes.writeGenotypes(genoFile, chrIt->name, toBase);
		writeInvariantSites(haplotypes.getHaplotypesFirstIndividual(), invariantSitesFile);
		logfile->write(" done!");

		//now simulate and write reads
		simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesFirstIndividual(), bamFile);

		//end of chromosome
		logfile->endIndent();
	}

	//close stuff
	bamFile.close();
	genoFile.close();

	//clear memory
	for(int i=0; i<4; ++i)
		delete[] mutTable[i];
	delete[] mutTable;
}



//--------------------------------------------------------------------
//Functions to simulate a pair of individuals according to distances
//--------------------------------------------------------------------
void TSimulatorGenotypeCombination::fillTables(std::vector<double> & phis, std::vector<float> & baseFreq){
	//file cumulative frequencies of cases (phis)
	double sum = 0.0;
	int genoCase = 0;
	for(std::vector<double>::iterator it = phis.begin(); it != phis.end(); ++it, ++genoCase){
		sum += *it;
		cumulGenoCaseFrequencies[genoCase] = sum;
	}
	cumulGenoCaseFrequencies[genoCase] = 1.0;
	if(fabs(sum - 1.0) > 0.0000000001)
		throw "Phis do not sum to 1.0! They sum to " + toString(sum) + ".";

	//prepare genotype frequency tables for each case
	cumulGenoCombinationFreq = new double*[9];
	genoTrans = new short**[9];

	//some variables
	int a,b,c,d;
	int index;

	//case 0: aa/aa
	//-----------------------------------------
	cumulGenoCombinationFreq[0] = new double[4];
	numGenotypeCombinations[0] = 4;
	genoTrans[0] = new short*[4];
	sum = 0.0;
	for(a=0; a<4; ++a){
		sum += baseFreq[a];
		cumulGenoCombinationFreq[0][a] = sum;
		genoTrans[0][a] = new short[4];
		genoTrans[0][a][0] = a;
		genoTrans[0][a][1] = a;
		genoTrans[0][a][2] = a;
		genoTrans[0][a][3] = a;
	}

	//cases 1 to 4: aa/ab, ab/aa, aa/bb, ab/ab
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	double* cumul = new double[12];
	index = 0;
	sum = 0.0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				sum += baseFreq[a] * baseFreq[b];
				cumul[index] = sum;
				++index;
			}
		}
	}
	//normalize
	for(index=0; index<12; ++index)
		cumul[index] /= sum;

	//now initialize
	for(int ca = 1; ca<5; ++ca){
		cumulGenoCombinationFreq[ca] = new double[12];
		numGenotypeCombinations[ca] = 12;
		genoTrans[ca] = new short*[12];
		for(index=0; index<12; ++index){
			genoTrans[ca][index] = new short[4];
			cumulGenoCombinationFreq[ca][index] = cumul[index];
		}
	}

	//assign genotype translations
	index = 0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				genoTrans[1][index][0] = a; genoTrans[1][index][1] = a; genoTrans[1][index][2] = a; genoTrans[1][index][3] = b;
				genoTrans[2][index][0] = a; genoTrans[2][index][1] = b; genoTrans[2][index][2] = a; genoTrans[2][index][3] = a;
				genoTrans[3][index][0] = a; genoTrans[3][index][1] = a; genoTrans[3][index][2] = b; genoTrans[3][index][3] = b;
				genoTrans[4][index][0] = a; genoTrans[4][index][1] = b; genoTrans[4][index][2] = a; genoTrans[4][index][3] = b;
				++index;
			}
		}
	}

	delete[] cumul;

	//cases 5 to 7: ab/ac, aa/bc, ab/cc
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	cumul = new double[24];
	index = 0;
	sum = 0.0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						sum += baseFreq[a] * baseFreq[b] * baseFreq[c];
						cumul[index] = sum;
						++index;
					}
				}
			}
		}
	}
	//normalize
	for(index=0; index<24; ++index)
		cumul[index] /= sum;

	//now initialize
	for(int ca = 5; ca<8; ++ca){
		cumulGenoCombinationFreq[ca] = new double[24];
		numGenotypeCombinations[ca] = 24;
		genoTrans[ca] = new short*[24];
		for(index=0; index<24; ++index){
			cumulGenoCombinationFreq[ca][index] = cumul[index];
			genoTrans[ca][index] = new short[4];
		}
	}

	//assign genotype translations
	index = 0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						genoTrans[5][index][0] = a; genoTrans[5][index][1] = b; genoTrans[5][index][2] = a; genoTrans[5][index][3] = c;
						genoTrans[6][index][0] = a; genoTrans[6][index][1] = a; genoTrans[6][index][2] = b; genoTrans[6][index][3] = c;
						genoTrans[7][index][0] = a; genoTrans[7][index][1] = b; genoTrans[7][index][2] = c; genoTrans[7][index][3] = c;
						++index;
					}
				}
			}
		}
	}

	delete[] cumul;

	//case 8: ab/cd
	//-----------------------------------------
	cumulGenoCombinationFreq[8] = new double[24];
	numGenotypeCombinations[8] = 24;
	genoTrans[8] = new short*[24];
	index = 0;
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						for(d=0; d<4; ++d){
							if(d!=a && d!=b && d!=c){
								cumulGenoCombinationFreq[8][index] = (double) (index+1.0) / 24.0;
								genoTrans[8][index] = new short[4];
								genoTrans[8][index][0] = a;
								genoTrans[8][index][1] = b;
								genoTrans[8][index][2] = c;
								genoTrans[8][index][3] = d;
								++index;
							}
						}
					}
				}
			}
		}
	}

	//prepare haplotype order table to randomly pick
	//----------------------------------------------
	orderLookup = new short*[4];
	orderLookup[0] = new short[4];
	orderLookup[0][0] = 0; orderLookup[0][1] = 1; orderLookup[0][2] = 2; orderLookup[0][3] = 3;
	orderLookup[1] = new short[4];
	orderLookup[1][0] = 0; orderLookup[1][1] = 1; orderLookup[1][2] = 3; orderLookup[1][3] = 2;
	orderLookup[2] = new short[4];
	orderLookup[2][0] = 1; orderLookup[2][1] = 0; orderLookup[2][2] = 2; orderLookup[2][3] = 3;
	orderLookup[3] = new short[4];
	orderLookup[3][0] = 1; orderLookup[3][1] = 0; orderLookup[3][2] = 3; orderLookup[3][3] = 2;

	//set as initialized
	tablesInitialized = true;
}

void TSimulatorGenotypeCombination::deleteTables(){
	if(tablesInitialized){
		for(int i=0; i<9; ++i){
			delete[] cumulGenoCombinationFreq[i];
			for(int j=0; j<numGenotypeCombinations[i]; ++j)
				delete[] genoTrans[i][j];
			delete[] genoTrans[i];
		}
		delete[] cumulGenoCombinationFreq;
		delete[] genoTrans;

		//and lookup
		for(int i=0; i<4; ++i)
			delete[] orderLookup[i];
		delete[] orderLookup;

		tablesInitialized = false;
	}
}

void TSimulatorGenotypeCombination::simulateHaplotypes(short** haplotypesInd0, short** haplotypesInd1, short* ref, float referenceDivergence, long length, TRandomGenerator* randomGenerator){
	//prepare cumulative probability distribution for reference
	float cumulRef[4];
	cumulRef[0] = 1.0 - referenceDivergence;
	cumulRef[1] = cumulRef[0] + referenceDivergence / 3.0;
	cumulRef[2] = cumulRef[1] + referenceDivergence / 3.0;
	cumulRef[3] = 1.0;

	//run across loci
	int c, g, o, r;
	for(long l=0; l<length; ++l){
		//pick a case
		c = randomGenerator->pickOne(9, cumulGenoCaseFrequencies);

		//pick genotypes
		g = randomGenerator->pickOne(numGenotypeCombinations[c], cumulGenoCombinationFreq[c]);

		//pick order
		o = randomGenerator->pickOne(4);

		//assign to haplotypes
		haplotypesInd0[0][l] = genoTrans[c][g][orderLookup[o][0]];
		haplotypesInd0[1][l] = genoTrans[c][g][orderLookup[o][1]];
		haplotypesInd1[0][l] = genoTrans[c][g][orderLookup[o][2]];
		haplotypesInd1[1][l] = genoTrans[c][g][orderLookup[o][3]];

		//simulate reference
		if(c == 0){
			ref[l] = (genoTrans[c][g][0] + randomGenerator->pickOne(4, cumulRef)) % 4;
		} else {
			r = randomGenerator->pickOne(4);
			ref[l] = genoTrans[c][g][r];
		}
	}
}


void TSimulator::simulateIndividualPair(std::vector<double> & phis, double referenceDivergence, std::string outname){
	//open BAM files
	logfile->startIndent("Opening bam files for writing:");
	TSimulatorBamFile bamFile0(outname + "_ind0.bam", readGroupName, chromosomes, logfile);
	TSimulatorBamFile bamFile1(outname + "_ind1.bam", readGroupName, chromosomes, logfile);
	bamFileOpen = true;
	logfile->endIndent();

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	TSimulatorReference referenceObj(filename, toBase, logfile);

	//prepare haplotypes and
	TSimulatorHaplotypes haplotypes(2);

	//open file for true genotypes
	filename = outname + "_trueGenotypes.txt";
	std::ofstream genoFile(filename.c_str());

	//initialize genotype combination tables
	TSimulatorGenotypeCombination genoComb(phis, baseFreq);

	//simulate sequences
	int refId = 0;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId){
		logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

		//update reference storage and update haplotype lengths
		referenceObj.setChr(chrIt->name, chrIt->length);
		haplotypes.setLength(chrIt->length);

		//simulate genotypes according to their distributions
		logfile->listFlush("Simulating genotypes ...");
		genoComb.simulateHaplotypes(haplotypes.getHaplotypesOfIndividual(0), haplotypes.getHaplotypesOfIndividual(1), referenceObj.getPointerToRef(), referenceDivergence, chrIt->length, randomGenerator);
		haplotypes.writeGenotypes(genoFile, chrIt->name, toBase);
		logfile->write(" done!");

		//simulating reads
		logfile->startIndent("Simulating reads for individual 1:");
		simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesOfIndividual(0), bamFile0);
		logfile->endIndent();
		logfile->startIndent("Simulating reads for individual 2:");
		simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesOfIndividual(1), bamFile1);
		logfile->endIndent();

		logfile->endIndent();
	}

	//close stuff
	bamFile0.close();
	bamFile1.close();
	genoFile.close();
}


//--------------------------------------------------------------------
//Functions to simulate pooled data
//--------------------------------------------------------------------
//TODO: Need to switch to haplotype model

/*
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
