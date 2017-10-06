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
	header.ReadGroups.Add(readGroupName + "\tPU:UNKNOWN\tLB:UNKNOWN\tSM:Sim1\tCN:UNKNOWN\tPL:ILLUMINA");
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

TSimulator::TSimulator(TLog* Logfile, TRandomGenerator* RandomGenerator, TParameters & params){
	logfile = Logfile;
	randomGenerator = RandomGenerator;

	//set basic things
	bamFileOpen = false;
	qualToErroTableInitialized = false;
	pmdInitialized = false;
	readLengthDistInitialized = false;
	initializeQualToErrorTable();

	//read basic simulation settings
	logfile->startIndent("Reading simulation parameters:");

	//depth
	float depth = params.getParameterDoubleWithDefault("depth", 10.0);
	logfile->list("Will simulate to an average depth of " + toString(depth) + ".");
	setDepth(depth);

	//read length
	std::string tmp = params.getParameterStringWithDefault("readLength", "100");
	setReadLength(tmp);
	logfile->list(readLengthDist->getFunctionString());
	logfile->conclude("Implies an average read length of " + toString(readLengthDist->mean()));

	//base frequencies
	std::vector<float> freq;
	tmp = params.getParameterStringWithDefault("baseFreq", "0.25,0.25,0.25,0.25");
	fillVectorFromString(tmp, freq, ',');
	if(freq.size() != 4) throw "baseFreq vector must have size = 4!";
	setBaseFreq(freq);

	//reference divergence
	referenceDivergence = params.getParameterDoubleWithDefault("refDiv", 0.01);
	logfile->list("Will simulate data with reference divergence = " + toString(referenceDivergence) + ".");

	//quality distribution
	double mQ = params.getParameterDoubleWithDefault("meanQual", 30);
	double sdQ = params.getParameterDoubleWithDefault("sdQual", 10);
	int maxQ = params.getParameterDoubleWithDefault("maxQual", 500);
	setQualityDistribution(mQ, sdQ, maxQ);
	logfile->list("Will simulate normal distributed quality scores with mean = " + toString(meanQual) + " and sd = " + toString(sdQual) + ", capped at " + toString(maxQual) + ".");

	//quality transformation
	initializeQualityTransform(params);

	//PMD
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		pmdInitialized = true;
		pmdObject = new TPMD(params, logfile);
	} else pmdInitialized = false;


	//chromosomes
	initializeChromosomes(params, logfile);
	logfile->endIndent();

	//set default parameters
	bamAlignment.Name = "*";
	bamAlignment.MapQuality = 50;
	setReadGroupName("SimReadGroup");

	//helper tools
	toBase[0] = 'A'; toBase[1] = 'C'; toBase[2] = 'G'; toBase[3] = 'T';
};

void TSimulator::setQualityDistribution(double mean, double sd, int MaxQual){
	meanQual = mean + 33.0; //add 33 to mean quality to get in in char
	sdQual = sd;
	maxQual = MaxQual;

}

void TSimulator::initializeQualityTransform(TParameters & params){
	std::vector<std::string> string_vec;

	if(params.parameterExists("qualTransform")){
		params.fillParameterIntoVector("qualTransform", string_vec, ',');
		std::vector<double> beta;
		repeatIndexes(string_vec, beta);
		if(beta.size() != 24)
			throw "Wrong number of beta values for quality transformation (" + toString(beta.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";
		std::string s = concatenateString(beta, ",");
		logfile->list("Will transform qualities with beta = {" + s + "}");
		qualityTransformation = new TSimulatorRecalTransform(beta, readLengthDist);
	} else if(params.parameterExists("recal")){
		std::string filename = params.getParameterString("recal");
		logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");
		std::ifstream file(filename.c_str());
		if(!file) throw "Failed to open file '" + filename + "' for reading!";

		//tmp variables for reading
		std::string tmp;
		std::vector<std::string> vec;
		std::vector<double> beta;

		//skip header
		std::getline(file, tmp);

		//parse file to read details for each read group
		std::getline(file, tmp);
		fillVectorFromString(tmp, vec, "\t");
		logfile->done();

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 25) throw "Found " + toString(vec.size()) + " instead of 25 columns in '" + filename + "' on line 2!";
			for(int i=1; i < 25; ++i) beta.push_back(stringToFloat(vec[i]));
			logfile->done();
			if(beta.size() != 24)
				throw "Wrong number of beta values for quality transformation (" + toString(beta.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";
			std::string s = concatenateString(beta, ",");
			logfile->list("Will transform qualities with beta = {" + s + "}");
			//setQualityTransformation(beta);
			qualityTransformation = new TSimulatorRecalTransform(beta, readLengthDist);
		}
	} else {
		qualityTransformation = new TSimulatorQuality(readLengthDist);
		logfile->list("Will write quality scores into BAM without distortion");
	}
}

void TSimulator::initializeChromosomes(TParameters & params, TLog* logfile){
	std::vector<std::string> string_vec;
	std::vector<long> chrLength;
	params.fillParameterIntoVectorWithDefault("chrLength", string_vec, ',', "1000000");
	repeatIndexes(string_vec, chrLength);
	std::vector<int> ploidy;
	if(params.parameterExists("ploidy")){
		params.fillParameterIntoVector("ploidy", string_vec, ',');
		repeatIndexes(string_vec, ploidy);
	} else {
		for(size_t i=0; i<chrLength.size(); ++i)
			ploidy.push_back(2);
	}
	if(ploidy.size() != chrLength.size())
		throw "List of chromosome lengths and ploidies differ in length!";
	std::vector<bool> haploid;
	for(std::vector<int>::iterator it=ploidy.begin(); it!=ploidy.end(); ++it){
		if(*it == 1) haploid.push_back(true);
		else if(*it == 2) haploid.push_back(false);
		else throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!";
	}

	if(chrLength.size() < 1)
		throw "Issue understanding length of chromosomes!";
	if(chrLength.size() == 1){
		int numChr = params.getParameterIntWithDefault("numChr", 1);
		std::string text = "Will simulate " + toString(numChr) ;
		if(haploid[0]) text += " haploid";
		else text += " diploid";
		text += " chromosome(s) of length " + toString(chrLength[0]) + " each.";
		logfile->list(text);
		initializeChromosomes(numChr, chrLength[0], haploid[0]);
	} else {
		logfile->startIndent("Will simulate " + toString(chrLength.size()) + " chromosome(s) of the following length:");
		std::vector<bool>::iterator hIt=haploid.begin();
		std::string text;
		for(std::vector<long>::iterator it=chrLength.begin(); it!=chrLength.end(); ++it, ++hIt){
			text = toString(*it) + " (";
			if(*hIt) text += "haploid)";
			else text += "diploid)";
			logfile->list(text);
		}
		initializeChromosomes(chrLength, haploid);
		logfile->endIndent();
	}
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

void TSimulator::setReadLength(std::string s){
//	if(qualTransformationInitialized)
//		throw "TSimulator: Can not change read length after quality transformation was initialized!";

	if(readLengthDistInitialized)
		delete readLengthDist;

	//check if it is a specific function
	size_t pos = s.find("(");
	std::string tmp;

	if(pos != std::string::npos){
		//is a function
		std::string type = s.substr(0, pos);
		s.erase(0, pos);
		if(type == "gamma")
			readLengthDist = new TSimulatorReadLengthGamma(randomGenerator, s);
		else if(type == "gammaMode")
			readLengthDist = new TSimulatorReadLengthGammaMode(randomGenerator, s);
		else throw "Unknown read length distribution '" + type + "'!";
	} else {
		//fixed length
		readLengthDist = new TSimulatorReadLength(randomGenerator, s);
	}
	readLengthDistInitialized = true;

	//check simulation effort
	if(readLengthDist->probAcceptance() < 0.9)
		logfile->warning("The chosen read length distribution will only result in " + toString(readLengthDist->probAcceptance()) + " of draws being accepted.");
}

void TSimulator::setDepth(float depth){
	seqDepth = depth;
}

void TSimulator::setBaseFreq(std::vector<float> & freq){
	float sum = 0.0;
	for(int i=0; i<4; ++i){
		baseFreq[i] = freq[i];
		sum += baseFreq[i];
	}
	for(int i=0; i<4; ++i){
		baseFreq[i] /= sum;
	}
	cumulBaseFreq[0] = baseFreq[0];
	cumulBaseFreq[1] = cumulBaseFreq[0] + baseFreq[1];
	cumulBaseFreq[2] = cumulBaseFreq[1] + baseFreq[2];
	cumulBaseFreq[3] = 1.0;

	logfile->list("Simulating with base frequencies A:" + toString(baseFreq[0]) + " C:" + toString(baseFreq[1])+ " G:" + toString(baseFreq[1])+ " T:" + toString(baseFreq[1]));
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
	if(4.0 * beta[1] * constant > beta[0] * beta[0]) throw "beta[0]^2 cannot be smaller than 4beta[1](position + context constants)";
	if(beta[1] == 0.0){
		q = -constant / beta[0];
	} else {
		tmp = sqrt(beta[0] * beta[0] - 4.0 * beta[1] * constant);
		q = (-tmp - beta[0]) / 2.0 / beta[1];
		//if(q < 0) q = (-tmp - beta[0]) / 2.0 / beta[1];
	}

	tmp = exp(q);
	if(tmp == 0) throw "choose different quality transformation parameters! tmp == 0";
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
void TSimulator::simulateReadsFromHaplotypes(std::vector<TSimulatorChromosome>::iterator & thisChr, short** haplotypes, TSimulatorBamFile & bamFile, std::string extraProgressText){
	//Initialize probabilities to simulate reads
	long numReads = thisChr->length * seqDepth / readLengthDist->mean();
	long chrLengthForStart = thisChr->length - readLengthDist->max();
	double probReadPerSite = 1.0 / (double) chrLengthForStart;
	long numReadsSimulated = 0;
	int numReadsHere;
	int r;

	//prepare bam alignment
	bamAlignment.RefID = thisChr->refID;

	//initialize progress reporting
	int prog;
	int oldProg = 0;
	std::string progressString = "Simulating about " + toString(numReads) + " reads" + extraProgressText + " ...";
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
	readLengthContainer rl;

	//pick a fragment and read length
	readLengthDist->sample(rl);

	bamAlignment.Length = rl.readLength;
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', rl.readLength));

	//simulate a read starting here
	bamAlignment.Position = pos;
	bamAlignment.QueryBases = "";
	bamAlignment.Qualities = "";

	//choose haplotype
	previousBase = 4; //means N
	for(p=0; p<rl.readLength; ++p){
		//get true nucleotide
		base = haplotype[p + pos];

		//apply PMD
		if(pmdInitialized){
			if(base == 1 ){ //means is C
				if(randomGenerator->getRand() < pmdObject->getProbCT(p))
					base = 3; //means T
			} else if(base == 2){ //means is G
				if(randomGenerator->getRand() < pmdObject->getProbGA(rl.fragmentLength - p - 1))
					base = 0; //means A
			}
		}

		//sample quality and add error
		qual = sampleQuality();
		if(randomGenerator->getRand() < qualToErroTable[qual])
			base = (base + randomGenerator->pickOne(3) + 1) % 4;

		//add to bam alignment
		//int returnQual(int & qual, int & pos, TGenotypeMap & genoMap, int & previousBase, int & base);
		std::cout << "original qual " << qual << std::endl;
		int transQual = qualityTransformation->returnQual(qual, p, genoMap.getContext(previousBase, base));
		std::cout << "transqual " << transQual << std::endl;
		if(transQual > maxQual)	bamAlignment.Qualities += (char) maxQual;
		else bamAlignment.Qualities += (char) transQual;
		previousBase = base;
		bamAlignment.QueryBases += toBase[base];


		/*if(qualTransformationInitialized){
			int transQual = transformQuality(qual, p, genoMap.getContext(previousBase, base));
			//cap at highest possible illumina score
			if(transQual > maxQual)	bamAlignment.Qualities += (char) maxQual;
			else bamAlignment.Qualities += (char) transQual;
			previousBase = base;
		} else {
			if(qual > maxQual) bamAlignment.Qualities += (char) maxQual;
			else bamAlignment.Qualities += (char) qual;
		}
		bamAlignment.QueryBases += toBase[base];*/
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
		mutTable[i][1] += mutTable[i][0];
		mutTable[i][2] += mutTable[i][1];
		mutTable[i][3] = 1.0;
	}
}

void TSimulator::simulateDiploidHaplotypesCurChromosome(short** haplotypes, float** & mutTable, short* ref){
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

void TSimulator::writeBEDFiles(short** haplotypes, gz::ogzstream & invariantSitesFile, gz::ogzstream & variantSitesFile){
	//0-based
	for(int l=0; l<chrIt->length; ++l){
		if(haplotypes[0][l] == haplotypes[1][l]){
			invariantSitesFile << chrIt->name << "\t" << l << "\t" << l+1 << "\t" << toBase[haplotypes[0][l]] << "\t" << toBase[haplotypes[1][l]] << "\n";
		} else variantSitesFile << chrIt->name << "\t" << l << "\t" << l+1 << "\t" << toBase[haplotypes[0][l]] << "\t" << toBase[haplotypes[1][l]] << "\n";

	}
}

void TSimulator::simulateSingleIndividual(double theta, std::string & outname){
	std::vector<double> thetas;
	for(unsigned int i=0; i<chromosomes.size(); ++i)
		thetas.push_back(theta);
	simulateSingleIndividual(thetas, outname);
}

void TSimulator::simulateSingleIndividual(std::vector<double> theta, std::string & outname){
	//one thet aper chromosome
	if(theta.size() != chromosomes.size())
		throw "Number of theta values provided does not match number of chromosomes to simulate!";

	//open BAM file
	TSimulatorBamFile bamFile(outname + ".bam", readGroupName, chromosomes, logfile);
	bamFileOpen = true;

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	TSimulatorReference referenceObj(filename, toBase, logfile);

	//prepare haplotypes and
	TSimulatorHaplotypes haplotypes(1);

	//open file for true genotypes
	filename = outname + "_trueGenotypes.txt.gz";
	gz::ogzstream genoFile(filename.c_str());

	//open file for invariant positions
	filename = outname + "_invariantSites.bed.gz";
	gz::ogzstream invariantSitesFile(filename.c_str());

	//open file for variant positions
	filename = outname + "_variantSites.bed.gz";
	gz::ogzstream variantSitesFile(filename.c_str());


	//prepare mutation table
	float** mutTable;
	mutTable = new float*[4];
	for(int i=0; i<4; ++i)
		mutTable[i] = new float[4];

	//simulate sequences
	int refId = 0;
	double oldTheta = -1.0;
	std::vector<double>::iterator thetaIt = theta.begin();
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId, ++thetaIt){
		logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

		//create mutation table
		if(*thetaIt != oldTheta){
			fillMutationTable(mutTable, *thetaIt);
			oldTheta = *thetaIt;
		}

		//update reference storage and update haplotype lengths
		referenceObj.setChr(chrIt->name, chrIt->length);
		haplotypes.setLength(chrIt->length);

		//simulate genotypes
		logfile->listFlush("Simulating genotypes ...");
		simulateDiploidHaplotypesCurChromosome(haplotypes.getHaplotypesFirstIndividual(), mutTable, referenceObj.getPointerToRef());
		haplotypes.writeGenotypes(genoFile, chrIt->name, toBase);
		writeBEDFiles(haplotypes.getHaplotypesFirstIndividual(), invariantSitesFile, variantSitesFile);
		logfile->write(" done!");

		//now simulate and write reads
		simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesFirstIndividual(), bamFile, "");

		//end of chromosome
		logfile->endIndent();
	}

	//close stuff
	bamFile.close();
	bamFileOpen = false;
	genoFile.close();

	//clear memory
	for(int i=0; i<4; ++i)
		delete[] mutTable[i];
	delete[] mutTable;
}

//-------------------------------------------------------------------------------------
//Functions to simulate a pair of individuals according to a specific genetic distance
//-------------------------------------------------------------------------------------
void TSimulatorGenotypeCombination::fillTables(std::vector<double> & phis, float* baseFreq){
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
	double* cumul = new double[24];

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

	//cases 1 to 3: aa/ab, ab/aa, aa/bb
	//-----------------------------------------
	//build normalized cumulative vector for these cases
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
				++index;
			}
		}
	}

	//cases 4: ab/ab
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	index = 0;
	sum = 0.0;
	genoTrans[4] = new short*[6];
	for(a=0; a<3; ++a){
		for(b=a+1; b<4; ++b){
			sum += baseFreq[a] * baseFreq[b];
			cumul[index] = sum;
			genoTrans[4][index] = new short[4];
			genoTrans[4][index][0] = a; genoTrans[4][index][1] = b; genoTrans[4][index][2] = a; genoTrans[4][index][3] = b;
			++index;
		}
	}
	//normalize
	for(index=0; index<6; ++index)
		cumul[index] /= sum;

	//now initialize
	cumulGenoCombinationFreq[4] = new double[6];
	numGenotypeCombinations[4] = 6;
	for(index=0; index<6; ++index){
		cumulGenoCombinationFreq[4][index] = cumul[index];

	}

	//case 5: ab/ac
	//-----------------------------------------
	//build normalized cumulative vector for these cases
	index = 0;
	sum = 0.0;
	genoTrans[5] = new short*[24];
	for(a=0; a<4; ++a){
		for(b=0; b<4; ++b){
			if(a!=b){
				for(c=0; c<4; ++c){
					if(c!=a && c!=b){
						sum += baseFreq[a] * baseFreq[b] * baseFreq[c];
						cumul[index] = sum;
						genoTrans[5][index] = new short[4];
						genoTrans[5][index][0] = a; genoTrans[5][index][1] = b; genoTrans[5][index][2] = a; genoTrans[5][index][3] = c;
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
	cumulGenoCombinationFreq[5] = new double[24];
	numGenotypeCombinations[5] = 24;
	for(index=0; index<24; ++index){
		cumulGenoCombinationFreq[5][index] = cumul[index];
	}

	//cases 6 and 7: aa/bc, ab/cc
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

	//clean up
	delete[] cumul;

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

void TSimulator::simulateIndividualPair(std::vector<double> & phis, std::string & outname){
	if(phis.size() != 9)
		throw "Wrong number of phi! Required are nine values for genotype combinations 00/00, 00/01, 01/00, 00/11, 01/01, 01/02, 00/12, 01/22, 01/23";

	//normalize phis
	double sum = 0.0;
	for(std::vector<double>::iterator it=phis.begin(); it!=phis.end(); ++it)
		sum += *it;
	if(sum != 1.0){
		logfile->list("Normalizing phi to sum to one (currently summing to " + toString(sum) + ").");
		for(std::vector<double>::iterator it=phis.begin(); it!=phis.end(); ++it)
			*it /= sum;
	}
	logfile->list("Used phi are: " + concatenateString(phis, ", "));

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
	gz::ogzstream genoFile(filename.c_str());

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
		//TODO: add functionality for haploid chromosomes!
		logfile->listFlush("Simulating genotypes ...");
		genoComb.simulateHaplotypes(haplotypes.getHaplotypesOfIndividual(0), haplotypes.getHaplotypesOfIndividual(1), referenceObj.getPointerToRef(), referenceDivergence, chrIt->length, randomGenerator);
		haplotypes.writeGenotypes(genoFile, chrIt->name, toBase);
		logfile->write(" done!");

		//simulating reads
		simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesOfIndividual(0), bamFile0, " for individual 1");
		simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesOfIndividual(1), bamFile1, " for individual 2");

		logfile->endIndent();
	}

	//close stuff
	bamFile0.close();
	bamFile1.close();
	bamFileOpen = false;
	genoFile.close();
}

//--------------------------------------------------------
//Functions to simulate multiple individuals bases on SFS
//--------------------------------------------------------
void TSimulator::fillMutationTable(float** & mutTable){
	//table to pick derived allele
	double sum;
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			mutTable[i][j] = baseFreq[i] * baseFreq[j];
		}
		mutTable[i][i] = 0.0;

		//normalize within row
		sum = 0.0;
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
}

static inline int is_odd(int x){ return x % 2 != 0; }

void TSimulator::simulateHaplotypes(TSimulatorHaplotypes & haplotypes, SFS* sfs, float** & mutTable, short* ref){
	//prepare cumulative probability distribution for reference
	float cumulRef[4];
	cumulRef[0] = 1.0 - referenceDivergence;
	cumulRef[1] = cumulRef[0] + referenceDivergence / 3.0;
	cumulRef[2] = cumulRef[1] + referenceDivergence / 3.0;
	cumulRef[3] = 1.0;

	//variables
	short ancestral, derived;
	int alleleCount;
	int i, j;
	int numInd = haplotypes.size();
	int numNeeded;

	//now simulate haplotypes
	//TODO: add functionality for haploid chromosomes!
	if(chrIt->haploid){
		for(int l=0; l<chrIt->length; ++l){
			//pick alleles
			ancestral = randomGenerator->pickOne(4, cumulBaseFreq);
			derived = randomGenerator->pickOne(4, mutTable[ancestral]);

			//pick derived allele frequency
			alleleCount = sfs->getRandomAlleleCount(randomGenerator);

			//pick haplotypes that are derived
			numNeeded = alleleCount;
			for(i=0; i<numInd; ++i){
				if(randomGenerator->getRand() < (double) numNeeded / (double) (numInd - i)){
					haplotypes(i,0,l) = derived;
					--numNeeded;
					if(numNeeded == 0){
						for(j=i+1; j<numInd; ++j){
							haplotypes(i,0,l) = ancestral;
							haplotypes(i,1,l) = ancestral;
						}
						break;
					}
				} else
					haplotypes(i,0,l) = ancestral;

				//make homozygous
				haplotypes(i,1,l) = haplotypes(i,0,l);
			}

			//decide on reference sequence
			if(alleleCount > 0){
				if(randomGenerator->getRand() < (double) alleleCount / (double) numInd)
					ref[l] = derived;
				else
					ref[l] = ancestral;
			} else
				ref[l] = (ancestral + randomGenerator->pickOne(4, cumulRef)) % 4;
		}
	} else {
		int numHaplotypes = 2 * numInd;

		std::string f = "alleleCounts.txt";
		std::ofstream oo(f.c_str());

		for(int l=0; l<chrIt->length; ++l){
			//pick alleles
			ancestral = randomGenerator->pickOne(4, cumulBaseFreq);
			derived = randomGenerator->pickOne(4, mutTable[ancestral]);

			//pick derived allele frequency
			alleleCount = sfs->getRandomAlleleCount(randomGenerator);
			oo << alleleCount << "\n";

			//pick haplotypes that are derived
			if(alleleCount == 0){
				for(i=0; i<numInd; ++i){
					haplotypes(i,0,l) = ancestral;
					haplotypes(i,1,l) = ancestral;
				}

				//decide on reference sequence
				ref[l] = (ancestral + randomGenerator->pickOne(4, cumulRef)) % 4;
			} else {
				numNeeded = alleleCount;
				for(i=0; i<numHaplotypes; ++i){
					double prob = (double) numNeeded / (double) (numHaplotypes - i);
					if(randomGenerator->getRand() < prob){
						haplotypes(i / 2, is_odd(i), l) = derived;
						--numNeeded;
						if(numNeeded == 0){
							for(j=i+1; j<numHaplotypes; ++j)
								haplotypes(j / 2, is_odd(j), l) = ancestral;
							break;
						}
					} else
						haplotypes(i / 2, is_odd(i), l) = ancestral;
				}

				//decide on reference sequence
				if(randomGenerator->getRand() < (double) alleleCount / (double) numHaplotypes)
					ref[l] = derived;
				else
					ref[l] = ancestral;
			}
		}

		oo.close();
	}
}

void TSimulator::simulatePopulationFromSFS(double theta, int numIndividuals, std::string & outname){
	std::vector<double> thetas;
	for(unsigned int i=0; i<chromosomes.size(); ++i)
		thetas.push_back(theta);

	simulatePopulationFromSFS(thetas, numIndividuals, outname);
}

void TSimulator::simulatePopulationFromSFS(std::vector<double> & thetas, int numIndividuals, std::string & outname){
	if(thetas.size() != chromosomes.size())
		throw "Number of theta values does not match number of chromosomes!";

	//generate SFS for each chromosome
	std::vector<SFS*> sfs;
	logfile->listFlush("Generating SFS from provided theta values ...");
	int chr = 1;
	std::string filename;
	std::vector<double>::iterator thetaIt = thetas.begin();
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++thetaIt, ++chr){
		sfs.push_back(new SFS((2 - chrIt->haploid) * numIndividuals, (float) *thetaIt));

		//save true SFS
		filename = outname + "_trueSFS_chr" + toString(chr) + ".txt";
		(*sfs.rbegin())->writeToFile(filename);
	}
	logfile->done();
	logfile->conclude("True SFS written to '" + outname + "_trueSFS_chr*.txt'.");

	//Now simulate
	simulatePopulationFromSFS(sfs, numIndividuals, outname);

	//deleting SFS
	for(std::vector<SFS*>::iterator it=sfs.begin(); it!=sfs.end(); ++it)
		delete *it;
}

void TSimulator::simulatePopulationFromSFS(std::vector<std::string> & sfsFileNames, bool folded, int numIndividuals, std::string & outname){
	if(sfsFileNames.size() != chromosomes.size())
		throw "Number of SFS files does not match number of chromosomes!";

	//read the SFS of each chromosome from the corresponding file
	std::vector<SFS*> sfs;
	std::vector<std::string>::iterator it = sfsFileNames.begin();
	int nChr;
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++it){
		logfile->listFlush("Reading the sfs of chromosome '" + chrIt->name + "' from file '" + *it + "' ...");
		if(folded) sfs.push_back(new SFSfolded(*it));
		else sfs.push_back(new SFS(*it));
		logfile->done();

		nChr = (2-chrIt->haploid) * numIndividuals;
		if((*sfs.rbegin())->numChromosomes != nChr)
			throw "SFS does not match sample size! It contains data for " + toString((*sfs.rbegin())->numChromosomes) + " instead of " + toString(nChr) + " chromosomes.";
	}

	//Now simulate
	simulatePopulationFromSFS(sfs, numIndividuals, outname);

	//deleting SFS
	for(std::vector<SFS*>::iterator it=sfs.begin(); it!=sfs.end(); ++it)
		delete *it;
}

void TSimulator::simulatePopulationFromSFS(std::vector<SFS*> sfs, int numIndividuals, std::string & outname){
	//one SFS per chromosome! Each SFS needs to have the same number of
	if(sfs.size() != chromosomes.size())
			throw "Number of theta values provided does not match number of chromosomes to simulate!";

	//open bam files
	logfile->startIndent("Opening BAM files:");
	TSimulatorBamFile** bamFiles = new TSimulatorBamFile*[numIndividuals];
	for(int i=0; i<numIndividuals; ++i)
		bamFiles[i] = new TSimulatorBamFile(outname + "_ind" + toString(i+1) + ".bam", readGroupName, chromosomes, logfile);
	bamFileOpen = true;
	logfile->endIndent();

	//open FASTA file for reference sequences
	std::string filename = outname + ".fasta";
	TSimulatorReference referenceObj(filename, toBase, logfile);

	//prepare haplotypes and
	TSimulatorHaplotypes haplotypes(numIndividuals);

	//open file for true genotypes
	filename = outname + "_trueGenotypes.txt";
	gz::ogzstream genoFile(filename.c_str());

	//prepare mutation table
	float** mutTable;
	mutTable = new float*[4];
	for(int i=0; i<4; ++i)
		mutTable[i] = new float[4];
	fillMutationTable(mutTable);

	//simulate sequences
	int refId = 0;
	std::vector<SFS*>::iterator sfsIt = sfs.begin();
	for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId, ++sfsIt){
		logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

		//update reference storage and update haplotype lengths
		referenceObj.setChr(chrIt->name, chrIt->length);
		haplotypes.setLength(chrIt->length);

		//simulate genotypes
		logfile->listFlush("Simulating genotypes ...");
		simulateHaplotypes(haplotypes, *sfsIt, mutTable, referenceObj.getPointerToRef());
		haplotypes.writeGenotypes(genoFile, chrIt->name, toBase);
		logfile->write(" done!");

		//now simulate and write reads
		logfile->startIndent("Simulating reads:");
		for(int i=0; i<numIndividuals; ++i)
			simulateReadsFromHaplotypes(chrIt, haplotypes.getHaplotypesOfIndividual(i), *bamFiles[i], " for individual " + toString(i+1));
		logfile->endIndent();

		//end of chromosome
		logfile->endIndent();
	}

	//close stuff
	logfile->startIndent("Indexing BAM files:");
	for(int i=0; i<numIndividuals; ++i){
		delete bamFiles[i];
	}
	delete[] bamFiles;
	bamFileOpen = false;
	logfile->endIndent();
	genoFile.close();

	//clear memory
	for(int i=0; i<4; ++i)
		delete[] mutTable[i];
	delete[] mutTable;
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
