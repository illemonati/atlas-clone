/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */


#include "TGenome.h"

//-------------------------------------------------------
//TGenome
//-------------------------------------------------------
TGenome::TGenome(TLog* Logfile, TParameters & params){
	logfile = Logfile;
	randomGeneratorInitialized = false;

	//read parameters
	filename = params.getParameterString("bam");
	windowSize = params.getParameterDouble("window");
	numWindowsOnChr = 0;
	//if(windowSize < 1000) throw "Window size should be at least 1Kb!";
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 0.95);

	//outputname
	outputName = params.getParameterStringWithDefault("out", "");
	if(outputName == ""){
		//guess from filename
		outputName = filename;
		outputName = extractBeforeLast(outputName, ".");
	}

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;
	curStart = -1;
	curEnd = -1;
	oldPos = -1;

	//open BAM file
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";

	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + filename + "'!";

	//read header
	bamHeader = bamReader.GetHeader();
	readGroups.fill(bamHeader);
	chrIterator = bamHeader.Sequences.End();

	//initialize post mortem damage
	initializePostMortemDamage(params);

	//initialize recalibration
	initializeRecalibration(params);

	//check if we mask sites
	if(params.parameterExists("mask")){
		doMasking = true;
		std::string maskFile = params.getParameterString("mask");
		logfile->startIndent("Will mask all sites listed in BED file '" + maskFile + "':");
		logfile->listFlush("Reading file ...");
		mask = new TBedReader(maskFile, windowSize);
		logfile->write(" done!");
		logfile->endIndent();
		//mask->print();
	} else doMasking = false;

	//for debugging: work on one window only
	limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	limitChr = params.getParameterIntWithDefault("limitChr", 1000000);
};

void TGenome::jumpToEnd(){
	chrIterator = bamHeader.Sequences.End();
}

void TGenome::restartChromosome(TWindowPair & windowPair){
	chrIterator = bamHeader.Sequences.Begin();
	chrNumber = 0;

	moveChromosome(windowPair);
}

bool TGenome::iterateChromosome(TWindowPair & windowPair){
	if(chrIterator == bamHeader.Sequences.End()){
		chrIterator = bamHeader.Sequences.Begin();
		chrNumber = 0;
	} else {
		logfile->endNumbering();
		//move to next
		++chrIterator;
		++chrNumber;
	}

	//did we reach end?
	if(chrIterator == bamHeader.Sequences.End() || chrNumber >= limitChr){
		curEnd = 0;
		chrIterator = bamHeader.Sequences.End();
		return false;
	}

	moveChromosome(windowPair);
	return true;
}

void TGenome::moveChromosome(TWindowPair & windowPair){
	//jump reader
	bamReader.Jump(chrNumber, 0);

	//restart windows
	chrLength = stringToLong(chrIterator->Length);
	numWindowsOnChr = ceil(chrLength / (double) windowSize);

	curStart = 0;
	curEnd = 0;
	oldPos = -1;
	windowNumber = 0;
	int nextEnd = windowSize;
	if(nextEnd > chrLength) nextEnd = chrLength;
	windowPair.nextPointer->move(0, nextEnd);

	//advance mask
	if(doMasking) mask->setChr(chrIterator->Name);

	//write progress
	logfile->startNumbering("Parsing chromosome '" + chrIterator->Name + "':");
}

bool TGenome::iterateWindow(TWindowPair & windowPair){
	if(curEnd > 0) logfile->endIndent();

	//swap window pairs
	windowPair.swap();

	//move to next region
	curStart = curEnd;
	curEnd = windowPair.curPointer->end;
	if(curStart >= chrLength || windowNumber>= limitWindows) return false;

	//move next
	long nextEnd = curEnd + windowSize;
	if(nextEnd > chrLength) nextEnd = chrLength + 1;
	windowPair.nextPointer->move(curEnd, nextEnd);

	++windowNumber;

	//report
	logfile->number("Window [" + toString(curStart) + ", " + toString(curEnd) + ") of " + toString(numWindowsOnChr) + " on '" + chrIterator->Name + "':");
	logfile->addIndent();

	return true;
};

bool TGenome::readData(TWindowPair & windowPair){
	logfile->listFlush("Reading data ...");

	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//parse through reads
	int numReads = 0;
	while(bamReader.GetNextAlignment(bamAlignement) && bamAlignement.RefID==chrNumber){
		//std::cout << "REF ID = " << bamAlignement.RefID << "\tpos = " << bamAlignement.Position << std::endl;
		//only take those reads that pass QC
		if(!bamAlignement.IsFailedQC() && !bamAlignement.IsDuplicate() && bamAlignement.Position >= curStart){
			//check if bam file is sorted
			if(bamAlignement.Position < oldPos){
				throw "BAM file must be sorted by position!";
			}
			oldPos = bamAlignement.Position;

			//check if still within current window and add to window
			if(bamAlignement.Position >= curEnd){
				//nextWindow->addFromRead(bamAlignement);
				break;
			} else {
				++numReads;
				if(windowPair.curPointer->addFromRead(bamAlignement, pmdObjects, &readGroups)){
					//add also to next window in case reads overhangs current window -> function returns true
					windowPair.nextPointer->addFromRead(bamAlignement, pmdObjects, &readGroups);
				}
			}
		}
	}

	//apply mask
	if(doMasking) windowPair.curPointer->applyMask(mask);

	//calc coverage
	windowPair.curPointer->calcCoverage();

	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");
	logfile->conclude("read data from " + toString(numReads) + " reads.");
	logfile->conclude("coverage is " + toString(windowPair.curPointer->coverage));
	logfile->conclude(toString(windowPair.curPointer->fractionsitesCoverageAtLeastTwo * 100) + "% of all sites are covered at least twice");
	logfile->conclude(toString(windowPair.curPointer->fractionSitesNoData * 100) + "% of all sites have no data");

	return true;
};

void TGenome::initializePostMortemDamage(TParameters & params){
	logfile->startIndent("Initializing Post Mortem Damage (PMD):");
	//create an array of TPMD objects for each read group
	pmdObjects = new TPMD[readGroups.numGroups];

	//now fill them!
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		logfile->list("Initializing one PMD function for all read groups.");
		//all read groups have the same pmd
		if(params.parameterExists("pmd")){
			std::string pmdString = params.getParameterString("pmd");
			logfile->list("Initializing PMD for both C->T and G->A with function '" + pmdString +"'.");
			for(int i=0; i<readGroups.numGroups; ++i){
				pmdObjects[i].initializeFunction(pmdString, pmdGA);
				pmdObjects[i].initializeFunction(pmdString, pmdCT);
			}
			logfile->conclude(pmdObjects[0].getFunctionString(pmdCT));
		} else {
			if(!params.parameterExists("pmdCT")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdCT' has to be provided!";
			std::string pmdStringCT = params.getParameterString("pmdCT");
			logfile->list("Initializing post mortem C->T damage with function '" + pmdStringCT +"'.");
			for(int i=0; i<readGroups.numGroups; ++i){
				pmdObjects[i].initializeFunction(pmdStringCT, pmdCT);
			}
			logfile->conclude(pmdObjects[0].getFunctionString(pmdCT));

			//second G->A
			if(!params.parameterExists("pmdGA")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdGA' has to be provided!";
			std::string pmdStringGA = params.getParameterString("pmdGA");
			logfile->list("Initializing post mortem G->A damage with function '" + pmdStringGA +"'.");
			for(int i=0; i<readGroups.numGroups; ++i){
				pmdObjects[i].initializeFunction(pmdStringGA, pmdGA);
			}
			logfile->conclude(pmdObjects[0].getFunctionString(pmdGA));
		}
	} else if(params.parameterExists("pmdFile")){
		//read from file for each read group
		std::string filename = params.getParameterString("pmdFile");
		logfile->list("Reading PMD from  from file '" + filename + "'.");
		std::ifstream file(filename.c_str());
		if(!file) throw "Failed to open PMD file '" + filename + "'!";

		//parse file that has structure: readGroup PMD(CT) PMD(GA)
		int lineNum = 0;
		std::string line;
		std::vector<std::string> vec;
		int readGroupId;
		while(file.good() && !file.eof()){
			++lineNum;
			//skip empty lines or those that start with //
			std::getline(file, line);
			line = extractBefore(line, "//");
			trimString(line);
			if(!line.empty()){
				fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
				if(vec.size() != 3) throw "Found " + toString(vec.size()) + " instead of 3 columns in '" + filename + "' on line " + toString(lineNum) + "!";
				//get read group
				if(readGroups.readGroupExists(vec[0])){ //ignore if it does not exist
					readGroupId = readGroups.find(vec[0]);
					//initialize functions
					pmdObjects[readGroupId].initializeFunction(vec[1], pmdCT);
					logfile->conclude("For read group '" + vec[0] + "', C->T: " + pmdObjects[readGroupId].getFunctionString(pmdCT));
					pmdObjects[readGroupId].initializeFunction(vec[2], pmdGA);
					logfile->conclude("For read group '" + vec[0] + "', G->A: " + pmdObjects[readGroupId].getFunctionString(pmdGA));
				}
			}
		}

		//close file
		file.close();

		//test if we have a function for all read groups
		for(int i=0; i<readGroups.numGroups; ++i){
			if(!pmdObjects[i].functionInitialized(pmdCT)) throw "PMD C->T for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
			if(!pmdObjects[i].functionInitialized(pmdGA)) throw "PMD G->A for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
		}
	} else {
		//no post mortem damage
		logfile->list("Assuming there is no PMD in the data.");
		std::string pmdString = "none";
		for(int i=0; i<readGroups.numGroups; ++i){
			pmdObjects[i].initializeFunction(pmdString, pmdGA);
			pmdObjects[i].initializeFunction(pmdString, pmdCT);
		}
	}
	logfile->endIndent();
}

void TGenome::initializeRecalibration(TParameters & params){
	if(params.parameterExists("BQSRQuality")){
		recalObject = new TRecalibrationBQSR(&bamHeader, params, logfile);
	} else {
		recalObject = new TRecalibration();
	}
}

void TGenome::openThetaOutputFile(std::ofstream & out){
	std::string filename = outputName + "_theta_estimates.txt";
	logfile->list("Writing theta estimates to '" + filename + "'");
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	out << std::setprecision(9) << "Chr\t";
	out << "start\tend\tcoverage\tmissing\ttwoOrMore\tpi(A)\tpi(C)\tpi(G)\tpi(T)\ttheta_MLE\ttheta_C95_l\ttheta_C95_u\tLL";
	out << "\n";
}

void TGenome::initializeRandomGenerator(TParameters & params){
	logfile->listFlush("Initializing random generator ...");

	if(params.parameterExists("fixedSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator=new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
	randomGeneratorInitialized = true;
}

void TGenome::estimateTheta(TParameters & params){
	//EM params
	EMParameters EMParams(params, logfile);

	//open output
	std::ofstream out; openThetaOutputFile(out);

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//check if we have data -> can be extended to ensure
			if(windows.cur->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				//estimate Theta
				out << chrIterator->Name << "\t";
				windows.cur->estimateTheta(EMParams, recalObject, out, logfile);
			}
		}
	}

	//clean up
	out.close();
}

void TGenome::calcLikelihoodSurfaces(TParameters & params){
	//read params
	int steps = params.getParameterIntWithDefault("steps", 100);
	int numWindows = params.getParameterIntWithDefault("numWindows", 1);

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	int windowsCalculated = 0;
	std::string filename;
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//check if we have data -> can be extended to ensure
			if(windows.cur->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				//open file
				std::ofstream out;
				filename = outputName + chrIterator->Name + "_" + toString(windows.cur->start) + "_LLsurface.txt";
				out.open(filename.c_str());
				if(!out) throw "Failed to open output file '" + outputName + "'!";

				//calc surface
				logfile->listFlush("Calculating likelihood surface ...");
				windows.cur->calcLikelihoodSurface(recalObject, out, steps);
				logfile->write(" done!");

				//close output
				out.close();

				//check if we break
				++windowsCalculated;
				if(windowsCalculated >= numWindows) break;
			}
		}
	}
}

void TGenome::callMLEGenotypes(TParameters & params){
	//initialize random generator
	initializeRandomGenerator(params);

	//do we print sites with no data?
	bool printIfNoData = params.parameterExists("printAll");
	if(printIfNoData) logfile->list("Will print all sites, even those without data");

	//open FASTA reference
	BamTools::Fasta reference;
	bool printRefBase = false;
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		logfile->list("Adding reference base from '" + fastaFile + "'.");
		std::string fastaIndex = fastaFile + ".fai";
		reference.Open(fastaFile, fastaIndex);
		printRefBase = true;
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream out;
	if(params.parameterExists("vcf")){
		if(!printRefBase) throw "Can not print VCF file without reference!";
		writeVCF = true;

		//open file
		filename = outputName + "_MLEGenotypes.vcf.gz";
		logfile->list("Writing estimates of allele presence in VCF format to '" + filename + "'");
		out.open(filename.c_str());
		if(!out) throw "Failed to open output file '" + filename + "'!";

		//write header
		out << "##fileformat=VCFv4.2\n";
		out << "##source=estimHet\n";
		out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		out << "##INFO=<ID=GG,Number=10,Type=Integer,Description=\"Phred-scaled genotype likelihoods for all genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT and TT\">\n";
		out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		out << "##FORMAT=<ID=PL,Number=1,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n";
		out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		filename = outputName + "_MLEGenotypes.txt.gz";
		logfile->list("Writing estimates of allele presence to '" + filename + "'");
		out.open(filename.c_str());
		if(!out) throw "Failed to open output file '" + filename + "'!";

		//write header
		out << "chr\tpos";
		if(printRefBase) out << "\tRef";
		out << "\tcoverage\tP(A|D)\tP(C|D)\tP(G|D)\tP(T|D)\tMAP\tQ\n";
	}

	//write header
	out << std::setprecision(3);
	out << "chr\tpos";
	if(printRefBase) out << "\tRef";
	out << "\tcoverage\tL(AA)\tL(AC)\tL(AG)\tL(AT)\tL(CC)\tL(CG)\tL(CT)\tL(GG)\tL(GT)\tL(TT)\tMLE\tQ\n";

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//add reference data
			if(printRefBase) windows.cur->addReferenceBaseToSites(reference, chrNumber);

			//call genotypes
			logfile->listFlush("Calling MLE genotypes ...");
			windows.cur->callMLEGenotype(recalObject, *randomGenerator, out, chrIterator->Name, printIfNoData, printRefBase, writeVCF);
			logfile->write(" done!");
		}
	}

	//clean up
	out.close();
}


void TGenome::callBayesianGenotypes(TParameters & params){
	//initialize random generator
	initializeRandomGenerator(params);

	//do we estimate theta or is it given?
	double theta;
	bool estimateTheta;
	EMParameters* EMParams;
	std::ofstream out;
	if(params.parameterExists("theta")){
		estimateTheta = false;
		theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
	} else {
		estimateTheta = true;
		//read EM params
		EMParams = new EMParameters(params, logfile);
		openThetaOutputFile(out);
	}

	//do we print sites with no data?


	//open FASTA reference
	BamTools::Fasta reference;
	bool printRefBase = false;
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		logfile->list("Adding reference base from '" + fastaFile + "'.");
		std::string fastaIndex = fastaFile + ".fai";
		reference.Open(fastaFile, fastaIndex);
		printRefBase = true;
	}

	//limit to a set of sites?
	bool limitToSitesWithKnownAlleles;
	bool printIfNoData;
	TSiteSubset* subset;
	if(params.parameterExists("sites")){
		subset = new TSiteSubset(params.getParameterString("sites"), windowSize);
		limitToSitesWithKnownAlleles = true;
		printIfNoData = true;
	} else {
		printIfNoData = params.parameterExists("printAll");
		if(printIfNoData) logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream output;
	if(params.parameterExists("vcf")){
		if(!printRefBase) throw "Can not print VCF file without reference!";
		writeVCF = true;

		//open file
		filename = outputName + "_BayesianGenotypes.vcf.gz";
		logfile->list("Writing Bayesian genotypes in VCF format to '" + filename + "'");
		output.open(filename.c_str());
		if(!output) throw "Failed to open output file '" + filename + "'!";

		//write header
		output << "##fileformat=VCFv4.2\n";
		output << "##source=estimHet\n";
		output << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		output << "##INFO=<ID=PP,Number=10,Type=Integer,Description=\"Phred-scaled posterior probabilities of all genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT and TT\">\n";
		output << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		output << "##FORMAT=<ID=GP,Number=1,Type=Integer,Description=\"Genotype posterior probabilities (phred-scaled)\">\n";
		output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		filename = outputName + "_BayesianGenotypes.txt.gz";
		logfile->list("Writing Bayesian genotypes to '" + filename + "'");
		output.open(filename.c_str());
		if(!output) throw "Failed to open output file '" + filename + "'!";

		//write header
		output << "chr\tpos";
		if(printRefBase) output << "\tRef";
		output << "\tcoverage\tP(AA|D)\tP(AC|D)\tP(AG|D)\tP(AT|D)\tP(CC|D)\tP(CG|D)\tP(CT|D)\tP(GG|D)\tP(GT|D)\tP(TT|D)\tMAP\tQ\n";
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//check if we have data -> can be extended to ensure
			if(windows.cur->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				out << chrIterator->Name << "\t";

				//set Theta
				if(estimateTheta){
					windows.cur->estimateTheta((*EMParams), recalObject, out, logfile);
				} else {
					windows.cur->calculateEmissionProbabilities(recalObject);
					windows.cur->estimateBaseFrequencies();
					windows.cur->setTheta(theta);
				}

				//call allele presence
				logfile->listFlush("Calling Bayesian Genotypes ...");
				if(limitToSitesWithKnownAlleles){
					windows.cur->callBayesianGenotypeKnownAlleles(subset, *randomGenerator, output, chrIterator->Name, printIfNoData, writeVCF);
				} else {
					if(printRefBase) windows.cur->addReferenceBaseToSites(reference, chrNumber);
					windows.cur->callBayesianGenotype(*randomGenerator, output, chrIterator->Name, printIfNoData, printRefBase, writeVCF);
				}
				logfile->write(" done!");
			}
		}
	}

	//clean up
	out.close();
	if(estimateTheta) delete EMParams;
}

void TGenome::callAllelePresence(TParameters & params){
	//initialize random generator
	initializeRandomGenerator(params);

	//do we estimate theta or is it given?
	double theta;
	bool estimateTheta;
	EMParameters* EMParams;
	if(params.parameterExists("theta")){
		estimateTheta = false;
		theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
	} else {
		estimateTheta = true;
		//read EM params
		EMParams = new EMParameters(params, logfile);
	}

	//do we print sites with no data?
	bool printIfNoData = params.parameterExists("printAll");
	if(printIfNoData) logfile->list("Will print all sites, even those without data");

	//open FASTA reference
	BamTools::Fasta reference;
	bool printRefBase = false;
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		logfile->list("Adding reference base from '" + fastaFile + "'.");
		std::string fastaIndex = fastaFile + ".fai";
		reference.Open(fastaFile, fastaIndex);
		printRefBase = true;
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream outAllelePresence;
	if(params.parameterExists("vcf")){
		if(!printRefBase) throw "Can not print VCF file without reference!";
		writeVCF = true;

		//open file
		filename = outputName + "_AllelePresence.vcf.gz";
		logfile->list("Writing estimates of allele presence in VCF format to '" + filename + "'");
		outAllelePresence.open(filename.c_str());
		if(!outAllelePresence) throw "Failed to open output file '" + filename + "'!";

		//write header
		outAllelePresence << "##fileformat=VCFv4.2\n";
		outAllelePresence << "##source=estimHet\n";
		outAllelePresence << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		outAllelePresence << "##INFO=<ID=PP,Number=10,Type=Integer,Description=\"Phred-scaled posterior probabilities of allele presence in the order A, C, G and T\">\n";
		outAllelePresence << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		outAllelePresence << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		filename = outputName + "_AllelePresence.txt.gz";
		logfile->list("Writing estimates of allele presence to '" + filename + "'");
		outAllelePresence.open(filename.c_str());
		if(!outAllelePresence) throw "Failed to open output file '" + filename + "'!";

		//write header
		outAllelePresence << "chr\tpos";
		if(printRefBase) outAllelePresence << "\tRef";
		outAllelePresence << "\tcoverage\tP(A|D)\tP(C|D)\tP(G|D)\tP(T|D)\tMAP\tQ\n";
	}

	//open output for theta
	std::ofstream out; openThetaOutputFile(out);

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//check if we have data -> can be extended to ensure
			if(windows.cur->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				out << chrIterator->Name << "\t";

				//set Theta
				if(estimateTheta){
					windows.cur->estimateTheta((*EMParams), recalObject, out, logfile);
				} else {
					windows.cur->calculateEmissionProbabilities(recalObject);
					windows.cur->estimateBaseFrequencies();
					windows.cur->setTheta(theta);
				}

				//add reference data
				if(printRefBase) windows.cur->addReferenceBaseToSites(reference, chrNumber);

				//call allele presence
				logfile->listFlush("Calling allele presence ...");
				windows.cur->callAllelePresence(*randomGenerator, outAllelePresence, chrIterator->Name, printIfNoData, printRefBase, writeVCF);
				logfile->write(" done!");
			}
		}
	}

	//clean up
	out.close();
	if(estimateTheta) delete EMParams;
}

void TGenome::printPileup(){
	//prepare windows
	TWindowPairDiploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_pileup.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	TGenotypeMap genoMap;
	out << "Chr\tposition\tcoverage\tbases";
	for(int i=0; i<10; ++i)
		out << "\tPem(" << genoMap.getGenotypeString(i) << ")";
	out << "\n";

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//print pileup
			windows.cur->printPileup(recalObject, out, chrIterator->Name);
		}
	}

	//clean up
	out.close();
}


void TGenome::estimateErrorCalibrationEM(TParameters & params){
	//Initialize calibration parameters.
	//We assume a linear model log(error) = eta = beta_0 + beta_1 * quality + beta_2 * position in reads
	//                                            + gamma_A * Ind(d = A) + gamma_C * Ind(d = C) + gamma_G * Ind(d = G) + gamma_T * Ind(d = T)

	/*

	//read EM parameters
	int numEMIterations = params.getParameterIntWithDefault("iterations", 10);
	logfile->list("Will perform at max " + toString(numEMIterations) + " EM iterations.");
	double maxEpsilon = params.getParameterDoubleWithDefault("maxEps", 0.000001);
	logfile->list("Will stop EM when deltaLL < " + toString(maxEpsilon));
	double NewtonRalphsonNumIterations = params.getParameterIntWithDefault("NRiterations", 10);
	logfile->list("Will conduct at max " + toString(NewtonRalphsonNumIterations) + " Newton-Ralphson iterations");
	double NewtonRalphsonMaxF = params.getParameterDoubleWithDefault("maxF", 0.00001);
	logfile->list("Will stop Newton-Ralphson when F < " + toString(NewtonRalphsonMaxF));

	//create recalibration object
	TRecalibrationEM recalObject(&params, logfile);

	//prepare windows
	TWindowPairHaploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_calibration.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";
	out << "iteration";
	recalObject.writeHeader(out);

	//run EM iterations
	for(int i=0; i < numEMIterations; ++i){
		logfile->startIndent("EM Iteration " + toString(i+1) + ":");

		//prepare recal object for next iteration
		recalObject.initEMStep();

		//run Newton-Ralphosn iteration
		for(int j=0; j<NewtonRalphsonNumIterations; ++j){
			logfile->startIndent("Calculating Jacobian for Newton-Ralphson iteration " + toString(j+1) + ":");
			recalObject.initNetwonRalphsonStep();

			while(iterateChromosome(windows)){
				while(iterateWindow(windows)){
					//read data for current window
					readData(windows);
					windows.cur->estimateBaseFrequencies();
					windows.cur->addToJacobian(&recalObject);
				}
			}

			//clean up memory
			windows.clear();

			//perform parameter update
			logfile->endIndent();
			logfile->listFlush("Updating parameters ...");
			recalObject.runNewtonRalphson();
			logfile->write(" done!");

			//check if we break
			if(recalObject.maxF < NewtonRalphsonMaxF){
				logfile->list("Stopping Newton-Ralphson since maxF < " + toString(NewtonRalphsonMaxF));
				break;
			}
		}

		//save current parameters
		recalObject.saveParams();

		//Report to file
		recalObject.writeParams(out);
		logfile->endIndent();

		//check if we break EM

	}
	logfile->endIndent();

	//clean up
	out.close();

	*/
}

void TGenome::fillSequence(std::vector<double> & vec, std::string & str){
	//it is either a number, or a sequence min-max:num steps
	std::string::size_type posDash = str.find_first_of('-', 1);
	if(posDash != std::string::npos){
		std::string::size_type posColon = str.find_first_of(':');
		if(posColon != std::string::npos){
			//fill sequence
			double min = stringToDoubleCheck(str.substr(0, posDash));
			double max = stringToDoubleCheck(str.substr(posDash+1, posColon-posDash-1));
			int numSteps = stringToIntCheck(str.substr(posColon+1));
			double step = (max - min) / (double) (numSteps - 1.0);
			for(int i=0; i<numSteps; ++i) vec.push_back(min + (double) i * step);
		} else throw "Unable to understand sequence '" + str + "'!";
	} else {
		//should be a number
		vec.push_back(stringToDoubleCheck(str));
	}
}

void TGenome::calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params){
	/*
	//read vectors of betas to test
	logfile->startIndent("Will calculate likelihood on a grid with these marginal ranges:");
	logfile->startIndent("Will calculate LL for these parameter values:");

	//beta0
	std::vector<double> beta0;
	std::string tmp = params.getParameterString("beta0");
	fillSequence(beta0, tmp);
	logfile->list("beta0 = " + concatenateString(beta0, ", "));

	//beta1
	std::vector<double> beta1;
	tmp = params.getParameterString("beta1");
	fillSequence(beta1, tmp);
	logfile->list("beta1 = " + concatenateString(beta1, ", "));

	//beta2
	std::vector<double> beta2;
	tmp = params.getParameterString("beta2");
	fillSequence(beta2, tmp);
	logfile->list("beta2 = " + concatenateString(beta2, ", "));
	logfile->endIndent();

	//create recalibration object
	TRecalibrationEM recalObject(logfile);

	//prepare windows
	TWindowPairHaploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_calibration_LL_surface.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";
	recalObject.writeHeader(out);
	out << "\tLL\n";

	//run loop over all combinations
	double* recalParams = new double[3];
	for(std::vector<double>::iterator itBeta0=beta0.begin(); itBeta0!=beta0.end(); ++itBeta0){
		recalParams[0] = *itBeta0;
		for(std::vector<double>::iterator itBeta1=beta1.begin(); itBeta1!=beta1.end(); ++itBeta1){
			recalParams[1] = *itBeta1;
			for(std::vector<double>::iterator itBeta2=beta2.begin(); itBeta2!=beta2.end(); ++itBeta2){
				recalParams[2] = *itBeta2;
				logfile->startIndent("Calculating LL for beta0 = " + toString(*itBeta0) + ", beta1 = " + toString(*itBeta1) + ", beta2 = " + toString(*itBeta2) + ":");

				//reset
				recalObject.resetLikelihood();

				//set parameters
				recalObject.setParams(recalParams);

				//loop over all windows
				while(iterateChromosome(windows)){
					while(iterateWindow(windows)){
						//read data for current window
						readData(windows);
						//calc LL
						windows.cur->estimateBaseFrequencies();
						windows.cur->addToLikelihoodRecalibration(&recalObject);
					}
				}

				//clean up memory
				windows.clear();

				//output
				logfile->conclude("LL = " + toString(recalObject.logLikelihood));
				recalObject.writeParams(out);
				out << "\t" << recalObject.logLikelihood << std::endl;
				logfile->endIndent();
			}
		}
	}

	//clean up
	out.close();
	*/
}

void TGenome::BQSR(TParameters & params){
	//prepare windows
	TWindowPairHaploid windows;

	//open FASTA reference
	std::string fastaFile = params.getParameterString("fasta");
	std::string fastaIndex = fastaFile + ".fai";
	BamTools::Fasta reference;
	reference.Open(fastaFile, fastaIndex);

	//create BQSR object
	TRecalibrationBQSR bqsr(&bamHeader, params, logfile);
	if(bqsr.allConverged()){
		logfile->list("No need to estimate any BQSR cells. Aborting Program.");
		return;
	}

	//tmp variables
	bool hasConverged = false;
	int loopNumber = 0;

	//check if we use first chromosome for initial convergence
	if(params.parameterExists("preConverge")){
		int numPreConvLoops = params.getParameterInt("preConverge");
		logfile->startIndent("Only using first chromosome to get initial estimates for " + toString(numPreConvLoops) + " loops :");

		//do we limit number of windows?
		bool limitWindows = params.parameterExists("limitWindows");
		int maxNumWindows = 0;
		if(limitWindows)
			maxNumWindows = params.getParameterInt("limitWindows");

		//run until it converges
		while(!hasConverged && loopNumber < numPreConvLoops){
			++loopNumber;
			logfile->startIndent("Running pre-recalibration loop " + toString(loopNumber) + ":");

			//jump to first chromosome
			if(loopNumber == 1)	iterateChromosome(windows);
			else restartChromosome(windows);

			//iterate over all windows
			while(iterateWindow(windows)){
				//read data for current window
				readData(windows);

				//add reference data
				windows.cur->addReferenceBaseToSites(reference, chrNumber);

				//add the base to BQSR
				windows.cur->addSitesToBQSR(bqsr, logfile);

				logfile->list("All done for this window!");

				//check if we break
				if(limitWindows && windowNumber == maxNumWindows) break;
			}
			logfile->endIndent();

			//clean up memory and restart from chr 1
			windows.clear();

			//estimate epsilon
			hasConverged = bqsr.estimateEpsilon();

			//write results to file
			bqsr.writeToFile(outputName + "_preConvergence_Loop" + toString(loopNumber));
			logfile->endIndent();
		}

		//reset counters and such
		hasConverged = false;
		loopNumber = 0;
		logfile->endIndent();
		jumpToEnd();
	}

	//loop over bam until BQSR converges
	while(!hasConverged){
		++loopNumber;
		logfile->startIndent("Running recalibration loop " + toString(loopNumber) + ":");
		//loop over all windows

		while(iterateChromosome(windows)){
			while(iterateWindow(windows)){
				//read data for current window
				readData(windows);

				//add reference data
				windows.cur->addReferenceBaseToSites(reference, chrNumber);

				//add the base to BQSR
				windows.cur->addSitesToBQSR(bqsr, logfile);

				logfile->list("All done for this window!!");
			}
		}


		//clean up memory
		windows.clear();

		//estimate epsilon
		hasConverged = bqsr.estimateEpsilon();

		//write results to file
		bqsr.writeToFile(outputName + "_Loop" + toString(loopNumber));

		logfile->endIndent();
	}

	//write results to file
	bqsr.writeToFile(outputName);
}


void TGenome::printQualityTransformation(TParameters & params){
	//prepare windows
	TWindowPairHaploid windows;

	//create table to store counts
	int maxQ = params.getParameterIntWithDefault("maxQ", 100);
	TQualityTransformTable QT(maxQ);

	//preapre output
	std::ofstream out;
	std::string filename;

	//loop over windows and add to table
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//add the base to BQSR
			windows.cur->addSitesToQualityTransformTable(recalObject, QT, logfile);
		}
		//print out after each chr
		filename = outputName + "_qualityTransformation_upToChr_" + toString(chrNumber) + ".txt";
		out.open(filename.c_str());
		if(!out) throw "Failed to open output file '" + outputName + "'!";
		QT.printTable(out);
		out.close();
	}

	//print final table
	filename = outputName + "_qualityTransformation.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";
	QT.printTable(out);
	out.close();

	//clean up memory
	windows.clear();
}
