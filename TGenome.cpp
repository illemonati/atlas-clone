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
	initializeRandomGenerator(params);

	//read parameters
	filename = params.getParameterString("bam");
	windowSize = params.getParameterDoubleWithDefault("window", 100000);
	numWindowsOnChr = 0;
	//if(windowSize < 1000) throw "Window size should be at least 1Kb!";
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 1.0);

	//outputname
	outputName = params.getParameterStringWithDefault("out", "");
	if(outputName == ""){
		//guess from filename
		outputName = filename;
		outputName = extractBeforeLast(outputName, ".");
	}
	logfile->list("Writing output files with prefix '" + outputName + "'.");

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;
	curStart = -1;
	curEnd = -1;
	oldPos = -1;
	oldAlignementMustBeConsidered = false;

	//open BAM file
	logfile->list("Reading data from BAM file '" + filename + "'.");
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";

	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + filename + "'!";

	//read header
	bamHeader = bamReader.GetHeader();
	readGroups.fill(bamHeader);
	chrIterator = bamHeader.Sequences.End();

	//open FASTA reference
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		std::string fastaIndex = fastaFile + ".fai";
		logfile->list("Reading reference sequence from '" + fastaFile + "'");
		if(!reference.Open(fastaFile, fastaIndex)) throw "Failed to open FASTA file '" + fastaFile + "'! Is index file present?";
		fastaReference = true;
	} else fastaReference = false;

	//initialize post mortem damage
	initializePostMortemDamage(params);
	doRecalibration = false;
	recalObjectInitialized = false;

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

	if(params.parameterExists("maskCpG")){
		if(!fastaReference) throw "Cannot mask CpG sites without reference!";
		doCpGMasking = true;
		std::string maskFile = params.getParameterString("maskCpG");
	} else doCpGMasking = false;

	//limit chrs and / or windows
	useChromosome = new bool[bamHeader.Sequences.Size()];
	if(params.parameterExists("chr")){
		logfile->startIndent("Will limit analysis to the following chromosomes:");

		//set all chromosomes to false
		for(int i=0; i<bamHeader.Sequences.Size(); ++i)
				useChromosome[i] = false;

		//parse chromosome names
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("chr"), vec, ',');
		int num;
		for(std::vector<std::string>::iterator it=vec.begin(); it!=vec.end(); ++it){
			//find chromosome
			num = 0;
			for(chrIterator = bamHeader.Sequences.Begin(); chrIterator != bamHeader.Sequences.End(); ++chrIterator, ++num){
				if(chrIterator->Name == *it){
					useChromosome[num] = true;
					logfile->list(*it);
					break;
				}

			}
			if(chrIterator == bamHeader.Sequences.End()) throw "Chromosome '" + *it + "' is not present in the bam header!";
		}
		chrIterator = bamHeader.Sequences.End();
		limitChr = 1000000;
		logfile->endIndent();
	} else {
		limitChr = params.getParameterIntWithDefault("limitChr", 1000000);
		if(params.parameterExists("limitWindows")) logfile->list("Will limit analysis to the first " + toString(limitChr) + " chromosomes.");
		for(int i=0; i<bamHeader.Sequences.Size(); ++i)
			useChromosome[i] = true;
	}
	limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows")) logfile->list("Will limit analysis to the first " + toString(limitWindows) + " windows per chromosome.");

};

void TGenome::jumpToEnd(){
	chrIterator = bamHeader.Sequences.End();
	chrNumber = -1;
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

	//do we use this chromosome? if not, move on!
	while(chrIterator != bamHeader.Sequences.End() && !useChromosome[chrNumber]){
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
	oldAlignementMustBeConsidered = false;

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

	//jump reader
	//bamReader.Jump(chrNumber, curStart);

	//report
	logfile->number("Window [" + toString(curStart) + ", " + toString(curEnd) + ") of " + toString(numWindowsOnChr) + " on '" + chrIterator->Name + "':");
	logfile->addIndent();

	return true;
};

bool TGenome::addAlignementToWindows(BamTools::BamAlignment & alignement, TWindowPair & windowPair){
	//std::cout << "REF ID = " << bamAlignement.RefID << "\tpos = " << bamAlignement.Position << std::endl;
	//only take those reads that pass QC
	if(!alignement.IsFailedQC() && !alignement.IsDuplicate() && alignement.Position >= curStart){
		//check if bam file is sorted
		if(alignement.Position < oldPos){
			throw "BAM file must be sorted by position!";
		}
		oldPos = alignement.Position;

		//check if still within current window and add to window
		if(alignement.Position >= curEnd) return false;
		else {
			if(windowPair.curPointer->addFromRead(alignement, pmdObjects, &readGroups)){
				//add also to next window in case reads overhangs current window -> function returns true
				windowPair.nextPointer->addFromRead(alignement, pmdObjects, &readGroups);
			}
		}
	}
	return true; //continue
}

bool TGenome::readData(TWindowPair & windowPair){
	logfile->listFlush("Reading data ...");

	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//parse through reads
	if(oldAlignementMustBeConsidered){
		oldAlignementMustBeConsidered = false;
		if(!addAlignementToWindows(bamAlignement, windowPair)){
			//next read is for a later window
			oldAlignementMustBeConsidered = true;
			gettimeofday(&end, NULL);
			logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");
			logfile->conclude("No data in this window.");
			return false; //still only in next window
		}
	}

	while(bamReader.GetNextAlignment(bamAlignement) && bamAlignement.RefID==chrNumber){
		if(!addAlignementToWindows(bamAlignement, windowPair)){
			//read is beyond window and should be reconsidered
			oldAlignementMustBeConsidered = true;
			break;
		}
	}

	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");

	if(windowPair.curPointer->numReadsInWindow > 0){
		//apply masks
		if(doMasking){
			logfile->listFlush("Masking sites ...");
			windowPair.curPointer->applyMask(mask);
			logfile->write(" done!");
		}
		if(doCpGMasking){
			logfile->listFlush("Masking CpG sites ...");
			windowPair.curPointer->maskCpG(reference, chrNumber);
			logfile->write(" done!");
		}

		//calc coverage
		windowPair.curPointer->calcCoverage();

		//report
		logfile->conclude("read data from " + toString(windowPair.curPointer->numReadsInWindow) + " reads.");
		logfile->conclude("coverage is " + toString(windowPair.curPointer->coverage));
		logfile->conclude(toString(windowPair.curPointer->fractionsitesCoverageAtLeastTwo * 100) + "% of all sites are covered at least twice");
		logfile->conclude(toString(windowPair.curPointer->fractionSitesNoData * 100) + "% of all sites have no data");
		return true;
	} else {
		logfile->conclude("No data in this window.");
		return false;
	}
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
	if(params.parameterExists("recal")){
		recalObject = new TRecalibrationEM(&bamHeader, params, logfile);
		doRecalibration = true;
	} else if(params.parameterExists("BQSRQuality")){
		recalObject = new TRecalibrationBQSR(&bamHeader, params, logfile);
		doRecalibration = true;
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
		doRecalibration = false;
		recalObject = new TRecalibration();
	}
	recalObjectInitialized = true;

	//check if estimation is required, in which case throw an error!
	if(recalObject->requiresEstimation()) throw "Can not use provided recalibration: estimation is required!";
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
	//initialize recalibration
	initializeRecalibration(params);

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
			if(readData(windows)){
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
	}

	//clean up
	out.close();
}

void TGenome::calcLikelihoodSurfaces(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

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
			if(readData(windows)){
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
}

void TGenome::callMLEGenotypes(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = true;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile);
		limitToSitesWithKnownAlleles = true;
	} else {
		printIfNoData = params.parameterExists("printAll");
		if(printIfNoData) logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream out;
	std::string outputFileName;
	if(params.parameterExists("vcf")){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;

		//open file
		outputFileName = outputName + "_MLEGenotypes.vcf.gz";
		logfile->list("Writing MLE genotypes in VCF format to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		out << "##fileformat=VCFv4.2\n";
		out << "##source=estimHet\n";
		out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		if(!limitToSitesWithKnownAlleles) out << "##INFO=<ID=GG,Number=10,Type=Integer,Description=\"Phred-scaled relative likelihoods of all genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT and TT\">\n";
		out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		out << "##FORMAT=<ID=PL,Number=1,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n";
		out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		outputFileName = outputName + "_MLEGenotypes.txt.gz";
		logfile->list("Writing MLE genotypes to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		out << "chr\tpos";
		if(limitToSitesWithKnownAlleles){
			out << "\tRef\tAlt\tcoverage\tL(RR)\tL(RA)\tL(AA)\tMLE\tQ\n";
		} else {
			if(fastaReference) out << "\tRef";
			out << "\tcoverage\tL(AA)\tL(AC)\tL(AG)\tL(AT)\tL(CC)\tL(CG)\tL(CT)\tL(GG)\tL(GT)\tL(TT)\tMLE\tQ\n";
		}
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				//read data for current window
				if(readData(windows)){
					//call genotypes
					logfile->listFlush("Calling MLE genotypes ...");
					if(limitToSitesWithKnownAlleles){
						windows.cur->addReferenceBaseToSites(subset);
						windows.cur->callMLEGenotypeKnownAlleles(recalObject, subset, *randomGenerator, out, chrIterator->Name, writeVCF);
					} else {
						if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
						windows.cur->callMLEGenotype(recalObject, *randomGenerator, out, chrIterator->Name, printIfNoData, fastaReference, writeVCF);
					}
					logfile->write(" done!");
				}
			}
		}
	}

	//clean up
	out.close();
}

void TGenome::callBayesianGenotypes(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	double theta;
	bool estimateTheta;
	EMParameters* EMParams = NULL;
	std::ofstream outTheta;
	if(params.parameterExists("theta")){
		estimateTheta = false;
		theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
	} else {
		estimateTheta = true;
		//read EM params
		EMParams = new EMParameters(params, logfile);
		openThetaOutputFile(outTheta);
	}

	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = true;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile);
		limitToSitesWithKnownAlleles = true;
	} else {
		printIfNoData = params.parameterExists("printAll");
		if(printIfNoData) logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream output;
	std::string outputFileName;
	if(params.parameterExists("vcf")){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;
		//open file
		outputFileName = outputName + "_BayesianGenotypes.vcf.gz";
		logfile->list("Writing Bayesian genotypes in VCF format to '" + outputFileName + "'");
		output.open(outputFileName.c_str());
		if(!output) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		output << "##fileformat=VCFv4.2\n";
		output << "##source=estimHet\n";
		output << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		if(!limitToSitesWithKnownAlleles) output << "##INFO=<ID=PP,Number=10,Type=Integer,Description=\"Phred-scaled posterior probabilities of all genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT and TT\">\n";
		output << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		output << "##FORMAT=<ID=GP,Number=1,Type=Integer,Description=\"Genotype posterior probabilities (phred-scaled)\">\n";
		output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		outputFileName = outputName + "_BayesianGenotypes.txt.gz";
		logfile->list("Writing Bayesian genotypes to '" + outputFileName + "'");
		output.open(outputFileName.c_str());
		if(!output) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		output << "chr\tpos";
		if(fastaReference) output << "\tRef";
		if(limitToSitesWithKnownAlleles) output << "\talt\tcoverage\tP(RR|D)\tP(RA|D)\tP(AA|D)\tMAP\tQ\n";
		else output << "\tcoverage\tP(AA|D)\tP(AC|D)\tP(AG|D)\tP(AT|D)\tP(CC|D)\tP(CG|D)\tP(CT|D)\tP(GG|D)\tP(GT|D)\tP(TT|D)\tMAP\tQ\n";
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			//read data for current window
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				if(readData(windows)){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} else {
						//set Theta
						if(estimateTheta){
							outTheta << chrIterator->Name << "\t";
							windows.cur->estimateTheta((*EMParams), recalObject, outTheta, logfile);
						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							windows.cur->setTheta(theta);
						}

						//call Bayesian genotypes
						logfile->listFlush("Calling Bayesian Genotypes ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callBayesianGenotypeKnownAlleles(subset, *randomGenerator, output, chrIterator->Name, writeVCF);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callBayesianGenotype(*randomGenerator, output, chrIterator->Name, printIfNoData, fastaReference, writeVCF);
						}
						logfile->write(" done!");
					}
				}
			}
		}
	}

	//clean up
	if(estimateTheta){
		outTheta.close();
		delete EMParams;
	}
	if(limitToSitesWithKnownAlleles) delete subset;
}

void TGenome::callAllelePresence(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	double theta;
	bool estimateTheta;
	EMParameters* EMParams = NULL;
	std::ofstream outTheta;
	if(params.parameterExists("theta")){
		estimateTheta = false;
		theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
	} else {
		estimateTheta = true;
		openThetaOutputFile(outTheta);
		//read EM params
		EMParams = new EMParameters(params, logfile);
	}


	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = true;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile);
		limitToSitesWithKnownAlleles = true;
	} else {
		printIfNoData = params.parameterExists("printAll");
		if(printIfNoData) logfile->list("Will print all sites, even those without data");
	}


	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream outAllelePresence;
	std::string outputFileName;
	if(params.parameterExists("vcf")){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;

		//open file
		outputFileName = outputName + "_AllelePresence.vcf.gz";
		logfile->list("Writing estimates of allele presence in VCF format to '" + outputFileName + "'");
		outAllelePresence.open(outputFileName.c_str());
		if(!outAllelePresence) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		outAllelePresence << "##fileformat=VCFv4.2\n";
		outAllelePresence << "##source=estimHet\n";
		outAllelePresence << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		outAllelePresence << "##INFO=<ID=PP,Number=10,Type=Integer,Description=\"Phred-scaled posterior probabilities of allele presence in the order A, C, G and T\">\n";
		outAllelePresence << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		outAllelePresence << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		outputFileName = outputName + "_AllelePresence.txt.gz";
		logfile->list("Writing estimates of allele presence to '" + outputFileName + "'");
		outAllelePresence.open(outputFileName.c_str());
		if(!outAllelePresence) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		outAllelePresence << "chr\tpos";
		if(limitToSitesWithKnownAlleles) outAllelePresence << "\tRef\tAlt\tcoverage\tP(Ref|D)\tP(Alt|D)\tMAP\tQ\n";
		else {
			if(fastaReference) outAllelePresence << "\tRef";
			outAllelePresence << "\tcoverage\tP(A|D)\tP(C|D)\tP(G|D)\tP(T|D)\tMAP\tQ\n";
		}
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				//read data for current window
				if(readData(windows)){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} else {
						//set Theta
						if(estimateTheta){
							outTheta << chrIterator->Name << "\t";
							windows.cur->estimateTheta((*EMParams), recalObject, outTheta, logfile);
						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							windows.cur->setTheta(theta);
						}

						//call allele presence
						logfile->listFlush("Calling allele presence ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callAllelePresenceKnwonAlleles(subset, *randomGenerator, outAllelePresence, chrIterator->Name, writeVCF);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callAllelePresence(*randomGenerator, outAllelePresence, chrIterator->Name, printIfNoData, fastaReference, writeVCF);
						}
						logfile->write(" done!");
					}
				}
			} else logfile->list("No positions in this window.");
		}
	}

	//clean up
	if(estimateTheta){
		outTheta.close();
		delete EMParams;
	}
	if(limitToSitesWithKnownAlleles) delete subset;
}

void TGenome::printPileup(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

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
			if(readData(windows)){
				//print pileup
				windows.cur->printPileup(recalObject, out, chrIterator->Name);
			}
		}
	}

	//clean up
	out.close();
}


void TGenome::estimateErrorCalibrationEM(TParameters & params){
	//create recalibration object
	TRecalibrationEM recalObjectEM(&bamHeader, params, logfile);
	if(!recalObjectEM.estimatetionRequired){
		logfile->list("No need to estimate anything. Aborting Program.");
		return;
	}

	//prepare windows
	TWindowPairHaploid windows;

	//add sites to EM object
	logfile->startIndent("Reading data from windows:");
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				windows.cur->addToRecalibrationEM(recalObjectEM);
			}
		}
	}
	//clean up memory
	windows.clear();
	logfile->endIndent();

	//run EM iterations
	recalObjectEM.runEM(outputName);
}

/*
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
}*/

void TGenome::calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params){
	//create recalibration object
	TRecalibrationEM recalObjectEM(&bamHeader, params, logfile);

	//prepare windows
	TWindowPairHaploid windows;

	//add sites to EM object
	logfile->startIndent("Reading data from windows:");
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				windows.cur->addToRecalibrationEM(recalObjectEM);
			}
		}
	}
	//clean up memory
	windows.clear();
	logfile->endIndent();

	//calc likelihood surface
	int numMarginalGridPoint = params.getParameterIntWithDefault("numGridPoints", 51);
	recalObjectEM.calcLikelihoodSurface(outputName + "_LLsurface.txt", numMarginalGridPoint);
}

void TGenome::BQSR(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not run BQSR recalibration without a provided FASTA reference!";

	//prepare windows
	TWindowPairHaploid windows;

	//create BQSR object
	TRecalibrationBQSR bqsr(&bamHeader, params, logfile);
	if(bqsr.allConverged()){
		logfile->list("No need to estimate any BQSR cells. Aborting Program.");
		return;
	}

	//read in max number of loops allowed
	int maxNumLoops = params.getParameterIntWithDefault("maxNumLoops", 500);

	//tmp variables
	bool hasConverged = false;
	int loopNumber = 0;

	//do we print temporary tables?
	bool printTmpTables = params.parameterExists("writeTmpTables");
	if(printTmpTables) logfile->list("Temporary BQSR tables will be written to disk.");

	/*
	//check if we use first chromosome for initial convergence
	if(params.parameterExists("preConverge")){
		int numPreConvLoops = params.getParameterInt("preConverge");
		int numChrPreConv = params.getParameterIntWithDefault("limitChrPreConv", 1) - 1;
		int numWindowsPreConv = params.getParameterIntWithDefault("limitWindowsPreConv", 100);
		logfile->startIndent("Running initial BQSR (pre-convergence):");
		logfile->list("Limiting pre-convergence to the first " + toString(numChrPreConv) + " chromosomes.");
		logfile->list("Limiting pre-convergence to the first " + toString(numWindowsPreConv) + " windows on each chromosome.");

		//run until it converges
		while(!hasConverged && loopNumber < numPreConvLoops){
			++loopNumber;
			logfile->startIndent("Running pre-recalibration loop " + toString(loopNumber) + ":");

			//jump to first chromosome
			while(chrNumber < numChrPreConv && iterateChromosome(windows)){
				//iterate over all windows
				while(windowNumber < numWindowsPreConv && iterateWindow(windows)){
					//read data for current window
					if(readData(windows)){
						//add reference data
						windows.cur->addReferenceBaseToSites(reference, chrNumber);

						//add the base to BQSR
						windows.cur->addSitesToBQSR(bqsr, logfile);
					}
				}
				logfile->endIndent();
			}
			logfile->endIndent();

			//clean up memory and restart from chr 1
			windows.clear();
			jumpToEnd();

			//estimate epsilon
			hasConverged = bqsr.estimateEpsilon();

			//write results to file
			bqsr.writeToFile(outputName + "_preConvergence_Loop" + toString(loopNumber));
			logfile->endIndent();
		}

		//reset counters and such
		bqsr.reopenEstimation();
		hasConverged = false;
		loopNumber = 0;
		logfile->endIndent();
		jumpToEnd();
	} */

	//loop over bam until BQSR converges
	while(!hasConverged){
		++loopNumber;
		logfile->startIndent("Running recalibration loop " + toString(loopNumber) + ":");
		//loop over all windows

		if(!bqsr.dataHasBeenStored()){
			logfile->startIndent("Reading data from BAM file:");
			while(iterateChromosome(windows)){
				while(iterateWindow(windows)){
					//read data for current window
					if(readData(windows)){
						//add reference data
						windows.cur->addReferenceBaseToSites(reference, chrNumber);

						//add the base to BQSR
						windows.cur->addSitesToBQSR(bqsr, logfile);
					}
				}
			}
			//clean up memory
			windows.clear();
			logfile->endIndent();
		}

		//estimate epsilon
		hasConverged = bqsr.estimateEpsilon(outputName);

		//write results to file
		if(printTmpTables) bqsr.writeCurrentTmpTable(outputName + "_Loop" + toString(loopNumber));

		logfile->endIndent();

		//check if max num loops is reached
		if(loopNumber >= maxNumLoops){
			logfile->write("Reached maximum number of loops (" + toString(maxNumLoops) + "), but BQSR has not yet converged!");
			break;
		}
	}
}


void TGenome::printQualityTransformation(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//prepare windows
	TWindowPairHaploid windows;

	//create table to store counts
	int maxQ = params.getParameterIntWithDefault("maxQ", 100);
	TQualityTransformTable QT(maxQ);

	//prepare output
	std::ofstream out;
	std::string filename;

	//loop over windows and add to table
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				windows.cur->addSitesToQualityTransformTable(recalObject, QT, logfile);
			}
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

void TGenome::recalibrateBamFile(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_recalibrated.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//create array of TBase
	TBase* basePointer;

	//other temp variables
	char base, quality;
	BaseContext context;
	TGenotypeMap genoMap;
	int posInRead, revPosInRead;
	double pmdCT, pmdGA;
	int len;
	std::string qual;
	std::string readGroup;
	int readGroupId;
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignement)){
		++counter;
		len = bamAlignement.Length;
		qual.clear();

		//get readgroup info
		bamAlignement.GetTag("RG", readGroup);
		readGroupId = readGroups.find(readGroup);

		//parse into bases
		for(int pos = 0; pos < len; ++pos){
			base = bamAlignement.QueryBases.at(pos);
			if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
				quality = bamAlignement.Qualities.at(pos);
				if((int) quality > 33){
					if(bamAlignement.IsReverseStrand()){
						if(pos == (len - 1)) context = genoMap.getContext('N', base);
						else context = genoMap.getContext(bamAlignement.QueryBases.at(pos + 1), base);
						posInRead = len - pos - 1;
						revPosInRead = pos;
						pmdCT = pmdObjects[readGroupId].getProbCT(len - pos);
						pmdGA = pmdObjects[readGroupId].getProbGA(pos + 1);
					} else {
						if(pos == 0) context = genoMap.getContext('N', base);
						else context = genoMap.getContext(bamAlignement.QueryBases.at(pos - 1), base);
						posInRead = pos;
						revPosInRead = len - pos - 1;
						pmdCT = pmdObjects[readGroupId].getProbCT(pos + 1);
						pmdGA = pmdObjects[readGroupId].getProbGA(len - pos);
					}

					//create base object
					if(base == 'A') basePointer = new TBaseDiploidA(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);
					else if(base == 'C') basePointer = new TBaseDiploidC(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);
					else if(base == 'G') basePointer = new TBaseDiploidG(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);
					else basePointer = new TBaseDiploidT(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);

					//get new quality
					qual += recalObject->getQualityAsChar(basePointer);

					//delete base
					delete basePointer;
				} else qual += quality;
			} else qual += quality;
		}

		//update and write
		bamAlignement.Qualities = qual;
		bamWriter.SaveAlignment(bamAlignement);

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAm file!");
	logfile->removeIndent();
}

void TGenome::splitSingleEndReadGroups(TParameters & params){
	//read read groups and their expected lengths
	std::string filename = params.getParameterString("readGroups");
	logfile->listFlush("Reading single end read groups from file '" + filename + "' ...");
	std::map<int, TReadGroupMaxLength> singleEndRG;
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "!";

	//parse file
	int lineNum = 0;
	std::vector<std::string> vec;
	int len;
	int readGroupId;
	int truncatedReadGroupId;
	std::string readGroup;
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		if(!vec.empty()){
			if(vec.size() != 2) throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + filename + "'!";
			readGroupId = readGroups.find(vec[0]);
			len = stringToInt(vec[1]);
			if(len < 1) throw "Max length of read group '" + vec[0] + "' is < 1!";
			//add a new readgroup for the truncated reads to the header
			readGroup = vec[0] + "_truncated";
			bamHeader.ReadGroups.Add(readGroup);
			readGroups.fill(bamHeader);
			truncatedReadGroupId = readGroups.find(readGroup);

			//ad to map
			singleEndRG.emplace(readGroupId, TReadGroupMaxLength(len, truncatedReadGroupId, readGroup));
		}
	}
	logfile->write(" done!");
	logfile->conclude("read " + toString(singleEndRG.size()) + " read groups to be splitted.");

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_splitRG.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignement)){
		++counter;

		//get read group info
		bamAlignement.GetTag("RG", readGroup);
		readGroupId = readGroups.find(readGroup);

		//check if this RG needs to be parse
		singleEndRGIT = singleEndRG.find(readGroupId);
		if(singleEndRGIT != singleEndRG.end()){
			//check length
			if(bamAlignement.Length < singleEndRGIT->second.maxLen)
				bamAlignement.EditTag("RG", "Z", singleEndRGIT->second.truncatedReadGroup);
			else if(bamAlignement.Length > singleEndRGIT->second.maxLen) throw "Length of read in read group '" + readGroup + "' is > max length provided!";
		}

		//write
		bamWriter.SaveAlignment(bamAlignement);

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}


void TGenome::mergeReadGroups(TParameters & params){
	//read read groups and their expected lengths
	std::string filename = params.getParameterString("readGroups");
	logfile->listFlush("Reading read groups to be merged from file '" + filename + "' ...");
	std::vector< std::vector<std::string> > readGroupsToMerge;
	std::vector< std::vector<std::string> >::reverse_iterator rIt;
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "!";

	//construct new read groups in new header object
	BamTools::SamHeader newHeader(bamHeader);
	newHeader.ReadGroups.Clear();

	//parse file and construct new read groups in new header object
	int lineNum = 0;
	std::vector<std::string> vec;
	std::string readGroup;
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		if(!vec.empty()){
			if(vec.size() < 2) throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + filename + "'!";
			//add to new header
			newHeader.ReadGroups.Add(vec[0]);
			//others are those to be merged: find read group in header and store int
			readGroupsToMerge.push_back(std::vector<std::string>());
			rIt = readGroupsToMerge.rbegin();
			for(unsigned int i=1; i<vec.size(); ++i){
				rIt->push_back(vec[i]);
			}
		}
	}
	TReadGroups newReadGroupObject;
	newReadGroupObject.fill(newHeader);
	logfile->write(" done!");

	//report and construct map
	int* readGroupMap = new int[bamHeader.ReadGroups.Size()];
	for(int i=0; i<bamHeader.ReadGroups.Size(); ++i) readGroupMap[i] = -1;
	logfile->startIndent("The following read groups will be merged:");
	std::vector< std::vector<std::string> >::iterator mergeIt = readGroupsToMerge.begin();
	int oldId;
	for(int rg = 0; rg < newReadGroupObject.size(); ++rg, ++mergeIt){
		logfile->startIndent("New read group '" + newReadGroupObject.getName(rg) + "' will contain read groups:");
		for(std::vector<std::string>::iterator it = mergeIt->begin(); it != mergeIt->end(); ++it){
			logfile->list(*it);
			oldId = readGroups.find(*it);
			if(readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
			readGroupMap[oldId] = rg;
		}
		logfile->endIndent();
	}
	logfile->endIndent();

	//now add read groups that will not be merged
	bool printed = false;
	std::string name;
	for(int i = 0; i < readGroups.size(); ++i){
		//check if it is mapped, otherwise add
		if(readGroupMap[i] < 0){
			if(!printed){
				logfile->startIndent("The following read groups will be kept as is:");
				printed = true;
			}
			name = readGroups.getName(i);
			logfile->list(name);
			newHeader.ReadGroups.Add(name);
			newReadGroupObject.fill(newHeader);
			readGroupMap[i] = newReadGroupObject.find(name);
		}
	}
	if(printed) logfile->endIndent();
	else logfile->list("All existing read groups will be merged into a new read group.");

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_mergedRG.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, newHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignement)){
		++counter;

		//get read group info
		bamAlignement.GetTag("RG", readGroup);
		oldId = readGroups.find(readGroup);

		//save as new RG
		bamAlignement.EditTag("RG", "Z", newReadGroupObject.getName(readGroupMap[oldId]));
		bamWriter.SaveAlignment(bamAlignement);

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::addReadToPMD(TWindowDiploid* window, TGenotypeMap & genoMap, std::string & ref, TPMDTables & pmdTables){
	//Extract Read Group Info
	std::string readGroup;
	bamAlignement.GetTag("RG", readGroup);
	int readGroupId = readGroups.find(readGroup);
	int length = bamAlignement.AlignedBases.size();
	char base;
	int quality;
	Base readBase, refBase;

	//add to PMD
	//distinguish between cases
	int internalPos = bamAlignement.Position - window->start;
	if(bamAlignement.IsProperPair()){
		throw "Not yet done for paired end!";
	} else {
		//single end
		if(bamAlignement.IsReverseStrand()){
			//single end & reverse
			//forward position = len - pos - 1
			//reverse position = pos
			//FLIP BASES!
			for(int pos = 0; pos < length; ++pos, ++internalPos){
				base = bamAlignement.AlignedBases.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip ann other
					quality = bamAlignement.AlignedQualities.at(pos);
					if((int) quality > 32){ //skip if quality dies not make sense
						readBase = genoMap.flipBase(base);
						//std::cout << " " << internalPos << "," << ref[internalPos] << std::flush;
						refBase = genoMap.flipBase(ref[internalPos]);

						if(refBase == N){
							std::cout << "ISSUE WITH REF at pos = " << internalPos << " = " << internalPos + window->start << std::endl;
						}

						pmdTables.addForward(readGroupId, length - pos - 1, refBase, readBase);
						pmdTables.addReverse(readGroupId, pos, refBase, readBase);
					}
				}
			}
		} else {
			//single end & forward
			//forward position = pos
			//reverse position = len - pos -1
			for(int pos = 0; pos < length; ++pos, ++internalPos){
				base = bamAlignement.AlignedBases.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
					quality = bamAlignement.AlignedQualities.at(pos);
					if((int) quality > 32){ //skip if quality dies not make sense
						readBase = genoMap.getBase(base);
						refBase = genoMap.getBase(ref[internalPos]);
						pmdTables.addForward(readGroupId, pos, refBase, readBase);
						pmdTables.addReverse(readGroupId, length - pos - 1, refBase, readBase);
					}
				}
			}
		}
	}
}

void TGenome::estimatePMD(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//prepare windows
	TWindowPairDiploid windows;

	//prepare PMD table
	int maxLength = params.getParameterIntWithDefault("length", 100);
	logfile->list("Estimating PMD at the first " + toString(maxLength) + " positions.");
	TPMDTables pmdTables(&readGroups, maxLength);

	//measure runtime
	struct timeval start, end;

	//tmp variables
	int fastaEnd;
	std::string ref;
	TGenotypeMap genoMap;
	long numreadsAdded;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			gettimeofday(&start, NULL);
			logfile->listFlush("Adding reads to PMD tables ...");
			numreadsAdded = 0;

			//get fasta reference
			fastaEnd = windows.cur->end + 500; //note that end is last position + 1
			reference.GetSequence(chrNumber, windows.cur->start, fastaEnd, ref);

			//still have old alignment to condider?
			if(oldAlignementMustBeConsidered && bamAlignement.Position >= windows.cur->end){
					logfile->write(" done!");
					logfile->conclude("No data in this window.");
			} else {
				if(oldAlignementMustBeConsidered){
					//add this read to window
					oldAlignementMustBeConsidered = false;
					++numreadsAdded;
					addReadToPMD(windows.cur, genoMap, ref, pmdTables);
				}

				//parse additional reads
				while(bamReader.GetNextAlignment(bamAlignement) && bamAlignement.RefID==chrNumber){
					//check if this read is beyond window
					if(bamAlignement.Position >= windows.cur->end){
						oldAlignementMustBeConsidered = true;
						break;
					}
					addReadToPMD(windows.cur, genoMap, ref, pmdTables);
				}

				//report
				gettimeofday(&end, NULL);
				logfile->write(" done (added " + toString(numreadsAdded) + " reads in " + toString(end.tv_sec  - start.tv_sec) + "s)!");
			}
		}
	}

	//print tables and data
	std::string filename = outputName + "_PMD_Table.txt";
	logfile->listFlush("Writing PMD table to '" + filename + "' ...");
	pmdTables.writeTable(filename);
	logfile->write(" done!");
	filename = outputName + "_PMD_Table_counts.txt";
	logfile->listFlush("Writing PMD table of counts to '" + filename + "' ...");
	pmdTables.writeTableWithCounts(filename);
	logfile->write(" done!");
	filename = outputName + "_PMD_input_Table.txt";
	logfile->listFlush("Writing PMD input file to '" + filename + "' ...");
	pmdTables.writePMDFile(filename);
	logfile->write(" done!");
}

void TGenome::mergePairedEndReads(TParameters & params){
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_mergedReads.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;

	//create storage for reads until their mate was found
	std::vector< std::pair<BamTools::BamAlignment*, bool> > alignmentStorage;
	std::vector< std::pair<BamTools::BamAlignment*, bool> >::iterator it;
	BamTools::BamAlignment* alignmentPointer;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;
	int curChr = -1;


    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignement)){
		++counter;

		//if on new chromosome, empty storage
		if(curChr != bamAlignement.RefID){
			for(it = alignmentStorage.begin(); it != alignmentStorage.end(); ++it){
				bamWriter.SaveAlignment(*(it->first));
				delete it->first;
			}
			alignmentStorage.clear();
			curChr = bamAlignement.RefID;
		}

		//add alignment to storage
		if(bamAlignement.IsPaired()){
			if(bamAlignement.IsFirstMate()) alignmentStorage.push_back(std::pair<BamTools::BamAlignment*, bool>(new BamTools::BamAlignment(bamAlignement), false));
			else if(bamAlignement.IsSecondMate()){
				//find first mate -> should be in storage
				for(it=alignmentStorage.begin(); it!=alignmentStorage.end(); ++it){
					if(it->first->Name == bamAlignement.Name){
						//check if this read accepts mate
						if(it->second) throw "First read of '" + bamAlignement.Name + "' is not paired or has already been merged!";

						//merge reads
						alignmentPointer = it->first;
						if(bamAlignement.Position > alignmentPointer->Position + alignmentPointer->Length){
							//reads do not overlap -> add Ns in between
							int numN = bamAlignement.Position - alignmentPointer->Position + alignmentPointer->AlignedBases.size() - 1;
							for(int i=0; i<numN; ++i){
								alignmentPointer->AlignedBases += 'N';
								alignmentPointer->AlignedQualities += '!';
							}
							alignmentPointer->AlignedBases += bamAlignement.AlignedBases;
							alignmentPointer->AlignedQualities += bamAlignement.AlignedQualities;

						} else if(bamAlignement.Position < alignmentPointer->Position) throw "Second ate of read '" + bamAlignement.Name + "' ha spos < pos of first mate!";
						else {
							//reads do overlap
							std::string alignment;
							std::string quality;
							int firstOverlap = bamAlignement.Position - alignmentPointer->Position;
							int lastOverlapPlusOne = alignmentPointer->AlignedBases.size() + 1;

							//first copy from first mate
							alignment = alignmentPointer->AlignedBases.substr(0, firstOverlap);
							quality = alignmentPointer->AlignedQualities.substr(0, firstOverlap);

							//decide which alignment has higher quality in overlap
							for(int i=firstOverlap; i<lastOverlapPlusOne; ++i){
								//decide which quality is higher
								if(alignmentPointer->AlignedQualities.at(i) > bamAlignement.AlignedQualities.at(i - firstOverlap)){
									alignment += alignmentPointer->AlignedBases.at(i);
									quality += alignmentPointer->AlignedQualities.at(i);
								} else {
									alignment += bamAlignement.AlignedBases.at(i - firstOverlap);
									quality += bamAlignement.AlignedQualities.at(i - firstOverlap);
								}
							}

							//add rest from second
							if(alignmentPointer->Position + alignmentPointer->AlignedBases.size() < bamAlignement.Position + bamAlignement.AlignedBases.size()){
								alignment += bamAlignement.AlignedBases.substr(lastOverlapPlusOne - firstOverlap);
								quality += bamAlignement.AlignedQualities.at(lastOverlapPlusOne - firstOverlap);
							}

							//set
							alignmentPointer->AlignedBases += alignment;
							alignmentPointer->AlignedQualities += quality;
						}

						//update
						alignmentPointer->QueryBases = alignmentPointer->AlignedBases;
						alignmentPointer->Qualities = alignmentPointer->AlignedQualities;
						alignmentPointer->Length = alignmentPointer->AlignedBases.size();
						alignmentPointer->CigarData.clear();
						alignmentPointer->CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MATCH_CHAR, alignmentPointer->Length));
						alignmentPointer->SetIsFirstMate(false);
						alignmentPointer->SetIsPaired(false);
						alignmentPointer->SetIsProperPair(false);
						alignmentPointer->SetIsSecondMate(false);
						it->second = true;

						//write if is first in vector
						if(it == alignmentStorage.begin()){
							//write all that are OK
							for(; it != alignmentStorage.end(); ++it){
								if(it->second){
									bamWriter.SaveAlignment(*(it->first));
									delete it->first;
								} else {
									//first that can not be written -> erease all previous ones!
									alignmentStorage.erase(alignmentStorage.begin(), it);
									break;
								}
							}

						}
						break;
					}
				}
				if(it == alignmentStorage.end()) throw "One read of '" + bamAlignement.Name + "' is second mate, but first one has not been read!";
			} else throw "One read of '" + bamAlignement.Name + "' is paired, but neither first nor second mate!";
		} else {
			//read is not paired: add to storage or write
			if(alignmentStorage.empty()) bamWriter.SaveAlignment(bamAlignement);
			else alignmentStorage.push_back(std::pair<BamTools::BamAlignment*, bool>(new BamTools::BamAlignment(bamAlignement), true));
		}

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}



void TGenome::generatePSMCInput(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//read in parameters required
	double theta = params.getParameterDoubleWithDefault("theta", 0.001);
	logfile->list("Using theta = " + toString(theta));
	double confidence = params.getParameterDoubleWithDefault("confidence", 0.99);
	logfile->list("Calling heterozygosity state with confidence > " + toString(confidence));
	int blockSize = params.getParameterIntWithDefault("block", 100);
	//make sure window size is a multiple of block length!
	if(windowSize % blockSize != 0) throw "Window size is not a multiple of block size!";

	//open output file
	std::ofstream output;
	std::string outputFileName = outputName + ".psmcfa";
	logfile->list("Writing PSMC input file to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int nCharOnLine = 0;

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		if(nCharOnLine > 0) output << '\n';
		output << '>' << chrIterator->Name << '\n';
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			//set Theta
			logfile->listFlush("Calculating emission probabilities ...");
			windows.cur->calculateEmissionProbabilities(recalObject);
			windows.cur->estimateBaseFrequencies();
			windows.cur->setTheta(theta);
			logfile->write(" done!");

			//create PSMC input
			logfile->listFlush("Estimating heterozygosity status ...");
			if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
			windows.cur->generatePSMCInput(blockSize, confidence, output, nCharOnLine);
			logfile->write(" done!");
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
}

void TGenome::downSampleBamFile(TParameters & params){
	//read downsampling rate
	std::string prob = params.getParameterString("prob");
	//check if prob is a vector of multiple probabilities
	std::vector<double> downSampleProbVector;
	if(!stringContainsOnly(prob, "-0123456789.,")) throw "Wrong format on probability list: use floating point numbers delimited by commas (e.g. 0.1,0.2,0.5).";
	fillVectorFromString(prob, downSampleProbVector, ',');

	//check if probs are between 0 and 1, save in array and print them
	std::vector<double>::iterator it;
	int numProbs = downSampleProbVector.size();
	double* downSampleProb = new double[numProbs];
	logfile->listFlush("Will accept reads with probabilities");
	bool first = true;
	int i=0;
	for(it=downSampleProbVector.begin(); it!=downSampleProbVector.end(); ++it, ++i){
		if(first) first = false;
		else logfile->flush(",");
		logfile->flush(" " + toString(*it));
		if(*it <= 0.0 || *it >= 1.0) throw "All probabilities have to be between >0 and <1!";
		downSampleProb[i] = *it;
	}

	//open bam files for writing
	BamTools::BamWriter* bamWriter = new BamTools::BamWriter[numProbs];
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->startIndent("Writing results to the following files:");
	for(i=0; i<numProbs; ++i){
		//construct and print filename
		filename = outputName + "_downsampled" + toString(downSampleProb[i]) + ".bam";
		logfile->list(filename);

		//open file
		if(!bamWriter[i].Open(filename, bamHeader, references))
				throw "Failed to open BAM file '" + filename + "'!";
	}
	logfile->endIndent();

	//other temp variables
	long counter = 0;
	double r;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignement)){
		++counter;

		//accept read or not?
		r = randomGenerator->getRand();
		for(i=0; i<numProbs; ++i){
			if(r < downSampleProb[i]){
				bamWriter[i].SaveAlignment(bamAlignement);
			}
		}

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer and clean up memory
	for(i=0; i<numProbs; ++i){
		bamWriter[i].Close();
	}
	delete[] downSampleProb;
	delete[] bamWriter;

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::estimateApproximateCoverage(TParameters & params){	//get genome length
	double totLength = 0.0;
	for(chrIterator = bamHeader.Sequences.Begin(); chrIterator!=bamHeader.Sequences.End(); ++chrIterator)
		totLength += stringToLong(chrIterator->Length);

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;

	//other temp variables
	long counter = 0;
	double toNumAlignedBases = 0.0;

    //now parse through bam file and sum number of aligned bases
	while (bamReader.GetNextAlignment(bamAlignement)){
		++counter;

		//accept read or not?
		toNumAlignedBases += bamAlignement.AlignedBases.length();

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

	//report approximate coverage
	logfile->list("Approximate coverage was estimated at " + toString(toNumAlignedBases/totLength));
}





















