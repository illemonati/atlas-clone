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

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;
	curStart = -1;
	curEnd = -1;
	oldPos = -1;
	oldAlignementMustBeConsidered = false;

	//open BAM file
	filename = params.getParameterString("bam");
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

	maxReadLength = params.getParameterIntWithDefault("maxReadLength", 1000);
	logfile->list("Will only consider reads up to " + toString(maxReadLength) + " bp.");

	//initialize alignment parser
	alignmentParser.init(&readGroups, maxReadLength, logfile);

	//read window parameter
	if(!params.parameterExists("window") && params.parameterExists("windows")) logfile->warning("Argument 'windows' specified, but unknown. Did you mean 'window'?");
	std::string tmp = params.getParameterStringWithDefault("window", "1000000");
	//check if it is a number
	if(stringContainsOnly(tmp, "1234567890.Ee-+")){
		windowsPredefined = false;
		windowSize = stringToInt(tmp);
		logfile->list("Setting window size to " + toString(windowSize));
		if(windowSize < maxReadLength) throw "Window size " + tmp + " out of range! Windows must be at least as large as the max read length (" + toString(maxReadLength) + " bp)!";
	} else {
		windowsPredefined = true;
		logfile->listFlush("Limiting analysis to windows defined in '" + tmp + "'...");
		predefinedWindows = new TBed(tmp);
		logfile->done();
		logfile->conclude("read " + toString(predefinedWindows->size()) + " on " + toString(predefinedWindows->getNumChromosomes()) + " chromosomes");
	}
	numWindowsOnChr = 0;

	//outputname
	outputName = params.getParameterStringWithDefault("out", "");
	if(outputName == ""){
		//guess from filename
		outputName = filename;
		outputName = extractBeforeLast(outputName, ".");
	}
	logfile->list("Writing output files with prefix '" + outputName + "'.");

	//open FASTA reference
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		std::string fastaIndex = fastaFile + ".fai";
		logfile->list("Reading reference sequence from '" + fastaFile + "'");
		if(!reference.Open(fastaFile, fastaIndex)) throw "Failed to open FASTA file '" + fastaFile + "'! Is index file present?";
		fastaReference = true;
		alignmentParser.addReference(&reference);
	} else fastaReference = false;

	//initialize post mortem damage
	hasPMD = false;
	initializePostMortemDamage(params);
	doRecalibration = false;
	recalObjectInitialized = false;

	//check if we mask sites
	if(params.parameterExists("mask")){
		if(windowsPredefined) throw "Masking is currently not implemented if windows are predefined from a BED file.";
		if(params.parameterExists("sites")) throw "Masking is currently not implemented if variant positions are also specified with \"sites\"";
		if(params.parameterExists("regions")) throw "Cannot use mask and regions at the same time";
		doMasking = true;
		std::string maskFile = params.getParameterString("mask");
		logfile->startIndent("Will mask all sites listed in BED file '" + maskFile + "':");
		logfile->listFlush("Reading file ...");
		mask = new TBedReader(maskFile, windowSize, bamHeader.Sequences, logfile);
		logfile->done();
		logfile->endIndent();
		//mask->print();
	} else doMasking = false;

	if(params.parameterExists("maskCpG")){
		if(!fastaReference) throw "Cannot mask CpG sites without reference!";
		doCpGMasking = true;
		std::string maskFile = params.getParameterString("maskCpG");
		logfile->list("Will mask all CpG sites");
	} else doCpGMasking = false;

	if(params.parameterExists("regions")){
		if(windowsPredefined) throw "Regions is currently not implemented if windows are predefined from a BED file.";
		if(params.parameterExists("sites")) throw "Regions is currently not implemented if variant positions are also specified with \"sites\"";
		considerRegions = true;
		std::string regionsFile = params.getParameterString("regions");
		logfile->startIndent("Will limit analysis to all regions listed in BED file '" + regionsFile + "':");
		logfile->listFlush("Reading file ...");
		mask = new TBedReader(regionsFile, windowSize, bamHeader.Sequences, logfile);
		logfile->done();
		logfile->endIndent();
	} else considerRegions = false;

	//filters
	if(params.parameterExists("minDepth") || params.parameterExists("maxDepth")){
		applyDepthFilter = true;
		int tmpInt;
		tmpInt = params.getParameterIntWithDefault("minDepth", 0);
		if(tmpInt < 0) throw "minDepth must be >= 0!";
		minDepth = tmpInt;
		tmpInt = params.getParameterIntWithDefault("maxDepth", 1000000);
		if(tmpInt < minDepth) throw "maxDepth must be >= minDepth!";
		maxDepth = tmpInt;
		logfile->list("Will filter out sites with sequencing depth < " + toString(minDepth) + " or > " + toString(maxDepth));
	} else {
		applyDepthFilter = false;
		minDepth = 0;
		maxDepth = 1000000;
	}

	//quality filters
	minPhredInt = params.getParameterIntWithDefault("minQual", 1);
	if(minPhredInt < 0) throw "minQual must be >= 0!";
	maxPhredInt = params.getParameterIntWithDefault("maxQual", 93);
	if(maxPhredInt < minPhredInt) throw "maxQual must be >= minQual!";
	alignmentParser.setQualityFilters(minPhredInt+33, maxPhredInt+33);
	logfile->list("Will filter out bases with quality outside the range [" + toString(minPhredInt) + ", " + toString(maxPhredInt) + "]");

	//quality filters for printing
	minOutQual = params.getParameterIntWithDefault("minOutQual", 1) + 33;
	if(minOutQual < 0) throw "minOutQual must be >= 0!";
	maxOutQual = params.getParameterIntWithDefault("maxOutQual", 93) + 33;
	if(maxOutQual < minOutQual) throw "maxOutQual must be >= minOutQual!";
	alignmentParser.setQualityRangeForPrinting(minOutQual, maxOutQual);
	logfile->list("Will print qualities truncated to [" + toString(minOutQual) + ", " + toString(maxOutQual) + "]");

	//trimming ends
	if(params.parameterExists("trim3") || params.parameterExists("trim5")){
		int trim3 = params.getParameterIntWithDefault("trim3", 0);
		if(trim3 < 0) throw "trimming distance trim3 must be >= 0!";
		int trim5 = params.getParameterIntWithDefault("trim5", 0);
		if(trim5 < 0) throw "trimming distance trim5 must be >= 0!";
		if(trim3>0 || trim5>0){
			alignmentParser.setReadTrimming(trim3, trim5);
			logfile->list("Will trim first " + toString(trim3) + " and " + toString(trim5) + " bases from the 3' and 5' end, respectively.");
		}
	}

	//other filters
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 1.0);
	if(maxMissing > 1.0) throw "maxMissing must be smaller or equal to 1.0!";

	maxRefN = params.getParameterDoubleWithDefault("maxRefN", 1.0);
	if(maxRefN > 1.0) throw "maxRefN must be smaller or equal to 1.0!";
	if(maxRefN < 1.0 && fastaReference == false) throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided.";

	if(params.parameterExists("keepDuplicates")){
		alignmentParser.keepDuplicates();
		logfile->list("Will keep duplicate reads.");
	}

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
		if(params.parameterExists("limitChr")) logfile->list("Will limit analysis to the first " + toString(limitChr) + " chromosomes.");
		for(int i=0; i<bamHeader.Sequences.Size(); ++i)
			useChromosome[i] = true;
	}
	limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows")) logfile->list("Will limit analysis to the first " + toString(limitWindows) + " windows per chromosome.");

	//limit readGroups
	if(params.parameterExists("readGroup")){
		limitReadGroups = true;
		fillVectorFromString(params.getParameterString("readGroup"), readGroupsInUse, ',');
		logfile->startIndent("Will limit analysis to the following read groups:");
		for(int i=0; i < readGroups.numGroups; i++){
			if(std::find(readGroupsInUse.begin(), readGroupsInUse.end(), readGroups.getName(i)) != readGroupsInUse.end()){
				readGroups.inUse[i] = true;
				logfile->list(readGroups.getName(i));
			} else readGroups.inUse[i] = false;
		}
		logfile->endIndent();
	} else limitReadGroups = false;

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
	chrLength = stringToLong(chrIterator->Length);

	//restart windows
	curStart = 0;
	curEnd = 0;
	oldPos = -1;
	windowNumber = 0;
	if(windowsPredefined){
		predefinedWindows->setChr(chrIterator->Name);
		numWindowsOnChr = predefinedWindows->getNumWindowsOnCurChr();
		int nextEnd = predefinedWindows->curWindowEnd();
		if(nextEnd > chrLength) nextEnd = chrLength + 1;
		else windowPair.nextPointer->move(chrIterator->Name, chrNumber, predefinedWindows->curWindowStart(), nextEnd);
	} else {
		numWindowsOnChr = ceil(chrLength / (double) windowSize);
		int nextEnd = windowSize;
		if(nextEnd > chrLength) nextEnd = chrLength; //!!! removed +1 because we are zero-based. Chedk if true!
		windowPair.nextPointer->move(chrIterator->Name, chrNumber, 0, nextEnd);
	}

	//advance mask
	if(doMasking || considerRegions) mask->setChr(chrIterator->Name);

	//write progress
	logfile->startNumbering("Parsing chromosome '" + chrIterator->Name + "':");
}

bool TGenome::iterateWindow(TWindowPair & windowPair){
	if(curEnd > 0) logfile->endIndent();

	//swap window pairs
	windowPair.swap();

	//move to next region
	curStart = windowPair.curPointer->start;
	curEnd = windowPair.curPointer->end;
	if(curStart >= chrLength || windowNumber >= limitWindows) return false;

	//move next
	if(windowsPredefined){
		if(numWindowsOnChr < 1){
			logfile->conclude("No windows on this chromosome.");
			return false;
		}
		//jump reader if large gap to previous window
		//TODO:: check if this does not mean we miss reads starting prior to the window but extending into it.
		if(windowPair.nextPointer->start - windowPair.curPointer->end > maxReadLength){
			if(curStart-maxReadLength < 0)
				bamReader.Jump(chrNumber, 0);
			else
				bamReader.Jump(chrNumber, curStart-maxReadLength);
		}

		//now move coordinates of next window
		if(predefinedWindows->nextWindow()){
			int nextEnd = predefinedWindows->curWindowEnd();
			if(nextEnd > chrLength) nextEnd = chrLength;
			windowPair.nextPointer->move(chrIterator->Name, chrNumber, predefinedWindows->curWindowStart(), nextEnd);
		} else {
			windowPair.nextPointer->move(chrIterator->Name, chrNumber, chrLength, chrLength+1);
		}
	} else {
		long nextEnd = curEnd + windowSize;
		if(nextEnd > chrLength) nextEnd = chrLength;
		windowPair.nextPointer->move(chrIterator->Name, chrNumber, curEnd, nextEnd);
	}

	++windowNumber;

	//report
	logfile->number("Window [" + toString(curStart) + ", " + toString(curEnd) + "] of " + toString(numWindowsOnChr) + " on '" + chrIterator->Name + "':");
	logfile->addIndent();

	return true;
};

bool TGenome::addAlignementToWindows(TAlignmentParser & alignment, TWindowPair & windowPair){
	//check if bam file is sorted
	if(alignment.position < oldPos)
		throw "BAM file must be sorted by position!";
	oldPos = alignment.position;

	//and add
	int readEnd = alignment.position + alignment.lastAlignedPos;
	if(alignment.position >= curStart || readEnd >= curStart){ //
		//check if still within current window and add to window
		if(alignment.position >= curEnd) return false;
		else {
			if(windowPair.addToCur(alignmentParser, pmdObjects)){
				//add also to next window in case reads overhangs current window -> function returns true
				windowPair.addToNext(alignmentParser, pmdObjects);
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
		if(!addAlignementToWindows(alignmentParser, windowPair)){
			//next read is for a later window
			oldAlignementMustBeConsidered = true;
			gettimeofday(&end, NULL);
			logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");
			logfile->conclude("No data in this window.");
			return false; //still only in next window
		}
	}

	while(alignmentParser.readAlignment(bamReader) && alignmentParser.chrNumber==chrNumber){
		if(alignmentParser.passedFilters && readGroups.readGroupInUse(alignmentParser.readGroupId)){
			//now parse alignment
			alignmentParser.parse();

			//and add to windows
			if(!addAlignementToWindows(alignmentParser, windowPair)){
				//read is beyond window and should be reconsidered
				oldAlignementMustBeConsidered = true;
				break;
			}
		}
	}

	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");

	if(windowPair.curPointer->numReadsInWindow > 0){
		//apply masks and filters
		if(doMasking){
			logfile->listFlush("Masking sites ...");
			windowPair.curPointer->applyMask(mask, considerRegions);
			logfile->done();
		} else if(considerRegions){
			logfile->listFlush("Masking sites outside regions ...");
			windowPair.curPointer->applyMask(mask, considerRegions);
			logfile->done();
		} else if(doCpGMasking){
			logfile->listFlush("Masking CpG sites ...");
			windowPair.curPointer->maskCpG(reference, chrNumber);
			logfile->done();
		} if(applyDepthFilter){
			windowPair.curPointer->applyDepthFilter(minDepth, maxDepth);
		} if(maxRefN < 1.0 && fastaReference == true){
			windowPair.curPointer->addReferenceBaseToSites(reference, chrNumber);
			windowPair.curPointer->calcFracN();
		}

		//calc sequencing depth
		windowPair.curPointer->calcDepth();

		//report
		logfile->conclude("read data from " + toString(windowPair.curPointer->numReadsInWindow) + " reads.");
		logfile->conclude("sequencing depth is " + toString(windowPair.curPointer->depth));
		logfile->conclude(toString(windowPair.curPointer->fractionsitesDepthAtLeastTwo * 100) + "% of all sites are covered at least twice");
		logfile->conclude(toString(windowPair.curPointer->fractionSitesNoData * 100) + "% of all sites have no data");
		if(windowPair.curPointer->fractionSitesNoData > maxMissing){
			logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			return false;
		}
		if(maxRefN < 1.0 && fastaReference == true){
			logfile->conclude(toString(windowPair.curPointer->fractionRefIsN * 100) + "% of all reference bases are 'N'");
			if(windowPair.curPointer->fractionRefIsN > maxRefN){
				logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
				return false;
			}
		}
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
		//all read groups have the same pmd
		logfile->list("Initializing one PMD function for all read groups.");
		pmdObjects[0].initialize(params, logfile);
		for(int i=1; i<readGroups.numGroups; ++i)
			pmdObjects[i].initialize(pmdObjects[0]);
		hasPMD = true;
	} else if(params.parameterExists("pmdFile")){
		//read from file for each read group
		std::string filename = params.getParameterString("pmdFile");
		logfile->list("Reading PMD from file '" + filename + "'.");
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
				//	logfile->conclude("For read group '" + vec[0] + "', C->T: " + pmdObjects[readGroupId].getFunctionString(pmdCT));
					pmdObjects[readGroupId].initializeFunction(vec[2], pmdGA);
				//	logfile->conclude("For read group '" + vec[0] + "', G->A: " + pmdObjects[readGroupId].getFunctionString(pmdGA));
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
		hasPMD = true;
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
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);
	if(params.parameterExists("recal")){
		std::string filename = params.getParameterString("recal");
		recalObject = new TRecalibrationEM(&bamHeader, filename, params, logfile, readGroupMap);
		doRecalibration = true;
	} else if(params.parameterExists("BQSRQuality")){
		recalObject = new TRecalibrationBQSR(&bamHeader, params, logfile, readGroupMap);
		doRecalibration = true;
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
		doRecalibration = false;
		recalObject = new TRecalibration(readGroupMap);
	}
	recalObjectInitialized = true;

	//check if estimation is required, in which case throw an error!
	if(recalObject->requiresEstimation()) throw "Can not use provided recalibration: estimation is required!";
}

void TGenome::initializeRandomGenerator(TParameters & params){
	logfile->listFlush("Initializing random generator ...");

	if(params.parameterExists("fixedSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else {
		randomGenerator=new TRandomGenerator();
	}
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
	randomGeneratorInitialized = true;
}

void TGenome::indexBamFile(std::string & filename){
	logfile->listFlush("Creating index of BAM file '" + filename + "' ...");
	BamTools::BamReader reader;
	if(!reader.Open(filename))
		throw "Failed to open BAM file '" + filename + "' for indexing!";

	// create index for BAM file
	reader.CreateIndex(BamTools::BamIndex::STANDARD);

	//close BAM file
	reader.Close();
	logfile->done();
}

//-----------------------------------------------------
//Functions for theta estimation
//-----------------------------------------------------
void TGenome::estimateTheta(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//Theta estimator
	TThetaEstimator thetaEstimator(params, logfile);

	//open output
	std::string filename = outputName + "_theta_estimates.txt.gz";
	TThetaOutputFile thetaOut(filename, &thetaEstimator, logfile);

	//check for which segements theta is to be estimated
	if(params.parameterExists("thetaGenomeWide") || considerRegions){
		if(params.parameterExists("thetaGenomeWide"))
			logfile->startIndent("Estimating theta genome-wide:");
		else logfile->startIndent("Estimating theta at specific sites:");

		//HACK!!
		bool onlyBootstrap = params.parameterExists("onlyBootstrap");

		int numBootstraps = 0;
		if(params.parameterExists("bootstraps")){
			int numBootstraps = params.getParameterInt("bootstraps");
			estimateThetaGenomeWide(thetaEstimator, thetaOut, onlyBootstrap, numBootstraps);
			bootstrapTetaEstimation(numBootstraps, thetaEstimator);
		} else
			estimateThetaGenomeWide(thetaEstimator, thetaOut, onlyBootstrap, numBootstraps);

		logfile->endIndent();
		if(params.parameterExists("bootstraps")){
			int numBootstraps = params.getParameterInt("bootstraps");
			bootstrapTetaEstimation(numBootstraps, thetaEstimator);
		}
	} else
		estimateThetaWindows(thetaEstimator, thetaOut);

	//clean up
	thetaOut.close();
}

void TGenome::estimateThetaWindows(TThetaEstimator & thetaEstimator, TThetaOutputFile & out){
	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			if(readData(windows)){
				if(windows.cur->fractionSitesNoData > maxMissing){
					logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
				} if(windows.cur->fractionRefIsN > maxRefN){
					logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
				} else {
					logfile->startIndent("Estimating Theta:");

					//measure runtime
					struct timeval startTime, endTime;
					gettimeofday(&startTime, NULL);

					//adding sites to estimator
					logfile->listFlush("Calculating emission probabilities ...");
					thetaEstimator.clear();
					windows.cur->addSitesToThetaEstimator(recalObject, thetaEstimator.pointerToDataContainer());
					logfile->done();

					//estimate Theta
					if(thetaEstimator.estimateTheta())
						out.writeWindow(chrIterator->Name, windows.cur->start, windows.cur->end, &thetaEstimator);

					//clear theta estimator
					thetaEstimator.clear();

					//finish
					gettimeofday(&endTime, NULL);
					logfile->list("Total computation time for this window was ", endTime.tv_sec  - startTime.tv_sec, "s");
					logfile->endIndent();

				}
			} else logfile->list("No relevant positions -> skipping this window.");
		}
	}
}

void TGenome::estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool onlyReadData, int numBootstraps){
	if(considerRegions)
		logfile->startIndent("Estimating theta at specific sites:");

	//prepare windows
	TWindowPairDiploid windows;
	thetaEstimator.clear();

	//add sites to estimator
	logfile->startIndent("Adding sites to data structure:");
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			if(readData(windows)){
				//adding sites to estimator
				logfile->listFlush("Calculating emission probabilities ...");
				try{
					windows.cur->addSitesToThetaEstimator(recalObject, thetaEstimator.pointerToDataContainer());
				} catch(...){
					throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size, selecting fewer regions or limiting to sites with a minimal depth (>=2 recommended).";
				}
				logfile->done();
			}
		}
	}
	logfile->endIndent();

	//estimate Theta
	//HACK!!!!
	if(!onlyReadData){
		logfile->startIndent("Estimate theta based on a total of " + toString(thetaEstimator.sizeWithData()) + " sites:");
		thetaEstimator.estimateTheta();
	}

	if(considerRegions)
		out.writeWindow("regions", "-", "-", &thetaEstimator);
	else
		out.writeWindow("genome-wide", "-", "-", &thetaEstimator);
	if(numBootstraps == 0) thetaEstimator.clear();
}

void TGenome::bootstrapTetaEstimation(int numBootstraps, TThetaEstimator & thetaEstimator){
	if(numBootstraps < 1) throw "Number of bootstraps must be > 1!";
	logfile->startIndent("Generating " + toString(numBootstraps) + " bootstrap estimates of theta:");

	//measure runtime
	struct timeval startTime, endTime;
	gettimeofday(&startTime, NULL);

	//open output file
	gz::ogzstream bootstrapOut;
	std::string bootstrapFilename = outputName + "_theta_bootstraps.txt.gz";
	logfile->list("Writing theta bootstraps to '" + bootstrapFilename + "'");
	bootstrapOut.open(bootstrapFilename.c_str());
	if(!bootstrapOut) throw "Failed to open output file '" + bootstrapFilename + "'!";

	//write header
	bootstrapOut << std::setprecision(9) << "Bootstrap";
	thetaEstimator.writeHeader(bootstrapOut);
	bootstrapOut << "\n";

	//loop over bootstraps
	for(int s=0; s<numBootstraps; ++s){
		logfile->startIndent("Bootstrap " + toString(s+1) + " of " + toString(numBootstraps) + ":");

		//run bootstrap
		bootstrapOut << s+1;
		thetaEstimator.bootstrapTheta(*randomGenerator, bootstrapOut);

		logfile->endIndent();
	}

	//finish
	gettimeofday(&endTime, NULL);
	logfile->list("Total computation time for theta bootstrapping was ", round((endTime.tv_sec  - startTime.tv_sec) / 6.0)/10.0, " min");
	logfile->endIndent();
}

void TGenome::calcThetaLikelihoodSurfaces(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//read params
	int steps = params.getParameterIntWithDefault("steps", 100);

	//prepare windows
	TWindowPairDiploid windows;

	//Theta estimator
	TThetaEstimator estimator(logfile);

	//iterate through windows
	std::string filename;
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				//check if we have data -> can be extended to ensure
				logfile->startIndent("Calculating likelihood surface for Theta:");

				//measure runtime
				struct timeval startTime, endTime;
				gettimeofday(&startTime, NULL);

				//adding sites to estimator
				logfile->listFlush("Calculating emission probabilities ...");
				windows.cur->addSitesToThetaEstimator(recalObject, estimator.pointerToDataContainer());
				logfile->done();

				//open file
				gz::ogzstream out;
				filename = outputName + chrIterator->Name + "_" + toString(windows.cur->start) + "_LLsurface.txt.gz";
				out.open(filename.c_str());
				if(!out) throw "Failed to open output file '" + outputName + "'!";


				//estimate Theta
				logfile->listFlush("Calculating likelihood surface ...");
				estimator.calcLikelihoodSurface(out, steps);
				logfile->done();

				//clear theta estimator
				estimator.clear();

				//finish
				out.close();
				gettimeofday(&endTime, NULL);
				logfile->list("Total computation time for this window was ", endTime.tv_sec  - startTime.tv_sec, "s");
				logfile->endIndent();
			}
		}
	}
};

void TGenome::estimateThetaRatio(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//Theta estimator
	TThetaEstimatorRatio thetaEstimatorRatio(params, logfile);

	//read the two regions to be used
	logfile->startIndent("Reading regions:");

	//region 1
	std::string regionsFile1 = params.getParameterString("regions1");
	logfile->listFlush("Reading regions 1 from BED file '" + regionsFile1 + "' ...");
	TBedReader region1(regionsFile1, windowSize, bamHeader.Sequences, logfile);
	logfile->done();

	//region 2
	std::string regionsFile2 = params.getParameterString("regions2");
	logfile->listFlush("Reading regions 2 from BED file '" + regionsFile2 + "' ...");
	TBedReader region2(regionsFile2, windowSize, bamHeader.Sequences, logfile);
	logfile->done();
	logfile->endIndent();

	//prepare windows
	TWindowPairDiploid windows;

	//add sites to estimator
	logfile->startIndent("Adding sites to data structures:");
	while(iterateChromosome(windows)){
		region1.setChr(chrIterator->Name);
		region2.setChr(chrIterator->Name);

		while(iterateWindow(windows)){
			if(readData(windows)){
				//adding sites to estimator
				logfile->listFlush("Calculating emission probabilities ...");
				windows.cur->calculateEmissionProbabilities(recalObject);
				try{
					windows.cur->addSitesToThetaEstimator(thetaEstimatorRatio.pointerToDataContainer(), region1);
					windows.cur->addSitesToThetaEstimator(thetaEstimatorRatio.pointerToDataContainer2(), region2);
				} catch(...){
					throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size, selecting fewer regions or limiting to sites with a minimal depth (>=2 recommended).";
				}
				logfile->done();
			}
		}
	}
	logfile->endIndent();

	//estimate Theta ratio
	thetaEstimatorRatio.estimateRatio(*randomGenerator, outputName);

	//clean up
}


//------------------------------------------
//Callers (NEW)
//------------------------------------------
TGenotypePrior* TGenome::initializeGenotypePrior(TParameters & params){
	TGenotypePrior* prior;
	logfile->startIndent("Initializing genotype prior:");
	//read prior from parameters
	std::string priorMethod = params.getParameterStringWithDefault("prior", "theta");
	if(priorMethod == "unif"){
		prior = new TGenotypePriorUniform();
		logfile->list("Will use a uniform prior with equal weights for all genotypes.");
	} else if(priorMethod == "theta"){
		if(params.parameterExists("fixedTheta")){
			double theta = params.getParameterDouble("fixedTheta");
			logfile->list("Will use a fixed theta = " + toString(theta));
			bool equalBaseFreq = params.parameterExists("equalBaseFreq");
			if(equalBaseFreq)
				logfile->list("Will use equal base frequencies.");
			else
				logfile->list("Will estimate base frequencies individually for each window.");
			prior = new TGenotypePriorFixedTheta(theta, equalBaseFreq, logfile);
		} else {
			logfile->list("Will use a prior based on theta and base frequencies estimated individually for each window.");
			std::string thetaOuputName = outputName + "_theta_estimates.txt.gz";
			if(params.parameterExists("defaultTheta")){
				double defaultTheta = params.getParameterDouble("defaultTheta");
				logfile->list("Will use a default theta of ", defaultTheta, " for windows with limited data.");
				prior = new TGenotypePriorTheta(params, thetaOuputName, defaultTheta, logfile);
			} else
				prior = new TGenotypePriorTheta(params, thetaOuputName, logfile);
		}
	} else throw "Unknown prior type '" + priorMethod + "'!";
	logfile->endIndent();

	return prior;
}

void TGenome::callGenotypesNew(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//initialize recalibration
	initializeRecalibration(params);

	//--------------------------
	//initialize caller
	//--------------------------
	logfile->startIndent("Initializing caller:");
	TCaller* caller;
	std::string method = params.getParameterStringWithDefault("method", "MLE");
	if(method == "randomBase"){
		caller = new TCallerRandomBase(randomGenerator);
	} else if(method == "majorityBase"){
		caller = new TCallerMajorityBase(randomGenerator);
	} else if(method == "allelePresence"){
		caller = new TCallerAllelePresence(randomGenerator);
	} else if(method == "MLE"){
		caller = new TCallerMLE(randomGenerator);
	} else if(method == "Bayesian"){
		caller = new TCallerBayes(randomGenerator);
	} else if(method == "gVCF"){
		throw "GVCF NOT YET IMPLEMENTED!";
		caller->printSitesWithNoData();
	} else throw "Unknown calling method '" + method + "'! Use randomBase, allelePresence, MLE, Bayesian or gVCF.";

	//read output settings
	if(params.parameterExists("infoFields"))
		caller->printInfoFields(params.getParameterString("infoFields"));
	if(params.parameterExists("formatFields"))
		caller->printGenotypeFields(params.getParameterString("formatFields"));
	if(params.parameterExists("printAll")) caller->setPrintSitesWithNoData(true);
	if(params.parameterExists("noAltIfHomoRef")) caller->setPrintAltIfHomoRef(false);
	if(params.parameterExists("noTriallelic")) caller->setAllowTriallelic(false);

	//report output settings
	caller->reportSettings(logfile);

	//prior setting
	TGenotypePrior* prior;
	if(caller->usesPrior()){
		prior = initializeGenotypePrior(params);
		caller->setPrior(prior->getPointerToPrior());
	} else prior = new TGenotypePrior();

	//open output file
	std::string sampleName = params.getParameterStringWithDefault("indName", outputName);
	caller->openVCF(outputName, sampleName);
	logfile->endIndent();

	//prepare windows
	//Allow for haploid windows for some callers?
	TWindowPairDiploid windows;

	//---------------------------------------------------------------------
	// Now call, either all sites or limiting to sites with known alleles.
	//---------------------------------------------------------------------
	if(params.parameterExists("sites")){
		//Limit to sites with known alleles
		logfile->startIndent("Will limit calls to sites with known alleles:");
		TSiteSubset subset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, false);
		logfile->endIndent();

		while(iterateChromosome(windows)){
			while(iterateWindow(windows)){
				//read data for current window
				if(readData(windows) || caller->printSitesWithNoData()){
					//update genotype prior
					prior->update(windows.cur, logfile);

					//now call using known alleles
					logfile->listFlush("Calling genotypes ...");
					windows.cur->callKnwonAlleles(*caller, *recalObject, reference, subset);
					logfile->done();
				}
			}
		}
	} else { //not limiting to sites with known alleles
		//Use all sites and identify alleles
		while(iterateChromosome(windows)){
			while(iterateWindow(windows)){
				//read data for current window
				if(readData(windows) || caller->printSitesWithNoData()){
					//update genotype prior
					prior->update(windows.cur, logfile);

					//now call
					logfile->listFlush("Calling genotypes ...");
					windows.cur->call(*caller, *recalObject, reference);
					logfile->done();
				}
			}
		}
	}

	//clean up
	delete caller;
	delete prior;
}


//------------------------------------------
//Callers
//------------------------------------------
void TGenome::callMLEGenotypes(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//set all booleans in the booleans class
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = false;
	bool noAltIfHomoRef = false;

	bool gVCF = false;
	bool writeVCF = false;
	bool beagle = false, printOnlyGL = false;


	std::string indName;
	TSiteSubset* subset = NULL;

	//only call at specific sites?
	if(params.parameterExists("sites")){
		bool invariantSites = false;

		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile, invariantSites);
		limitToSitesWithKnownAlleles = true;

	//if not, how much information should be printed?
	} else {
		if(params.parameterExists("printAll")){
			printIfNoData = true;
			logfile->list("Will print all sites, even those without data");
		}
		if(params.parameterExists("noAltIfHomoRef")){
			noAltIfHomoRef = true;
			logfile->list("Will not print alternative alleles when genotype is 0/0");
		}
		if(params.parameterExists("gVCF")){
			gVCF = true;
			if(!printIfNoData) throw "gVCF format includes calls for all sites. Use parameter \"printAll\".";
			if(noAltIfHomoRef) throw "gVCF format includes printing alternative alleles even if genotype is 0/0. Remove \"printIfNoData\".";
			if(!fastaReference) throw "Can not print VCF file without reference!";
			logfile->list("Will print output in gVCF format");
		}
	}

	if(params.parameterExists("beagle")){
		if(limitToSitesWithKnownAlleles == false) throw "Need sites file specifying major and minor alleles for beagle format!";
		beagle=true;
		logfile->list("Will print output in beagle format");
		printOnlyGL = params.parameterExists("printOnlyGL");
		if(printOnlyGL) logfile->list("Will print only genotype likelihoods");
		indName = params.getParameterStringWithDefault("indName", outputName);
	}

	if(params.parameterExists("vcf")){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;
	}

	if((writeVCF + gVCF + beagle) > 1) throw "More than one output format specified!";

	//open output file
	gz::ogzstream out;
	std::string outputFileName;

//	boolsForVCF boolaeans(writeVCF, noAltIfHomoRef, gVCF);



	if(params.parameterExists("vcf") || gVCF){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;
		std::string sName = params.getParameterStringWithDefault("sampleName", outputName);


		//open file
		if(gVCF) outputFileName = outputName + "_MLEGenotypes.gvcf.gz";
		else outputFileName = outputName + "_MLEGenotypes.vcf.gz";
		if(gVCF) logfile->list("Writing MLE genotypes in gVCF format to '" + outputFileName + "'");
		else logfile->list("Writing MLE genotypes in VCF format to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		bool onlyPhredGP = false;
		writeVcfHeader(&out, limitToSitesWithKnownAlleles, onlyPhredGP);

	} else if(beagle){
		//open file
		outputFileName = outputName + "_MLEGenotypes.beagle.gz";
		logfile->list("Writing MLE genotypes to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		if(printOnlyGL) out << indName << "\t" << indName << "\t" << indName << "\n";
		else out << "#marker\talleleA\talleleB\t" << indName << "\t" << indName << "\t" << indName << "\n";

	} else {
		//open file
		outputFileName = outputName + "_MLEGenotypes.txt.gz";
		logfile->list("Writing MLE genotypes to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		out << "chr\tpos\tRef";
		if(limitToSitesWithKnownAlleles) out << "\tAlt\tdepth\tL(RR)\tL(RA)\tL(AA)\tMLE\tQ\n";
		else out << "\tdepth\tL(AA)\tL(AC)\tL(AG)\tL(AT)\tL(CC)\tL(CG)\tL(CT)\tL(GG)\tL(GT)\tL(TT)\tMLE\tQ\n";
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles || beagle) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			if((!limitToSitesWithKnownAlleles && !beagle) || subset->hasPositionsInWindow(windows.cur->start)){
				//read data for current window
				if(readData(windows) || printIfNoData){
					//call genotypes
					logfile->listFlush("Calling MLE genotypes ...");
					if(limitToSitesWithKnownAlleles || beagle){
						windows.cur->addReferenceBaseToSites(subset);
						windows.cur->callMLEGenotypeKnownAlleles(recalObject, subset, *randomGenerator, out, chrIterator->Name, writeVCF, noAltIfHomoRef, beagle, printOnlyGL);
					} else {
						if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
						windows.cur->callMLEGenotype(recalObject, *randomGenerator, out, chrIterator->Name, printIfNoData, fastaReference, writeVCF, gVCF, noAltIfHomoRef);
					}
					logfile->done();
				}
			}
		}
	}
	//clean up
	out.close();
}

bool TGenome::initThetaEstimatorForCallers(TParameters & params, TThetaEstimator* & thetaEstimator){
	if(params.parameterExists("theta")){
		double theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
		thetaEstimator = new TThetaEstimator(logfile);
		thetaEstimator->setTheta(theta);
		return false;
	} else {
		//prepare theta estimator
		thetaEstimator = new TThetaEstimator(params, logfile);
		return true;
	}
}

void TGenome::writeVcfHeader(gz::ogzstream* output, bool limitToSitesWithKnownAlleles, bool onlyPhredGP){
	//write header
	(*output) << "##fileformat=VCFv4.2\n";
	(*output) << "##source=atlas\n";
	(*output) << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	if(!limitToSitesWithKnownAlleles)
		(*output) << "##INFO=<ID=PP,Number=10,Type=Integer,Description=\"Phred-scaled posterior probabilities of all genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT and TT\">\n";
	(*output) << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	if(onlyPhredGP)
		(*output) << "##FORMAT=<ID=GP,Number=G,Type=Integer,Description=\"Genotype posterior probabilities (phred-scaled)\">\n";
	else
		(*output) << "##FORMAT=<ID=GP,Number=G,Type=Integer,Description=\"Genotype posterior probabilities round(phred(1-Posterior Prob))\">\n";
	(*output) << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	(*output) << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
	(*output) << "##FORMAT=<ID=GG,Number=10,Type=Integer,Description=\"All phred-scaled normalized genotype likelihoods\">\n";
	(*output) << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled normalized genotype likelihoods\">\n";
	(*output) << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth\">\n";
	(*output) << "##FORMAT=<ID=AB,Number=1,Type=Float,Description=\"Allelic imbalance\">\n";
	(*output) << "##FORMAT=<ID=AI,Number=1,Type=Float,Description=\"Binomial probability of allelic imbalance if Hz site\">\n";
	(*output) << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
}

void TGenome::callBayesianGenotypes(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	TThetaEstimator* thetaEstimator;
	TThetaOutputFile* thetaOut = NULL;
	bool estimateTheta = initThetaEstimatorForCallers(params, thetaEstimator);
	if(estimateTheta){
		std::string filename = outputName + "_theta_estimates.txt.gz";
		thetaOut->open(filename, thetaEstimator, logfile);
	}

	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool noAltIfHomoRef = false;
	bool printIfNoData = true;
	bool printPP = false;
	bool onlyPhredGP = false;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(windowsPredefined) throw "Using site subsets is currently not implemented if windows are predefined from a BED file.";
		bool invariantSites = false;
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile, invariantSites);
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
		if(params.parameterExists("noAltIfHomoRef")){
			noAltIfHomoRef = true;
			logfile->list("Will not print alternative alleles when genotype is 0/0");
		}
		if(params.parameterExists("printPP")){
			printPP = true;
			logfile->list("Will print the phred scaled posterior probabilities for all genotypes");
		}
		if(params.parameterExists("onlyPhredGP")){
			onlyPhredGP = true;
			logfile->list("Will print GP in phred format");
		} else {
			logfile->list("Will print GP as phred(1-posterior probability)");
		}
		//open file
		outputFileName = outputName + "_BayesianGenotypes.vcf.gz";
		logfile->list("Writing Bayesian genotypes in VCF format to '" + outputFileName + "'");
		output.open(outputFileName.c_str());
		if(!output) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		writeVcfHeader(&output, limitToSitesWithKnownAlleles, onlyPhredGP);

	} else {
		//open file
		outputFileName = outputName + "_BayesianGenotypes.txt.gz";
		logfile->list("Writing Bayesian genotypes to '" + outputFileName + "'");
		output.open(outputFileName.c_str());
		if(!output) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		output << "chr\tpos";
		if(fastaReference) output << "\tRef";
		if(limitToSitesWithKnownAlleles) output << "\talt\tdepth\tP(RR|D)\tP(RA|D)\tP(AA|D)\tMAP\tQ\n";
		else output << "\tdepth\tP(AA|D)\tP(AC|D)\tP(AG|D)\tP(AT|D)\tP(CC|D)\tP(CG|D)\tP(CT|D)\tP(GG|D)\tP(GT|D)\tP(TT|D)\tMAP\tQ\n";
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			//read data for current window
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				if(readData(windows) || printIfNoData){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing && estimateTheta){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN && estimateTheta){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//estimate Theta?
						if(estimateTheta){
							//adding sites to estimator
							logfile->listFlush("Calculating emission probabilities ...");
							windows.cur->addSitesToThetaEstimator(recalObject, thetaEstimator->pointerToDataContainer());
							logfile->done();

							//estimate Theta
							thetaEstimator->estimateTheta();

							//write results to file
							thetaOut->writeWindow(chrIterator->Name, windows.cur->start, windows.cur->end, thetaEstimator);

							//clear theta estimator
							(*thetaEstimator).clear();

						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							thetaEstimator->setBaseFreq(windows.cur->baseFreq);
						}

						//call Bayesian genotypes
						logfile->listFlush("Calling Bayesian Genotypes ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callBayesianGenotypeKnownAlleles(subset, *thetaEstimator, *randomGenerator, output, chrIterator->Name, writeVCF);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callBayesianGenotype(*thetaEstimator, *randomGenerator, output, chrIterator->Name, printIfNoData, fastaReference, writeVCF, noAltIfHomoRef, printPP, onlyPhredGP);
						}
						logfile->done();
					}
				}
			}
		}
	}

	//clean up
	if(estimateTheta){
		delete thetaEstimator;
		delete thetaOut;
	}
	if(limitToSitesWithKnownAlleles) delete subset;
}

void TGenome::callAllelePresence(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	TThetaEstimator* thetaEstimator;
	TThetaOutputFile* thetaOut = NULL;
	bool estimateTheta = initThetaEstimatorForCallers(params, thetaEstimator);
	if(estimateTheta){
		std::string filename = outputName + "_theta_estimates.txt.gz";
		thetaOut->open(filename, thetaEstimator, logfile);
	}


	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = false;
	bool noAltIfHomoRef = false;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(windowsPredefined) throw "Using site subsets is currently not implemented if windows are predefined from a BED file.";
		bool invariantSites = false;
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile, invariantSites);
		limitToSitesWithKnownAlleles = true;
	} if(params.parameterExists("noAltIfHomoRef")){
		noAltIfHomoRef = true;
		logfile->list("Will not print alternative alleles when genotype is 0/0");
	} else {
		if(params.parameterExists("printAll")){
			printIfNoData = true;
			logfile->list("Will print all sites, even those without data");
		}
		if(params.parameterExists("noAltIfHomoRef")){
			noAltIfHomoRef = true;
			logfile->list("Will not print alternative alleles when genotype is 0/0");
		}
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
		outAllelePresence << "##source=ATLAS\n";
		outAllelePresence << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		outAllelePresence << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
		outAllelePresence << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		outAllelePresence << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		outAllelePresence << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
		if (limitToSitesWithKnownAlleles) outAllelePresence << "##FORMAT=<ID=PP,Number=2,Type=Integer,Description=\"Phred-scaled posterior probabilities of allele presence in the order A, C, G and T\">\n";
		else outAllelePresence << "##FORMAT=<ID=PP,Number=4,Type=Integer,Description=\"Phred-scaled posterior probabilities of allele presence in the order A, C, G and T\">\n";
		outAllelePresence << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		outputFileName = outputName + "_AllelePresence.txt.gz";
		logfile->list("Writing estimates of allele presence to '" + outputFileName + "'");
		outAllelePresence.open(outputFileName.c_str());
		if(!outAllelePresence) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		outAllelePresence << "chr\tpos";
		if(limitToSitesWithKnownAlleles) outAllelePresence << "\tRef\tAlt\tdepth\tP(Ref|D)\tP(Alt|D)\tMAP\tQ\n";
		else {
			if(fastaReference) outAllelePresence << "\tRef";
			outAllelePresence << "\tdepth\tP(A|D)\tP(C|D)\tP(G|D)\tP(T|D)\tMAP\tQ\n";
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
				//TODO: ask Dan if following line should actually be if(readData || printIfNoData)
				if(readData(windows) || printIfNoData){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing && estimateTheta){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN && estimateTheta){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//estimate Theta?
						if(estimateTheta){
							//adding sites to estimator
							logfile->listFlush("Calculating emission probabilities ...");
							(*thetaEstimator).clear();
							windows.cur->addSitesToThetaEstimator(recalObject, thetaEstimator->pointerToDataContainer());
							logfile->done();

							//estimate Theta
							thetaEstimator->estimateTheta();

							//write results to file
							thetaOut->writeWindow(chrIterator->Name, windows.cur->start, windows.cur->end, thetaEstimator);

							//clear theta estimator
							(*thetaEstimator).clear();

						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							thetaEstimator->setBaseFreq(windows.cur->baseFreq);
						}

						//call allele presence
						logfile->listFlush("Calling allele presence ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callAllelePresenceKnwonAlleles(subset, *thetaEstimator, *randomGenerator, outAllelePresence, chrIterator->Name, writeVCF, noAltIfHomoRef);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callAllelePresence(*thetaEstimator, *randomGenerator, outAllelePresence, chrIterator->Name, printIfNoData, fastaReference, writeVCF, noAltIfHomoRef);
						}
						logfile->done();
					}
				}
			} else logfile->list("No positions in this window.");
		}
	}

	//clean up
	if(estimateTheta){
		delete thetaOut;
		delete thetaEstimator;
	}
	if(limitToSitesWithKnownAlleles) delete subset;
}

void TGenome::randomBaseCaller(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	bool printIfNoData = false;
	if(params.parameterExists("printAll")){
		printIfNoData = true;
		logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	gz::ogzstream randomBases;
	std::string outputFileName;

	//open file
	outputFileName = outputName + "_randomBase.txt.gz";
	logfile->list("Writing random base calls to '" + outputFileName + "'");
	randomBases.open(outputFileName.c_str());
	if(!randomBases) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	randomBases << "chr\tpos\tref\tdepth\tpileup\trandom_base\n";

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows) || printIfNoData){
				//call random allele
				logfile->listFlush("Calling random base ...");
				if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
				windows.cur->callRandomBase(*randomGenerator, randomBases, chrIterator->Name, printIfNoData);
				logfile->done();

			} else logfile->list("No positions in this window.");
		}
	}
}

void TGenome::majorityBaseCaller(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	bool printIfNoData = false;
	if(params.parameterExists("printAll")){
		printIfNoData = true;
		logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	gz::ogzstream randomBases;
	std::string outputFileName;

	//open file
	outputFileName = outputName + "_majorityCalls.txt.gz";
	logfile->list("Writing random base calls to '" + outputFileName + "'");
	randomBases.open(outputFileName.c_str());
	if(!randomBases) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	randomBases << "chr\tpos\tref\tdepth\tpileup\trandom_base\n";

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows) || printIfNoData){
				//call random allele
				logfile->listFlush("Calling random base ...");
				if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
				windows.cur->majorityCall(*randomGenerator, randomBases, chrIterator->Name, printIfNoData);
				logfile->done();

			} else logfile->list("No positions in this window.");
		}
	}
}


//---------------------------------------------------
// I/O
//---------------------------------------------------
void TGenome::writeGLF(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//print all?
	bool printIfNoData = params.parameterExists("printAll");
	if(printIfNoData)
		logfile->list("Will print all sites, even those without data");

	//open GLF file
	std::string outputFileName = outputName + ".glf.gz";
	logfile->list("Will write data to GLF file '" + outputFileName + "'");
	TGlfWriter writer(outputFileName);

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		writer.newChromosome(chrIterator->Name, stringToLong(chrIterator->Length));
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				//write to GLF
				logfile->listFlush("Adding window to GLF file ...");
				windows.cur->calculateEmissionProbabilities(recalObject);
				windows.cur->addToGLF(writer, printIfNoData);
				logfile->done();
			}
		}
	}
	//clean up
	writer.close();
}

void TGenome::printPileup(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//prepare windows
	TWindowPairDiploid windows;

	bool printOnlySitesWithData = false;
	if(params.parameterExists("printOnlySitesWithData")){
		printOnlySitesWithData = true;
		logfile->list("Will print only sites with data");
	}

	//open output
	gz::ogzstream out;
	std::string filename = outputName + "_pileup.txt.gz";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	TGenotypeMap genoMap;
	out << "Chr\tposition\tref\tdepth\trefDepth\tbases";
	for(int i=0; i<10; ++i)
		out << "\tPem(" << genoMap.getGenotypeString(i) << ")";
	out << "\n";

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			if(readData(windows)){
				if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
				windows.cur->printPileup(recalObject, out, chrIterator->Name, printOnlySitesWithData);
			}
		}
	}

	//clean up
	out.close();
}


void TGenome::generatePSMCInput(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//read in parameters required
	double theta = params.getParameterDoubleWithDefault("theta", 0.001);
	logfile->list("Using theta = " + toString(theta));
	TThetaEstimator thetaEstimator(logfile);
	thetaEstimator.setTheta(theta);

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

			//set base frequencies
			logfile->listFlush("Calculating emission probabilities ...");
			windows.cur->calculateEmissionProbabilities(recalObject);
			windows.cur->estimateBaseFrequencies();
			thetaEstimator.setBaseFreq(windows.cur->baseFreq);
			logfile->done();

			//create PSMC input
			logfile->listFlush("Estimating heterozygosity status ...");
			if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
			windows.cur->generatePSMCInput(thetaEstimator, blockSize, confidence, output, nCharOnLine);
			logfile->done();
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
}

void TGenome::createDepthMask(TParameters & params){
	int minDepthForMask = params.getParameterInt("minDepthForMask");
	int maxDepthForMask = params.getParameterInt("maxDepthForMask");
	if(params.parameterExists("maxDepth") || params.parameterExists("minDepth")) throw "Cannot mask sites for sequencing depth while creating the mask!";

	std::ofstream output;
	std::string outputFileName = outputName + "_minDepth"+ toString(minDepthForMask) + "_maxDepth" + toString(maxDepthForMask) + "_depthMask.bed";
	logfile->list("Writing sites with depth < " + toString(minDepthForMask) + " and with depth > " + toString(maxDepthForMask) + " to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			logfile->listFlush("Writing sites to mask to output file ...");
			windows.cur->createDepthMask(minDepthForMask, maxDepthForMask, output, chrIterator->Name);
			logfile->done();
		}
	}
	output.close();
}

//---------------------------------------------------
//recalibration
//---------------------------------------------------
void TGenome::estimateErrorCalibrationEM(TParameters & params){
	//create recalibration object
	std::string filename;
	if(params.parameterExists("recal")) filename = params.getParameterString("recal");
	else filename = "empty";

	bool writeTmpTables = false;
	if(params.parameterExists("writeTmpTables")){
		writeTmpTables = true;
		logfile->list("Will write intermediate estimates of EM and Newton-Raphson to file.");
	}
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);
	TRecalibrationEM recalObjectEM(&bamHeader, filename, params, logfile, readGroupMap);
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
			if(readData(windows)) windows.cur->addToRecalibrationEM(recalObjectEM);
			else logfile->list("No positions in this window.");
		}
	}
	//clean up memory
	windows.clear();
	logfile->endIndent();

	//run EM iterations
	recalObjectEM.runEM(outputName, writeTmpTables);
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

void TGenome::calculateLikelihoodErrorCalibrationEM(TParameters & params){
	//create recalibration object
	std::string filename = params.getParameterString("recal");
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);
	TRecalibrationEM recalObjectEM(&bamHeader, filename, params, logfile, readGroupMap);

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
	//int numMarginalGridPoint = params.getParameterIntWithDefault("numGridPoints", 51);
	//recalObjectEM.calcLikelihoodSurface(outputName + "_LLsurface.txt", numMarginalGridPoint);

	logfile->list("LL = " + toString(recalObjectEM.calcLL()));
}

void TGenome::BQSR(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not run BQSR recalibration without a provided FASTA reference!";

	//prepare windows
	TWindowPairHaploid windows;

	//create BQSR object
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);
	TRecalibrationBQSR bqsr(&bamHeader, params, logfile, readGroupMap);
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

	//do we consider only specific sites?
	bool invariantSites = false;
	TSiteSubset* subset = NULL;

	if(params.parameterExists("sites")){
		subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		invariantSites = true;
	}

	//loop over bam until BQSR converges
	while(!hasConverged){
		++loopNumber;
		logfile->startIndent("Running recalibration loop " + toString(loopNumber) + ":");
		//loop over all windows

		if(!bqsr.dataHasBeenStored()){
			logfile->startIndent("Reading data from BAM file:");
			while(iterateChromosome(windows)){
				if(invariantSites) subset->setChr(chrIterator->Name);
				while(iterateWindow(windows)){
					if((invariantSites && subset->hasPositionsInWindow(windows.cur->start)) || !invariantSites){
						//read data for current window
						if(readData(windows)){
							//add reference data
							windows.cur->addReferenceBaseToSites(reference, chrNumber);

							//add the base to BQSR
							if(invariantSites) windows.cur->addSitesToBQSR(bqsr, subset, logfile);
							else windows.cur->addSitesToBQSR(bqsr, logfile);
						}
					} else logfile->list("No positions in this window.");
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
};

void TGenome::printQualityDistribution(TParameters & params){
	//Assemble quality distribution
	int maxQinPrintQualityDistribution = maxPhredInt;
//			params.getParameterIntWithDefault("maxQ", 100);

	maxQinPrintQualityDistribution += 33; //internally, quality is phred(error) + 33!
	logfile->list("Will assemble quality distribution up to a quality of " + toString(maxQinPrintQualityDistribution-33) + " (" + (char) maxQinPrintQualityDistribution + ").");

	//initialize tables: one overall, one per read group
	std::vector<TQualityTable> qualDist;
	for(int i=0; i<readGroups.numGroups; ++i)
		qualDist.emplace_back(maxQinPrintQualityDistribution);

	//vars
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
    long counter = 0;

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader)){
        if(useChromosome[alignmentParser.chrNumber] && alignmentParser.passedFilters && readGroups.readGroupInUse(alignmentParser.readGroupId)){
			++counter;

			//update and write (only if alignment qualities could be calculated)
			alignmentParser.addToQualityTable(qualDist[alignmentParser.readGroupId]);

			//report
			reportProgressParsingBamFile(counter, start);
        }
	}

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

	//print per read group table
	logfile->startIndent("Writing distributions:");
	std::string outFileName;
	for(int i=0; i<readGroups.numGroups; ++i){
		//open output file
		outFileName = outputName + "_" + readGroups.getName(i) + "_qualityDistribution.txt";
		logfile->listFlush("Writing distribution for read group '" + readGroups.getName(i) + "' to '" + outFileName + "' ...");
		qualDist[i].write(outFileName);
		logfile->done();
	}

	//print overall table
	outFileName = outputName + "_total_qualityDistribution.txt";
	logfile->listFlush("Writing total distribution to '" + outFileName + " ...");
	TQualityTable allQualDist(maxQinPrintQualityDistribution);
	for(int i=0; i<readGroups.numGroups; ++i)
		allQualDist.add(qualDist[i]);
	allQualDist.write(outFileName);
	logfile->done();
	logfile->endIndent();
}

void TGenome::printQualityTransformation(TParameters & params){
	//TODO: use TAlignmentParser. No need to use windows!
	//initialize recalibration
	//compare to a second recalibration definition?
	if(params.parameterExists("recal2") && !params.parameterExists("recal")) throw "use recal instead of recal2 for comparison of recalibrated qualities to original qualities!";
	recalObjectInitialized2 = false;
	if(params.parameterExists("recal")){
		std::string nameRecal = params.getParameterString("recal");
		TReadGroupMap readGroupMap(&bamHeader, params, logfile);
		recalObject = new TRecalibrationEM(&bamHeader, nameRecal, params, logfile, readGroupMap);
		if(params.parameterExists("recal2")){
			std::string nameRecal2 = params.getParameterString("recal2");
			recalObject2 = new TRecalibrationEM(&bamHeader, nameRecal2, params, logfile, readGroupMap);
			doRecalibration2 = true;
			recalObjectInitialized2 = true;
		} else if(params.parameterExists("BQSRQuality")){
			TReadGroupMap readGroupMap(&bamHeader, params, logfile);
			recalObject2 = new TRecalibrationBQSR(&bamHeader, params, logfile, readGroupMap);
			doRecalibration2 = true;
			recalObjectInitialized2 = true;
		}
	} else if(params.parameterExists("BQSRQuality")){
		TReadGroupMap readGroupMap(&bamHeader, params, logfile);
		recalObject = new TRecalibrationBQSR(&bamHeader, params, logfile, readGroupMap);
		doRecalibration = true;
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
		doRecalibration = false;
		TReadGroupMap readGroupMap(&bamHeader, params, logfile);
		recalObject = new TRecalibration(readGroupMap);
	}
	//recalObjectInitialized = true;

	//check if estimation is required, in which case throw an error!
	if(recalObject->requiresEstimation()) throw "Can not use provided recalibration: estimation is required!";

	//prepare windows
	TWindowPairHaploid windows;
	int maxQ = params.getParameterIntWithDefault("maxQ", 100);

	//create table to store counts
	std::vector<TQualityTransformTable*> QTtables;
	for(int i=0; i<readGroups.numGroups + 1; ++i){
		TQualityTransformTable* QT = new TQualityTransformTable(maxQ);
		QTtables.push_back(QT);
	}

	//prepare output
	std::ofstream out;
	std::string filename;

	//loop over windows and add to table
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				if(recalObjectInitialized2) windows.cur->addSitesToQualityTransformTable(recalObject, recalObject2, QTtables, logfile);
				else windows.cur->addSitesToQualityTransformTable(recalObject, QTtables, logfile);
			}
		}
	}

	//print final tables for read groups
	for(int i=0; i<readGroups.numGroups; ++i){
		filename = outputName + "_" + readGroups.getName(i) + "_qualityTransformation.txt";
		out.open(filename.c_str());
		if(!out) throw "Failed to open output file '" + filename + "'!";
		QTtables[i]->printTable(out);
		out.close();

		//clean up vector
		delete QTtables.at(i);
	}
	//print table for total data
	filename = outputName + "_total_qualityTransformation.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	QTtables.at(QTtables.size()-1)->printTable(out);
	out.close();
	delete QTtables.at(QTtables.size()-1);

	//clean up memory
	windows.clear();

}

void TGenome::reportProgressParsingBamFile(const long & counter, const struct timeval & start){
	if(counter % 1000000 == 0){
		static struct timeval end;
		gettimeofday(&end, NULL);
		float runtime = (end.tv_sec  - start.tv_sec)/60.0;
		logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	}
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

	//do we also account for PMD?
	bool withPMD = params.parameterExists("withPMD");
	if(!withPMD && hasPMD) logfile->list("Note: PMD will not be reflected in the quality scores (preferred option when using ATLAS). If you want the quality scores to reflect pmd, use \"withPMD\"!");
	else if(withPMD && hasPMD) logfile->list("Probability of PMD will be reflected in new quality scores");
	else if(withPMD && !hasPMD) throw "Probability of PMD is unknown. Provide PMD patterns or remove \"withPMD\"";
	if(withPMD && !fastaReference) throw "Cannot run recalBAM withPMD without reference!";

	//other tmp variables
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	if(withPMD){
		while(alignmentParser.readAlignment(bamReader)){
	        if(useChromosome[alignmentParser.chrNumber] && alignmentParser.passedFilters && readGroups.readGroupInUse(alignmentParser.readGroupId)){
				++counter;

				alignmentParser.recalibrate(*recalObject, pmdObjects);
				alignmentParser.save(bamWriter);

				reportProgressParsingBamFile(counter, start);
			}
        }
	} else {
		while(alignmentParser.readAlignment(bamReader)){
	        if(useChromosome[alignmentParser.chrNumber] && alignmentParser.passedFilters && readGroups.readGroupInUse(alignmentParser.readGroupId)){
				++counter;

				alignmentParser.recalibrate(*recalObject);
				alignmentParser.save(bamWriter);

				reportProgressParsingBamFile(counter, start);
	        }
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

	//create index of new bam file
	logfile->listFlush("Creating index of recalibrated BAM file '" + filename + "' ...");
	BamTools::BamReader reader;
	if(!reader.Open(filename))
		throw "Failed to open BAM file '" + filename + "' for indexing!";

	// create index for BAM file
	reader.CreateIndex(BamTools::BamIndex::STANDARD);

	//close BAM file
	reader.Close();
	logfile->done();
}

void TGenome::binQualityScores(TParameters & params){
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_binnedQualityScores.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	TGenotypeMap genoMap;
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader)){
		++counter;

		//update and write (only if alignment qualities could be calculated)
		alignmentParser.binQualityScores();
		alignmentParser.save(bamWriter);

		//report
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

//---------------------------------------------------
//BAM manipulation / statistics
//---------------------------------------------------
void TGenome::assessSoftClipping(TParameters & params){
	//build table ??

	//open output file
	std::string filename = outputName + "_clippingStats.txt.gz";
	gz::ogzstream out(filename.c_str());
	if(!out)
		throw "Failed to open file '" + filename + "' for writing!";
	out << "Read\tposition\tClippedLeft\tnNotClipped\tnClippedRight\n";

	//other temp variables
	std::vector<BamTools::CigarOp>::iterator it;
	int S_left, S_right, middle;
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
	gettimeofday(&start, NULL);

	//now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader)){
		alignmentParser.assessSoftClipping(S_left, middle, S_right);

		//report
		out << alignmentParser.bamAlignment.Name << "\t" << alignmentParser.bamAlignment.Position << "\t" << S_left << "\t" << middle << "\t" << S_right << "\n";

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close output file
	out.close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::assessOverlap(TParameters & params){
	//initialize table
	int* counts = new int[maxReadLength]();

	//other variables
	float numProperPairs = 0.0;
	int counter = 0;

	//open output file
	std::string filename = outputName + "_overlapStats.txt";
	std::ofstream out(filename.c_str());
	if(!out)
		throw "Failed to open file '" + filename + "' for writing!";
	out << "overlap\tcount\tproportion\n";

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

	//now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader)){
		int overlap = alignmentParser.measureOverlap();
		if(overlap >= 0){
			++counts[overlap];
			++numProperPairs;
		}
		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//write counts to table
	for(int i=0; i<maxReadLength; ++i){
        out << i << "\t" << counts[i] << "\t" << (float) counts[i] / numProperPairs << "\n";
	}

	out.close();
	delete[] counts;
}

void TGenome::splitSingleEndReadGroups(TParameters & params){
	//read read groups and their expected lengths
	std::string filename = params.getParameterString("readGroups");
	bool allowForLarger = params.parameterExists("allowForLarger");

	logfile->listFlush("Reading single end read groups from file '" + filename + "' ...");
	std::map<int, TReadGroupMaxLength> singleEndRG;
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "!";

	//parse file
	int lineNum = 0;
	std::vector<std::string> vec;
	std::vector<std::string>::iterator it;
	int len;
	int readGroupId;
	int truncatedReadGroupId;
	BamTools::SamReadGroupIterator trunc, orig;
	std::string readGroup;

	//parse read groups
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

			//copy original tags to truncated read groups
			trunc = bamHeader.ReadGroups.Begin()+truncatedReadGroupId;
			orig = bamHeader.ReadGroups.Begin()+readGroupId;
			trunc->Library = orig->Library;
			trunc->PlatformUnit = orig->PlatformUnit;
			trunc->PredictedInsertSize = orig->PredictedInsertSize;
			trunc->ProductionDate = orig->ProductionDate;
			trunc->Program = orig->Program;
			trunc->Sample = orig->Sample;
			trunc->SequencingCenter = orig->SequencingCenter;
			trunc->SequencingTechnology = orig->SequencingTechnology;

			//add to map
			singleEndRG.emplace(int(readGroupId), TReadGroupMaxLength(len, truncatedReadGroupId, readGroup));
		}
	}
	logfile->done();
	logfile->conclude("read " + toString(singleEndRG.size()) + " read groups to be split.");

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_splitRG.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;
 	BamTools::BamAlignment bamAlignment;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;


    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		//check if this RG needs to be parsed
		std::string tmp;
		bamAlignment.GetTag("RG", tmp);
		readGroupId = readGroups.find(tmp);
		singleEndRGIT = singleEndRG.find(readGroupId);
		if(singleEndRGIT != singleEndRG.end()){
			//check length
			if(bamAlignment.Length == singleEndRGIT->second.maxLen)
				bamAlignment.EditTag("RG", "Z", singleEndRGIT->second.truncatedReadGroup);
			else if(bamAlignment.Length > singleEndRGIT->second.maxLen){
				if(allowForLarger) bamAlignment.EditTag("RG", "Z", singleEndRGIT->second.truncatedReadGroup);
				else logfile->warning("Length of read " + bamAlignment.Name + " in read group '" + readGroup + "' is > max length provided! Ignoring read.");
			}
		}

		//write
		bamWriter.SaveAlignment(bamAlignment);

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//now generate bam index
	indexBamFile(filename);

	//report
	reportProgressParsingBamFile(counter, start);
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
	logfile->done();

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
 	BamTools::BamAlignment bamAlignment;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		//get read group info
		bamAlignment.GetTag("RG", readGroup);
		oldId = readGroups.find(readGroup);

		//save as new RG
		bamAlignment.EditTag("RG", "Z", newReadGroupObject.getName(readGroupMap[oldId]));
		bamWriter.SaveAlignment(bamAlignment);

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}


void TGenome::mergePairedEndReads(TParameters & params){
	//TODO: make merge function that accounts for indels!
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_mergedReads.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	if(hasPMD) logfile->warning("PMD is given but relevant for read merging.");

	//which read groups are paired-end?
	std::string pairedRG = params.getParameterStringWithDefault("pairedReadGroups", "all");
	bool* pairedReadGroups = new bool[readGroups.numGroups];
	bool allReadGroupsPaired;
	if(pairedRG == "all"){
		//all are used, initialize to true
		for(int i=0; i<readGroups.numGroups; ++i)
			pairedReadGroups[i] = true;
		allReadGroupsPaired = true;
		logfile->list("Will merge pairs in all read groups");
	} else {
		//initialize all to false
		for(int i=0; i<readGroups.numGroups; ++i)
			pairedReadGroups[i] = false;
		//change the paired to true
		std::vector<std::string> vec;
		fillVectorFromString(pairedRG, vec, ',');
		logfile->startIndent("Will only merge pairs in the following read groups:");
		for(unsigned int i=0; i<vec.size(); ++i){
			pairedReadGroups[readGroups.find(vec.at(i))] = true;
			logfile->list(vec.at(i));
		}
		logfile->endIndent();
		allReadGroupsPaired = false;
	}

	//should we omit some reads that did not pass validateSamFile?
	bool blacklistGiven = false;
	std::map <std::string, int> readsToOmit;

	if(params.parameterExists("blacklist")){
		blacklistGiven = true;
		//open blacklist file
		std::string blacklist = params.getParameterString("blacklist");
		logfile->listFlush("Reading reads to be omitted from '" + blacklist + "...");
		std::ifstream file(blacklist.c_str());
		if(!file) throw "Failed to open file '" + blacklist + "!";

		int lineNum = 0;
		std::vector<std::string> vec;

		//fill list of reads to omit
		while(file.good() && !file.eof()){
			++lineNum;
			fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
			if(!vec.empty()) readsToOmit.emplace(vec[0], 1);
		}
		logfile->write("done! Read " + toString(lineNum) + " read names");
	}


	//other temp variables
	long counter = 0;
 	BamTools::BamAlignment bamAlignment;

	//create storage for reads until their mate was found
	std::vector< std::pair<BamTools::BamAlignment*, bool> > alignmentStorage;
	std::vector< std::pair<BamTools::BamAlignment*, bool> >::iterator it;
	std::vector< std::pair<BamTools::BamAlignment*, bool> >::iterator IT;

	BamTools::BamAlignment* alignmentPointer;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
	int curChr = -1;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		if((readsToOmit.count(bamAlignment.Name) > 0))
			continue;
		else if(allReadGroupsPaired || pairedReadGroups[readGroups.find(bamAlignment)]){

			if(abs(bamAlignment.InsertSize) < bamAlignment.AlignedBases.size()){
				if(bamAlignment.IsProperPair()){
					logfile->warning("read " + bamAlignment.Name + " was filtered out because it was longer than the insert size (" + toString(abs(bamAlignment.InsertSize)) + "<" + toString(bamAlignment.AlignedBases.size()) + ")");
				//	std::cout << "aligned bases: "<< bamAlignment.AlignedBases << std::endl;
				}
				readsToOmit.emplace(bamAlignment.Name, 1);
				continue;

			} else if(!bamAlignment.IsProperPair() || !bamAlignment.IsPrimaryAlignment()|| bamAlignment.IsDuplicate() || bamAlignment.IsSupplementary() ||(bamAlignment.IsReverseStrand() && bamAlignment.IsMateReverseStrand()) || (!bamAlignment.IsReverseStrand() && !bamAlignment.IsMateReverseStrand())){
				continue;
			}

			else {
				//if on new chromosome, empty storage
				if(curChr != bamAlignment.RefID){
					for(it = alignmentStorage.begin(); it != alignmentStorage.end(); ++it){
						bamWriter.SaveAlignment(*(it->first));
						delete it->first;
					}
					alignmentStorage.clear();
					curChr = bamAlignment.RefID;
				}
				//add alignment to storage
				if(bamAlignment.IsProperPair()){
					if(!bamAlignment.IsReverseStrand()) {
						//mate on forward strand is always first in bam file
						alignmentStorage.push_back(std::pair<BamTools::BamAlignment*, bool>(new BamTools::BamAlignment(bamAlignment), false));
					}
					else if(bamAlignment.IsReverseStrand()){
						//find first mate -> should be in storage
						for(it=alignmentStorage.begin(); it!=alignmentStorage.end(); ++it){
							if(it->first->Name == bamAlignment.Name){
								//check if this read accepts mate
								if(it->second)
									throw "First read of '" + bamAlignment.Name + "' is not paired or has already been merged!";
								//merge reads
								alignmentPointer = it->first;
								if(bamAlignment.Position > alignmentPointer->Position + alignmentPointer->AlignedBases.size()){
									//reads do not overlap -> add Ns in between
									//int numN = bamAlignment.Position - (alignmentPointer->Position + alignmentPointer)->AlignedBases.size() - 1;
									int numN = abs(bamAlignment.InsertSize) - bamAlignment.AlignedBases.size() - alignmentPointer->AlignedBases.size();
									for(int i=0; i<numN; ++i){
										alignmentPointer->AlignedBases += 'N';
										alignmentPointer->AlignedQualities += '!';
									}
									alignmentPointer->AlignedBases += bamAlignment.AlignedBases;
									alignmentPointer->AlignedQualities += bamAlignment.AlignedQualities;
								} else if(bamAlignment.Position < alignmentPointer->Position) std::cout << "Second mate of read '" << bamAlignment.Name << "' has pos < pos of first mate!";
								else {
									//reads do overlap -> merge them
									std::string alignment;
									std::string quality;
									int firstOverlap = bamAlignment.Position - alignmentPointer->Position;
									int lastOverlapPlusOne = -1;
									if((alignmentPointer->Position + alignmentPointer->AlignedBases.size()) <= (bamAlignment.Position + bamAlignment.AlignedBases.size())) lastOverlapPlusOne = alignmentPointer->AlignedBases.size();
									else lastOverlapPlusOne = bamAlignment.AlignedBases.size();

									//first copy from first mate
									alignment = alignmentPointer->AlignedBases.substr(0, firstOverlap);
									quality = alignmentPointer->AlignedQualities.substr(0, firstOverlap);
									//decide which alignment has higher quality in overlap
									for(int i=firstOverlap; i<lastOverlapPlusOne; ++i){
										//decide which quality is higher
										if(alignmentPointer->AlignedQualities.at(i) > bamAlignment.AlignedQualities.at(i - firstOverlap)){
											alignment += alignmentPointer->AlignedBases.at(i);
											quality += alignmentPointer->AlignedQualities.at(i);
										} else {
											alignment += bamAlignment.AlignedBases.at(i - firstOverlap);
											quality += bamAlignment.AlignedQualities.at(i - firstOverlap);
										}
									}

									//add rest from second
									if(alignmentPointer->Position + alignmentPointer->AlignedBases.size() < bamAlignment.Position + bamAlignment.AlignedBases.size()){
										alignment += bamAlignment.AlignedBases.substr(lastOverlapPlusOne - firstOverlap);
										quality += bamAlignment.AlignedQualities.substr(lastOverlapPlusOne - firstOverlap);
									}

									if(alignment.length() != abs(bamAlignment.InsertSize) + 1  && alignment.length() != abs(bamAlignment.InsertSize)) logfile->warning("merged alignment length of reads " + bamAlignment.Name + " is not equal to original insert size (+1)!");

									//set
									alignmentPointer->AlignedBases = alignment;
									alignmentPointer->AlignedQualities = quality;
								}
								//update

								alignmentPointer->QueryBases = alignmentPointer->AlignedBases;
								alignmentPointer->Qualities = alignmentPointer->AlignedQualities;
								alignmentPointer->Length = alignmentPointer->AlignedBases.size();
								alignmentPointer->CigarData.clear();
								alignmentPointer->CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MATCH_CHAR, alignmentPointer->Length));
								alignmentPointer->SetIsFirstMate(false);
								alignmentPointer->SetIsSecondMate(false);
								alignmentPointer->SetIsPaired(false);
								alignmentPointer->SetIsProperPair(false);
								alignmentPointer->SetIsMateReverseStrand(false);
								alignmentPointer->SetIsReverseStrand(false); //the read that comes first in BAM is always fwd strand
								alignmentPointer->MateRefID = -1;
								alignmentPointer->MatePosition = -1;
								alignmentPointer->InsertSize = 0;

								it->second = true;
//								if(bamAlignment.Name == "HWI-ST558:352:C7T9BACXX:1:1201:7676:94794") std::cout << "aligment " << bamAlignment.Name << "was updated" << std::endl;

								//write if is first in vector
								if(it == alignmentStorage.begin()){

									//write all that are OK
									for(; it != alignmentStorage.end(); ++it){
										if(it->second){
											bamWriter.SaveAlignment(*(it->first)); //saves the alignment to the bam file
//											if(bamAlignment.Name == "HWI-ST558:352:C7T9BACXX:1:2101:1301:32927") std::cout << "aligment " << bamAlignment.Name << " was output" << std::endl;
											delete it->first;
										} else {
											//first that can not be written -> erase all previous ones!
											it = alignmentStorage.erase(alignmentStorage.begin(), it);
											break;
										}
									}
									if(it == alignmentStorage.end()){
										alignmentStorage.clear();
									}
								} break;
							}
						}

						if(!alignmentStorage.empty() && it == alignmentStorage.end()) throw "One read of '" + bamAlignment.Name + "' is reverse mate, but forward one has not been read!";
					} else{
						for(it=alignmentStorage.begin(); it != alignmentStorage.end(); ++it){
							delete it->first;
						}
						alignmentStorage.clear();
						throw "One read of '" + bamAlignment.Name + "' is paired, but neither first nor second mate!";
					}
				} else {
					//read is not paired: add to storage or write
					if(alignmentStorage.empty()) bamWriter.SaveAlignment(bamAlignment);
					else alignmentStorage.push_back(std::pair<BamTools::BamAlignment*, bool>(new BamTools::BamAlignment(bamAlignment), true));
				}
			}
		}

		//read is in single-end read group
		else{
			bamWriter.SaveAlignment(bamAlignment);
		}

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();
	delete[] pairedReadGroups;

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

	//create index of new bam file
	logfile->listFlush("Creating index of recalibrated BAM file '" + filename + "' ...");
	BamTools::BamReader reader;
	if(!reader.Open(filename))
		throw "Failed to open BAM file '" + filename + "' for indexing!";

	// create index for BAM file
	reader.CreateIndex(BamTools::BamIndex::STANDARD);

	//close BAM file
	reader.Close();
	logfile->done();
}

void TGenome::downSampleBamFile(TParameters & params){
	//read downsampling rate
	std::string prob = params.getParameterString("prob");
	int times = params.getParameterIntWithDefault("times", 1);
	//check if prob is a vector of multiple probabilities
	std::vector<double> downSampleProbVector;
	if(!stringContainsOnly(prob, "-0123456789.,")) throw "Wrong format on probability list: use floating point numbers delimited by commas (e.g. 0.1,0.2,0.5).";
	fillVectorFromString(prob, downSampleProbVector, ',');

	std::vector<double>::iterator it;
	int numProbs = downSampleProbVector.size();

	if(times > 1 && numProbs > 1) throw "Replicated downsampling is not implemented for more than one prob";
	if(times > 1) numProbs = times;

	//check if probs are between 0 and 1, save in array and print them
	double* downSampleProb = new double[numProbs];
	logfile->listFlush("Will accept reads with probabilities");
	bool first = true;
	int i=0;
	if(times == 1){
		for(it=downSampleProbVector.begin(); it!=downSampleProbVector.end(); ++it, ++i){
			if(first) first = false;
			else logfile->flush(",");
			logfile->flush(" " + toString(*it));
			if(*it <= 0.0 || *it > 1.0) throw "All probabilities have to be between  0 and 1!";
			if(*it == 1.0) logfile->warning("Probability of 1 will result in identical file");
			if(*it == 0.0) logfile->warning("Probability of 0 will result in empty file");
			downSampleProb[i] = *it;
		}
	} else {
		for(int i=0; i<times; ++i){
			downSampleProb[i] = *(downSampleProbVector.begin());
		}
	}
	logfile->newLine();

	//open bam files for writing
	BamTools::BamWriter* bamWriter = new BamTools::BamWriter[numProbs];
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->startIndent("Writing results to the following files:");
	if(times > 1){
		for(i=0; i<numProbs; ++i){
			//construct and print filename
			filename = outputName + "_downsampled_" + toString(downSampleProb[i]) + "_" + toString(i) + ".bam";
			logfile->list(filename);
			//open file
			if(!bamWriter[i].Open(filename, bamHeader, references))	throw "Failed to open BAM file '" + filename + "'!";
		}
	} else {
		for(i=0; i<numProbs; ++i){
			//construct and print filename
			filename = outputName + "_downsampled_" + toString(downSampleProb[i]) + ".bam";
			logfile->list(filename);
			//open file
			if(!bamWriter[i].Open(filename, bamHeader, references))	throw "Failed to open BAM file '" + filename + "'!";
		}
	}
	logfile->endIndent();

	//other temp variables
	long counter = 0;
	double r;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader)){
		//filters
        if(!readGroups.readGroupInUse(alignmentParser.readGroupId)) continue;
        if(!useChromosome[alignmentParser.chrNumber]) continue;
		++counter;

		//accept read or not?
		for(i=0; i<numProbs; ++i){
			r = randomGenerator->getRand(); //inside loop to avoid correlation when multiple probs
			if(r < downSampleProb[i])
				alignmentParser.save(bamWriter[i]);
		}

		//report
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer and clean up memory
	for(i=0; i<numProbs; ++i){
		bamWriter[i].Close();
	}
	delete[] downSampleProb;
	delete[] bamWriter;

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::downSampleReads(TParameters & params){
	double fraction = params.getParameterDoubleWithDefault("fraction", 0.1);
	logfile->list("Each base has a probability of " + toString(fraction)+ " of being masked.");


	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_downsampledReads" + toString(fraction) + ".bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader)){

		alignmentParser.downsampleAlignment(fraction, *randomGenerator);
		alignmentParser.save(bamWriter);

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::diagnoseBamFile(TParameters & params){
    //open output files
    std::ofstream outputDepth;
    std::string outputFileNameCov = outputName + "_approximateDepth.txt";
    logfile->list("Writing sequencing depth estimates to '" + outputFileNameCov + "'");
    outputDepth.open(outputFileNameCov.c_str());
    if(!outputDepth) throw "Failed to open output file '" + outputFileNameCov + "'!";

    std::ofstream outputMQ;
    std::string outputFileNameMQ = outputName + "_MQ.txt";
    logfile->list("Writing MQ histogram to '" + outputFileNameMQ + "'");
    outputMQ.open(outputFileNameMQ.c_str());
    if(!outputMQ) throw "Failed to open output file '" + outputFileNameMQ + "'!";

    std::ofstream outputReadLen;
    std::string outputFileNameRL = outputName + "_readLength.txt";
    logfile->list("Writing read length histogram to '" + outputFileNameRL + "'");
    outputReadLen.open(outputFileNameRL.c_str());
    if(!outputReadLen) throw "Failed to open output file '" + outputFileNameRL + "'!";

    std::ofstream fragmentStats;
    std::string outputFileNameFL = outputName + "_fragmentStats.txt";
    logfile->list("Writing fragment length mean and variance to '" + outputFileNameFL + "'");
    fragmentStats.open(outputFileNameFL.c_str());
    if(!fragmentStats) throw "Failed to open output file '" + outputFileNameFL + "'!";

    double totLength = 0.0;
    int chrNum = 0;
    for(chrIterator = bamHeader.Sequences.Begin(); chrIterator!=bamHeader.Sequences.End(); ++chrIterator, ++chrNum)
        if(useChromosome[chrNum]) totLength += stringToLong(chrIterator->Length);

    //prepare reporting
    logfile->startIndent("Parsing through BAM file:");
    struct timeval start;
    gettimeofday(&start, NULL);

    //other temp variables
    long counter = 0;
    int RGInd = -1;
    std::vector<double> cov;
    double totCov = 0.0;
    int numProperPairs = 0;
    long sumFragLen = 0;
    long sumSquaredFragLen = 0;

    long** MQ = new long*[readGroups.numGroups];
    long** RL = new long*[readGroups.numGroups];
    for(int i = 0; i < readGroups.numGroups; ++i){
    	cov.push_back(0);
    	MQ[i] = new long[100];
    	RL[i] = new long[500];
    	for(int j=0; j<100; ++j) MQ[i][j]=0;
    	for(int j=0; j<500; ++j) RL[i][j]=0;
    }

    //now parse through bam file and sum number of aligned bases
    //TODO: avoid getting properties from bamAlignment, use alignmentParser
    while(alignmentParser.readAlignment(bamReader)){
    	//filters
        if(!readGroups.readGroupInUse(alignmentParser.readGroupId)) continue;
        if(!useChromosome[alignmentParser.chrNumber]) continue;
        if(!alignmentParser.passedFilters) continue;
        if(alignmentParser.isProperPair){
        	if(!alignmentParser.isReverseStrand){
        		++numProperPairs;
        		sumFragLen += abs(alignmentParser.bamAlignment.InsertSize);
        		sumSquaredFragLen += (alignmentParser.bamAlignment.InsertSize * alignmentParser.bamAlignment.InsertSize);
        	}
        }

        RGInd = alignmentParser.readGroupId;
        totCov += alignmentParser.bamAlignment.AlignedBases.length();
        cov[RGInd] += alignmentParser.bamAlignment.AlignedBases.length();
        ++MQ[RGInd][alignmentParser.bamAlignment.MapQuality];
        ++RL[RGInd][alignmentParser.bamAlignment.Length];

        //report
        ++counter;
        reportProgressParsingBamFile(counter, start);
    }

    //report
    reportProgressParsingBamFile(counter, start);logfile->list("Reached end of BAM file!");
    logfile->removeIndent();
    logfile->list("Approximate sequencing depth was estimated at " + toString(totCov/totLength));

    logfile->listFlush("Writing to output files ...");

    //cov
    outputDepth << "RG\tApproximate_depth";
    outputDepth << "\nallReadGroups\t" << totCov/totLength;
    for(int r=0; r<readGroups.numGroups; ++r){
        outputDepth << "\n" << readGroups.getName(r) << "\t" << cov[r]/totLength;
    }
    outputDepth << "\n";


    //MQ
    long tot;
    outputMQ << "RG\tMapping_quality\tCount";
    for(int i=0; i<100; ++i){
    	tot = 0;
    	for(int r=0; r<readGroups.numGroups; ++r) tot += MQ[r][i];
		outputMQ << "\nallReadGroups\t" << i << "\t" << tot;

    }
    for(int r=0; r<readGroups.numGroups; ++r){
        for(int i=0; i<100; ++i){
            outputMQ << "\n" << readGroups.getName(r) << "\t" << i << "\t" << MQ[r][i];
        }
    }
    outputMQ << "\n";


    //RL
    outputReadLen << "RG\tRead_length\tCount";
    for(int i=0; i<500; ++i){
    	tot = 0;
    	for(int r=0; r<readGroups.numGroups; ++r) tot += RL[r][i];
		outputReadLen << "\nallReadGroups\t" << i << "\t" << tot;
    }
    for(int r=0; r<readGroups.numGroups; ++r){
        for(int i=0; i<500; ++i){
            outputReadLen << "\n" << readGroups.getName(r)<< "\t" << i << "\t" << RL[r][i];
        }
    }
    outputReadLen << "\n";


    //FL
    float mean = float(sumFragLen)/float(numProperPairs);
    float var = float(sumSquaredFragLen) / float(numProperPairs) - (mean*mean);
    fragmentStats << "mean: " << mean << "\n" << "variance: " << var << "\n";

    logfile->done();

    outputDepth.close();
    outputMQ.close();
    outputReadLen.close();
    fragmentStats.close();

    for(int i = 0; i < readGroups.numGroups; ++i){
    	delete MQ[i];
    	delete RL[i];
    }
    delete [] MQ;
    delete [] RL;
}

void TGenome::allelicDepth(TParameters & params){
	std::ofstream output;
	std::string outputFileName = outputName + "_allelicDepth.txt";
	logfile->list("Writing allelic imbalance table to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	//int maxCov = params.getParameterIntWithDefault("maxCov", 20);
	int maxAllelicDepth = params.getParameterInt("maxAllelicDepth");
	int size = maxAllelicDepth+1; // need 0 bin
	int nCharOnLine = 0;

	//prepare array
	long**** siteCounts = new long***[size];
	for(int i=0; i<size; ++i){
		siteCounts[i] = new long**[size];
		for(int j=0; j<size; ++j){
			siteCounts[i][j] = new long*[size];
			for(int k=0; k<size; ++k){
				siteCounts[i][j][k] = new long[size];
				for(int l=0; l<size; ++l){
					siteCounts[i][j][k][l] = 0;
				}
			}
		}
	}

	//write header
	output << "A\tC\tG\tT\tCounts\tDepth" << std::endl;

	//prepare windows
	TWindowPairDiploid windows;
	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			windows.cur->countAlleles(siteCounts, maxAllelicDepth);
			logfile->listFlush("Adding imbalance values to table ...");
			logfile->write(" done!");
		}
	}

	//write to file
	for(int i=0; i<(size); ++i){
		for(int j=0; j<(size); ++j){
			for(int k=0; k<(size); ++k){
				for(int l=0; l<(size); ++l){
					output << i << "\t" << j << "\t" << k << "\t" << l << "\t" << siteCounts[i][j][k][l] << "\t" << i + j + k + l;
					output << std::endl;
				}
			}
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
	delete[] siteCounts;
}

void TGenome::estimateApproximateDepthPerWindow(TParameters & params){
	//open output file
	std::ofstream output;
	std::string outputFileName = outputName + "_depthPerWindow.txt";
	logfile->list("Writing sequencing depth estimates to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int nCharOnLine = 0;

	//write header
	output << "chr\tstart\tend\tdepth" << std::endl;

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//write to file
			logfile->listFlush("Writing sequencing depth estimates to file ...");
			if(windows.cur->depth == -1.0) output << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t" << "0" << "\n";
			else output << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t" << windows.cur->depth << "\n";
			logfile->done();
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
}

void TGenome::estimateDepthPerSite(TParameters & params){
	std::ofstream output;
	std::ofstream outputNormalized;
	std::ofstream outputQuantiles;

	std::string outputFileName = outputName + "_depthPerSite.txt";
	std::string outputFileNameNormalized = outputName + "_normalizedCumulativeDepthPerSite.txt";
	std::string outputFileNameQuantiles = outputName + "_quantilesDepthPerSite.txt";


	logfile->list("Writing cumulative depth distribution to '" + outputFileName + "'");
	logfile->list("Writing normalized cumulative depth distribution to '" + outputFileNameNormalized + "'");
	logfile->list("Writing quantiles of cumulative depth distribution to '" + outputFileNameQuantiles + "'");


	output.open(outputFileName.c_str());
	outputNormalized.open(outputFileNameNormalized.c_str());
	outputQuantiles.open(outputFileNameQuantiles.c_str());


	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	if(!outputNormalized) throw "Failed to open output file '" + outputFileNameNormalized + "'!";
	if(!outputQuantiles) throw "Failed to open output file '" + outputFileNameQuantiles + "'!";


	int maxCov = params.getParameterIntWithDefault("maxDepth", 20);
	if(!maxCov) throw "No maximum depth specified!";
	int size = maxCov + 2; // need 0 bin and >maxCov bin
	int nCharOnLine = 0;
	double totalSites = 0.0;

	//prepare arrays
	long * siteDepth = new long[size];
	for(int i=0; i<size; ++i){
		siteDepth[i] = 0;
	}

	long * siteDepthNormalized = new long[size];
	for(int i=0; i<size; ++i){
		siteDepthNormalized[i] = 0;
	}

	//write headers
	output << "depth\tcounts" << std::endl;
	outputNormalized << "depth\tcounts" << std::endl;
	outputQuantiles << "percent\tquantile" << std::endl;

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			logfile->listFlush("Adding depth to table ...");
			windows.cur->calcDepthPerSite(siteDepth, maxCov);
			logfile->done();
		}
	}

	//write to files
	//read counts
	for(int i=0; i<(size-1); ++i){
		output << i << "\t" << siteDepth[i] << "\n";
		totalSites += (double) siteDepth[i];
	}
	totalSites += siteDepth[size - 1];
	output << ">" << maxCov << "\t" << siteDepth[size - 1] << std::endl;

	//normalized cumulative distribution and quantiles
	double cumul = 0.0;
	double norm = 0.0;
	float percentages[14] = {0.001, 0.005, 0.01, 0.025, 0.05, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999};
	int numberOfPercentages = 14;
	int quantiles[numberOfPercentages];
	std::fill_n(quantiles, numberOfPercentages, -1);
	for(int i=0; i<(size-1); ++i){
		cumul += (double) (siteDepth[i]);
		if(i > 0 && norm == 1.0) outputNormalized << i << "\t" << 0 << "\n";
		else {
			norm = cumul / totalSites;
			outputNormalized << i << "\t" << norm << "\n";
		}
		for(int p=0; p<numberOfPercentages; ++p){
			//if smallest quantiles = 0
			if(i == 0 && norm > percentages[p]) quantiles[p] = i;
			else if(quantiles[p] == -1 && norm >= percentages[p]){
				quantiles[p] = i;
			}
		}
		if(quantiles[numberOfPercentages] == -1 && i >= percentages[numberOfPercentages]) quantiles[numberOfPercentages] = i;
	}
	if(norm == 1.0) outputNormalized << ">" << maxCov << "\t" << 0 << "\n";
	else {
		cumul += (double) siteDepth[size - 1];
		norm = cumul / totalSites;
		outputNormalized << ">" << maxCov << "\t" << norm << std::endl;
	}

	for(int p=0; p<numberOfPercentages; ++p){
		if(quantiles[p] == -1) quantiles[p] = maxCov;
	}

	for(int i=0; i < numberOfPercentages; ++i){
		outputQuantiles << percentages[i] << "\t" << quantiles[i] << std::endl;
	}

	//clean up
	if(nCharOnLine > 0){
		output << '\n';
		outputNormalized << '\n';
	}
	output.close();
	outputNormalized.close();
	outputQuantiles.close();

	delete[] siteDepth;
	delete[] siteDepthNormalized;

}

void TGenome::writeDepthPerSite(TParameters & params){
	gz::ogzstream out;

	std::string outputFileName = outputName + "_depthPerSite.txt.gz";
	logfile->list("Writing per site depth to '" + outputFileName + "'");

	out.open(outputFileName.c_str());
	if(!out) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	out << "chr\tpos\tdepth" << std::endl;

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			logfile->listFlush("Writing depth per site ...");
			windows.cur->printDepthPerSite(out, chrIterator->Name);
			logfile->done();
		}
	}

	//clean up
	out.close();
}


//---------------------------------------------------
//PMD
//---------------------------------------------------
void TGenome::estimatePMD(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//prepare readGroup map
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);

	//prepare PMD table
	int maxLengthForInference = params.getParameterIntWithDefault("length", 50);
	logfile->list("Estimating PMD at the first " + toString(maxLengthForInference) + " positions.");
	TPMDTables pmdTables(&readGroups, maxLengthForInference, maxReadLength, readGroupMap);

	//measure progress and runtime
	struct timeval start;
	long numreadsAdded = 0;
	gettimeofday(&start, NULL);

	//iterate through BAM file
	while(alignmentParser.readAlignment(bamReader)){
        if(useChromosome[alignmentParser.chrNumber] && alignmentParser.passedFilters && readGroups.readGroupInUse(alignmentParser.readGroupId))
			alignmentParser.addToPMDTables(pmdTables);

		//report
		++numreadsAdded;
		reportProgressParsingBamFile(numreadsAdded, start);
	}

	//report
	reportProgressParsingBamFile(numreadsAdded, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

	//print tables and data
	std::string filename = outputName + "_PMD_Table.txt";
	logfile->listFlush("Writing PMD table to '" + filename + "' ...");
	pmdTables.writeTable(filename);
	logfile->done();
	filename = outputName + "_PMD_Table_counts.txt";
	logfile->listFlush("Writing PMD table of counts to '" + filename + "' ...");
	pmdTables.writeTableWithCounts(filename);
	logfile->done();
	filename = outputName + "_PMD_input_Empiric.txt";
	logfile->listFlush("Writing PMD input file to '" + filename + "' ...");
	pmdTables.writePMDFile(filename);
	logfile->done();

	//estimate exponential model
	filename = outputName + "_PMD_input_Exponential.txt";
	logfile->listFlush("Estimating PMD exponential models and writing them to '" + filename + "' ...");
	int numNRIterations = params.getParameterIntWithDefault("numNRIterations", 100);
	double eps = params.getParameterDoubleWithDefault("eps", 0.001);
	pmdTables.fitExponentialModel(numNRIterations, eps, filename, logfile);
	logfile->done();
}

void TGenome::estimateGTPattern(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//prepare readGroup map
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);

	//prepare PMD table
	int maxLengthForInference = params.getParameterIntWithDefault("length", 50);
	logfile->list("Estimating GT pattern at the first " + toString(maxLengthForInference) + " positions.");
	TPMDTables pmdTables(&readGroups, maxLengthForInference, maxReadLength, readGroupMap);

	//measure progress and runtime
	struct timeval start;
	long numreadsAdded = 0;
	gettimeofday(&start, NULL);

	//iterate through BAM file
	while(alignmentParser.readAlignment(bamReader)){
        if(useChromosome[alignmentParser.chrNumber] && alignmentParser.passedFilters && readGroups.readGroupInUse(alignmentParser.readGroupId))
			alignmentParser.addToPMDTables(pmdTables);

		//report
		++numreadsAdded;
		reportProgressParsingBamFile(numreadsAdded, start);
	}

	//report
	reportProgressParsingBamFile(numreadsAdded, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

	//print tables and data
	std::string filename = outputName + "_PMD_Table.txt";
	logfile->listFlush("Writing PMD table to '" + filename + "' ...");
	pmdTables.writeTable(filename);
	logfile->done();
	filename = outputName + "_PMD_Table_counts.txt";
	logfile->listFlush("Writing PMD table of counts to '" + filename + "' ...");
	pmdTables.writeTableWithCounts(filename);
	logfile->done();
	filename = outputName + "_GT_Empiric.txt";
	pmdTables.writeGTFile(filename);
	logfile->listFlush("Writing GT empiric pattern to '" + filename + "' ...");


	logfile->done();
}

void TGenome::runPMDS(TParameters & params){
	//parse bam file and calculate PMDS for each read (seeSkoglund et al. 2014)
	//write new bam file with PMDS score added
	//parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)

	initializeRecalibration(params);
	if(!fastaReference) throw "Cannot run PMDS without reference!";

	//get parameters
	double pi = params.getParameterDoubleWithDefault("pi", 0.001);
	logfile->list("Running PMDS with rate of polymorphism (pi) = " + toString(pi));
	double minPMDS = params.getParameterDoubleWithDefault("minPMDS", -10000);
	double maxPMDS = params.getParameterDoubleWithDefault("maxPMDS", 10000);
	logfile->list("Filtering out reads with " + toString(minPMDS) + " > PMDS > " + toString(maxPMDS));


	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
	gettimeofday(&start, NULL);
	float runtime;
	long counter = 0, counterF = 0;

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_PMDS.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//now parse through bam file and write alignments
	double PMDS;
	while(alignmentParser.readAlignment(bamReader)){
		++counter;

        if(useChromosome[alignmentParser.chrNumber] && alignmentParser.passedFilters && readGroups.readGroupInUse(alignmentParser.readGroupId)){
			//recalibrate quality scores
			alignmentParser.recalibrate(*recalObject);

			//calc PMD
			PMDS = alignmentParser.calculatePMDS(pi, pmdObjects);

			//update and write
			if(PMDS > minPMDS && PMDS < maxPMDS){
				alignmentParser.updateOptionalSamField("DS", PMDS);
				alignmentParser.save(bamWriter);
			} else ++counterF;
		} else alignmentParser.save(bamWriter);

		//report progress
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			logfile->list("Analyzed " + toString(counter) + " reads in " + toString((end.tv_sec  - start.tv_sec)/60.0) + " min and filtered out " + toString(counterF) + " of them!");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Analyzed " + toString(counter) + " reads in " + toString(runtime) + " min. and filtered out " + toString(counterF) + " of them!");

}

