/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"


//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
TFastaBuffer::TFastaBuffer(BamTools::Fasta* Reference){
	bufferSize = 10000000;
	reference = Reference;
	referenceSequence = "";
	curStart = -1;
	curChr = -1;
	curEnd = -1;
}


void TFastaBuffer::moveTo(const int & chr, const int32_t & pos){
	curChr = chr;
	curStart = pos;
	curEnd = pos + bufferSize;
	if(!reference->GetSequence(chr, curStart, curEnd, referenceSequence))
		throw "Problem reading " + toString(chr) + ":" + toString(curStart) + "-" + toString(curEnd) + " from fasta file!";
}

void TFastaBuffer::fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref){
	//move buffer, if necessary
	if(chr != curChr || end > curEnd || start < curStart)
		moveTo(chr, start);

	//now copy to string
	ref.assign(referenceSequence, start - curStart, end - start + 1);
}

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
TAlignmentParser::TAlignmentParser(){
	logfile = NULL;
	_keepDuplicates = false;
	parse = false;
	previousAlignmentPos = -1;
	previousAlignmentChr = -1;
	oldAlignment = NULL;
	oldAlignmentMustBeConsidered = false;

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;
//	curStart = -1;
//	curEnd = -1;

	//window parameters
	windowNumber = -1;
	windowSize = -1;
	numWindowsOnChr = -1;
	maxReadLength = -1;
	maxMissing = 1.0;
	maxRefN = 1.0;
	windowsPredefined = false;
	predefinedWindows = NULL;

	//masks
	doMasking = false;
	doCpGMasking = false;
	considerRegions = false;
	mask = NULL;

	//filters
	applyQualityFilter = false;
	minQual = 33;
	maxQual = 126;
	applyDepthFilter = false;
	minDepth = 0;
	maxDepth = 10000;
	minPhredInt = 0;
	maxPhredInt = 93;
	minQualForPrinting = 33;
	maxQualForPrinting = 126;
	trimReads = false;
	trimmingLength3Prime = 0;
	trimmingLength5Prime = 0;

	//limit chr and windows
	limitWindows = -1;
	limitChr = -1;
	useChromosome = NULL;

	//reference
	hasReference = false;
	fastaReference = NULL;
	fastaBuffer = NULL;

	//post mortem damage and recalibration
	hasPMD = false;
	pmdObjects = NULL;
	doRecalibration = false;
	doRecalibration2 = false;
	recalObjectInitialized = false;
	recalObjectInitialized2 = false;
	recalObject = NULL;
	recalObject2 = NULL;
};

TAlignmentParser::TAlignmentParser(TParameters & params, TLog* Logfile){
	TAlignmentParser();
	logfile->list("Will only consider reads up to " + toString(maxReadLength) + " bp.");

	init(params, Logfile);
};

void TAlignmentParser::init(TParameters & params, TLog* Logfile){
	logfile = Logfile;

	//open BAM file
	filename = params.getParameterString("bam");
	logfile->list("Reading data from BAM file '" + filename + "'.");
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";
	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + filename + "'!";


	//initialize bam stuff
	bamHeader = bamReader.GetHeader();
	chrIterator = bamHeader.Sequences.End();

	//---------------------
	//window parameters
	//---------------------
	maxReadLength = params.getParameterIntWithDefault("maxReadLength", 1000);

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

	//--------------------
	//limit chr and windows
	//--------------------

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

	//------------
	//masks
	//------------
	//normal mask
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

	//reverse masking
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

	//CpG mask
	if(params.parameterExists("maskCpG")){
		if(!fastaReference) throw "Cannot mask CpG sites without reference!";
		doCpGMasking = true;
		std::string maskFile = params.getParameterString("maskCpG");
		logfile->list("Will mask all CpG sites");
	} else doCpGMasking = false;

	//------------
	//filters
	//------------
	//depth filter
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
	setQualityFilters(minPhredInt+33, maxPhredInt+33);
	logfile->list("Will filter out bases with quality outside the range [" + toString(minPhredInt) + ", " + toString(maxPhredInt) + "]");

	//quality filters for printing
	int minOutQual = params.getParameterIntWithDefault("minOutQual", 1) + 33;
	if(minOutQual < 0) throw "minOutQual must be >= 0!";
	int maxOutQual = params.getParameterIntWithDefault("maxOutQual", 93) + 33;
	if(maxOutQual < minOutQual) throw "maxOutQual must be >= minOutQual!";
	setQualityRangeForPrinting(minOutQual, maxOutQual);
	logfile->list("Will print qualities truncated to [" + toString(minOutQual) + ", " + toString(maxOutQual) + "]");

	//filter for missing reference
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 1.0);
	if(maxMissing > 1.0) throw "maxMissing must be smaller or equal to 1.0!";

	maxRefN = params.getParameterDoubleWithDefault("maxRefN", 1.0);
	if(maxRefN > 1.0) throw "maxRefN must be smaller or equal to 1.0!";
	if(maxRefN < 1.0 && hasReference == false) throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided.";

	//-----------------
	//read groups
	//-----------------

	readGroups.fill(bamHeader);

	//limit readGroups
	if(params.parameterExists("readGroup")){
		readGroups.filterReadGroups(params.getParameterString("readGroup"));
		logfile->startIndent("Will limit analysis to the following read groups:");
		readGroups.printReadgroupsInUse(logfile);
		logfile->endIndent();
	}

	//------------
	//recal and pmd
	//------------

	initializePostMortemDamage(params);
	initializeRecalibration(params);

	//------------
	//other
	//------------

	if(params.parameterExists("keepDuplicates")){
		keepDuplicates();
		logfile->list("Will keep duplicate reads.");
	}
};

void TAlignmentParser::setQualityFilters(int MinQual, int MaxQual){
	applyQualityFilter = true;
	minQual = MinQual;
	maxQual = MaxQual;
};

void TAlignmentParser::setQualityRangeForPrinting(int minQual, int maxQual){
	minQualForPrinting = minQual;
	maxQualForPrinting = maxQual;
};

void TAlignmentParser::setReadTrimming(int trim3Prime, int trim5Prime){
	trimmingLength3Prime = trim3Prime;
	trimmingLength5Prime = trim5Prime;
	trimReads = true;
};

void TAlignmentParser::addReference(BamTools::Fasta* reference){
	hasReference = true;
	fastaReference = reference;
	fastaBuffer = new TFastaBuffer(reference);
};

void TAlignmentParser::fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment){
	if(!hasReference) //is this check really necessary?
		throw "No reference provided!";

	fastaBuffer->fill(alignment.chrNumber, alignment.position, alignment.position + alignment.alignedPos[alignment.length-1], referenceSequence);
};

//--------------
//move genome
//--------------
void TAlignmentParser::jumpToEnd(){
	chrIterator = bamHeader.Sequences.End();
	chrNumber = -1;
}

void TAlignmentParser::restartChromosome(TWindow & window){
	chrIterator = bamHeader.Sequences.Begin();
	chrNumber = 0;

	moveChromosome(window);
}

bool TAlignmentParser::iterateChromosome(TWindow & window){
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
		window.end = 0;
		chrIterator = bamHeader.Sequences.End();
		return false;
	}

	moveChromosome(window);
	return true;
}

void TAlignmentParser::moveChromosome(TWindow & window){
	std::cout << "moving chr" << std::endl;
	//jump reader
	bamReader.Jump(chrNumber, 0);
	oldAlignmentMustBeConsidered = false;
	chrLength = stringToLong(chrIterator->Length);

	//restart windows
//	curStart = 0;
//	curEnd = 0;
	previousAlignmentPos = -1;
	windowNumber = 0;
	if(windowsPredefined){
		predefinedWindows->setChr(chrIterator->Name);
		numWindowsOnChr = predefinedWindows->getNumWindowsOnCurChr();
		int nextEnd = predefinedWindows->curWindowEnd();
		if(nextEnd > chrLength) nextEnd = chrLength + 1;
		else window.move(predefinedWindows->curWindowStart(), nextEnd);
	} else {
		numWindowsOnChr = ceil(chrLength / (double) windowSize);
		int nextEnd = windowSize;
		//TODO:!!! removed +1 because we are zero-based. Check if true!
		if(nextEnd > chrLength) nextEnd = chrLength;
		window.move(0, nextEnd);
		std::cout << "moved window to " << window.start << " and " << window.end << std::endl;
	}

	//advance mask
	if(doMasking || considerRegions) mask->setChr(chrIterator->Name);

	//write progress
	logfile->startNumbering("Parsing chromosome '" + chrIterator->Name + "':");
}

bool TAlignmentParser::moveToNextWindow(TWindow & window){
	if(window.end > 0) logfile->endIndent();

	//move to next region
//	curStart = window.start;
//	curEnd = window.end;
	if(window.end >= chrLength || windowNumber >= limitWindows) return false;

	//move next
	if(windowsPredefined){
		if(numWindowsOnChr < 1){
			logfile->conclude("No windows on this chromosome.");
			return false;
		}
		//jump reader if large gap to previous window
		//TODO:: check if this does not mean we miss reads starting prior to the window but extending into it.
		if(window.start - window.end > maxReadLength)
			bamReader.Jump(chrNumber, window.start);

		//now move coordinates of next window
		if(predefinedWindows->nextWindow()){
			int nextEnd = predefinedWindows->curWindowEnd();
			if(nextEnd > chrLength) nextEnd = chrLength;
			window.move(predefinedWindows->curWindowStart(), nextEnd);
		} else {
			window.move(chrLength, chrLength+1);
		}
	} else {
		long nextEnd = window.end + windowSize;
		if(nextEnd > chrLength)
			nextEnd = chrLength;
		window.move(window.end, nextEnd);
	}

	++windowNumber;

	//report
	logfile->number("Window [" + toString(window.start) + ", " + toString(window.end) + "] of " + toString(numWindowsOnChr) + " on '" + chrIterator->Name + "':");
	logfile->addIndent();

	return true;
};

//------------------------------
//initialize PMD and recalibration
//------------------------------
void TAlignmentParser::initializePostMortemDamage(TParameters & params){
	logfile->startIndent("Initializing Post Mortem Damage (PMD):");
	//create an array of TPMD objects for each read group
	pmdObjects = new TPMD[readGroups.size()];

	//now fill them!
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		//all read groups have the same pmd
		logfile->list("Initializing one PMD function for all read groups.");
		pmdObjects[0].initialize(params, logfile);
		for(int i=1; i<readGroups.size(); ++i)
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
		for(int i=0; i<readGroups.size(); ++i){
			if(!pmdObjects[i].functionInitialized(pmdCT)) throw "PMD C->T for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
			if(!pmdObjects[i].functionInitialized(pmdGA)) throw "PMD G->A for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
		}
		hasPMD = true;
	} else {
		//no post mortem damage
		logfile->list("Assuming there is no PMD in the data.");
		std::string pmdString = "none";
		for(int i=0; i<readGroups.size(); ++i){
			pmdObjects[i].initializeFunction(pmdString, pmdGA);
			pmdObjects[i].initializeFunction(pmdString, pmdCT);
		}
	}
	logfile->endIndent();
}

void TAlignmentParser::initializeRecalibration(TParameters & params){
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

//------------------------------
//reading data
//------------------------------
bool TAlignmentParser::readAlignment(BamTools::BamReader & bamReader, TAlignment & alignment){
	//make sure container is empty
	alignment.clear();

	if(!bamReader.GetNextAlignment(bamAlignment))
		return false;

	//check read length
	if(bamAlignment.AlignedBases.size() > maxReadLength)
		throw "Alignment '" +  bamAlignment.Name + "' is longer than the max read length! Please change max read length to parse this data.";

	//fill alignment
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
	int readGroupId = readGroups.find(readGroup);
	alignment.fill(bamAlignment, readGroupId);

	//check if insert size is shorter than read, this means we are reading the adaptor sequence
	bool filtersPassed = true;
	if(bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) <= bamAlignment.AlignedBases.length()){
		logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
		filtersPassed = false;
	} else {
		//apply filters: read group in use and basic QC
		filtersPassed = readGroups.readGroupInUse(readGroupId)
						&& bamAlignment.IsMapped() && !bamAlignment.IsFailedQC()
						&& bamAlignment.IsPrimaryAlignment()
						&& (_keepDuplicates || !bamAlignment.IsDuplicate());
	}

	//set filter
	alignment.setFiltersPassed(filtersPassed);

	if(parse){
		std::cout << "########### parsing alignment!" << std::endl;
		alignment.parse(genoMap, qualityMap);

		//add missing information to bases
		alignment.fillReadGroupInfo(readGroupId);

		if(doRecalibration)
			alignment.recalibrate(*recalObject, qualityMap);
		if(hasPMD)
			alignment.fillPmdProbabilities(pmdObjects);
		if(hasReference)
			fillReferenceSequence(fastaBuffer, alignment);
		if(applyQualityFilter)
			alignment.filterForBaseQuality(minQual, maxQual);
	}

	return true;
}


//---------------------
//read data in windows
//---------------------
bool TAlignmentParser::readDataInWindows(TWindow & window, TReadGroups & readGroups){
	std::cout << "reading data in windows" << std::endl;
	setParsingToTrue();

	while(iterateChromosome(window)){
		//report
		logfile->number("Window [" + toString(window.start) + ", " + toString(window.end) + "] of " + toString(numWindowsOnChr) + " on '" + chrIterator->Name + "':");
		logfile->addIndent();
		for(int w=0; w<numWindowsOnChr; ++w){
			//measure runtime
			struct timeval start, end;
			gettimeofday(&start, NULL);

			logfile->listFlush("Reading data ...");
			while(readAlignmentsIntoWindow(window, readGroups)) continue;
			window.fillSites();
			if(hasReference) window.addReferenceBaseToSites(*fastaReference, previousAlignmentChr);

			gettimeofday(&end, NULL);
			logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");

			applyFilters(window);
			moveToNextWindow(window);
		}
/*		while(iterateWindow(window)){
			std::cout << "after iterateWindow start and end are " <<window.start << " and " << window.end << std::endl;
			while(readAndFilterWindow(window, readGroups))
			return true;
		}*/
	}
	return false;
}

bool TAlignmentParser::readAlignmentsIntoWindow(TWindow & window, TReadGroups & readGroups){
	//add old alignment if necessary
	if(oldAlignmentMustBeConsidered){
		oldAlignmentMustBeConsidered = false;
		if(!addToWindow((*oldAlignment), window)){
			logfile->conclude("No data in this window.");
			return false; //still only in next window
		} else {
			//update previous alignment info
			previousAlignmentPos = oldAlignment->position;
			previousAlignmentChr = oldAlignment->chrNumber;
		}
	}
	//get alignment to fill from window and make sure it's empty
	TAlignment* alignmentP = window.getNewAlignment(maxReadLength);
	alignmentP->clear();

	//pass it to reader
	if(readAlignment(bamReader, *alignmentP)){
		if(readGroups.readGroupInUse(alignmentP->readGroupId)){
			if(!addToWindow((*alignmentP), window)){
				//read is beyond window and should be added to next
				oldAlignment = alignmentP;
				oldAlignmentMustBeConsidered = true;
				return false;
			} else {
				//update previous alignment info
				previousAlignmentPos = alignmentP->position;
				previousAlignmentChr = alignmentP->chrNumber;
			}
		}
		return true;
	}
	return false;
};

bool TAlignmentParser::addToWindow(TAlignment& alignment, TWindow & window){
	//check if bam file is sorted
	if(alignment.getPosition() < previousAlignmentPos)
		throw "BAM file must be sorted by position!";

	//and add
	if((alignment.getPosition() >= window.start) && alignment.getPosition() <= window.end){
		//add alignment with complete bases to window
		window.addAlignment(&alignment);
		return true; //continue
	} else
		return false;
}

bool TAlignmentParser::applyFilters(TWindow & window){
	if(window.numReadsInWindow > 0){
		//apply masks and filters
		if(doMasking){
			logfile->listFlush("Masking sites ...");
			window.applyMask(mask, considerRegions);
			logfile->done();
		} else if(considerRegions){
			logfile->listFlush("Masking sites outside regions ...");
			window.applyMask(mask, considerRegions);
			logfile->done();
		} else if(doCpGMasking){
			logfile->listFlush("Masking CpG sites ...");
			window.maskCpG(*fastaReference, previousAlignmentChr);
			logfile->done();
		} if(applyDepthFilter){
			window.applyDepthFilter(minDepth, maxDepth);
		} if(maxRefN < 1.0 && hasReference == true){
			window.addReferenceBaseToSites(*fastaReference, previousAlignmentChr);
			window.calcFracN();
		}

		//calc sequencing depth
		window.calcDepth();

		//report
		logfile->conclude("read data from " + toString(window.numReadsInWindow) + " reads.");
		logfile->conclude("sequencing depth is " + toString(window.depth));
		logfile->conclude(toString(window.fractionsitesDepthAtLeastTwo * 100) + "% of all sites are covered at least twice");
		logfile->conclude(toString(window.fractionSitesNoData * 100) + "% of all sites have no data");
		if(window.fractionSitesNoData > maxMissing){
			logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			return false;
		}
		if(maxRefN < 1.0 && hasReference == true){
			logfile->conclude(toString(window.fractionRefIsN * 100) + "% of all reference bases are 'N'");
			if(window.fractionRefIsN > maxRefN){
				logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
				return false;
			}
		}
		return true;
	} else {
		logfile->conclude("No data in this window.");
		return false;
	}
}


