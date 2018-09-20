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
	oldAlignmentInitialized = false;
	oldAlignmentMustBeConsidered = false;

	curReadGroupID = -1;

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;

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
	minQual = 0;
	maxQual = 93;
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

TAlignmentParser::TAlignmentParser(int MaxReadLength, TParameters & params, TLog* Logfile){
	TAlignmentParser();
	logfile->list("Will only consider reads up to " + toString(maxReadLength) + " bp.");

	init(MaxReadLength, params, Logfile);
};

TAlignmentParser::~TAlignmentParser(){
		if(hasReference)
			delete fastaBuffer;
		if(doMasking)
			delete mask;
		if(windowsPredefined)
			delete predefinedWindows;
		if(useChromosome)
			delete[] useChromosome;
		if(recalObjectInitialized) delete recalObject;
		if(pmdObjects) delete[] pmdObjects;
		if(oldAlignmentInitialized)
			delete oldAlignment;
	}

void TAlignmentParser::init(int MaxReadLength, TParameters & params, TLog* Logfile){
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
	maxReadLength = MaxReadLength;
	oldAlignment = new TAlignment(maxReadLength);
	oldAlignmentInitialized = true;

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
		for(std::vector<std::string>::iterator it=vec.begin(); it!=vec.end(); ++it){
			//find chromosome
			int num = 0;
			for(BamTools::SamSequenceIterator chrIt = bamHeader.Sequences.Begin(); chrIt != bamHeader.Sequences.End(); ++chrIt, ++num){
				if(chrIt->Name == *it){
					useChromosome[num] = true;
					logfile->list(*it);
					break;
				}
				if(chrIt == bamHeader.Sequences.End()) throw "Chromosome '" + *it + "' is not present in the bam header!";
			}
		}
		logfile->endIndent();


	} else {
		if(params.parameterExists("limitChr")){
			//set all chromosomes to false
			int num = 0;
			for(BamTools::SamSequenceIterator chrIt = bamHeader.Sequences.Begin(); chrIt != bamHeader.Sequences.End(); ++chrIt, ++num)
				useChromosome[num] = false;

			limitChr = params.getParameterIntWithDefault("limitChr", 1000000);
			logfile->list("Will limit analysis to the first " + toString(limitChr) + " chromosomes.");

			num = 0;
			for(BamTools::SamSequenceIterator chrIt = bamHeader.Sequences.Begin(); chrIt != bamHeader.Sequences.End(); ++chrIt, ++num){
				if(num == limitChr)
					break;
				useChromosome[num] = true;
				logfile->list(chrIt->Name);
			}
		} else {
			for(int i=0; i<bamHeader.Sequences.Size(); ++i)
				useChromosome[i] = true;
		}
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

	fastaBuffer->fill(alignment.chrNumber, alignment.position, alignment.position + alignment.bases[alignment.length-1].alignedPos, referenceSequence);
};

std::string TAlignmentParser::chrNumberToName(int chrNumber){
	int counter = 0;
	for(BamTools::SamSequenceIterator chrIt=bamHeader.Sequences.Begin(); chrIt!=bamHeader.Sequences.End(); ++chrIt){
		if(counter == chrNumber)
			return chrIt->Name;
		++counter;

	}
	throw "chrNumber not in header";
}

//--------------
//move genome
//--------------
void TAlignmentParser::jumpToEnd(){
	chrIterator = bamHeader.Sequences.End();
	chrNumber = -1;
}

void TAlignmentParser::restartChromosomes(TWindow & window){
	chrIterator = bamHeader.Sequences.Begin();
	chrNumber = 0;

	moveChromosome(window);
}

void TAlignmentParser::moveChromosome(TWindow & window){
	//jump reader
	oldAlignmentMustBeConsidered = false;
	chrLength = stringToLong(chrIterator->Length);

	//restart windows
	previousAlignmentPos = -1;
	windowNumber = 0;

	if(windowsPredefined){
		//find next used chromosome with windows
		do {
			predefinedWindows->setChr(chrIterator->Name);
			numWindowsOnChr = predefinedWindows->getNumWindowsOnCurChr();
			if(numWindowsOnChr < 1 || useChromosome[chrNumber] == false ){
				if(useChromosome[chrNumber])
					logfile->conclude("No windows on chromosome " + chrIterator->Name + ".");
				++chrIterator;
				++chrNumber;
				predefinedWindows->setChr(chrIterator->Name);
				numWindowsOnChr = predefinedWindows->getNumWindowsOnCurChr();
			}
		} while((numWindowsOnChr < 1 || !useChromosome[chrNumber]) && chrIterator != bamHeader.Sequences.End());

		//now jump
		window.move(predefinedWindows->curWindowStart(), predefinedWindows->curWindowEnd(), chrNumber);
		bamReader.Jump(chrNumber, window.start);

	} else {
		while(!useChromosome[chrNumber]){
			++chrIterator;
			++chrNumber;
			chrLength = stringToLong(chrIterator->Length);
		}
		bamReader.Jump(chrNumber, 0);
		numWindowsOnChr = ceil(chrLength / (double) windowSize);
		int nextEnd = windowSize;
		//TODO:!!! removed +1 because we are zero-based. Check if true!
		if(nextEnd > chrLength){
			nextEnd = chrLength;
		}
		window.move(0, nextEnd, chrNumber);
	}


	if(chrIterator == bamHeader.Sequences.End()){
		return;
	}

	//advance mask
	if(doMasking || considerRegions) mask->setChr(chrIterator->Name);

	//write progress
	logfile->endIndent();
	logfile->startNumbering("Parsing chromosome '" + chrIterator->Name + "':");
}

bool TAlignmentParser::moveToNextWindowOnChr(TWindow & window){

	if(window.end > 0) logfile->endIndent();

	//move to next region
	++windowNumber;
	if(window.end >= chrLength || windowNumber >= limitWindows)
		return false;

	//move next
	long nextEnd = window.end + windowSize;
	if(nextEnd > chrLength)
		nextEnd = chrLength;
	window.move(window.end, nextEnd, chrNumber);

	return true;
};

bool TAlignmentParser::moveToNextPredefinedWindow(TWindow & window){

	if(window.end > 0) logfile->endIndent();

	++windowNumber;
	if(windowNumber >= limitWindows)
		return false;
	if(predefinedWindows->nextWindow()){
		window.move(predefinedWindows->curWindowStart(), predefinedWindows->curWindowEnd(), chrNumber);
		//should we jump or are we already close enough to next window
		if(abs(window.start - previousAlignmentPos) > maxReadLength){
			if(window.start - maxReadLength < 0)
				bamReader.Jump(chrNumber, 0);
			else
				bamReader.Jump(chrNumber, window.start - maxReadLength);
		}
		return true;
	} else
		return false;
}

bool TAlignmentParser::moveWindow(TWindow & window){
	//returns false when end of genome is reached
	if(windowsPredefined){
		//if at beginning of BAM file
		if(chrIterator == bamHeader.Sequences.End()){
			restartChromosomes(window);

			if(chrIterator == bamHeader.Sequences.End())
				throw "found no predefined windows in BED file!";

		} else {
			//now move coordinates of next window
			if(!moveToNextPredefinedWindow(window)){
				//no more windows left on chr
				++chrIterator;
				++chrNumber;

				if(chrIterator == bamHeader.Sequences.End())
					return false;

				moveChromosome(window);


				if(chrIterator == bamHeader.Sequences.End()){
	//				throw "after move chromosome, about to return false for move window";
					return false;
				}
				std::cout << "incrementing windows c" << std::endl;

				++windowNumber;
			}
		}

	} else {
		//if at beginning of BAM file
		if(chrIterator == bamHeader.Sequences.End()){
			restartChromosomes(window);
			std::cout << "####### at beginning of BAM file" << std::endl;
		}
		else {
			if(!moveToNextWindowOnChr(window)){
				//there is no window left on chr
				++chrIterator;
				++chrNumber;

				//do we use this chromosome? if not, move on!
				while(chrIterator != bamHeader.Sequences.End() && !useChromosome[chrNumber]){
					++chrIterator;
					++chrNumber;
				}

				//did we reach end?
				if(chrIterator == bamHeader.Sequences.End() || chrNumber >= limitChr){
					window.end = 0;
//					chrIterator = bamHeader.Sequences.End();
					return false;
				}

				moveChromosome(window);
			}
		}
	}

	//report
	logfile->number("Window [" + toString(window.start) + ", " + toString(window.end) + ") of " + toString(numWindowsOnChr) + " on '" + chrIterator->Name + "':");
	logfile->addIndent();
	return true;

}

//------------------------------
//reading data
//------------------------------
bool TAlignmentParser::readAlignment(){
	bool filtersPassed;
	do {
		if(!bamReader.GetNextAlignment(bamAlignment))
			return false;

		//check if bam file is sorted
		if(bamAlignment.RefID != previousAlignmentChr){
			previousAlignmentPos = -1;
			previousAlignmentChr = bamAlignment.RefID;
		}
		if(bamAlignment.Position < previousAlignmentPos)
			throw "BAM file must be sorted by position!";
		previousAlignmentPos = bamAlignment.Position;

		//check read length
		if(bamAlignment.AlignedBases.size() > maxReadLength)
			throw "Alignment '" +  bamAlignment.Name + "' is longer than the max read length! Please change max read length to parse this data.";

		//store read group ID
		std::string readGroup;
		bamAlignment.GetTag("RG", readGroup);
		curReadGroupID = readGroups.find(readGroup);

		//filter
		//TODO: add functionality to not filter at all (i.e. _keepAll switch)
		filtersPassed = true;
		//check if insert size is shorter than read, this means we are reading the adaptor sequence
		//TODO: should add insertions to bamAlignment.AlignedBases.length()
		if(bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) <= bamAlignment.AlignedBases.length()){;// + alignment.numInsertions)){
			logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
			filtersPassed = false;
		} else {
			//apply filters: read group in use and basic QC
			filtersPassed = readGroups.readGroupInUse(curReadGroupID)
							&& bamAlignment.IsMapped() && !bamAlignment.IsFailedQC()
							&& bamAlignment.IsPrimaryAlignment()
							&& !bamAlignment.IsSupplementary()
							&& (_keepDuplicates || !bamAlignment.IsDuplicate());
		}
	} while(!filtersPassed);

	return true;
}


void TAlignmentParser::fillAlignment(TAlignment & alignment){
	//make sure container is empty
	alignment.clear(); //necessary? Not do it in TWindow?

	//fill alignment
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
	int readGroupId = readGroups.find(readGroup);

	alignment.fill(bamAlignment, readGroupId);

	if(parse){
		//add all info from bamAlignment to bases
		alignment.parse(genoMap, qualMap);

		//add missing information to bases
		alignment.fillReadGroupInfo(readGroupId);
		alignment.fillPmdProbabilities(pmdObjects);

		if(doRecalibration)
			recalibrate(alignment);
		if(hasReference)
			fillReferenceSequence(fastaBuffer, alignment);
		if(applyQualityFilter)
			alignment.filterForBaseQuality(minQual, maxQual);
	}
}


bool TAlignmentParser::readNextAligment(TAlignment & alignment){
	//use this in TGenome for functionalities that don't need windows
	if(readAlignment()){
		fillAlignment(alignment);
		return true;
	}
	return false;
}

//---------------------
//read data in windows
//---------------------
bool TAlignmentParser::readDataInNextWindow(TWindow & window){
	setParsingToTrue();

	//move window
	if(!moveWindow(window)){
		return false;
	}

	//read data
	readAlignmentsIntoWindow(window);

	return true;
};


void TAlignmentParser::readAlignmentsIntoWindow(TWindow & window){
	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);
	logfile->listFlush("Reading data ...");

	//check if old alignment is to be used.
	if(oldAlignmentMustBeConsidered){
		if(bamAlignment.Position >= window.end)
			return;

		oldAlignmentMustBeConsidered = false;
		if(oldAlignment->lastPositionPlusOne > window.start)
			oldAlignment = window.swapUsedForEmptyAlignment(oldAlignment, maxReadLength);
	}

	//read alignments
	int counter = 0;
	while(readAlignment()){
		//fill alignment
		fillAlignment(*oldAlignment);

		//check if alignment is used in current window
		if(oldAlignment->position >= window.end || oldAlignment->chrNumber != window.chrNumber){
			oldAlignmentMustBeConsidered = true;
			break;
		}

		//check if alignment end is after window start
		if(oldAlignment->lastPositionPlusOne > window.start){
			oldAlignment = window.swapUsedForEmptyAlignment(oldAlignment, maxReadLength);
		}
	}

	//fill sites
	window.fillSites();
	if(hasReference) window.addReferenceBaseToSites(*fastaReference, previousAlignmentChr);

	//report
	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");

	//apply filters
	applyFilters(window);
};


void TAlignmentParser::applyFilters(TWindow & window){
	window.passedFilters = false;
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
			return;
		}
		if(maxRefN < 1.0 && hasReference == true){
			logfile->conclude(toString(window.fractionRefIsN * 100) + "% of all reference bases are 'N'");
			if(window.fractionRefIsN > maxRefN){
				logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
				return;
			}
		}
		window.passedFilters = true;
	} else {
		logfile->conclude("No data in this window.");
	}
}

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
	std::string task = params.getParameterString("task");
	if(task == "qualityTransformation" || task == "recalBAM"){
		//don't change base error rates in fillAlignment!
		doRecalibration = false;
		initializeRecalibrationForQualityTransformation(params);
	} else {
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
}

void TAlignmentParser::initializeRecalibrationForQualityTransformation(TParameters & params){
	//are we comparing to original quality or to another recalibration?
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

	//check if estimation is required, in which case throw an error!
	if(recalObject->requiresEstimation()) throw "Can not use provided recalibration: estimation is required!";
}

void TAlignmentParser::recalibrate(TAlignment & alignment){
	//make sure read is parsed and has reference
	if(!alignment.parsed) throw "Read was not parsed!";

	if(recalObject->recalibrationChangesQualities()){
		//recalibrate quality scores
		for(int d=0; d<alignment.length; ++d){
			if(alignment.bases[d].aligned){
				alignment.bases[d].errorRate = recalObject->getErrorRate(alignment.bases[d]);
			}
		}

		alignment.changed = true;
	} else alignment.changed = false;
	alignment.recalibrated = true;
};

void TAlignmentParser::recalibrateWithPMD(TAlignment & alignment){
	//make sure read is parsed and has reference
	if(!alignment.parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	for(int d=0; d<alignment.length; ++d){
//		int k = length - d - 1;
		if(alignment.bases[d].aligned){
			//recalibrate quality scores
			if(recalObject->recalibrationChangesQualities())
				alignment.bases[d].errorRate = recalObject->getErrorRate(alignment.bases[d]);

			if(alignment.bases[d].context > 20) throw "there is a invalid context in alignment " + alignment.alignmentName + " at position " + toString(d);

			//now add effect of PMD
			if(alignment.bases[d].base == T && referenceSequence[alignment.bases[d].alignedPos] == 'C')
				alignment.bases[d].errorRate = 1.0 - ((1.0 - alignment.bases[d].errorRate)*(1.0 - alignment.bases[d].PMD_CT)); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
			else if(alignment.bases[d].base == A && referenceSequence[alignment.bases[d].alignedPos] == 'G')
				alignment.bases[d].errorRate = 1.0 - ((1.0 - alignment.bases[d].errorRate)*(1.0 - alignment.bases[d].PMD_GA)); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
		} else {
			alignment.bases[d].errorRate = qualMap.qualityToErrorMap[alignment.qualityOriginal[d]];
		}
	}

	//set pointer to recalibrated scores
//	alignment.quality = alignment.qualityRecalibrated;
	alignment.recalibrated = true;
	alignment.changed = true;
};

void TAlignmentParser::addSitesToQualityTransformTable(TAlignment & alignment, TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile){
	for(int i=0; i<alignment.length; ++i){
		QTtables.at(alignment.readGroupId)->add(alignment.qualityOriginal[i], qualMap.errorToQuality(recalObject->getErrorRate(alignment.bases[i])));
		QTtables.at(QTtables.size() - 1)->add(alignment.qualityOriginal[i], qualMap.errorToQuality(recalObject->getErrorRate(alignment.bases[i])));
	}
}

void TAlignmentParser::addSitesToQualityTransformTable(TAlignment & alignment, TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile){
	for(int i=0; i<alignment.length; ++i){
		QTtables.at(alignment.readGroupId)->add(qualMap.errorToQuality(recalObject->getErrorRate(alignment.bases[i])), qualMap.errorToQuality(otherRecalObject->getErrorRate(alignment.bases[i])));
		QTtables.at(QTtables.size() - 1)->add(qualMap.errorToQuality(recalObject->getErrorRate(alignment.bases[i])), qualMap.errorToQuality(otherRecalObject->getErrorRate(alignment.bases[i])));

	}
}


