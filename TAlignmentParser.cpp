/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"



//------------------------------
//reading alignments
//------------------------------
void TAlignmentParser::_fillAlignment(BAM::TAlignment & alignment){
	//This functions parses an alignment into bases
	bamFile.fill(alignment);

	//add all info from bamAlignment to bases
	alignment.parse(genoMap, qualMap, genotypeLikelihoodCalculator.getSequencingErrorModels());

	if(trimReads)
		alignment.trimRead(trimmingLength3Prime, trimmingLength5Prime);
	if(applyQualityFilter)
		alignment.filterForBaseQualityAsPhredInt(minQualityAsPhredInt, maxQualityAsPhredInt);
	if(applyContextFilter)
		alignment.filterForContext(ignoreTheseContexts, genoMap);
	if(hasReference)
		alignment.addReference(*fastaBuffer);
};

bool TAlignmentParser::readNextAlignment(BAM::TAlignment & alignment){
	if(!bamFile.readNextAlignmentThatPassesFilters()){
		return false;
	}
	_fillAlignment(alignment);
	return true;
};

//-----------------------------------------------------
//TAlignmentParserWindows
//-----------------------------------------------------
TAlignmentParser::TAlignmentParser(){
	logfile = nullptr;
	bamFile = nullptr;
	genotypeLikelihoodCalculator = nullptr;
	_parse = false;
	oldAlignment = NULL;
	oldAlignmentInitialized = false;
	oldAlignmentMustBeConsidered = false;
	hasWindowIndent = false;

	curReadGroupID = -1;

	//window parameters
	windowNumber = -1;
	windowSize = 0;
	numWindowsOnChr = -1;
	maxMissing = 1.0;
	maxRefN = 1.0;
	windowsPredefined = false;
	predefinedWindows = NULL;
	sitesProvided = false;

	//masks
	doMasking = false;
	considerRegions = false;
	mask = NULL;

	//filters
	applyQualityFilter = false;
	minQualityAsPhredInt = 0;
	maxQualityAsPhredInt = 93;
	applyContextFilter = false;
	trimReads = false;
	trimmingLength3Prime = 0;
	trimmingLength5Prime = 0;
	applyDepthFilter = false;
	readUpToDepth = 10000;
	minDepth = 0;
	maxDepth = 10000;

	//blacklist
	_updateBlacklist = false;

	//limit chr and windows
	limitWindows = -1;
	skipWindows = 0;
	doLimitReads = false;
	limitReads = -1;

	//reference
	hasReference = false;
	fastaBuffer = NULL;
	chrChangedWindow = false;
};

TAlignmentParser::~TAlignmentParser(){
	if(doMasking)
		delete mask;
	if(windowsPredefined)
		delete predefinedWindows;
	if(subset)
		delete subset;
	if(oldAlignmentInitialized){
		delete oldAlignment;
	}
};


void TAlignmentParser::_setWindowParameters(TParameters & params){
	if(!params.parameterExists("window") && params.parameterExists("windows")) logfile->warning("Argument 'windows' specified, but unknown. Did you mean 'window'?");
	std::string tmp = params.getParameterStringWithDefault("window", "1000000");
	//check if it is a number
	if(stringContainsOnly(tmp, "1234567890.Ee-+")){
		windowsPredefined = false;
		windowSize = stringToInt(tmp);
		logfile->list("Setting window size to " + toString(windowSize) + ". (parameter 'window')");
		if(windowSize < bamFile->maxReadLength())
			throw "Window size " + tmp + " out of range! Windows must be at least as large as the max read length (" + toString(bamFile->maxReadLength()) + " bp)!";
	} else {
		windowsPredefined = true;
		logfile->listFlush("Limiting analysis to windows defined in '" + tmp + "'...");
		predefinedWindows = new TBed(tmp);
		logfile->done();
		logfile->conclude("Read " + toString(predefinedWindows->size()) + " of cumulative length " + toString(predefinedWindows->length()) + " bp on " + toString(predefinedWindows->getNumChromosomes()) + " chromosomes.");
	}
	numWindowsOnChr = 0;
};

void TAlignmentParser::_setWindowFilters(TParameters & params){
	//filter for missing reference
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 1.0);
	if(maxMissing > 1.0) throw "maxMissing must be smaller or equal to 1.0!";
	logfile->list("Will filter out windows with a missing data fraction > " + toString(maxMissing) + ". (parameter 'maxMissing')");

	maxRefN = params.getParameterDoubleWithDefault("maxRefN", 1.0);
	if(maxRefN > 1.0) throw "maxRefN must be smaller or equal to 1.0!";
	if(maxRefN < 1.0 && hasReference == false) throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided.";
	logfile->list("Will filter out windows with a fraction of 'N' in reference > " + toString(maxMissing) + ". (parameter 'maxRefN')");
};

void TAlignmentParser::_setSiteFilters(TParameters & params){
	//depth filter
	readUpToDepth = params.getParameterIntWithDefault("readUpToDepth", 1000);
	if(params.parameterExists("minDepth") || params.parameterExists("maxDepth")){
		applyDepthFilter = true;
		unsigned int tmpInt;
		tmpInt = params.getParameterIntWithDefault("minDepth", 0);
		if(tmpInt < 0)
			throw "minDepth must be >= 0!";
		minDepth = tmpInt;
		tmpInt = params.getParameterIntWithDefault("maxDepth", readUpToDepth);
		if(tmpInt < minDepth) throw "maxDepth must be >= minDepth!";
		maxDepth = tmpInt;
		readUpToDepth = maxDepth + 1;
		logfile->list("Will filter out sites with sequencing depth < " + toString(minDepth) + " or > " + toString(maxDepth) + ". (parameters 'minDepth', 'maxDepth')");
	} else {
		applyDepthFilter = false;
		minDepth = 0;
		maxDepth = readUpToDepth;
	}
	logfile->list("Will read data up to depth " + toString(readUpToDepth) + " and ignore additional bases. (parameter 'readUpToDepth')");
};


void TAlignmentParser::_setMasks(TParameters & params){
	//normal mask
	if(params.parameterExists("mask")){
		if(windowsPredefined) throw "Masking is currently not implemented if windows are predefined from a BED file.";
		if(params.parameterExists("alleles")) throw "Masking is currently not implemented if variant positions are also specified with 'sites'";
		if(params.parameterExists("regions")) throw "Cannot use mask and regions at the same time";
		doMasking = true;
		std::string maskFile = params.getParameterString("mask");

		//limit sites?
		int siteLimit = -1;
		if(params.parameterExists("siteLimit")){
			siteLimit = params.getParameterInt("siteLimit");
			if(siteLimit < 0)
				throw "site limit cannot be smaller than 0!";
			logfile->startIndent("Will mask the first " + toString(siteLimit) + " sites listed in BED file '" + maskFile + "':");
		} else {
			logfile->startIndent("Will mask all sites listed in BED file '" + maskFile + "':");
		}
		logfile->listFlush("Reading file ...");
		mask = new BAM::TBedReader(maskFile, windowSize, bamFile->chromosomes, siteLimit, logfile);
		logfile->done();
		logfile->endIndent();
		//mask->print();
	} else doMasking = false;

	//reverse masking
	if(params.parameterExists("regions")){
		if(windowsPredefined) throw "Regions is currently not implemented if windows are predefined from a BED file.";
		if(params.parameterExists("alleles")) throw "Regions is currently not implemented if variant positions are also specified with \"sites\"";
		considerRegions = true;
		std::string regionsFile = params.getParameterString("regions");

		//limitSites
		int siteLimit = -1;
		if(params.parameterExists("siteLimit")){
			siteLimit = params.getParameterInt("siteLimit");
			if(siteLimit < 0)
				throw "site limit cannot be smaller than 0!";
			logfile->startIndent("Will limit analysis to the first " + toString(siteLimit) + " sites listed in BED file '" + regionsFile + "':");
		} else {
			logfile->startIndent("Will limit analysis to all regions listed in BED file '" + regionsFile + "' (parameter 'regions'):");
		}
		logfile->listFlush("Reading file ...");
		mask = new BAM::TBedReader(regionsFile, windowSize, bamFile->chromosomes, siteLimit, logfile);
		logfile->done();
		logfile->endIndent();
	} else considerRegions = false;
};

void TAlignmentParser::_setParsingLimits(TParameters & params){
	//limit windows
	skipWindows = params.getParameterIntWithDefault("skipWindows", 0);
	if(skipWindows > 0) logfile->list("Will skip the first " + toString(skipWindows) + " windows per chromosome. (parameter 'skipWindows')");
	limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows"))
		logfile->list("Will limit analysis to the first " + toString(limitWindows) + " windows per chromosome. (parameter 'limitWindows')");
	if(limitWindows <= skipWindows)
		throw "limitWindows has to be larger than skipWindows!";
};


void TAlignmentParser::setUpdateBlacklistToTrue(){
	_updateBlacklist = true;
	logfile->list("Storing ignored reads in a blacklist.");
};

void TAlignmentParser::setWriteBlacklistToFileToTrue(const std::string filename){
	logfile->list("Writing ignored reads to '" + filename + "'.");
	blacklist.enableWriting(filename);
};



//--------------
//move genome
//--------------

void TAlignmentParser::_jumpToEnd(){
	bamFile->chromosomes.jumpToEnd();
};

void TAlignmentParser::_restartChromosomes(TWindow_base & window){
	bamFile->chromosomes.begin();

	_moveChromosome(window);
};

void TAlignmentParser::_moveChromosome(TWindow_base & window){
	//get reference to chromosome to declutter code
	BAM::TChromosomes& chromosomes = bamFile->chromosomes;

	//jump reader
	oldAlignmentMustBeConsidered = false;

	//restart windows
	if(hasWindowIndent){
		logfile->removeIndent();
		hasWindowIndent = false;
	}
	windowNumber = 0;

	if(windowsPredefined){
		//find next used chromosome with windows
		do {
			predefinedWindows->setChr(chromosomes.curName());
			numWindowsOnChr = predefinedWindows->getNumWindowsOnCurChr();
			if(numWindowsOnChr < 1 || chromosomes.curInUse() == false){
				if(chromosomes.curInUse())
					logfile->conclude("No windows on chromosome " + chromosomes.curName() + ".");
				chromosomes.next();
				if(chromosomes.end()){
					return;
				}

				predefinedWindows->setChr(chromosomes.curName());
				numWindowsOnChr = predefinedWindows->getNumWindowsOnCurChr();
			}
		} while((numWindowsOnChr < 1 || !chromosomes.curInUse()) && !chromosomes.end());

		//now jump
		window.move(predefinedWindows->curWindowStart(), predefinedWindows->curWindowEnd(), chromosomes.curIndex(), logfile);
		window._chrName = chromosomes.curName();
		bamFile->jump(chromosomes.curIndex(), window.startPos);

	} else {
		while(chromosomes.curInUse() == false || skipWindows * windowSize > chromosomes.curLength()){
			chromosomes.next();
		}
		numWindowsOnChr = ceil(chromosomes.curLength() / (double) windowSize);

		uint32_t curStart = skipWindows * windowSize;
		bamFile->jump(chromosomes.curIndex(), curStart);
		uint32_t nextEnd = curStart + windowSize;

		if(nextEnd > chromosomes.curLength()){
			nextEnd = chromosomes.curLength();
		}
		window.move(curStart, nextEnd, chromosomes.curIndex(), logfile);
		window._chrName = chromosomes.curName();
	}

	if(chromosomes.end())
		return;

	//advance mask
	if(doMasking || considerRegions) mask->setChr(chromosomes.curName());
	if(sitesProvided) subset->setChr(chromosomes.curName());

	//write progress
	if(chromosomes.curIndex() > 0) logfile->endIndent();
	logfile->startNumbering("Parsing chromosome '" + chromosomes.curName() + "':");
};

bool TAlignmentParser::_moveToNextWindowOnChr(TWindow_base & window){
	//get reference to chromosome to declutter code
	TChromosomes& chromosomes = bamFile->chromosomes;

	//if sites defined
	int counter = 0;
	do{
		//move possible?
		++windowNumber;
		++counter;
	} while(sitesProvided && !subset->hasPositionsInWindow(window.endPos) && window.endPos + window.length * counter < chromosomes.curLength());

	if(window.endPos >= chromosomes.curLength() || windowNumber >= limitWindows)
		return false;

	//calculate new end
	long nextEnd = window.endPos + windowSize;
	if(nextEnd > chromosomes.curLength())
		nextEnd = chromosomes.curLength();
	window.move(window.endPos, nextEnd, chromosomes.curIndex(), logfile);

	return true;
};

bool TAlignmentParser::_moveToNextPredefinedWindow(TWindow_base & window){
	//get reference to chromosome to declutter code
	TChromosomes& chromosomes = bamFile->chromosomes;

	++windowNumber;
	if(windowNumber >= limitWindows)
		return false;
	if(predefinedWindows->nextWindow()){
		window.move(predefinedWindows->curWindowStart(), predefinedWindows->curWindowEnd(), chromosomes.curIndex(), logfile);

		//should we jump or are we already close enough to next window
		if(bamFile->curPosition() > window.startPos || bamFile->curPosition() < window.startPos - maxReadLength){
			if(window.startPos < maxReadLength)
				bamFile->jump(chromosomes.curIndex(), 0);
			else{
				bamFile->jump(chromosomes.curIndex(), window.startPos - maxReadLength);
			}
		}
		return true;
	} else
		return false;
};

bool TAlignmentParser::_moveWindow(TWindow_base & window){
	//get reference to chromosome to declutter code
	TChromosomes& chromosomes = bamFile->chromosomes;

	//returns false when end of genome is reached
	if(windowsPredefined){
		//if at beginning of BAM file
		if(chromosomes.endPos()){
			_restartChromosomes(window);

			if(chromosomes.endPos())
				throw "found no predefined windows in BED file! Does file exist?";
			chrChangedWindow = true;

		} else {
			//now move coordinates of next window
			if(!_moveToNextPredefinedWindow(window)){
				//no more windows left on chr
				chromosomes.next();

				if(chromosomes.endPos()){
					if(hasWindowIndent){
						logfile->removeIndent();
						hasWindowIndent = false;
					}
					return false;
				}

				_moveChromosome(window);
				chrChangedWindow = true;

				if(chromosomes.endPos()){
					if(hasWindowIndent){
						logfile->removeIndent();
						hasWindowIndent = false;
					}
					return false;
				}
				++windowNumber;
			} else
				//was able to move to next window on chr
				chrChangedWindow = false;
		}

	} else {
		//if at beginning of BAM file
		if(chromosomes.endPos()){
			_restartChromosomes(window);
			chrChangedWindow = true;
		} else {
			if(!_moveToNextWindowOnChr(window)){
				//there is no window left on chr
				chromosomes.next();

				//do we use this chromosome? if not, move on!
				while(!chromosomes.endPos() && !chromosomes.curInUse()){
					chromosomes.next();
				}

				//did we reach end?
				if(chromosomes.endPos()){
					window.endPos = 0;
					if(hasWindowIndent){
						logfile->removeIndent();
						hasWindowIndent = false;
					}
					return false;
				}
				_moveChromosome(window);
				chrChangedWindow = true;
			} else {
				chrChangedWindow = false;
			}
		}
	}

	//report
	if(hasWindowIndent) logfile->removeIndent();
	logfile->number("Window [" + toString(window.startPos) + ", " + toString(window.endPos) + ") of " + toString(numWindowsOnChr) + " on '" + chromosomes.curName() + "':");
	logfile->addIndent();
	hasWindowIndent = true;

	return true;
};

//---------------------
//read data in windows
//---------------------
bool TAlignmentParser::readDataInNextWindow(TWindow & window){
	setParsingToTrue();

	//move window
	if(!_moveWindow(window)){
		return false;
	}

	//read data
	_readAlignmentsIntoWindow(window);
	return true;
};

void TAlignmentParser::_readAlignmentsIntoWindow(TWindow & window){
	//measure runtime
	logfile->listFlushTime("Reading data ...");

	//check if old alignment is to be used.
	if(oldAlignmentMustBeConsidered){
		if(bamFile->curPosition() >= window.endPos){
			logfile->warning("Old alignment is after window!");
			return;
		}

		oldAlignmentMustBeConsidered = false;
		if(oldAlignment->lastAlignedPositionWithRespectToRef >= window.startPos)
			oldAlignment = window.swapUsedForEmptyAlignment(oldAlignment, maxReadLength);
	}

	//read alignments
	int counter = 0;
	while(_readAlignment()){
		//fill alignment
		//return false if a post-parsing filer is not passed
		if(!_fillAlignment(*oldAlignment)){
			continue;
		}

		++counter;

		//check if alignment starts after current window end -> break
		if(oldAlignment->position >= window.endPos || oldAlignment->refID != window.chrNumber){
			oldAlignmentMustBeConsidered = true;
			break;
		}

		//check if alignment contains part of the window
		//if read continues outside of window, this is dealt with by window object
		if(oldAlignment->position >= window.startPos || oldAlignment->lastAlignedPositionWithRespectToRef >= window.startPos){
			oldAlignment = window.swapUsedForEmptyAlignment(oldAlignment, maxReadLength);
		}
	}

	//fill sites
	if(sitesProvided){
		window.fillSitesSubset(subset, readUpToDepth);
		window.addReferenceBaseToSites(subset);
	} else {
		window.fillSites(readUpToDepth);
		if(hasReference) window.addReferenceBaseToSites(*fastaReference);
	}

	//report
	logfile->doneTime();

	//apply filters
	_applyWindowFilters(window);
};

void TAlignmentParser::downsampleWindow(TWindow_base & destination, TWindow & source, const double downsamplingProb, TRandomGenerator* randomGenerator){
	//downsample
	logfile->listFlush("Downsampling reads with p = " + toString(downsamplingProb) + " ...");
	if(sitesProvided)
		destination.downsampleFromOther(source, subset, readUpToDepth, downsamplingProb, randomGenerator);
	else
		destination.downsampleFromOther(source, readUpToDepth, downsamplingProb, randomGenerator);
	logfile->done();

	//apply filters
	_applyWindowFilters(destination);
};

void TAlignmentParser::_applyWindowFilters(TWindow_base & window){
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
		} if(applyDepthFilter){
			window.applyDepthFilter(minDepth, maxDepth);
		} if(maxRefN < 1.0 && hasReference == true){
			window.calcFracN();
		}

		//calc sequencing depth
		window.calcDepth();

		//report
		logfile->conclude("read data from " + toString(window.numReadsInWindow) + " reads.");
		logfile->conclude("sequencing depth is " + toString(window.depth));
		logfile->conclude(toString(window.fractionDepthAtLeastTwo * 100) + "% of all sites are covered at least twice");
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
};


