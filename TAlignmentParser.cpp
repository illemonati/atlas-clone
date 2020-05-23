/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"



//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
TAlignmentParser::TAlignmentParser(){
	logfile = nullptr;
	bamFile = nullptr;
	genotypeLikelihoodCalculator = nullptr;
	_keepAll = false;
	_keepDuplicates = false;
	_filterSoftClips = false;
	_keepImproperPairs = false;
	_keepUnmappedReads = false;
	_keepFailedQC = false;
	_keepSecondary = false;
	_keepSupplementary = false;
	applyFragmentLengthFilter = false;
	minFragmentLength = -1;
	maxFragmentLength = -1;
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
	maxReadLength = 0;
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
	applyDepthFilter = false;
	readUpToDepth = 10000;
	minDepth = 0;
	maxDepth = 10000;
	applyQualityFilter = false;
	applyMQFilter = false;
	applyContextFilter = false;
	minMQ = 0;
	maxMQ = 10000;
	minQualityAsPhredInt = 0;
	maxQualityAsPhredInt = 93;
	trimReads = false;
	trimmingLength3Prime = 0;
	trimmingLength5Prime = 0;
	_keepReadsLongerThanInsertSize = false;
	useStrand[0] = true; useStrand[1] = true;
	useMate[0] = true; useMate[1] = true;

	//blacklist
	_updateBlacklist = false;
	_writeBlackList = false;

	//limit chr and windows
	limitWindows = -1;
	skipWindows = 0;
	doLimitReads = false;
	limitReads = -1;

	//reference
	hasReference = false;
	fastaReference = NULL;
	fastaBuffer = NULL;
	chrChangedWindow = false;
};

TAlignmentParser::TAlignmentParser(TParameters & params, TBamFile * BamFile, GenotypeLikelihoods::TGenotypeLikelihoodCalculator* GenotypeLikelihoodCalculator, TLog* Logfile){
	TAlignmentParser();

	init(params,BamFile, GenotypeLikelihoodCalculator, Logfile);
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

	if(_writeBlackList)
		ignoredReads.close();
};

void TAlignmentParser::init(TParameters & params, TBamFile * BamFile, GenotypeLikelihoods::TGenotypeLikelihoodCalculator* GenotypeLikelihoodCalculator, TLog* Logfile){
	logfile = Logfile;
	bamFile = BamFile;
	genotypeLikelihoodCalculator = GenotypeLikelihoodCalculator;

	//alignments
	oldAlignment = new TAlignment(bamFile->maxReadLength());
	oldAlignmentInitialized = true;

	//settings
	_setWindowParameters(params);
	_setParsingLimits(params);

	_setMasks(params);
	_setFilters(params);
	_setChrPloidy(params);
};

/*
void TAlignmentParser::_openBamFile(std::string filename, bool indexNotRequired){
	//open BAM file
	logfile->list("Reading data from BAM file '" + filename + "'.");
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";
	//load index file
	if(!bamReader.LocateIndex() && !indexNotRequired)
		throw "No index file found for BAM file '" + filename + "'!";

	//initialize bam stuff
	bamHeader = bamReader.GetHeader();

	//initialize chromosomes
	chromosomes = TChromosomes(&bamHeader);

	//get file size
	chromosomes.jumpToBeginningOfLastChr();
	bamReader.Jump(bamHeader.Sequences.Size() - 1, 0);
	BamTools::BamAlignment bamAlignment;
	bamReader.GetNextAlignment(bamAlignment);
	sizeOfBamFile = bamReader.tell();
	bamReader.Rewind();

	chromosomes.jumpToEnd();
};
*/

void TAlignmentParser::setOutName(std::string outputName){
	outname = outputName;
}

void TAlignmentParser::_setWindowParameters(TParameters & params){
	if(!params.parameterExists("window") && params.parameterExists("windows")) logfile->warning("Argument 'windows' specified, but unknown. Did you mean 'window'?");
	std::string tmp = params.getParameterStringWithDefault("window", "1000000");
	//check if it is a number
	if(stringContainsOnly(tmp, "1234567890.Ee-+")){
		windowsPredefined = false;
		windowSize = stringToInt(tmp);
		logfile->list("Setting window size to " + toString(windowSize) + ". (parameter 'window')");
		if(windowSize < maxReadLength)
			throw "Window size " + tmp + " out of range! Windows must be at least as large as the max read length (" + toString(maxReadLength) + " bp)!";
	} else {
		windowsPredefined = true;
		logfile->listFlush("Limiting analysis to windows defined in '" + tmp + "'...");
		predefinedWindows = new TBed(tmp);
		logfile->done();
		logfile->conclude("Read " + toString(predefinedWindows->size()) + " of cumulative length " + toString(predefinedWindows->length()) + " bp on " + toString(predefinedWindows->getNumChromosomes()) + " chromosomes.");
	}
	numWindowsOnChr = 0;
};

void TAlignmentParser::_setFilters(TParameters & params){
	_setWindowFilters(params);

	if(params.parameterExists("keepAllReads")){
		_setAlignmentFiltersToKeepAll();
		logfile->list("Will keep all reads. All filters turned off (parameter 'keepAll')");
	} else {
		_setAlignmentFilters(params);
	}

	_setSiteFilters(params);

	_setBaseFilters(params);
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

void TAlignmentParser::_setContextFilter(std::vector<std::string> contexts){
	applyContextFilter = true;
	for(unsigned int i=0; i < contexts.size(); i++){
		if(contexts[i].size() != 2)
			throw "Context " + contexts[i] + " does not have a length of 2! (parameter 'maskContext')";
		BaseContext co = genoMap.getContext(contexts[i][0], contexts[i][1]);
		ignoreTheseContexts.emplace(co, 1);
	}
	logfile->startIndent("Will mask the following contexts (parameter 'maskContext'):");
	std::map<BaseContext, int>::iterator it;
	for(it = ignoreTheseContexts.begin(); it != ignoreTheseContexts.end(); ++it){
		logfile->list(genoMap.getContextString(it->first));
	}
	logfile->endIndent();
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

void TAlignmentParser::_setBaseFilters(TParameters & params){
	//quality filters
	_setQualityFilters(params.getParameterIntWithDefault("minQual", 1), params.getParameterIntWithDefault("maxQual", 93));

	//quality filters for printing
	_setQualityRangeForPrinting(params.getParameterIntWithDefault("minOutQual", 0), params.getParameterIntWithDefault("maxOutQual", 93));

	//context filter
	if(params.parameterExists("ignoreContexts")){
		std::vector<std::string> contexts;
		fillVectorFromString(params.getParameterString("ignoreContexts"), contexts, ',');
		_setContextFilter(contexts);
	}
};

void TAlignmentParser::_setQualityFilters(int MinPhredInt, int MaxPhredInt){
	minQualityAsPhredInt = MinPhredInt;
	maxQualityAsPhredInt = MaxPhredInt;
	if(minQualityAsPhredInt < 0) throw "minQual must be >= 0!";
	if(maxQualityAsPhredInt < minQualityAsPhredInt) throw "maxQual must be >= minQual!";
	applyQualityFilter = true;

	logfile->list("Will filter out bases with quality outside the range [" + toString(minQualityAsPhredInt) + ", " + toString(maxQualityAsPhredInt) + "] (parameters 'minQual', 'maxQual')");
};

void TAlignmentParser::_setQualityRangeForPrinting(int minQual, int maxQual){
	logfile->list("Will print qualities truncated to [" + toString(minQual) + ", " + toString(maxQual) + "] (parameters 'minOutQual', 'maxOutQual')");
	if(minQual < 0 || minQual > 255) throw "minOutQual " + toString(minQual) + " is outside accepted range [0, 255]!";
	if(maxQual < 0 || maxQual > 255) throw "maxOutQual " + toString(minQual) + " is outside accepted range [0, 255]!";
	if(maxQual < minQual) throw "maxOutQual must be >= minOutQual!";

	//set in quality map
	qualMap.setQualityLimits(minQual, maxQual);
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
		mask = new TBedReader(maskFile, windowSize, bamFile->chromosomes, siteLimit, logfile);
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
		mask = new TBedReader(regionsFile, windowSize, bamFile->chromosomes, siteLimit, logfile);
		logfile->done();
		logfile->endIndent();
	} else considerRegions = false;
};

void TAlignmentParser::_setFragmentLengthFilter(TParameters & params){
	int MinFragmentLength = params.getParameterIntWithDefault("minFragmentLength", 0);
	int MaxFragmentLength = params.getParameterIntWithDefault("maxFragmentLength", 255);
	if(MinFragmentLength < 0 || MinFragmentLength > 255)
		throw "minFragmentLength '" + toString(MinFragmentLength) + "' is outside the accepted range [0,255]!";
	if(MinFragmentLength < 0 || MinFragmentLength > 255)
		throw "maxFragmentLength '" + toString(MinFragmentLength) + "' is outside the accepted range [0,255]!";
	if(MinFragmentLength > MaxFragmentLength)
		throw "minFragmentLength must be <= maxFragmentLength!";

	minFragmentLength = MinFragmentLength;
	maxFragmentLength = MaxFragmentLength;
	applyFragmentLengthFilter = true;
	setParsingToTrue(); //filter requires parsing!
	logfile->list("Will filter out reads with fragment length outside the range [" + toString(MinFragmentLength) + ", " + toString(MaxFragmentLength) + "]. (parameters 'minFragmentLength', 'maxFragmentLength')");
};

void TAlignmentParser::_setParsingLimits(TParameters & params){
	//Limit chromosomes
	if(params.parameterExists("chr")){
		logfile->startIndent("Will limit analysis to the following chromosomes (parameter 'chr':");
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("chr"), vec, ',');
		bamFile->chromosomes.useSpecifiedChr(vec);
		bamFile->chromosomes.writeUsedChromosomes(logfile);
	} else {
		if(params.parameterExists("limitChr")){
			std::string limitName = params.getParameterString("limitChr");
			logfile->list("Will limit analysis to all chromosomes up to and including " + limitName + ". (parameter 'limitChr')");
			bamFile->chromosomes.limitChr(limitName);
		}
	}

	//limit windows
	skipWindows = params.getParameterIntWithDefault("skipWindows", 0);
	if(skipWindows > 0) logfile->list("Will skip the first " + toString(skipWindows) + " windows per chromosome. (parameter 'skipWindows')");
	limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows"))
		logfile->list("Will limit analysis to the first " + toString(limitWindows) + " windows per chromosome. (parameter 'limitWindows')");
	if(limitWindows <= skipWindows)
		throw "limitWindows has to be larger than skipWindows!";

	//limit reads
	if(params.parameterExists("limitReads")){
		doLimitReads = true;
		limitReads = params.getParameterLong("limitReads");
		logfile->list("Will limit analysis to " + toString(limitReads) + " reads.");
	}
};

void TAlignmentParser::_setReadTrimming(TParameters & params){
	//trimming ends
	if(params.parameterExists("trim3") || params.parameterExists("trim5")){
		trimmingLength3Prime = params.getParameterIntWithDefault("trim3", 0);
		if(trimmingLength3Prime < 0) throw "trimming distance trim3 must be >= 0!";
		trimmingLength5Prime = params.getParameterIntWithDefault("trim5", 0);
		if(trimmingLength5Prime < 0) throw "trimming distance trim5 must be >= 0!";
		if(trimmingLength3Prime > 0 || trimmingLength5Prime > 0){
			logfile->list("Will trim first " + toString(trimmingLength3Prime) + " and " + toString(trimmingLength5Prime) + " bases from the 3' and 5' end, respectively. (parameters 'trim3', 'trim5')");
		}
	}
	trimReads = true;
};


void TAlignmentParser::_setChrPloidy(TParameters & params){
	logfile->list("Chromosomes with no further specifications are assumed to be diploid (parameters 'ploidy' or 'haploid' to change ploidy).");
	if(params.parameterExists("ploidy")){
		std::string ploidyFileName = params.getParameterString("ploidy");
		logfile->list("Reading ploidy specification per chromosome from file '" + ploidyFileName + "'");
		std::ifstream ploidyFile(ploidyFileName.c_str());
		if(!ploidyFile)
			throw "Failed to open file '" + ploidyFileName + "'!";
		bamFile->chromosomes.specifyPloidy(ploidyFile, logfile);
	}

	if(params.parameterExists("haploid")){
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("haploid"), vec, ',');
		bamFile->chromosomes.setToHaploid(vec, logfile);
	}
};

void TAlignmentParser::addReference(TFastaBuffer* FastaBuffer){
	hasReference = true;
	fastaBuffer = FastaBuffer;
};

void TAlignmentParser::setUpdateBlacklistToTrue(){
	_updateBlacklist = true;
	logfile->list("Storing ignored reads in a blacklist.");
};

void TAlignmentParser::setWriteBlacklistToFileToTrue(){
	std::string filename = outname + "_ignoredReads.txt.gz";
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
	TChromosomes& chromosomes = bamFile->chromosomes;

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
		window.chrName = chromosomes.curName();
		bamFile->jump(chromosomes.curIndex(), window.start);

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
		window.chrName = chromosomes.curName();
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
	} while(sitesProvided && !subset->hasPositionsInWindow(window.end) && window.end + window.length * counter < chromosomes.curLength());

	if(window.end >= chromosomes.curLength() || windowNumber >= limitWindows)
		return false;

	//calculate new end
	long nextEnd = window.end + windowSize;
	if(nextEnd > chromosomes.curLength())
		nextEnd = chromosomes.curLength();
	window.move(window.end, nextEnd, chromosomes.curIndex(), logfile);

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
		if(bamFile->curPosition() > window.start || bamFile->curPosition() < window.start - maxReadLength){
			if(window.start < maxReadLength)
				bamFile->jump(chromosomes.curIndex(), 0);
			else{
				bamFile->jump(chromosomes.curIndex(), window.start - maxReadLength);
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
		if(chromosomes.end()){
			_restartChromosomes(window);

			if(chromosomes.end())
				throw "found no predefined windows in BED file! Does file exist?";
			chrChangedWindow = true;

		} else {
			//now move coordinates of next window
			if(!_moveToNextPredefinedWindow(window)){
				//no more windows left on chr
				chromosomes.next();

				if(chromosomes.end()){
					if(hasWindowIndent){
						logfile->removeIndent();
						hasWindowIndent = false;
					}
					return false;
				}

				_moveChromosome(window);
				chrChangedWindow = true;

				if(chromosomes.end()){
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
		if(chromosomes.end()){
			_restartChromosomes(window);
			chrChangedWindow = true;
		} else {
			if(!_moveToNextWindowOnChr(window)){
				//there is no window left on chr
				chromosomes.next();

				//do we use this chromosome? if not, move on!
				while(!chromosomes.end() && !chromosomes.curInUse()){
					chromosomes.next();
				}

				//did we reach end?
				if(chromosomes.end()){
					window.end = 0;
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
	logfile->number("Window [" + toString(window.start) + ", " + toString(window.end) + ") of " + toString(numWindowsOnChr) + " on '" + chromosomes.curName() + "':");
	logfile->addIndent();
	hasWindowIndent = true;

	return true;
};

//------------------------------
//reading alignments
//------------------------------
bool TAlignmentParser::_readAlignment(){
	bool filtersPassed = false;
	do {
		if(!bamFile->readNextAlignment()){
			return false;
		}
		if(doLimitReads && bamFile->numAlignmentsRead() > limitReads){
			return false;
		}

		//filter
		/* check if insert size is shorter than read-insertions+deletions=alignedBases.length() -> this means we are reading the adaptor sequence.
		Insert size is determined by mapping -> insertions are not in ref and should not count. If we don't add deletions, adapter at end could be sequenced but we still keep read
		(deletions in aligned bases are represented as dashes) */
		if(bamFile->curIsLongerThanInsertSize()){
			if(_updateBlacklist){
				addToBlacklist(bamFile->_curBamAlignment, "longer than insert size (TLEN)");
			} else {
				logfile->warning("The following alignment is longer than its insert size: " + bamFile->_curBamAlignment.Name);
			}
			filtersPassed = false;
		} else {
			//apply filters: read group in use and basic QC
			filtersPassed = bamFile->curPassedQC();
			if(_updateBlacklist && !filtersPassed){
				addToBlacklist(bamFile->_curBamAlignment, "did not pass QC filters");
			}
		}
	} while(!filtersPassed);

	return true;
};

bool TAlignmentParser::_fillAlignment(TAlignment & alignment){
	//This functions parses an alignment into bases
	bamFile->fill(alignment);

	if(_parse){
		//add all info from bamAlignment to bases
		alignment.parse(genoMap, qualMap, genotypeLikelihoodCalculator->getSequencingErrorModels());

		//check for fragment length filter
		if(applyFragmentLengthFilter && (alignment.fragmentLength < minFragmentLength || alignment.fragmentLength > maxFragmentLength))
			return false;

		if(trimReads)
			alignment._trimRead(trimmingLength3Prime, trimmingLength5Prime);
		if(applyQualityFilter)
			alignment.filterForBaseQualityAsPhredInt(minQualityAsPhredInt, maxQualityAsPhredInt);
		if(applyContextFilter)
			alignment.filterForContext(ignoreTheseContexts);
		if(hasReference)
			alignment.addReference(fastaBuffer);
	}

	return true;
};

//------------------------
//read data in alignments
//------------------------

bool TAlignmentParser::readNextAlignment(TAlignment & alignment){
	//use this in TGenome for functionalities that don't need windows
	if(_readAlignment()){
		return _fillAlignment(alignment);
	}
	return false;
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
		if(bamFile->curPosition() >= window.end){
			logfile->warning("Old alignment is after window!");
			return;
		}

		oldAlignmentMustBeConsidered = false;
		if(oldAlignment->lastAlignedPositionWithRespectToRef >= window.start)
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
		if(oldAlignment->position >= window.end || oldAlignment->refID != window.chrNumber){
			oldAlignmentMustBeConsidered = true;
			break;
		}

		//check if alignment contains part of the window
		//if read continues outside of window, this is dealt with by window object
		if(oldAlignment->position >= window.start || oldAlignment->lastAlignedPositionWithRespectToRef >= window.start){
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

void TAlignmentParser::addSitesToQualityTransformTable(TAlignment & alignment, TQualityTransformTables & QTtables){
	alignment.addSitesToQualityTransformTable(QTtables);
};

void TAlignmentParser::addSitesToQualityTransformTable(TAlignment & alignment, GenotypeLikelihoods::TSequencingErrorModels & otherSeqErrors, TQualityTransformTables & QTtables){
	alignment.addSitesToQualityTransformTable(otherSeqErrors, QTtables);
};

void TAlignmentParser::_adaptQualityWhenMerging(TBase & bestBase, TBase & worstBase, const bool & adaptQuality){
	if(adaptQuality){
		double likelihood[4];
		double sum = 0.0;
		for(int i=0; i<4; ++i){
			if(bestBase.data.base == i){
				likelihood[i] = 1.0 - bestBase.errorRate;
			} else {
				likelihood[i] = bestBase.errorRate / 3.0;
			}

			if(worstBase.data.base == i){
				likelihood[i] *= 1.0 - worstBase.errorRate;
			} else {
				likelihood[i] *= worstBase.errorRate / 3.0;
			}
			sum += likelihood[i];

		}
		bestBase.errorRate = 1.0 - likelihood[bestBase.data.base] / sum;
	} else {
		worstBase.data.recalibratedQualityAsPhredInt = 0.0;
//		worstBase.base = N;
	}
};

void TAlignmentParser::mergeAlignedBasesOneRead(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator){
	//deletions and insertions are kept as is. these positions are not compared
	if(fwdAlignment->lastAlignedPositionWithRespectToRef >= revAlignment->position){
		fwdAlignment->setAlignmentHasChanged();
		revAlignment->setAlignmentHasChanged();

		//which read to keep
		bool keepFwd = true;
		if(randomGenerator->getRand() < 0.5)
			keepFwd = false;

		//reads overlap -> check if there are bases overlapping same position in ref
		//alignedPos is with respect to read
		int fwdP = 0;
		int revP = 0;
		while(fwdP <= fwdAlignment->lastAlignedPos && revP <= revAlignment->lastAlignedPos){
			if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] == revAlignment->position + revAlignment->alignedPosition[revP]){
				//bases overlap same position in ref -> choose at random which to keep
				if(keepFwd){
					_adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					_adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] < revAlignment->position + revAlignment->alignedPosition[revP]){
				++fwdP;
			} else {
				++revP;
			}
		}
	}
}

void TAlignmentParser::mergeAlignedBasesBamReadsRandom(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator){
	//deletions and insertions are kept as is. these positions are not compared
	if(fwdAlignment->lastAlignedPositionWithRespectToRef >= revAlignment->position){
		fwdAlignment->setAlignmentHasChanged();
		revAlignment->setAlignmentHasChanged();

		//reads overlap -> check if there are bases overlapping same position in ref
		//alignedPos is with respect to read
		int fwdP = 0;
		int revP = 0;
		while(fwdP <= fwdAlignment->lastAlignedPos && revP <= revAlignment->lastAlignedPos){
			if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] == revAlignment->position + revAlignment->alignedPosition[revP]){
				//bases overlap same position in ref -> choose at random which to keep
				if(randomGenerator->getRand() < 0.5){
					_adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					_adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] < revAlignment->position + revAlignment->alignedPosition[revP]){
				++fwdP;
			} else {
				++revP;
			}
		}
	}
};

void TAlignmentParser::mergeAlignedBasesBamReads(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality){
	//deletions and insertions are kept as is. these positions are not compared
	if(fwdAlignment->lastAlignedPositionWithRespectToRef >= revAlignment->position){
		fwdAlignment->setAlignmentHasChanged();
		revAlignment->setAlignmentHasChanged();

		//reads overlap -> check if there are bases overlapping same position in ref
		//alignedPos is with respect to read
		int fwdP = 0;
		int revP = 0;
		while(fwdP <= fwdAlignment->lastAlignedPos && revP <= revAlignment->lastAlignedPos){
			if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] == revAlignment->position + revAlignment->alignedPosition[revP]){
				//bases overlap same position in ref -> keep the one with higher quality
				if(fwdAlignment->bases[fwdP].recalibratedQualityAsPhredInt > revAlignment->bases[revP].recalibratedQualityAsPhredInt){
					_adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					_adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] < revAlignment->position + revAlignment->alignedPosition[revP]){
				++fwdP;
			} else {
				++revP;
			}
		}
	}
};

//-----------------------------------------------------
// TBamProgressReporter
//-----------------------------------------------------
TBamProgressReporter::TBamProgressReporter(int Frequency, TAlignmentParser* Parser, TLog* Logfile){
	_init(Frequency, Parser, Logfile);
};

TBamProgressReporter::TBamProgressReporter(TAlignmentParser* Parser, TLog* Logfile){
	_init(1000000, Parser, Logfile);
};

void TBamProgressReporter::_init(int Frequency, TAlignmentParser* Parser, TLog* Logfile){
	progressFrequency = Frequency;
	parser = Parser;
	logfile = Logfile;
	timer.start();
	lastProgressPrinted = 0;

	logfile->startIndent("Parsing through BAM file:");
};

void TBamProgressReporter::_printProgress(){
	std::string percentOfFile = to_string_with_precision(parser->getPositionInFile() * 100, 2);
	std::string millionReads = to_string_with_precision((double) parser->getNumAlignmentsRead() / 1000000.0, 1);
	logfile->list("Parsed " + millionReads + " million reads (est. " + percentOfFile + "%) in " + timer.minutes() + " min.");
};

void TBamProgressReporter::printProgress(){
	if(parser->getNumAlignmentsRead() - lastProgressPrinted >= progressFrequency){
		_printProgress();
		lastProgressPrinted = parser->getNumAlignmentsRead();
	}
};

void TBamProgressReporter::printEnd(){
	logfile->list("Reached end of BAM file.");
	std::string millionReads = to_string_with_precision((double) parser->getNumAlignmentsRead() / 1000000.0, 1);
	logfile->conclude("Parsed a total of " + millionReads + " million reads in " + timer.minutes() + " min.");
	logfile->endIndent();
};

void TBamProgressReporter::printEndNoEndIndent(){
	logfile->list("Reached end of BAM file.");
	std::string millionReads = to_string_with_precision((double) parser->getNumAlignmentsRead() / 1000000.0, 1);
	logfile->conclude("Parsed a total of " + millionReads + " million reads in " + timer.minutes() + " min.");
};



