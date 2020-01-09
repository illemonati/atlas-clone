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
	bufferSize = 100000;
	reference = Reference;
	referenceSequence = "";
	curStart = -1;
	curChr = -1;
	curEnd = -1;
};

void TFastaBuffer::moveTo(const int & chr, const int32_t & pos){
	curChr = chr;
	curStart = pos;
	curEnd = pos + bufferSize;
	if(!reference->GetSequence(chr, curStart, curEnd, referenceSequence))
		throw "Problem reading " + toString(chr) + ":" + toString(curStart) + "-" + toString(curEnd) + " from fasta file!";
};

void TFastaBuffer::fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref){
	//move buffer, if necessary
	if(chr != curChr || end > curEnd || start < curStart){
		if(end - start + 1 > bufferSize){
			bufferSize = end - start + 1;
		}
		moveTo(chr, start);
	}

	//now copy to string
	ref.assign(referenceSequence, start - curStart, end - start + 1);
};

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
TAlignmentParser::TAlignmentParser(){
	logfile = NULL;
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
	previousAlignmentPos = -1;
	previousAlignmentChr = -1;
	oldAlignment = NULL;
	oldAlignmentInitialized = false;
	oldAlignmentMustBeConsidered = false;

	totalNumberAlignmentsRead = 0;
	sizeOfBamFile = 0;

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
	minQual = 33;
	maxQual = 126;
	minPhredInt = 0;
	maxPhredInt = 93;
	minQualForPrinting = 33;
	maxQualForPrinting = 126;
	trimReads = false;
	trimmingLength3Prime = 0;
	trimmingLength5Prime = 0;
	applyFragmentLengthLongerThanInsertSizeFilter = false;
	useStrand[0] = true; useStrand[1] = true;
	useMate[0] = true; useMate[1] = true;

	//blacklist
	_updateBlacklist = false;
	_writeBlackList = false;

	//limit chr and windows
	limitWindows = -1;
	skipWindows = 0;
	indexOfLimitChr = -1;

	//reference
	hasReference = false;
	fastaReference = NULL;
	fastaBuffer = NULL;
	chrChangedAlignment = false;
	chrChangedWindow = false;

	//post mortem damage and recalibration
	hasPMD = false;
	pmdObjects = NULL;
	doRecalibration = false;
	recalObjectInitialized = false;
	recalObject = NULL;
};

TAlignmentParser::TAlignmentParser(int MaxReadLength, TParameters & params, TLog* Logfile){
	TAlignmentParser();

	init(MaxReadLength, params, Logfile);
};

TAlignmentParser::~TAlignmentParser(){
	if(hasReference){
		delete fastaBuffer;
	}
	if(doMasking)
		delete mask;
	if(windowsPredefined)
		delete predefinedWindows;
	if(subset)
		delete subset;
	if(recalObjectInitialized)
		delete recalObject;
	if(pmdObjects)
		delete[] pmdObjects;
	if(oldAlignmentInitialized){
		delete oldAlignment;
	}

	if(_writeBlackList)
		ignoredReads.close();
};

void TAlignmentParser::init(int MaxReadLength, TParameters & params, TLog* Logfile){
	logfile = Logfile;

	//BAM file
	filename = params.getParameterString("bam");
	openBamFile(filename);

	//alignments
	maxReadLength = MaxReadLength;
	oldAlignment = new TAlignment(maxReadLength);
	oldAlignmentInitialized = true;

	//initialize
	initializeReadGroups(params);
	initializePostMortemDamage(params);
	initializeRecalibration(params);

	//settings
	setWindowParameters(params);
	setChrAndWindowLimits(params);
	setMasks(params);
	setFilters(params);
	setChrPloidy(params);
};

void TAlignmentParser::openBamFile(std::string filename){
	//open BAM file
	logfile->list("Reading data from BAM file '" + filename + "'.");
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";
	//load index file
	if(!bamReader.LocateIndex())
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

void TAlignmentParser::setOutName(std::string outputName){
	outname = outputName;
}

void TAlignmentParser::setWindowParameters(TParameters & params){
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
		logfile->conclude("read " + toString(predefinedWindows->size()) + " on " + toString(predefinedWindows->getNumChromosomes()) + " chromosomes");
	}
	numWindowsOnChr = 0;
};

void TAlignmentParser::setFilters(TParameters & params){
	//depth filter
	readUpToDepth = params.getParameterIntWithDefault("readUpToDepth", 1000);
	if(params.parameterExists("minDepth") || params.parameterExists("maxDepth")){
		applyDepthFilter = true;
		unsigned int tmpInt;
		tmpInt = params.getParameterIntWithDefault("minDepth", 0);
		if(tmpInt < 0)
			throw "minDepth must be >= 0!";
		minDepth = tmpInt;
		tmpInt = params.getParameterIntWithDefault("maxDepth", 1000000);
		if(tmpInt < minDepth) throw "maxDepth must be >= minDepth!";
		maxDepth = tmpInt;
		readUpToDepth = maxDepth + 1;
		logfile->list("Will filter out sites with sequencing depth < " + toString(minDepth) + " or > " + toString(maxDepth) + ". (parameters 'minDepth', 'maxDepth')");
	} else {
		applyDepthFilter = false;
		minDepth = 0;
		maxDepth = 1000000;
	}
	logfile->list("Will read data up to depth " + toString(readUpToDepth) + " and ignore additional bases. (parameter 'readUpToDepth')");

	//Mapping quality filters
	if(params.parameterExists("minMQ") || params.parameterExists("maxMQ")){
		applyMQFilter = true;
		minMQ = params.getParameterIntWithDefault("minMQ", 0);
		if(minMQ < 0)
			throw "minMQ must be >= 0!";
		maxMQ = params.getParameterIntWithDefault("maxMQ", 1000000);
		if(maxMQ < minMQ)
			throw "maxMQ must be larger than minMQ";
		setMappingQualityFilters(minMQ, maxMQ);
		logfile->list("Will filter out reads with mapping quality outside the range [" + toString(minMQ) + ", " + toString(maxMQ) + "]. (parameters 'minMQ', 'maxMQ')");
	}

	//quality filters
	minPhredInt = params.getParameterIntWithDefault("minQual", 1);
	if(minPhredInt < 0) throw "minQual must be >= 0!";
	maxPhredInt = params.getParameterIntWithDefault("maxQual", 93);
	if(maxPhredInt < minPhredInt) throw "maxQual must be >= minQual!";
	setQualityFilters(minPhredInt, maxPhredInt);
	logfile->list("Will filter out bases with quality outside the range [" + toString(minPhredInt) + ", " + toString(maxPhredInt) + "] (parameters 'minQual', 'maxQual')");

	//quality filters for printing
	int minOutQual = params.getParameterIntWithDefault("minOutQual", 0) + 33;
	if(minOutQual < 0) throw "minOutQual must be >= 0!";
	int maxOutQual = params.getParameterIntWithDefault("maxOutQual", 93) + 33;
	if(maxOutQual < minOutQual) throw "maxOutQual must be >= minOutQual!";
	setQualityRangeForPrinting(minOutQual, maxOutQual);
	logfile->list("Will print qualities truncated to [" + toString(minOutQual) + ", " + toString(maxOutQual) + "] (parameters 'minOutQual', 'maxOutQual')");

	//context filter
	if(params.parameterExists("ignoreContexts")){
		std::vector<std::string> contexts;
		fillVectorFromString(params.getParameterString("ignoreContexts"), contexts, ',');
		setContextFilter(contexts);
	}

	//filter for missing reference
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 1.0);
	if(maxMissing > 1.0) throw "maxMissing must be smaller or equal to 1.0!";
	logfile->list("Will filter out windows with a missing data fraction > " + toString(maxMissing) + ". (parameter 'maxMissing')");

	maxRefN = params.getParameterDoubleWithDefault("maxRefN", 1.0);
	if(maxRefN > 1.0) throw "maxRefN must be smaller or equal to 1.0!";
	if(maxRefN < 1.0 && hasReference == false) throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided.";
	logfile->list("Will filter out windows with a fraction of 'N' in reference > " + toString(maxMissing) + ". (parameter 'maxRefN')");

	//duplicates
	if(params.parameterExists("keepDuplicates")){
		keepDuplicates();
		logfile->list("Will keep duplicate reads. (parameter 'keepDuplicates')");
	}

	//soft clips
	if(params.parameterExists("filterSoftClips")){
		filterSoftClips();
		logfile->list("Will filter out soft clipped reads. (parameter 'filterSoftClips')");
	}

	//improper pairs
	if(params.parameterExists("keepImproperPairs")){
		keepImproperPairs();
		logfile->list("Will keep improper pairs. (parameter 'keepImproperPairs')");
	}

	//unmapped reads
	if(params.parameterExists("keepUnmappedReads")){
		keepUnmappedReads();
		logfile->list("Will keep unmapped reads. (parameter 'keepUnmappedReads')");
	}

	//failed QC
	if(params.parameterExists("keepFailedQC")){
		keepFailedQC();
		logfile->list("Will keep reads that failed QC. (parameter 'keepFailedQC')");
	}

	//secondary reads
	if(params.parameterExists("keepSecondaryReads")){
		keepSecondaryReads();
		logfile->list("Will keep secondary reads. (parameter 'keepSecondaryReads')");
	}

	//supplementary reads
	if(params.parameterExists("keepSupplementaryReads")){
		keepSupplementaryReads();
		logfile->list("Will keep supplementary reads. (parameter 'keepSupplementaryReads')");
	}

	//fragment length
	if(params.parameterExists("keepReadsLongerThanFragment")){
		setApplyFragmentLengthLongerThanInsertSizeFilter(false);
		logfile->list("Will keep reads that are longer than the fragment size. (parameter 'keepReadsLongerThanFragment')");
	} else
		setApplyFragmentLengthLongerThanInsertSizeFilter(true);

	//strand
	if(params.parameterExists("keepOnlyFwd")){
		useStrand[1] = false;
		logfile->list("Will keep only forward mapping reads. (parameter 'keepOnlyFwd')");
	}
	else if(params.parameterExists("keepOnlyRev")){
		useStrand[0] = false;
		logfile->list("Will keep only reverse mapping reads. (parameter 'keepOnlyRev')");
	}

	//mate
	if(params.parameterExists("keepOnlyFirst")){
		useMate[1] = false;
		logfile->list("Will keep only the first mates. (parameter 'keepOnlyFirst')");
	}
	else if(params.parameterExists("keepOnlySecond")){
		useMate[0] = false;
		logfile->list("Will keep only the second mates. (parameter 'keepOnlySecond')");
	}
}

void TAlignmentParser::setQualityFilters(int MinPhredInt, int MaxPhredInt){
	applyQualityFilter = true;
	minPhredInt = MinPhredInt;
	maxPhredInt = MaxPhredInt;
	minQual = qualMap.phredIntToQuality(minPhredInt);
	maxQual = qualMap.phredIntToQuality(maxPhredInt);
};

void TAlignmentParser::setMappingQualityFilters(int MinMQ, int MaxMQ){
	applyMQFilter = true;
	minMQ = MinMQ;
	maxMQ = MaxMQ;
};

void TAlignmentParser::setQualityRangeForPrinting(int minQual, int maxQual){
	minQualForPrinting = minQual;
	maxQualForPrinting = maxQual;
};

void TAlignmentParser::setContextFilter(std::vector<std::string> contexts){
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

void TAlignmentParser::setMasks(TParameters & params){
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
		mask = new TBedReader(maskFile, windowSize, bamHeader.Sequences, siteLimit, logfile);
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
		mask = new TBedReader(regionsFile, windowSize, bamHeader.Sequences, siteLimit, logfile);
		logfile->done();
		logfile->endIndent();
	} else considerRegions = false;
}

void TAlignmentParser::setReadTrimming(int trim3Prime, int trim5Prime){
	if(trim3Prime < 0) throw "Trimming at 3 prime end must be >= 0!";
	if(trim5Prime < 0) throw "Trimming at 3 prime end must be >= 0!";
	trimmingLength3Prime = trim3Prime;
	trimmingLength5Prime = trim5Prime;
	trimReads = true;
};

void TAlignmentParser::setApplyFragmentLengthLongerThanInsertSizeFilter(bool filterYesNo){
	applyFragmentLengthLongerThanInsertSizeFilter = filterYesNo;
};

void TAlignmentParser::setFragmentLengthFilter(int MinFragmentLength, int MaxFragmentLength){
	if(MinFragmentLength < 0)
		throw "Min fragment length must be >= 0!";
	if(maxFragmentLength < 0)
		throw "Max fragment length must be >= 0!";
	if(MinFragmentLength > MaxFragmentLength)
		throw "Min fragment length must be <= max fragment length!";

	minFragmentLength = MinFragmentLength;
	maxFragmentLength = MaxFragmentLength;
	applyFragmentLengthFilter = true;
	setParsingToTrue(); //filter requires parsing!
};

void TAlignmentParser::initializeReadGroups(TParameters & params){
	readGroups.fill(bamHeader);

	//limit readGroups
	if(params.parameterExists("readGroup")){
		readGroups.filterReadGroups(params.getParameterString("readGroup"));
		logfile->startIndent("Will limit analysis to the following read groups:");
		readGroups.printReadgroupsInUse(logfile);
		logfile->endIndent();
	}
};

void TAlignmentParser::setChrAndWindowLimits(TParameters & params){
	if(params.parameterExists("chr")){
		//parse chromosome names
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("chr"), vec, ',');
		chromosomes.useSpecifiedChr(vec, logfile);
	} else {
		if(params.parameterExists("limitChr")){
			std::string limitName = params.getParameterString("limitChr");
			logfile->list("Will limit analysis to all chromosomes up to and including " + limitName + ". (parameter 'limitChr')");
			chromosomes.limitChr(limitName);
			indexOfLimitChr = chromosomes.getIndexFromName(limitName) + 1;
		}
	}

	skipWindows = params.getParameterIntWithDefault("skipWindows", 0);
	if(skipWindows > 0) logfile->list("Will skip the first " + toString(skipWindows) + " windows per chromosome. (parameter 'skipWindows')");
	limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows"))
		logfile->list("Will limit analysis to the first " + toString(limitWindows) + " windows per chromosome. (parameter 'limitWindows')");
	if(limitWindows <= skipWindows)
		throw "limitWindows has to be larger than skipWindows!";
};

void TAlignmentParser::setChrPloidy(TParameters & params){
	logfile->list("Chromosomes with no further specifications are assumed to be diploid (parameters 'ploidy' or 'haploid' to change ploidy).");
	if(params.parameterExists("ploidy")){
		std::string ploidyFileName = params.getParameterString("ploidy");
		logfile->list("Reading ploidy specification per chromosome from file '" + ploidyFileName + "'");
		std::ifstream ploidyFile(ploidyFileName.c_str());
		if(!ploidyFile)
			throw "Failed to open file '" + ploidyFileName + "'!";
		chromosomes.specifyPloidy(ploidyFile, logfile);
	}

	if(params.parameterExists("haploid")){
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("haploid"), vec, ',');
		chromosomes.setToHaploid(vec, logfile);
	}
};

void TAlignmentParser::addReference(BamTools::Fasta* reference){
	hasReference = true;
	fastaReference = reference;
	fastaBuffer = new TFastaBuffer(reference);
};

void TAlignmentParser::fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment){
	if(!hasReference) //is this check really necessary?
		throw "No reference provided!";
	std::string referenceSequence;
	fastaBuffer->fill(alignment.chrNumber, alignment.position, alignment.position + alignment.bases[alignment.length-1].alignedPos, referenceSequence);
	alignment.referenceSequence = referenceSequence;
	alignment.hasReference = true;
};

std::string TAlignmentParser::chrNumberToName(int chrNumber){
	int counter = 0;
	for(BamTools::SamSequenceIterator chrIt=bamHeader.Sequences.Begin(); chrIt!=bamHeader.Sequences.End(); ++chrIt, ++counter){
		if(counter == chrNumber)
			return chrIt->Name;
	}
	throw "chrNumber not in header";
};

int TAlignmentParser::chrNumberToLength(int chrNumber){
	int counter = 0;
	for(BamTools::SamSequenceIterator chrIt=bamHeader.Sequences.Begin(); chrIt!=bamHeader.Sequences.End(); ++chrIt, ++counter){
		if(counter == chrNumber)
			return stringToInt(chrIt->Length);
	}
	throw "chrNumber not in header";
};

long TAlignmentParser::calcReferenceLength(){
    return chromosomes.referenceLength();
};

std::string TAlignmentParser::getCurChrName(){
	return chromosomes.curName();
};

long TAlignmentParser::getCurChrLength(){
	return chromosomes.curLength();
};

int TAlignmentParser::getCurChrPloidy(){
	return chromosomes.curPloidy();
};

//--------------
//move genome
//--------------
void TAlignmentParser::jumpToEnd(){
	chromosomes.jumpToEnd();
};

void TAlignmentParser::restartChromosomes(TWindow_base & window){
	chromosomes.begin();

	moveChromosome(window);
};

void TAlignmentParser::moveChromosome(TWindow_base & window){
	//jump reader
	oldAlignmentMustBeConsidered = false;

	//restart windows
	previousAlignmentPos = -1;
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
		window.move(predefinedWindows->curWindowStart(), predefinedWindows->curWindowEnd(), chromosomes.curIndex());
		window.chrName = chromosomes.curName();
		bamReader.Jump(chromosomes.curIndex(), window.start);

	} else {
		while(chromosomes.curInUse() == false || skipWindows * windowSize > chromosomes.curLength()){
//			++chrIterator;
//			++chrNumber;
//			chrLength = stringToLong(chrIterator->Length);
			chromosomes.next();
		}
		numWindowsOnChr = ceil(chromosomes.curLength() / (double) windowSize);

		int curStart = skipWindows * windowSize;
		bamReader.Jump(chromosomes.curIndex(), curStart);
		int nextEnd = curStart + windowSize;

		if(nextEnd > chromosomes.curLength()){
			nextEnd = chromosomes.curLength();
		}
		window.move(curStart, nextEnd, chromosomes.curIndex());
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

bool TAlignmentParser::moveToNextWindowOnChr(TWindow_base & window){

	if(window.end > 0) logfile->endIndent();

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
	window.move(window.end, nextEnd, chromosomes.curIndex());

	return true;
};

bool TAlignmentParser::moveToNextPredefinedWindow(TWindow_base & window){

	if(window.end > 0) logfile->endIndent();

	++windowNumber;
	if(windowNumber >= limitWindows)
		return false;
	if(predefinedWindows->nextWindow()){
		window.move(predefinedWindows->curWindowStart(), predefinedWindows->curWindowEnd(), chromosomes.curIndex());
		//should we jump or are we already close enough to next window
		if(abs(window.start - previousAlignmentPos) > maxReadLength){
			previousAlignmentPos = -1;
			if(window.start - maxReadLength < 0)
				bamReader.Jump(chromosomes.curIndex(), 0);
			else{
				bamReader.Jump(chromosomes.curIndex(), window.start - maxReadLength);
			}
		}
		return true;
	} else
		return false;
};

bool TAlignmentParser::moveWindow(TWindow_base & window){
	//returns false when end of genome is reached
	if(windowsPredefined){
		//if at beginning of BAM file
		if(chromosomes.end()){
			restartChromosomes(window);

			if(chromosomes.end())
				throw "found no predefined windows in BED file! Does file exist?";
			chrChangedWindow = true;

		} else {
			//now move coordinates of next window
			if(!moveToNextPredefinedWindow(window)){
				//no more windows left on chr
				chromosomes.next();

				if(chromosomes.end())
					return false;

				moveChromosome(window);
				chrChangedWindow = true;

				if(chromosomes.end())
					return false;
				++windowNumber;
			} else
				//was able to move to next window on chr
				chrChangedWindow = false;
		}

	} else {
		//if at beginning of BAM file
		if(chromosomes.end()){
			restartChromosomes(window);
			chrChangedWindow = true;
		} else {
			if(!moveToNextWindowOnChr(window)){
				//there is no window left on chr
				chromosomes.next();

				//do we use this chromosome? if not, move on!
				while(!chromosomes.end() && !chromosomes.curInUse()){
					chromosomes.next();
				}

				//did we reach end?
				if(chromosomes.end() || (indexOfLimitChr != -1 && chromosomes.curIndex() >= indexOfLimitChr)){
					window.end = 0;
					return false;
				}
				moveChromosome(window);
				chrChangedWindow = true;
			} else {
				chrChangedWindow = false;
			}
		}
	}

	//report
	logfile->number("Window [" + toString(window.start) + ", " + toString(window.end) + ") of " + toString(numWindowsOnChr) + " on '" + chromosomes.curName() + "':");
	logfile->addIndent();
	return true;
};

//------------------------------
//reading alignments
//------------------------------
bool TAlignmentParser::readAlignment(){
	bool filtersPassed = false;
	do {
		if(!bamReader.GetNextAlignment(bamAlignment)){
			return false;
		}
		++totalNumberAlignmentsRead;

		//check if chromosome changed
		if(bamAlignment.RefID != previousAlignmentChr){
			previousAlignmentPos = -1;
			previousAlignmentChr = bamAlignment.RefID;
//			chrNumber = previousAlignmentChr;
			chrChangedAlignment = true;
		} else
			chrChangedAlignment = false;

		if(bamAlignment.Position < previousAlignmentPos)
			throw "BAM file must be sorted by position! Alignment '" + bamAlignment.Name + "' is at position " + toString(bamAlignment.Position) + ", which is before the position of the previous alignment (" + toString(previousAlignmentPos) + ")";
		previousAlignmentPos = bamAlignment.Position;

		//check read length
		/* check if insert size is shorter than read-insertions+deletions=alignedBases.length() -> this means we are reading the adaptor sequence.
		Insert size is determined by mapping -> insertions are not in ref and should not count. If we don't add deletions, adapter at end could be sequenced but we still keep read
		(deletions in aligned bases are represented as dashes) */
		if(bamAlignment.AlignedBases.size() > maxReadLength)
			throw "Alignment '" +  bamAlignment.Name + "' is longer than the max read length " + toString(maxReadLength) + "! Please change max read length to parse this data.";

		//store read group ID
		std::string readGroup;
		bamAlignment.GetTag("RG", readGroup);
		curReadGroupID = readGroups.find(readGroup);

		//filter
		//TODO: add functionality to not filter at all (i.e. _keepAll switch)
		filtersPassed = true;
		/* check if insert size is shorter than read-insertions+deletions=alignedBases.length() -> this means we are reading the adaptor sequence.
		Insert size is determined by mapping -> insertions are not in ref and should not count. If we don't add deletions, adapter at end could be sequenced but we still keep read
		(deletions in aligned bases are represented as dashes) */
		if(bamAlignment.IsPaired() && applyFragmentLengthLongerThanInsertSizeFilter && abs(bamAlignment.InsertSize) < bamAlignment.AlignedBases.length()){
			if(_updateBlacklist){
				addToBlacklist(bamAlignment, "longer than insert size (TLEN)");
			} else {
				logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
			}
			filtersPassed = false;
		} else {
			//apply filters: read group in use and basic QC
			filtersPassed = applyFilters();
			if(_updateBlacklist && !filtersPassed){
				addToBlacklist(bamAlignment, "did not pass parser filters");
			}
		}
	} while(!filtersPassed);

	return true;
};

bool TAlignmentParser::applyFilters(){
	std::vector<int> clipSizes;
	std::vector<int> readPositions;
	std::vector<int> genomePositions;

	bool filtersPassed = readGroups.readGroupInUse(curReadGroupID)
					&& (_keepImproperPairs || (!bamAlignment.IsPaired() || bamAlignment.IsProperPair()))
					&& (_keepUnmappedReads || bamAlignment.IsMapped())
					&& (_keepFailedQC || !bamAlignment.IsFailedQC())
					&& (_keepSecondary || bamAlignment.IsPrimaryAlignment())
					&& (_keepSupplementary || !bamAlignment.IsSupplementary())
					&& chromosomes.inUse(bamAlignment.RefID)
					&& (_keepDuplicates || !bamAlignment.IsDuplicate())
					&& (!_filterSoftClips || !bamAlignment.GetSoftClips(clipSizes, readPositions, genomePositions))
					&& useStrand[bamAlignment.IsReverseStrand()]
					&& useMate[bamAlignment.IsSecondMate()]
					&& (!applyMQFilter || (bamAlignment.MapQuality >= minMQ && bamAlignment.MapQuality <= maxMQ));

	return filtersPassed;
};

bool TAlignmentParser::fillAlignment(TAlignment & alignment){
	//make sure container is empty
	alignment.clear();

	//fill alignment
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
	int readGroupId = readGroups.find(readGroup);

	alignment.fill(bamAlignment, readGroupId);

	if(_parse){
		//add all info from bamAlignment to bases
		alignment.parse(genoMap, qualMap);

		//check for fragment length filter
		if(applyFragmentLengthFilter && (alignment.fragmentLength < minFragmentLength || alignment.fragmentLength > maxFragmentLength))
			return false;

		//add missing information to bases
		alignment.fillReadGroupInfo(readGroupId);
		alignment.fillPmdProbabilities(pmdObjects);

		if(trimReads)
			alignment.trimRead(trimmingLength3Prime, trimmingLength5Prime);
		if(applyQualityFilter)
			alignment.filterForBaseQuality(minQual, maxQual);
		if(applyContextFilter)
			alignment.filterForContext(ignoreTheseContexts);
		if(doRecalibration)
			recalibrate(alignment);
		if(hasReference)
			fillReferenceSequence(fastaBuffer, alignment);
	}

	return true;
};

//------------------------
//read data in alignments
//------------------------

bool TAlignmentParser::readNextAlignment(TAlignment & alignment){
	//use this in TGenome for functionalities that don't need windows
	if(readAlignment()){
		return fillAlignment(alignment);
	}
	return false;
};

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
		if(bamAlignment.Position >= window.end){
			return;
		}

		oldAlignmentMustBeConsidered = false;
		if(oldAlignment->lastAlignedPositionWithRespectToRef >= window.start)
			oldAlignment = window.swapUsedForEmptyAlignment(oldAlignment, maxReadLength);
	}

	//read alignments
	int counter = 0;
	while(readAlignment()){
		//fill alignment
		//return false if a post-parsing filer is not passed
		if(!fillAlignment(*oldAlignment)){
			continue;
		}

		++counter;

		//check if alignment starts after current window end -> break
		if(oldAlignment->position >= window.end || oldAlignment->chrNumber != window.chrNumber){
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
	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");

	//apply filters
	applyWindowFilters(window);
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
	applyWindowFilters(destination);
};

void TAlignmentParser::applyWindowFilters(TWindow_base & window){
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

//------------------------------
//initialize PMD and recalibration
//------------------------------
PMDType TAlignmentParser::getEnumPMDType(std::string pmdType){
	if(pmdType == "CT")
		return pmdCT;
	else if(pmdType == "GA")
		return pmdGA;
	else if(pmdType == "GT")
		return pmdGT;
	else if(pmdType == "CA")
		return pmdCA;
	else {
		throw "unknown pmdType: " + pmdType + "!";
	}
};

void TAlignmentParser::initializePostMortemDamage(TParameters & params){
	logfile->startIndent("Initializing Post Mortem Damage (PMD):");
	//create an array of TPMD objects for each read group
	pmdObjects = new TPMD[readGroups.size()];

	//now fill them!
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		//all read groups have the same pmd
		logfile->list("Initializing one PMD function for all read groups.");
		pmdObjects[0].initialize(params, logfile);
		for(size_t i=1; i<readGroups.size(); ++i)
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
				if(readGroups.readGroupExists(vec[0])){ //ignore if it does not exist
					//get read group and PMD type
					readGroupId = readGroups.find(vec[0]);
					if(!params.parameterExists("oldPMDFormat")){
						PMDType pmdType = getEnumPMDType(vec[1]);
						//initialize functions
						pmdObjects[readGroupId].initializeFunction(vec[2], pmdType);
					} else {
						//initialize functions
						pmdObjects[readGroupId].initializeFunction(vec[1], pmdCT);
					//	logfile->conclude("For read group '" + vec[0] + "', C->T: " + pmdObjects[readGroupId].getFunctionString(pmdCT));
						pmdObjects[readGroupId].initializeFunction(vec[2], pmdGA);
					//	logfile->conclude("For read group '" + vec[0] + "', G->A: " + pmdObjects[readGroupId].getFunctionString(pmdGA));
					}
				}
			}
		}

		//close file
		file.close();

		//test if we have a function for all read groups
		for(size_t i=0; i<readGroups.size(); ++i){
			if(!pmdObjects[i].functionInitialized(pmdCT)) throw "PMD C->T for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
			if(!pmdObjects[i].functionInitialized(pmdGA)) throw "PMD G->A for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
		}
		hasPMD = true;
	} else {
		//no post mortem damage
		logfile->list("Assuming there is no PMD in the data.");
		std::string pmdString = "none";
		for(size_t i=0; i<readGroups.size(); ++i){
			pmdObjects[i].initializeFunction(pmdString, pmdGA);
			pmdObjects[i].initializeFunction(pmdString, pmdCT);
		}
	}
	logfile->endIndent();
};

void TAlignmentParser::initializeRecalibration(TParameters & params){
	std::string task = params.getParameterString("task");
	if(params.parameterExists("BQSRQuality")){
		//initialize quality effect
		TRecalibrationBQSR* bqsr = new TRecalibrationBQSR(params.getParameterString("BQSRQuality"), logfile, &readGroups);

		//Do we consider the effect of the position in read (cycle)?
		if(params.parameterExists("BQSRPosition")){
			bqsr->addPositionEffect(params.getParameterString("BQSRPosition"));
		}

		if(params.parameterExists("BQSRPositionReverse")){
			bqsr->addPositionReverseEffect(params.getParameterString("BQSRPositionReverse"));
		}

		//Do we consider the context (dinucleotide)?
		if(params.parameterExists("BQSRContext")){
			bqsr->addContextEffect(params.getParameterString("BQSRContext"));
		}

		doRecalibration = true;

	} else if(params.parameterExists("recal")){
		recalObject = new TRecalibrationEM(params.getParameterString("recal"), &readGroups, logfile);
		doRecalibration = true;
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
		doRecalibration = false;
		recalObject = new TRecalibration();
	}
	recalObjectInitialized = true;

	//check if estimation is required, in which case throw an error!
	if(recalObject->_requiresEstimation()) throw "Can not use provided recalibration: estimation is required!";
};

void TAlignmentParser::recalibrate(TAlignment & alignment){
	//make sure read is parsed and has reference
	if(!alignment.parsed) throw "Read was not parsed!";

	if(recalObject->recalibrationChangesQualities()){
		//recalibrate quality scores
		for(int d=0; d<alignment.length; ++d){
			if(alignment.bases[d].aligned && alignment.bases[d].base != N){
				alignment.bases[d].errorRate = recalObject->getErrorRate(alignment.bases[d]);
			}
		}

		alignment.changed = true;
	} else alignment.changed = false;
	alignment.recalibrated = true;
};

void TAlignmentParser::addSitesToQualityTransformTable(TAlignment & alignment, TQualityTransformTables & QTtables){
	for(int i=0; i<alignment.length; ++i){
		if(alignment.bases[i].base != N){
			int newQual = qualMap.errorToQuality(alignment.bases[i].errorRate);
			QTtables.add(alignment.readGroupId, alignment.qualityOriginal[i], newQual);
		}
	}
};

void TAlignmentParser::addSitesToQualityTransformTable(TAlignment & alignment, TRecalibration* otherRecalObject, TQualityTransformTables & QTtables){
	for(int i=0; i<alignment.length; ++i){
		if(alignment.bases[i].base != N){
			int firstQual = qualMap.errorToQuality(alignment.bases[i].errorRate);
			double tmp = alignment.bases[i].errorRate;
			alignment.bases[i].errorRate = qualMap.qualityToError(alignment.qualityOriginal[i]);
			int secondQual = otherRecalObject->getQuality(alignment.bases[i]);
			alignment.bases[i].errorRate = tmp;
			QTtables.add(alignment.readGroupId, firstQual, secondQual);
		}
	}
};

void TAlignmentParser::adaptQualityWhenMerging(TBase & bestBase, TBase & worstBase, const bool & adaptQuality){
	if(adaptQuality){
		double likelihood[4];
		double sum = 0.0;
		for(int i=0; i<4; ++i){
			if(bestBase.base == i){
				likelihood[i] = 1.0 - bestBase.errorRate;
			} else {
				likelihood[i] = bestBase.errorRate / 3.0;
			}

			if(worstBase.base == i){
				likelihood[i] *= 1.0 - worstBase.errorRate;
			} else {
				likelihood[i] *= worstBase.errorRate / 3.0;
			}
			sum += likelihood[i];

		}
		bestBase.errorRate = 1.0 - likelihood[bestBase.base] / sum;
	} else {
		worstBase.errorRate = 1.0;
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
			if(fwdAlignment->position + fwdAlignment->bases[fwdP].alignedPos == revAlignment->position + revAlignment->bases[revP].alignedPos){
				//bases overlap same position in ref -> choose at random which to keep
				if(keepFwd){
					adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->bases[fwdP].alignedPos < revAlignment->position + revAlignment->bases[revP].alignedPos){
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
			if(fwdAlignment->position + fwdAlignment->bases[fwdP].alignedPos == revAlignment->position + revAlignment->bases[revP].alignedPos){
				//bases overlap same position in ref -> choose at random which to keep
				if(randomGenerator->getRand() < 0.5){
					adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->bases[fwdP].alignedPos < revAlignment->position + revAlignment->bases[revP].alignedPos){
				++fwdP;
			} else {
				++revP;
			}
		}
	}
}


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
			if(fwdAlignment->position + fwdAlignment->bases[fwdP].alignedPos == revAlignment->position + revAlignment->bases[revP].alignedPos){
				//bases overlap same position in ref -> keep the one with higher quality
				if(fwdAlignment->bases[fwdP].errorRate < revAlignment->bases[revP].errorRate){
					adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->bases[fwdP].alignedPos < revAlignment->position + revAlignment->bases[revP].alignedPos){
				++fwdP;
			} else {
				++revP;
			}
		}
	}
}

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
	gettimeofday(&start, NULL);
	lastProgressPrinted = 0;

	logfile->startIndent("Parsing through BAM file:");
};

std::string TBamProgressReporter::_getRunTime(){
	gettimeofday(&end, NULL);
	return to_string_with_precision((end.tv_sec  - start.tv_sec)/60.0, 2);
};

void TBamProgressReporter::_printProgress(){
	std::string percentOfFile = to_string_with_precision(parser->getPositionInFile() * 100, 2);
	std::string millionReads = to_string_with_precision((double) parser->getNumAlignmentsRead() / 1000000.0, 1);
	logfile->list("Parsed " + millionReads + " million reads (est. " + percentOfFile + "%) in " + _getRunTime() + " min.");
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
	logfile->conclude("Parsed a total of " + millionReads + " million reads in " + _getRunTime() + " min.");
	logfile->endIndent();
};


