/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */

#include "TGenome.h"

namespace GenomeTasks{

//---------------------------------------------------------------
// TGenome_basic
// A base class without filters and genotype likelihoods
//---------------------------------------------------------------
TGenome_basic::TGenome_basic(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator){
	_logfile = Logfile;
	_params = &Params;
	_randomGenerator = RandomGenerator;

	//open bam file
	//TODO: deal with index in better way: let tasks decide if they need an inex or not
	_bamFile.open(Params.getParameterString("bam"), Params.parameterExists("indexNotRequired"), _logfile);
	_bamFile.setLimits(Params, _logfile);

	//outputname
	_outputName = Params.getParameterStringWithDefault("out", "");
	if(_outputName == ""){
		//guess from BAM filename.
		_outputName = _bamFile.filename();
		_outputName = extractBeforeLast(_outputName, ".");
	}
	_logfile->list("Writing output files with prefix '" + _outputName + "'. (parameter 'out')");
};

void TGenome_basic::_openBamForWriting(const std::string filename, BAM::TOutputBamFile & outBam){
	_logfile->list("Writing alignments to new BAM to file '" + filename + "'.");
	outBam.open(filename, _bamFile);
	if(_params->parameterExists("writeBinnedQualities")){
		_logfile->list("When writing alignments, quality scores will be Illumina-binned. (parameter 'writeBinnedQualities').");
		outBam.binQualityScoresLikeIllumina();
	} else if(_params->parameterExists("minOutQual") || _params->parameterExists("minOutQual")){
				int MinPhredInt = _params->getParameterIntWithDefault("minOutQual", 0);
				int MaxPhredInt = _params->getParameterIntWithDefault("maxOutQual", 93);

				if(MinPhredInt < 0 || MinPhredInt > 255) throw "minOutQual " + toString(MinPhredInt) + " is outside accepted range [0, 255]!";
				if(MaxPhredInt < 0 || MaxPhredInt > 255) throw "maxOutQual " + toString(MaxPhredInt) + " is outside accepted range [0, 255]!";
				if(MaxPhredInt < MinPhredInt) throw "maxOutQual must be >= minOutQual!";

				_logfile->list("Will print qualities truncated to [" + toString(MinPhredInt) + ", " + toString(MaxPhredInt) + "] (parameters 'minOutQual', 'maxOutQual')");

				//set in quality map
				_qualMap.setQualityLimits(MinPhredInt, MaxPhredInt);
	} else {
		_logfile->list("Will use the full range of quality scores when writing alignments. (use 'writeBinnedQualities' to bin, 'minOutQual' and 'maxOutQual' to constrain).");
	}
};

//---------------------------------------------------------------
// TGenome_filtered
// A base class without genotype likelihoods but BAM filters enabled
//---------------------------------------------------------------
TGenome_filtered::TGenome_filtered(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_basic(Params, Logfile, RandomGenerator){
	_bamFile.setFilters(Params, _logfile);
};

void TGenome_filtered::_traverseBAMPassedQC(){
	//parse through bam file
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignmentThatPassesFilters()){
		//handle alignment by derived classes
		_handleAlignment();

		//report
		_bamFile.printProgress();
	}

	//report
	_bamFile.printEndWithSummary();
};


//---------------------------------------------------------------
// TGenome_parsed
// A base class with BAM filters and recalibration
//---------------------------------------------------------------
TGenome_parsed::TGenome_parsed(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_filtered(Params, Logfile, RandomGenerator){
	//initialize genotype likelihoods
	_genotypeLikelihoodCalculator.init(Params, &_bamFile.readGroups, _logfile);

	//set parsing filters
	_setReadTrimming(Params);
	_setQualityFilter(Params);
	_setContextFilter(Params);
	_setQualityRangeForPrinting(Params);
};

void TGenome_parsed::_openReference(bool required = false){
	if(!_reference.hasReference()){
		if(_params->parameterExists("fasta")){
			std::string fastaFile = _params->getParameterString("fasta");
			_logfile->list("Reading reference sequence from '" + fastaFile + "'");
			_reference.initialize(fastaFile);
		} else {
			if(required){
				throw "No reference provided! Use parameter 'fasta' to provide a reference fasta file.";
			}
		}
	}
};

void TGenome_parsed::_setReadTrimming(TParameters & params){
	//trimming ends
	if(params.parameterExists("trim3") || params.parameterExists("trim5")){
		_trimmingLength3Prime = params.getParameterIntWithDefault("trim3", 0);
		if(_trimmingLength3Prime < 0) throw "trimming distance trim3 must be >= 0!";
		_trimmingLength5Prime = params.getParameterIntWithDefault("trim5", 0);
		if(_trimmingLength5Prime < 0) throw "trimming distance trim5 must be >= 0!";
		if(_trimmingLength3Prime > 0 || _trimmingLength5Prime > 0){
			_logfile->list("Will trim first " + toString(_trimmingLength3Prime) + " and " + toString(_trimmingLength5Prime) + " bases from the 3' and 5' end, respectively. (parameters 'trim3', 'trim5')");
		}
	}
	_trimReads = true;
};

void TGenome_parsed::_setQualityFilter(TParameters & params){
	if(_qualityFilter.set(params, _logfile)){
		_applyQualityFilter = true;
	} else {
		_applyQualityFilter = false;
	}
};

void TGenome_parsed::_setContextFilter(TParameters & params){
	if(params.parameterExists("ignoreContexts")){
		std::vector<std::string> contexts;
		fillVectorFromString(params.getParameterString("ignoreContexts"), contexts, ',');
		_logfile->startIndent("Will mask the following contexts (parameter 'maskContext'):");
		for(auto& c : contexts){
			if(c.size() != 2 || !_genoMap.isValidBase(c[0]) || !_genoMap.isValidBase(c[1])){
				throw "Context " + c + " does not consist of two bases! (parameter 'maskContext')";
			}
			//ave context
			BaseContext co = _genoMap.getContext(c[0], c[1]);
			_ignoreTheseContexts.emplace(co, 1);
			_logfile->list(_genoMap.getContextString(co));
		}
		_logfile->endIndent();
		_applyContextFilter = true;
	} else {
		_applyContextFilter = false;
	}
};

void TGenome_parsed::_parseAlignment(BAM::TAlignment & alignment){
	//parse
	alignment.parse(_genoMap, _qualMap, _genotypeLikelihoodCalculator.getSequencingErrorModels());

	//apply filters
	if(_trimReads){
		alignment.trimRead(_trimmingLength3Prime, _trimmingLength5Prime);
	}
	if(_applyQualityFilter){
		alignment.filterForBaseQuality(_qualityFilter);
	}
	if(_applyContextFilter){
		alignment.filterForContext(_ignoreTheseContexts, _genoMap);
	}
	if(_reference.hasReference()){
		alignment.addReference(_reference);
	}
};

void TGenome_parsed::_traverseBAMPassedQC(){
	//parse through bam file
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignmentThatPassesFilters(_alignment)){
		//parse
		_parseAlignment(_alignment);

		//handle alignment by derived classes
		_handleAlignment();

		//report
		_bamFile.printProgress();
	}

	//report
	_bamFile.printEndWithSummary();
};

//---------------------------------------------------------------
// TGenome_windows
// A base class to traverse a BAM file in windows
//---------------------------------------------------------------
TGenome_windows::TGenome_windows(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):
		TGenome_parsed(Params, Logfile, RandomGenerator),
		_chromosomes(_bamFile.chromosomes){
	_setWindowParameters(Params);
	_setParsingLimits(Params);
	_setWindowFilters(Params);
	_setMasks(Params);
	_setSiteFilters(Params);
};

void TGenome_windows::_setWindowParameters(TParameters & params){
	if(!params.parameterExists("window") && params.parameterExists("windows")){
		_logfile->warning("Argument 'windows' specified, but unknown. Did you mean 'window'?");
	}
	std::string tmp = params.getParameterStringWithDefault("window", "1000000");

	//check if it is a number
	if(stringIsProbablyANumber(tmp)){
		_windowsPredefined = false;
		_windowSize = stringToInt(tmp);
		_logfile->list("Setting window size to " + toString(_windowSize) + ". (parameter 'window')");
		if(_windowSize < _bamFile.maxReadLength())
			throw "Window size " + tmp + " out of range! Windows must be at least as large as the max read length (" + toString(_bamFile.maxReadLength()) + " bp). (use parameter 'maxReadLength' to change)!";
	} else {
		_windowsPredefined = true;
		_logfile->listFlush("Limiting analysis to windows defined in '" + tmp + "'...");
		_predefinedWindows = new BAM::TBed(tmp);
		_logfile->done();
		_logfile->conclude("Read " + toString(_predefinedWindows->size()) + " of cumulative length " + toString(_predefinedWindows->length()) + " bp on " + toString(_predefinedWindows->getNumChromosomes()) + " chromosomes.");
	}
	_numWindowsOnChr = 0;
};

void TGenome_windows::_setParsingLimits(TParameters & params){
	//limit windows
	_skipWindows = params.getParameterIntWithDefault("skipWindows", 0);
	if(_skipWindows > 0) _logfile->list("Will skip the first " + toString(_skipWindows) + " windows per chromosome. (parameter 'skipWindows')");
	_limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows"))
		_logfile->list("Will limit analysis to the first " + toString(_limitWindows) + " windows per chromosome. (parameter 'limitWindows')");
	if(_limitWindows <= _skipWindows)
		throw "limitWindows has to be larger than skipWindows!";
};

void TGenome_windows::_setWindowFilters(TParameters & params){
	//filter for missing reference
	_maxMissing = params.getParameterDoubleWithDefault("maxMissing", 1.0);
	if(_maxMissing > 1.0) throw "maxMissing must be smaller or equal to 1.0!";
	_logfile->list("Will filter out windows with a missing data fraction > " + toString(_maxMissing) + ". (parameter 'maxMissing')");

	_maxRefN = params.getParameterDoubleWithDefault("maxRefN", 1.0);
	if(_maxRefN < 0.0 || _maxRefN > 1.0) throw "maxRefN must be within interval [0,1]!";
	_openReference();
	if(_maxRefN < 1.0 && !_reference.hasReference()) throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided! (use 'fasta' to provide a reference)";
	_logfile->list("Will filter out windows with a fraction of 'N' in reference > " + toString(_maxMissing) + ". (parameter 'maxRefN')");
};

void TGenome_windows::_setSiteFilters(TParameters & params){
	//depth filter
	_readUpToDepth = params.getParameterIntWithDefault("readUpToDepth", 100);
	_logfile->list("Will read data up to depth " + toString(_readUpToDepth) + " and ignore additional bases. (parameter 'readUpToDepth')");

	if(params.parameterExists("minDepth") || params.parameterExists("maxDepth")){
		_applyDepthFilter = true;
		unsigned int tmpInt;
		tmpInt = params.getParameterIntWithDefault("minDepth", 0);
		if(tmpInt < 0)
			throw "minDepth must be >= 0!";
		_minDepth = tmpInt;
		tmpInt = params.getParameterIntWithDefault("maxDepth", _readUpToDepth);
		if(tmpInt < _minDepth) throw "maxDepth must be >= minDepth!";
		_maxDepth = tmpInt;
		_readUpToDepth = _maxDepth + 1;
		_logfile->list("Will filter out sites with sequencing depth < " + toString(_minDepth) + " or > " + toString(_maxDepth) + ". (parameters 'minDepth', 'maxDepth')");
	} else {
		_applyDepthFilter = false;
		_minDepth = 0;
		_maxDepth = _readUpToDepth;
		_logfile->list("Will keep sites regardless of depth. (use 'minDepth' or 'maxDepth' to filter)");
	}

	//CpG filter
	if(params.parameterExists("filterCpG")){
		_filterCpG = true;
		_logfile->list("Will filter out CpG sites. (parameter 'filterCpG')");
		_openReference(true);
	} else {
		_filterCpG = false;
		_logfile->list("Will keep CpG sites. (use 'filterCpG' to remove)");
	}
};

void TGenome_windows::_setMasks(TParameters & params){
	//normal mask
	if(params.parameterExists("mask")){
		if(_windowsPredefined) throw "Masking is currently not implemented if windows are predefined from a BED file.";
		if(params.parameterExists("regions")) throw "Cannot use mask and regions at the same time.";
		_doMasking = true;
		std::string maskFile = params.getParameterString("mask");

		//limit sites?
		int siteLimit = -1;
		if(params.parameterExists("siteLimit")){
			siteLimit = params.getParameterInt("siteLimit");
			if(siteLimit < 0)
				throw "site limit cannot be smaller than 0!";
			_logfile->startIndent("Will mask the first " + toString(siteLimit) + " sites listed in BED file '" + maskFile + "':");
		} else {
			_logfile->startIndent("Will mask all sites listed in BED file '" + maskFile + "':");
		}
		_logfile->listFlush("Reading file ...");
		_mask = new BAM::TBedReaderWindows(maskFile, _windowSize, _chromosomes, siteLimit, _logfile);
		_logfile->done();
		_logfile->endIndent();
		//mask->print();
	} else _doMasking = false;

	//reverse masking
	if(params.parameterExists("regions")){
		if(_windowsPredefined) throw "Regions is currently not implemented if windows are predefined from a BED file.";
		_considerRegions = true;
		std::string regionsFile = params.getParameterString("regions");

		//limitSites
		int siteLimit = -1;
		if(params.parameterExists("siteLimit")){
			siteLimit = params.getParameterInt("siteLimit");
			if(siteLimit < 0)
				throw "site limit cannot be smaller than 0!";
			_logfile->startIndent("Will limit analysis to the first " + toString(siteLimit) + " sites listed in BED file '" + regionsFile + "':");
		} else {
			_logfile->startIndent("Will limit analysis to all regions listed in BED file '" + regionsFile + "' (parameter 'regions'):");
		}
		_logfile->listFlush("Reading file ...");
		_mask = new BAM::TBedReaderWindows(regionsFile, _windowSize, _chromosomes, siteLimit, _logfile);
		_logfile->done();
		_logfile->endIndent();
	} else _considerRegions = false;
};

void TGenome_windows::_openSiteSubset(const std::string paramName){
	if(_subset){
		throw "Site subset already initialized!";
	}

	std::string filename = _params->getParameterString(paramName, true);

	if(_considerRegions) throw "Site subsets (parameter '" + paramName + "') and regions (parameter 'regions') can not be used at the same time!";
	if(_doMasking) throw "Site subsets (parameter '" + paramName + "') and masks (parameter 'mask') can not be used at the same time!";

	_subset = std::make_unique<GenotypeLikelihoods::TSiteSubset>(filename, _windowSize, _logfile, false, _reference, _bamFile.chromosomes);
};

void TGenome_windows::_jumpToEnd(){
	_chromosomes.jumpToEnd();
};

void TGenome_windows::_restartChromosomes(GenotypeLikelihoods::TWindow_base & window){
	_chromosomes.begin();

	_moveChromosome(window);
};

void TGenome_windows::_moveChromosome(GenotypeLikelihoods::TWindow_base & window){
	//jump reader
	_oldAlignmentMustBeConsidered = false;

	//restart windows
	if(_hasWindowIndent){
		_logfile->removeIndent();
		_hasWindowIndent = false;
	}
	_windowNumber = 0;

	if(_windowsPredefined){
		//find next used chromosome with windows
		do {
			_predefinedWindows->setChr(_chromosomes.curName());
			_numWindowsOnChr = _predefinedWindows->getNumWindowsOnCurChr();
			if(_numWindowsOnChr < 1 || _chromosomes.curInUse() == false){
				if(_chromosomes.curInUse())
					_logfile->conclude("No windows on chromosome " + _chromosomes.curName() + ".");
				_chromosomes.next();
				if(_chromosomes.end()){
					return;
				}

				_predefinedWindows->setChr(_chromosomes.curName());
				_numWindowsOnChr = _predefinedWindows->getNumWindowsOnCurChr();
			}
		} while((_numWindowsOnChr < 1 || !_chromosomes.curInUse()) && !_chromosomes.end());

		//now jump
		window.move(_predefinedWindows->curWindowStart(), _predefinedWindows->curWindowEnd(), _chromosomes.curRefID(), _logfile);
		window._chrName = _chromosomes.curName();
		_bamFile.jump(_chromosomes.curRefID(), window.startPos);

	} else {
		while(_chromosomes.curInUse() == false || _skipWindows * _windowSize > _chromosomes.curLength()){
			_chromosomes.next();
		}
		_numWindowsOnChr = ceil(_chromosomes.curLength() / (double) _windowSize);

		uint32_t curStart = _skipWindows * _windowSize;
		_bamFile.jump(_chromosomes.curRefID(), curStart);
		uint32_t nextEnd = curStart + _windowSize;

		if(nextEnd > _chromosomes.curLength()){
			nextEnd = _chromosomes.curLength();
		}
		window.move(curStart, nextEnd, _chromosomes.curRefID(), _logfile);
		window._chrName = _chromosomes.curName();
	}

	if(_chromosomes.end())
		return;

	//advance mask
	if(_doMasking || _considerRegions) _mask->setChr(_chromosomes.curName());
	if(_subset) _subset->setChr(_chromosomes.curName());

	//write progress
	if(_chromosomes.curRefID() > 0) _logfile->endIndent();
	_logfile->startNumbering("Parsing chromosome '" + _chromosomes.curName() + "':");
};

bool TGenome_windows::_moveToNextWindowOnChr(GenotypeLikelihoods::TWindow_base & window){
	//if sites defined
	int counter = 0;
	do{
		//move possible?
		++_windowNumber;
		++counter;
	} while(_subset && !_subset->hasPositionsInWindow(window.endPos) && window.endPos + window.length * counter < _chromosomes.curLength());

	if(window.endPos >= _chromosomes.curLength() || _windowNumber >= _limitWindows)
		return false;

	//calculate new end
	long nextEnd = window.endPos + _windowSize;
	if(nextEnd > _chromosomes.curLength())
		nextEnd = _chromosomes.curLength();
	window.move(window.endPos, nextEnd, _chromosomes.curRefID(), _logfile);

	return true;
};

bool TGenome_windows::_moveToNextPredefinedWindow(GenotypeLikelihoods::TWindow_base & window){
	++_windowNumber;
	if(_windowNumber >= _limitWindows)
		return false;
	if(_predefinedWindows->nextWindow()){
		window.move(_predefinedWindows->curWindowStart(), _predefinedWindows->curWindowEnd(), _chromosomes.curRefID(), _logfile);

		//should we jump or are we already close enough to next window
		if(_bamFile.curPosition() > window.startPos || _bamFile.curPosition() < window.startPos - _bamFile.maxReadLength()){
			if(window.startPos < _bamFile.maxReadLength())
				_bamFile.jump(_chromosomes.curRefID(), 0);
			else{
				_bamFile.jump(_chromosomes.curRefID(), window.startPos - _bamFile.maxReadLength());
			}
		}
		return true;
	} else
		return false;
};

bool TGenome_windows::_moveWindow(GenotypeLikelihoods::TWindow_base & window){
	//returns false when end of genome is reached
	if(_windowsPredefined){
		//if at beginning of BAM file
		if(_chromosomes.end()){
			_restartChromosomes(window);

			if(_chromosomes.end())
				throw "found no predefined windows in BED file! Does file exist?";
			_chrChangedWindow = true;

		} else {
			//now move coordinates of next window
			if(!_moveToNextPredefinedWindow(window)){
				//no more windows left on chr
				_chromosomes.next();

				if(_chromosomes.end()){
					if(_hasWindowIndent){
						_logfile->removeIndent();
						_hasWindowIndent = false;
					}
					return false;
				}

				_moveChromosome(window);
				_chrChangedWindow = true;

				if(_chromosomes.end()){
					if(_hasWindowIndent){
						_logfile->removeIndent();
						_hasWindowIndent = false;
					}
					return false;
				}
				++_windowNumber;
			} else
				//was able to move to next window on chr
				_chrChangedWindow = false;
		}

	} else {
		//if at beginning of BAM file
		if(_chromosomes.end()){
			_restartChromosomes(window);
			_chrChangedWindow = true;
		} else {
			if(!_moveToNextWindowOnChr(window)){
				//there is no window left on chr
				_chromosomes.next();

				//do we use this chromosome? if not, move on!
				while(!_chromosomes.end() && !_chromosomes.curInUse()){
					_chromosomes.next();
				}

				//did we reach end?
				if(_chromosomes.end()){
					window.endPos = 0;
					if(_hasWindowIndent){
						_logfile->removeIndent();
						_hasWindowIndent = false;
					}
					return false;
				}
				_moveChromosome(window);
				_chrChangedWindow = true;
			} else {
				_chrChangedWindow = false;
			}
		}
	}

	//report
	if(_hasWindowIndent) _logfile->removeIndent();
	_logfile->number("Window [" + toString(window.startPos) + ", " + toString(window.endPos) + ") of " + toString(_numWindowsOnChr) + " on '" + _chromosomes.curName() + "':");
	_logfile->addIndent();
	_hasWindowIndent = true;

	return true;
};

//---------------------
//read data in windows
//---------------------
bool TGenome_windows::_readDataInNextWindow(GenotypeLikelihoods::TWindow & window){
	_windowTimer.start();

	//move window
	if(!_moveWindow(window)){
		return false;
	}

	//read data
	_readAlignmentsIntoWindow(window);
	return true;
};

void TGenome_windows::_readAlignmentsIntoWindow(GenotypeLikelihoods::TWindow & window){
	//measure runtime
	_logfile->listFlushTime("Reading data ...");

	//check if old alignment is to be used.
	if(_oldAlignmentMustBeConsidered){
		if(_bamFile.curPosition() >= window.endPos){
			_logfile->warning("Old alignment is after window!");
			return;
		}

		_oldAlignmentMustBeConsidered = false;
		if(_oldAlignment->lastAlignedPositionWithRespectToRef() >= window.startPos)
			_oldAlignment = window.swapUsedForEmptyAlignment(_oldAlignment);
	}

	//read alignments
	int counter = 0;
	while(_readAlignment()){
		//fill alignment
		//return false if a post-parsing filer is not passed
		if(!_fillAlignment(*_oldAlignment)){
			continue;
		}

		++counter;

		//check if alignment starts after current window end -> break
		if(_oldAlignment->position >= window.endPos || _oldAlignment->refID != window.refId()){
			_oldAlignmentMustBeConsidered = true;
			break;
		}

		//check if alignment contains part of the window
		//if read continues outside of window, this is dealt with by window object
		if(_oldAlignment->position >= window.startPos || _oldAlignment->lastAlignedPositionWithRespectToRef >= window.startPos){
			_oldAlignment = window.swapUsedForEmptyAlignment(_oldAlignment);
		}
	}

	//fill sites
	if(_subset){
		window.fillSitesSubset(*_subset, _readUpToDepth);
		window.addReferenceBaseToSites(*_subset);
	} else {
		window.fillSites(_readUpToDepth);
		window.addReferenceBaseToSites(_reference);
	}

	//report
	_logfile->doneTime();

	//apply filters
	_applyWindowFilters(window);
};

void TGenome_windows::_applyWindowFilters(GenotypeLikelihoods::TWindow_base & window){
	window._passedFilters = false;
	if(window._numReadsInWindow > 0){
		//apply masks and filters
		if(_doMasking){
			_logfile->listFlush("Masking sites ...");
			window.applyMask(_mask, _considerRegions);
			_logfile->done();
		} else if(_considerRegions){
			_logfile->listFlush("Masking sites outside regions ...");
			window.applyMask(_mask, _considerRegions);
			_logfile->done();
		}

		//filter sites
		if(_applyDepthFilter){
			window.applyDepthFilter(_minDepth, _maxDepth);
		}
		if(_filterCpG){
			window.maskCpG(_reference);
		}

		//report
		window.dataSummary(_logfile);

		//filter window
		window.filter(_maxMissing, _maxRefN, _logfile);
	} else {
		_logfile->conclude("No data in this window.");
	}
};

void TGenome_windows::_traverseBAMWindows(){
	//iterate through windows
	while(_readDataInNextWindow(_window)){
		if(_window.passedFilters()){
			if(_subset){
				_subset->setChr(_chromosomes.curName());
			}

			//do stuff in derived classes
			_handleWindow();

			//report end of window calculations
			_logfile->list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
		}
	}
};

}; //end namespace
