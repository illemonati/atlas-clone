/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */

#include "TGenome.h"


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

void TGenome_parsed::_traverseBAMPassedQC(){
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
		_predefinedWindows = new TBed(tmp);
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
		_mask = new BAM::TBedReader(maskFile, _windowSize, _chromosomes, siteLimit, _logfile);
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
		_mask = new BAM::TBedReader(regionsFile, _windowSize, _chromosomes, siteLimit, _logfile);
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

	_subset = std::make_unique<TSiteSubset>(filename, _windowSize, _logfile, false, _reference, _bamFile.chromosomes);
};

void TGenome_windows::_jumpToEnd(){
	_chromosomes.jumpToEnd();
};

void TGenome_windows::_restartChromosomes(TWindow_base & window){
	_chromosomes.begin();

	_moveChromosome(window);
};

void TGenome_windows::_moveChromosome(TWindow_base & window){
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

bool TGenome_windows::_moveToNextWindowOnChr(TWindow_base & window){
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

bool TGenome_windows::_moveToNextPredefinedWindow(TWindow_base & window){
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

bool TGenome_windows::_moveWindow(TWindow_base & window){
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
bool TGenome_windows::_readDataInNextWindow(TWindow & window){
	_windowTimer.start();

	//move window
	if(!_moveWindow(window)){
		return false;
	}

	//read data
	_readAlignmentsIntoWindow(window);
	return true;
};

void TGenome_windows::_readAlignmentsIntoWindow(TWindow & window){
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
		if(_oldAlignment->position >= window.endPos || _oldAlignment->refID != window.chrNumber){
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

void TGenome_windows::_applyWindowFilters(TWindow_base & window){
	window.passedFilters = false;
	if(window.numReadsInWindow > 0){
		//apply masks and filters
		if(_doMasking){
			_logfile->listFlush("Masking sites ...");
			window.applyMask(_mask, _considerRegions);
			_logfile->done();
		} else if(_considerRegions){
			_logfile->listFlush("Masking sites outside regions ...");
			window.applyMask(_mask, _considerRegions);
			_logfile->done();
		} if(_applyDepthFilter){
			window.applyDepthFilter(_minDepth, _maxDepth);
		} if(_maxRefN < 1.0 && _reference.hasReference() == true){
			window.calcFracN();
		}

		//calc sequencing depth
		window.calcDepth();

		//report
		_logfile->conclude("read data from " + toString(window.numReadsInWindow) + " reads.");
		_logfile->conclude("sequencing depth is " + toString(window.depth));
		_logfile->conclude(toString(window.fractionDepthAtLeastTwo * 100) + "% of all sites are covered at least twice");
		_logfile->conclude(toString(window.fractionSitesNoData * 100) + "% of all sites have no data");
		if(window.fractionSitesNoData > _maxMissing){
			_logfile->conclude("Level of missing data > threshold of " + toString(_maxMissing) + " -> skipping this window");
			return;
		}
		if(_maxRefN < 1.0 && _reference.hasReference() == true){
			_logfile->conclude(toString(window.fractionRefIsN * 100) + "% of all reference bases are 'N'");
			if(window.fractionRefIsN > _maxRefN){
				_logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(_maxRefN) + " -> skipping this window");
				return;
			}
		}
		window.passedFilters = true;
	} else {
		_logfile->conclude("No data in this window.");
	}
};

void TGenome_windows::_traverseBAMWindows(){
	//iterate through windows
	while(_readDataInNextWindow(_window)){
		if(_window.passedFilters){
			if(_subset){
				_subset->setChr(_chromosomes.curName());
			}

			//do stuff in derived classes
			_handleWindow();

			//report end of window calculations
			_logfile->list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
		} else _logfile->list("No relevant positions -> skipping this window.");
	}
};



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// OLD STUFF BELOW
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//-----------------------------------------------------
//Functions for theta estimation
//-----------------------------------------------------
bool TGenomeWindows::estimateTheta(TThetaEstimator & thetaEstimator, TWindow_base & window){
	//adding sites to estimator
	_logfile->listFlush("Calculating emission probabilities ...");
	thetaEstimator.clear();
	window.addSitesToThetaEstimator(thetaEstimator.pointerToDataContainer(), genotypeLikelihoodCalculator);
	_logfile->done();

	//estimate Theta
	return thetaEstimator.estimateTheta();
};

void TGenomeWindows::estimateTheta(TParameters & params){
	//Theta estimator
	TThetaEstimator thetaEstimator(params, _logfile, _randomGenerator);

	//open output
	std::string filename = _outputName + "_theta_estimates.txt.gz";
	TThetaOutputFile thetaOut(&thetaEstimator, filename, _logfile);

	//check for which segements theta is to be estimated
	if(params.parameterExists("thetaGenomeWide") || alignmentParser._considerRegions){
		if(params.parameterExists("thetaGenomeWide"))
			_logfile->startIndent("Estimating theta genome-wide:");
		else _logfile->startIndent("Estimating theta at specific sites:");

		//HACK!!
		bool onlyBootstrap = params.parameterExists("onlyBootstrap");

		int numBootstraps = 0;
		if(params.parameterExists("bootstraps")){
			int numBootstraps = params.getParameterInt("bootstraps");
			estimateThetaGenomeWide(thetaEstimator, thetaOut, onlyBootstrap, numBootstraps);
			bootstrapTetaEstimation(numBootstraps, thetaEstimator);
		} else
			estimateThetaGenomeWide(thetaEstimator, thetaOut, onlyBootstrap, numBootstraps);

		_logfile->endIndent();
	} else {
		estimateThetaWindows(thetaEstimator, thetaOut, params.parameterExists("printAll"));
	}

	//clean up
	thetaOut.close();
};

void TGenomeWindows::estimateThetaWindows(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool printAll){
	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters || printAll){
			TTimer timer;

			_logfile->startIndent("Estimating Theta:");
			if(estimateTheta(thetaEstimator, window) || printAll){
				out.write(window._chrName, window.startPos, window.endPos);
			}
			_logfile->endIndent();

			//finish
			_logfile->list("Total computation time for this window was ", timer.seconds(), "s");
		} else _logfile->list("No relevant positions -> skipping this window.");
	}
};

void TGenomeWindows::estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool onlyReadData, int numBootstraps){
	if(alignmentParser._considerRegions)
		_logfile->startIndent("Estimating theta at specific sites:");

	//prepare windows
	TWindow window;
	thetaEstimator.clear();

	//add sites to estimator
	_logfile->startIndent("Adding sites to data structure:");
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters){
			//adding sites to estimator
			_logfile->listFlushTime("Calculating emission probabilities ...");
			try{
				window.addSitesToThetaEstimator(thetaEstimator.pointerToDataContainer(), genotypeLikelihoodCalculator);
			} catch(...){
				throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size, selecting fewer regions or limiting to sites with a minimal depth (>=2 recommended).";
			}
			_logfile->doneTime();
		}
	}
	_logfile->endIndent();

	//estimate Theta
	//HACK!!!!
	if(!onlyReadData){
		_logfile->startIndent("Estimate theta based on a total of " + toString(thetaEstimator.sizeWithData()) + " sites:");
		thetaEstimator.estimateTheta();
	}

	if(alignmentParser._considerRegions)
		out.write("regions", "-", "-");
	else
		out.write("genome-wide", "-", "-");
	if(numBootstraps == 0)
		thetaEstimator.clear();
};

void TGenomeWindows::bootstrapTetaEstimation(int numBootstraps, TThetaEstimator & thetaEstimator){
	if(numBootstraps < 1) throw "Number of bootstraps must be > 1!";
	_logfile->startIndent("Generating " + toString(numBootstraps) + " bootstrap estimates of theta:");

	//measure runtime
	TTimer timer;

	//open output file
	TOutputFileZipped bootstrapOut;
	std::string bootstrapFilename = _outputName + "_theta_bootstraps.txt.gz";
	_logfile->list("Writing theta bootstraps to '" + bootstrapFilename + "'");
	bootstrapOut.open(bootstrapFilename.c_str());

	//write header
	bootstrapOut.setPrecision(9);
	std::vector<std::string> header = {"Bootstrap"};
	thetaEstimator.addToHeader(header);
	bootstrapOut.writeHeader(header);

	//loop over bootstraps
	for(int s=0; s<numBootstraps; ++s){
		_logfile->startIndent("Bootstrap " + toString(s+1) + " of " + toString(numBootstraps) + ":");

		//run bootstrap
		bootstrapOut << s+1;
		thetaEstimator.bootstrapTheta(&bootstrapOut);

		_logfile->endIndent();
	}

	//finish
	_logfile->list("Total computation time for theta bootstrapping was ", timer.minutes());
	_logfile->endIndent();
};

void TGenomeWindows::calcThetaLikelihoodSurfaces(TParameters & params){
	//read params
	int steps = params.getParameterIntWithDefault("steps", 100);

	//prepare windows
	TWindow window;

	//Theta estimator
	TThetaEstimator estimator(_logfile, _randomGenerator);

	//iterate through windows
	std::string filename;
	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//read data for current window
		if(window.passedFilters){
			//check if we have data -> can be extended to ensure
			_logfile->startIndent("Calculating likelihood surface for Theta:");

			//measure runtime
			TTimer timer;

			//adding sites to estimator
			_logfile->listFlush("Calculating emission probabilities ...");
			window.addSitesToThetaEstimator(estimator.pointerToDataContainer(), genotypeLikelihoodCalculator);
			_logfile->done();

			//open file
			gz::ogzstream out;
			filename = _outputName + window._chrName + "_" + toString(window.startPos) + "_LLsurface.txt";
			out.open(filename.c_str());
			if(!out) throw "Failed to open output file '" + _outputName + "'!";

			//estimate Theta
			_logfile->listFlush("Calculating likelihood surface ...");
			estimator.calcLikelihoodSurface(out, steps);
			_logfile->done();

			//clear theta estimator
			estimator.clear();

			//finish
			out.close();
			_logfile->list("Total computation time for this window was ", timer.seconds(), "s");
			_logfile->endIndent();
		}
	}
};

void TGenomeWindows::estimateThetaRatio(TParameters & params){
	//Theta estimator
	TThetaEstimatorRatio thetaEstimatorRatio(params, _logfile, _randomGenerator);

	//read the two regions to be used
	_logfile->startIndent("Reading regions:");

	int windowSize = alignmentParser.getWindowSize();

	//region 1
	std::string regionsFile1 = params.getParameterString("regions1");
	int siteLimit = -1;
	if(params.parameterExists("siteLimit1")){
		siteLimit = params.getParameterInt("siteLimit1");
		if(siteLimit < 0)
			throw "site limit cannot be smaller than 0!";
		_logfile->startIndent("Reading first " + toString(siteLimit) + " sites from regions 1 BED file '" + regionsFile1 + "':");
	} else {
		_logfile->listFlush("Reading regions 1 from BED file '" + regionsFile1 + "' ...");
	}
	TBedReader region1(regionsFile1, windowSize, bamFile.chromosomes, siteLimit,_logfile);
	_logfile->done();

	//region 2
	siteLimit = -1;
	std::string regionsFile2 = params.getParameterString("regions2");
	if(params.parameterExists("siteLimit2")){
		siteLimit = params.getParameterInt("siteLimit2");
		if(siteLimit < 0)
			throw "site limit cannot be smaller than 0!";
		_logfile->startIndent("Reading first " + toString(siteLimit) + " sites from regions 1 BED file '" + regionsFile1 + "':");
	} else {
		_logfile->listFlush("Reading regions 2 from BED file '" + regionsFile2 + "' ...");
	}
	TBedReader region2(regionsFile2, windowSize, bamFile.chromosomes, siteLimit, _logfile);
	_logfile->done();
	_logfile->endIndent();

	//prepare windows
	TWindow window;

	//add sites to estimator
	_logfile->startIndent("Adding sites to data structures:");
	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		region1.setChr(window._chrName);
		region2.setChr(window._chrName);
		if(window.passedFilters){
			//adding sites to estimator
			_logfile->listFlushTime("Calculating emission probabilities for sites within defined regions ...");
			try{
				window.addSitesToThetaEstimator(thetaEstimatorRatio.pointerToDataContainer(), genotypeLikelihoodCalculator, region1);
				window.addSitesToThetaEstimator(thetaEstimatorRatio.pointerToDataContainer2(), genotypeLikelihoodCalculator, region2);
			} catch(...){
				throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size, selecting fewer regions or limiting to sites with a minimal depth (>=2 recommended).";
			}
			_logfile->doneTime();
		}
	}
	_logfile->endIndent();

	//estimate Theta ratio
	thetaEstimatorRatio.estimateRatio(_outputName);
};

void TGenomeWindows::performDownsamplingThetaQC(TParameters & params){
	//parse downsampling parameters
	std::vector<double> downSampleProbVector;
	fillVectorOfDownsamplingProbabilities(params.getParameterStringWithDefault("prob", "1.0,0.5,0.2,0.1,0.05,0.02,0.01"), downSampleProbVector);

	//report probabilities
	_logfile->list("Will estimate theta after downsampling reads with probabilities (parameter 'prob'): " + concatenateString(downSampleProbVector, ", "));

	//check if full data is to be used (i.e. if prob = 1.0 is specified)
	bool printFull = false;
	if(*downSampleProbVector.begin() == 1.0){
		printFull = true;
		downSampleProbVector.erase(downSampleProbVector.begin());
	}

	//create window and estimator for full data
	TThetaEstimator estimator_full(params, _logfile, _randomGenerator);
	TWindow window_full;

	//create windows and estimators for downsamples
	std::vector<TThetaEstimator> estimators;
	std::vector<TWindow_base*> windows;
	if(downSampleProbVector.size() > 0){
		for(size_t i=0; i<downSampleProbVector.size(); ++i){
			estimators.emplace_back(estimator_full);
			windows.push_back(new TWindow_base());
		}
	}

	//open output
	std::string filename = _outputName + "_theta_estimates.txt.gz";
	TThetaOutputFile thetaOut;
	if(printFull){
		thetaOut.addEstimator(&estimator_full, "p1.0_");
	}
	for(size_t i=0; i<downSampleProbVector.size(); ++i){
		//assemble prefix without lagging zeros
		std::string prefix = toString(downSampleProbVector[i]);
		int pos = prefix.size() - 1;
		while(prefix[pos] == '0') pos--;
		prefix.erase(pos + 1, prefix.size() - pos - 1);
		prefix = "p" + prefix + "_";

		//now add estimator ot output file
		thetaOut.addEstimator(&estimators[i], prefix);
	}
	thetaOut.open(filename, _logfile);

	//messure runtime
	TTimer timer;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window_full)){
		if(window_full.passedFilters){
			timer.startPos();

			//estimate on full data
			if(printFull){
				_logfile->startIndent("Estimating Theta on full data:");
				estimateTheta(estimator_full, window_full);
				_logfile->endIndent();
			}

			if(windows.size() > 0){
				for(size_t i = 0; i<downSampleProbVector.size(); ++i){
					_logfile->startIndent("Estimating Theta on downsampled data (p = " + toString(downSampleProbVector[i]) + "):");

					//downsample
					alignmentParser.downsampleWindow(*windows[i], window_full, downSampleProbVector[i], _randomGenerator);

					//and estimate
					estimateTheta(estimators[i], *windows[i]);
					_logfile->endIndent();
				}
			}

			//write output
			thetaOut.write(window_full._chrName, window_full.startPos, window_full.endPos);

			//finish
			_logfile->list("Total computation time for this window was ", timer.seconds(), "s");
			_logfile->endIndent();
		} else _logfile->list("No relevant positions -> skipping this window.");
	}

	//clean up
	thetaOut.close();
};

//---------------------------------------------------
// I/O
//---------------------------------------------------

void TGenomeWindows::generatePSMCInput(TParameters & params){
	//read in parameters required
	double theta = params.getParameterDoubleWithDefault("theta", 0.001);
	_logfile->list("Using theta = " + toString(theta));
	TThetaEstimator thetaEstimator(_logfile, _randomGenerator);
	thetaEstimator.setTheta(theta);

	double confidence = params.getParameterDoubleWithDefault("confidence", 0.99);
	_logfile->list("Calling heterozygosity state with confidence > " + toString(confidence) + ". (parameter 'confidence')");
	int blockSize = params.getParameterIntWithDefault("block", 100);
	//make sure window size is a multiple of block length!
	if(alignmentParser.getWindowSize() % blockSize != 0) throw "Window size is not a multiple of block size!";

	//open output file
	std::ofstream output;
	std::string outputFileName = _outputName + ".psmcfa";
	_logfile->list("Writing PSMC input file to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int nCharOnLine = 0;

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(alignmentParser._chrChangedWindow){
			if(nCharOnLine > 0) output << '\n';
				output << '>' << bamFile.chromosomes.curName() << '\n';
		}
		if(window.passedFilters){
			//create PSMC input
			_logfile->listFlushTime("Estimating heterozygosity status ...");
			window.generatePSMCInput(thetaEstimator, genotypeLikelihoodCalculator, blockSize, confidence, output, nCharOnLine);
			_logfile->doneTime();
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
}

void TGenomeWindows::createDepthMask(TParameters & params){
	int minDepthForMask = params.getParameterInt("minDepthForMask");
	int maxDepthForMask = params.getParameterInt("maxDepthForMask");
	if(params.parameterExists("maxDepth") || params.parameterExists("minDepth"))
		throw "Cannot mask sites for sequencing depth while creating the mask!";

	std::ofstream output;
	std::string outputFileName = _outputName + "_minDepth"+ toString(minDepthForMask) + "_maxDepth" + toString(maxDepthForMask) + "_depthMask.bed";
	_logfile->list("Writing sites with depth < " + toString(minDepthForMask) + " and with depth > " + toString(maxDepthForMask) + " to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	//prepare windows
	TWindow window;

	//iterate through windows
	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//read data for current window
		if(window.passedFilters){
			_logfile->listFlush("Writing sites to mask to output file ...");
			window.createDepthMask(minDepthForMask, maxDepthForMask, output, window._chrName);
			_logfile->done();
		}
	}
	output.close();
}

//---------------------------------------------------
// Recalibration
//---------------------------------------------------
void TGenomeWindows::estimateErrorCalibrationEM(TParameters & params){
	if(alignmentParser.qualitiesScoresAreRecalibrated() && !params.parameterExists("rerecalibrate"))
		throw "Can not estimate recalibration: quality scores are already recalibrated while reading! (Use argument 'rerecalibrate' to overwrite this error)";

	//initialize maps
	std::string poolReadGroupsFile = params.getParameterString("poolReadGroups", false);
	TReadGroupMap readGroupMap(&(alignmentParser.readGroups), poolReadGroupsFile, _logfile);
	TQualityMap qualityMap;

	//create recalibration object
	bool writeTmpTables = false;
	if(params.parameterExists("writeTmpTables")){
		writeTmpTables = true;
		_logfile->list("Will write intermediate estimates of EM and Newton-Raphson to file.");
	}
	GenotypeLikelihoods::TRecalibrationEMEstimator recalObjectEM(params, &alignmentParser.readGroups, _logfile, &readGroupMap);

	//prepare windows
	TWindow window;

	//add sites to EM object
	if(params.parameterExists("alleles")){
		//Limit to sites with known alleles
		_logfile->startIndent("Will limit analysis to homozygous sites with known genotypes:");
		int windowSize = alignmentParser.getWindowSize();
		TSiteSubset subset(params.getParameterString("alleles"), windowSize, _logfile, true);

		//now parse through windows
		while(alignmentParser.readDataInNextWindow(window)){
			//make sure subset is on the correct chromosome
			subset.setChr(window._chrName);

			//read data for current window
			if(window.passedFilters && subset.hasPositionsInWindow(window.startPos))
				window.addToRecalibrationEM(recalObjectEM, subset, qualityMap);
			else _logfile->list("No positions in this window.");
		}
		_logfile->endIndent();

		//clean up memory
		window.clear();

		//now estimate
		recalObjectEM.performEstimationKnownGenotypes(_outputName, writeTmpTables);

	} else {
		_logfile->startIndent("Reading data from windows:");
		while(alignmentParser.readDataInNextWindow(window)){
			//read data for current window
			if(window.passedFilters)
				window.addToRecalibrationEM(recalObjectEM, qualityMap);
			else _logfile->list("No positions in this window.");
		}
		_logfile->endIndent();
		_logfile->endIndent();

		//clean up memory
		window.clear();

		//now estimate
		recalObjectEM.performEstimation(_outputName, writeTmpTables);
	}
};

void TGenomeWindows::calculateLikelihoodErrorCalibrationEM(TParameters & params){
	//create recalibration object
	TReadGroupMap readGroupMap(&alignmentParser.readGroups, params.getParameterString("poolReadGroups", false), _logfile);
	GenotypeLikelihoods::TRecalibrationEMEstimator recalObjectEM(params, &alignmentParser.readGroups, _logfile, &readGroupMap);
	recalObjectEM.initializeFromFile(params.getParameterString("recalForLL"));

	//prepare windows
	TWindow window;

	//other tmp variables
	TQualityMap qualMap;

	//add sites to EM object
	_logfile->startIndent("Reading data from windows:");
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters)
			window.addToRecalibrationEM(recalObjectEM, qualMap);
	}
	//clean up memory
	window.clear();
	_logfile->endIndent();

	_logfile->list("LL = " + toString(recalObjectEM.calcLL()));
};



//---------------------------------------------------
//BAM manipulation / statistics
//---------------------------------------------------

void TGenomeWindows::writeNonConservedBed(TParameters & params){
	if(!alignmentParser.hasReference){
		throw "Must provide reference to create nonRef mask";
	}

	//prepare windows
	TWindow window;

	//prepare output file
	std::ofstream output;
	std::string outputFileName = _outputName + "_nonRefMask.bed";
	_logfile->list("Writing mask to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			_logfile->listFlush("Writing invariant site coordinates to BED file ...");
			window.writeNonConservedBed(output);
			_logfile->write(" done!");
		}
	}
};

void TGenomeWindows::estimateApproximateDepthPerWindow(TParameters & params){
	//open output file
	std::ofstream output;
	std::string outputFileName = _outputName + "_depthPerWindow.txt";
	_logfile->list("Writing sequencing depth estimates to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int nCharOnLine = 0;

	//write header
	output << "chr\tstart\tend\tdepth" << std::endl;

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			//write to file
			_logfile->listFlush("Writing sequencing depth estimates to file ...");
			if(window.depth == -1.0) output << alignmentParser.getCurChrName() << "\t" << window.startPos << "\t" << window.endPos << "\t" << "0" << "\n";
			else output << alignmentParser.getCurChrName() << "\t" << window.startPos << "\t" << window.endPos << "\t" << window.depth << "\n";
			_logfile->done();
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
};

void TGenomeWindows::estimateDepthPerSite(TParameters & params){
	//initialize count object
	int maxDepth = params.getParameterIntWithDefault("maxDepth", 20);
	TDistributionOfCounts counts(maxDepth, "depth");

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			_logfile->listFlush("Adding depth to table ...");
			window.countDepthPerSite(counts);
			_logfile->done();
		}
	}

	//write
	std::string filename = _outputName + "_distDepthPerSite.txt";
	_logfile->listFlush("Writing depth distribution to '" + filename + "' ...");
	counts.writeCounts(filename);
	_logfile->done();

	filename = _outputName + "_cumulativeDepthPerSite.txt";
	_logfile->listFlush("Writing normalized cumulative depth distribution to '" + filename + "' ...");
	counts.writeNormalizedCumulativeCounts(filename);
	_logfile->done();

	filename = _outputName + "_quantilesDepthPerSite.txt";
	_logfile->listFlush("Writing quantiles of depth distribution '" + filename + "' ...");
	counts.writeQuantiles(filename);
	_logfile->done();
};

void TGenomeWindows::writeDepthPerSite(TParameters & params){
	gz::ogzstream out;

	std::string outputFileName = _outputName + "_depthPerSite.txt.gz";
	_logfile->list("Writing per site depth to '" + outputFileName + "'");

	out.open(outputFileName.c_str());
	if(!out) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	out << "chr\tpos\tdepth" << std::endl;

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			_logfile->listFlush("Writing depth per site ...");
			window.printDepthPerSite(out, alignmentParser.getCurChrName());
			_logfile->done();
		}
	}

	//clean up
	out.close();
};



void TGenomeWindows::printMateInformationPerSite(TParameters & params){
	//open output file
	std::string outputFileName = _outputName + "_mateInformation.txt.gz";
	_logfile->list("Writing mate information to file '" + outputFileName + "'.");
	TOutputFileZipped out(outputFileName);
	out.writeHeader({"chr", "pos", "depth", "numA", "numC", "numG", "numT", "numAlleles", "numFirstMate", "numSecondMate", "numFwd", "numRev"});

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			_logfile->listFlush("Writing mate info per site ...");
			window.printMateInformationPerSite(out, alignmentParser.getCurChrName());
			_logfile->done();
		}
	}

	//clean up
	out.close();
};



void TGenomeWindows::testGenotypeLikelihoods(TParameters & params){
	//create
	GenotypeLikelihoods::TGenotypeLikelihoodCalculator glcalc(params, &alignmentParser.readGroups, _logfile);

	//prepare windows
	TWindow window;
	GenotypeLikelihoods::TGenotypeLikelihoods genolik;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//for(int i=0; i<window.length; ++i){
		for(int i=0; i<3; ++i){

			std::cout << std::endl;
			std::cout << window._chrName << "\t" << window.startPos + i + 1 << "\t" << window.sites[i].getBases() << std::endl;

			//window.sites[i].calcEmissionProbabilities();
			glcalc.calculateGenotypeLikelihoods(window.sites[i].bases, genolik);


			for(uint8_t g=0; g<10; ++g){
				//std::cout << "\t" <<  alignmentParser.genoMap.getGenotypeString(g) << ": " << window.sites[i].emissionProbabilities[g] << ", " << genolik.at(g);
				std::cout << "\t" <<  alignmentParser.genoMap.getGenotypeString(g) << ": "  << ", " << genolik.at(g);
			}
			std::cout << std::endl;
		}
	}
};

