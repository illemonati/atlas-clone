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
TGenome_basic::TGenome_basic(){
    _logfile = nullptr;
    _params = nullptr;
    _randomGenerator = nullptr;
};

TGenome_basic::TGenome_basic(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator){
    _initialize(Params, Logfile, RandomGenerator);
};

void TGenome_basic::_initialize(TParameters &Params, TLog *Logfile, TRandomGenerator *RandomGenerator) {
    _logfile = Logfile;
    _params = &Params;
    _randomGenerator = RandomGenerator;

    //open bam file
    //TODO: deal with index in better way: let tasks decide if they need an index or not
    _bamFile.open(Params.getParameter<std::string>("bam"), Params.parameterExists("indexNotRequired"), _logfile);
    _bamFile.setLimits(Params, _logfile);

    //outputname
    if(Params.parameterExists("out")){
        _outputName = Params.getParameter<std::string>("out");
        _logfile->list("Writing output files with prefix '" + _outputName + "'. (parameter 'out')");
    } else {
        //guess from BAM filename.
        _outputName = _bamFile.filename();
        _outputName = extractBeforeLast(_outputName, ".");
        _logfile->list("Writing output files with prefix '" + _outputName + "'. (specify with 'out')");
    }
};

void TGenome_basic::_openBamForWriting(const std::string filename, BAM::TOutputBamFile & outBam){
	outBam.open(filename, _bamFile);
};

//---------------------------------------------------------------
// TGenome_filtered
// A base class without genotype likelihoods but BAM filters enabled
//---------------------------------------------------------------
TGenome_filtered::TGenome_filtered():TGenome_basic(){};

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
	_genotypeLikelihoodCalculator.init(Params, &_bamFile.readGroupsMutable(), _logfile);

	//set parsing filters
	_setReadTrimming(Params);
	_qualityFilter.set(Params, _logfile);
	_contextFilter.set(Params, _logfile);
};

void TGenome_parsed::_openReference(bool required){
	if(!_reference.hasReference()){
		if(_params->parameterExists("fasta")){
			std::string fastaFile = _params->getParameter<std::string>("fasta");
			_logfile->list("Reading reference sequence from '" + fastaFile + "'. (parameter fasta)");
			_reference.initialize(fastaFile);
		} else {
			if(required){
				throw "No reference provided! (Use parameter fasta to provide a reference)";
			}
		}
	}
};

void TGenome_parsed::_setReadTrimming(TParameters & params){
	//trimming ends
	if(params.parameterExists("trim3") || params.parameterExists("trim5")){
		_trimmingLength3Prime = params.getParameterWithDefault<int>("trim3", 0);
		if(_trimmingLength3Prime < 0) throw "trimming distance trim3 must be >= 0!";
		_trimmingLength5Prime = params.getParameterWithDefault<int>("trim5", 0);
		if(_trimmingLength5Prime < 0) throw "trimming distance trim5 must be >= 0!";
		if(_trimmingLength3Prime > 0 || _trimmingLength5Prime > 0){
			_logfile->list("Will trim first " + toString(_trimmingLength3Prime) + " and " + toString(_trimmingLength5Prime) + " bases from the 3' and 5' end, respectively. (parameters 'trim3', 'trim5')");
		}
		_trimReads = true;
	} else {
		_trimmingLength3Prime = 0;
		_trimmingLength5Prime = 0;
		_trimReads = false;
	}
};

void TGenome_parsed::_parseAlignment(BAM::TAlignment & alignment){
	//parse
	alignment.parse(_genotypeLikelihoodCalculator.getSequencingErrorModels());

	//apply filters
	if(_trimReads){
		alignment.trimRead(_trimmingLength3Prime, _trimmingLength5Prime);
	}

	alignment.filter(_qualityFilter);
	alignment.filter(_contextFilter);

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
		_chromosomes(_bamFile.chromosomes()){
	//reading parameters regarding windows
	_logfile->startIndent("Parsing window settings:");
	_setWindowParameters(Params);
	_setParsingLimits(Params);
	_setWindowFilters(Params);
	_setMasks(Params);
	_setSiteFilters(Params);
	_logfile->endIndent();

        _curAlignment = new BAM::TAlignment;
};

TGenome_windows::~TGenome_windows(){
	delete _curAlignment;
};

void TGenome_windows::_setWindowParameters(TParameters & params){
	if(!params.parameterExists("window") && params.parameterExists("windows")){
		_logfile->warning("Argument 'windows' specified, but unknown. Did you mean 'window'?");
	}
	std::string tmp = params.getParameterWithDefault<std::string>("window", "1000000");

	//check if it is a number
	if(stringIsProbablyANumber(tmp)){
		_windowSize = convertString<int>(tmp);
		_logfile->list("Setting window size to " + toString(_windowSize) + ". (parameter 'window')");
		if(_windowSize < _bamFile.maxReadLength()){
			throw "Window size " + tmp + " out of range! Windows must be at least as large as the max read length (" + toString(_bamFile.maxReadLength()) + " bp). (use parameter 'maxReadLength' to change)!";
		}
	} else {
		_logfile->listFlush("Limiting analysis to windows defined in BED file '" + tmp + "' (parameter window) ...");
		_predefinedWindows.add(tmp, _chromosomes);
		_logfile->done();
		_logfile->conclude("Read " + toString(_predefinedWindows.size()) + " of cumulative length " + toString(_predefinedWindows.length()) + " bp on " + toString(_predefinedWindows.numChromosomesWithWindows()) + " chromosomes.");
	}
	_numWindowsOnChr = 0;
};

void TGenome_windows::_setParsingLimits(TParameters & params){
	//limit windows
	_skipWindows = params.getParameterWithDefault<int>("skipWindows", 0);
	if(_skipWindows > 0) _logfile->list("Will skip the first " + toString(_skipWindows) + " windows per chromosome. (parameter 'skipWindows')");
	_limitWindows = params.getParameterWithDefault<long>("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows"))
		_logfile->list("Will limit analysis to the first " + toString(_limitWindows) + " windows per chromosome. (parameter 'limitWindows')");
	if(_limitWindows <= _skipWindows)
		throw "limitWindows has to be larger than skipWindows!";
};

void TGenome_windows::_setWindowFilters(TParameters & params){
	//filter for missing reference
	_maxMissing = params.getParameterWithDefault<double>("maxMissing", 1.0);
	if(_maxMissing < 0.0 || _maxMissing > 1.0) throw "maxMissing must be within [0, 1]!";
	if(_maxMissing < 1.0){
		_logfile->list("Will filter out windows with a missing data fraction > " + toString(_maxMissing) + ". (parameter 'maxMissing')");
	} else {
		_logfile->list("Will keep windows regardless of missingness. (use 'maxMissing' to filter)");
	}

	_maxRefN = params.getParameterWithDefault<double>("maxRefN", 1.0);
	if(_maxRefN < 0.0 || _maxRefN > 1.0) throw "maxRefN must be within interval [0,1]!";
	_openReference();
	if(_maxRefN < 1.0 && !_reference.hasReference()) throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided! (use 'fasta' to provide a reference)";
	_logfile->list("Will filter out windows with a fraction of 'N' in reference > " + toString(_maxMissing) + ". (parameter 'maxRefN')");
};

void TGenome_windows::_setSiteFilters(TParameters & params){
	//depth filter
	_readUpToDepth = params.getParameterWithDefault<int>("readUpToDepth", 100);
	_logfile->list("Will read data up to depth " + toString(_readUpToDepth) + " and ignore additional bases. (parameter 'readUpToDepth')");

	if(params.parameterExists("minDepth") || params.parameterExists("maxDepth")){
		_applyDepthFilter = true;
		unsigned int tmpInt;
		tmpInt = params.getParameterWithDefault<int>("minDepth", 0);
		if(tmpInt < 0)
			throw "minDepth must be >= 0!";
		_minDepth = tmpInt;
		tmpInt = params.getParameterWithDefault<int>("maxDepth", _readUpToDepth);
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

	//downsample?
	_downsampleDepth = params.getParameterWithDefault<int>("downsample", 0);
	if(_downsampleDepth > 0){
		_logfile->list("Will downsample sites to a depth <= ", _downsampleDepth, ". (parameter 'downsample')");
		if(_downsampleDepth >= _maxDepth){
			_logfile->warning("Downsample depth is >= max depth: no downsampling will occur.");
		}
		subsamplePicker = std::make_unique<TSubsamplePicker>(_randomGenerator, 30);
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
	if(params.parameterExists("mask") || params.parameterExists("regions")){
		_doMasking = true;
		std::string filename;
		if(params.parameterExists("mask")){
			//mask
			if(params.parameterExists("regions")) throw "Cannot use mask and regions at the same time.";
			filename = params.getParameter<std::string>("mask");
			_logfile->startIndent("Will mask all sites listed in BED file '" + filename + "':");
			_considerRegions = false;
		} else {
			//regions
			filename = params.getParameter<std::string>("regions");
			_logfile->startIndent("Will limit analysis to sites listed in BED file '" + filename + "' (parameter 'regions'):");
			_considerRegions = true;
		}

		//read file
		_logfile->listFlush("Reading file ...");
		_mask.add(filename, _chromosomes);
		_logfile->done();
		_logfile->conclude("Read "+ toString(_mask.size()) + " sites on " + toString(_mask.numChromosomesWithWindows()) + " chromosomes.");
		_logfile->endIndent();
	} else {
		_doMasking = false;
		_considerRegions = false;
	}
};

void TGenome_windows::_openSiteSubset(const std::string paramName){
	if(_subset){
		throw "Site subset already initialized!";
	}

	std::string filename = _params->getParameter<std::string>(paramName, true);

	if(_considerRegions) throw "Site subsets (parameter '" + paramName + "') and regions (parameter 'regions') can not be used at the same time!";
	if(_doMasking) throw "Site subsets (parameter '" + paramName + "') and masks (parameter 'mask') can not be used at the same time!";

	_subset = std::make_unique<GenotypeLikelihoods::TSiteSubset>(filename, _bamFile.chromosomes(), _logfile, false, _reference);
};

void TGenome_windows::_setCountersBeginningOfChromosome(){
	_chrChangedWindow = true;
	_windowNumber = 1;
};

bool TGenome_windows::_incrementWindow(GenotypeLikelihoods::TWindow_base & window){
	//move to next
	window += _windowSize;
	++_windowNumber;

	//Move to next chromosome if 1) we are at begininning of BAM (_curchromosome at end), 2) we are beyond _curChromosome or 3) reached window limit
	if(_curChromosome == _chromosomes.cend() || window.from() >= _curChromosome->chrEnd || _windowNumber > _limitWindows){
		//move to next chromosome
		if(_curChromosome == _chromosomes.cend()){ // beginning of chromosome
			_curChromosome = _chromosomes.cbegin();
		} else {
			++_curChromosome;
		}

		//advance if this chromosome is not used
		while(_curChromosome != _chromosomes.cend() && (!_curChromosome->inUse || _skipWindows * _windowSize > _curChromosome->length)){
			++_curChromosome;
		}

		if(_curChromosome == _chromosomes.cend()){
			return false;
		}

		_setCountersBeginningOfChromosome();
		_numWindowsOnChr = ceil(_curChromosome->length / (double) _windowSize);

		//move window to beginning of chromosome
		BAM::TGenomePosition newFrom = _curChromosome->chrStart + _skipWindows * _windowSize;
		window.move(newFrom, _windowSize , _curChromosome->name);
	} else {
		_chrChangedWindow = false;
	}

	return true;
};

bool TGenome_windows::_moveToNextWindow(GenotypeLikelihoods::TWindow_base & window){
	//move to next
	if(!_incrementWindow(window)){
		return false;
	}

	//if sites defined, check if there are sites
	if(_subset){
		while(!_subset->hasPositionsInWindow(window)){
			if(!_incrementWindow(window)){
				return false;
			}
		}
	}

	//make sure window does not go beyond chromosome end
	if(window.to() > _curChromosome->chrEnd){
		window.resize(_curChromosome->chrEnd - window.from());
	}

	return true;
};

bool TGenome_windows::_incrementPredefinedWindow(){
	uint32_t oldRefID = _curPredefinedWindow->refID();

	++_curPredefinedWindow;
	if(_curPredefinedWindow == _predefinedWindows.end()){
		return false;
	}

	//check if chromosome changed
	if(_curPredefinedWindow->refID() != oldRefID){
		_setCountersBeginningOfChromosome();
	} else {
		++_windowNumber;
		_chrChangedWindow = false;
	}

	return true;
};

bool TGenome_windows::_moveToNextPredefinedWindow(GenotypeLikelihoods::TWindow_base & window){
	//if at beginning of BAM file: restart
	if(_curChromosome == _chromosomes.cend()){
		_curPredefinedWindow = _predefinedWindows.begin();
		_curChromosome = _chromosomes.cbegin(_curPredefinedWindow->refID());
		_setCountersBeginningOfChromosome();
	} else {
		_incrementPredefinedWindow();
	}

	//go to next accepted window
	while(_curPredefinedWindow != _predefinedWindows.end() && (!_chromosomes.inUse(_curPredefinedWindow->refID()) || _windowNumber < _skipWindows || _windowNumber >= _limitWindows)){
		_incrementPredefinedWindow();
	}

	//return false if at end
	if(_curPredefinedWindow == _predefinedWindows.end()){
		return false;
	}

	//make sure we are on the right chromosome
	if(_curChromosome->refID() != _curPredefinedWindow->refID()){
		_curChromosome = _chromosomes.cbegin(_curPredefinedWindow->refID());

		//update num windows per chromosome
		_numWindowsOnChr = _predefinedWindows.numWindowsOnChr(_curChromosome->refID());
	}

	//move coordinates of window
	window.move(*_curPredefinedWindow, _curChromosome->name);

	//move in BAM: should we jump or are we already close enough to next window
	if(_bamFile.refID() == window.refID()){
		//same chromosome: jump only if we are far away
		if(_bamFile.curPosition() > window.from() || _bamFile.curPosition() < window.from() - _bamFile.maxReadLength()){
			_bamFile.jump(window.from() - _bamFile.maxReadLength());
			_curAlignment->clear();
		}
	} else {
		//different chromosome: jump
		_bamFile.jump(window.from() - _bamFile.maxReadLength());
		_curAlignment->clear();
    }

	//return true as we continue reading
	return true;
};

bool TGenome_windows::_moveWindow(GenotypeLikelihoods::TWindow_base & window){
	//returns false when end of genome is reached
	if(_predefinedWindows.empty()){
		//no predefined windows: regular traversing
		if(!_moveToNextWindow(window)){
			return false;
		}
	} else {
		//there are predefined windows
		if(!_moveToNextPredefinedWindow(window)){
			return false;
		}
	}

	//report chromosome
	if(_hasWindowIndent) _logfile->removeIndent();
	if(_chrChangedWindow){
		if(_curChromosome->refID() > 0) _logfile->endIndent();
			_logfile->startNumbering("Parsing chromosome '" + _curChromosome->name + "':");
	}

	//report window
	_logfile->number("Window [" + toString(window.from().position()+1) + ", " + toString(window.to().position()) + "] of " + toString(_numWindowsOnChr) + " on '" + _curChromosome->name + "':");
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
		//reached end
		if(_hasWindowIndent){
			_logfile->removeIndent();
			_hasWindowIndent = false;
		}
		return false;
	}

	//read data
	_readAlignmentsIntoWindow(window);
	return true;
};

bool TGenome_windows::_readAndParseAlignment(BAM::TAlignment & alignment){
    if(!_bamFile.readNextAlignmentThatPassesFilters(alignment)){
    	return false;
    }

    _parseAlignment(alignment);
    return true;
};

void TGenome_windows::_readAlignmentsIntoWindow(GenotypeLikelihoods::TWindow & window){
	//measure runtime
	_logfile->listFlushTime("Reading data ...");

	//make sure oldAligment is set
	if(_curAlignment->isEmpty()){
	    _readAndParseAlignment(*_curAlignment);
	}

	while(*(_curAlignment) < window.to()){
        //check if alignment contains part of the window
		//if read continues outside of window, this is dealt with by window object
        if(_curAlignment->lastAlignedPositionWithRespectToRef() >= window.from()){
            _curAlignment = window.swapUsedForEmptyAlignment(_curAlignment);
		}

		//read next alignment
        if (!_readAndParseAlignment(*_curAlignment))
            break;
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
	//apply site-specific filters
	if(window.numReadsInWindow() > 0){
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
		if(_downsampleDepth > 0){
			window.downsample(_downsampleDepth, *subsamplePicker);
		};
	}

	//apply filters on window
	window.filter(_maxMissing, _maxRefN, _logfile);

	//report
	if(window.numReadsInWindow() > 0){
		window.dataSummary(_logfile);
	} else {
		_logfile->conclude("No data in this window.");
	}
};

void TGenome_windows::_traverseBAMWindows(){
	_logfile->startIndent("Traversing BAM file in windows");

	//initializing
	_hasWindowIndent = false;
	_curChromosome = _chromosomes.cend(); //set chromosome to end to trigger restart.
	_curAlignment->clear();

	//iterate through windows
	while(_readDataInNextWindow(_window)){
		if(_window.passedFilters()){
			//do stuff in derived classes
			_handleWindow();

			//report end of window calculations
			_logfile->list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
		}
	}

	_logfile->endIndent();

	//report
	_bamFile.printEndWithSummary();
};

}; //end namespace
