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

void TGenome_basic::mergeReadGroups(TParameters & params){
	//get reference to declustter coe
	BAM::TReadGroups& readGroups = _bamFile.readGroups;

	//read read groups to be merged
	std::string filename = params.getParameterString("readGroups");
	_logfile->startIndent("Reading read groups to be merged from file '" + filename + "':");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "!";

	//create map oldId -> new Id. Fill with identity.
	uint16_t readGroupMap[readGroups.size()];
	for(size_t i=0; i<readGroups.size(); ++i){
		readGroupMap[i] = i;
	}

	//parse file and construct new read groups in new header object
	int lineNum = 0;
	std::vector<std::string> vec;
	std::set<std::string> readGroupsMerged;
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		if(!vec.empty()){
			if(vec.size() < 2) throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + filename + "'!";

			//create new read group
			uint16_t newId = readGroups.add(vec[0]).id;
			_logfile->startIndent("The following read groups will be merged into '" + vec[0] + "':");

			for(size_t i=1; i<vec.size(); ++i){
				//check for duplicates
				if(!readGroupsMerged.emplace(vec[i]).second){
					throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
				}

				uint16_t oldId = readGroups.getId(vec[i]);

				//set not to write to header
				readGroups.removeFromHeader(oldId);

				//update map
				readGroupMap[oldId] = newId;

				//report
				_logfile->list(vec[i]);
			}
			_logfile->endIndent();
		}
	}

	//report unaffected read groups
	std::vector<std::string> unaffectedReadGroups;
	for(size_t i=0; i<readGroups.size(); ++i){
		if(readGroupMap[i] ==  i){
			unaffectedReadGroups.emplace_back(readGroups.getName(i));
		}
	}

	if(unaffectedReadGroups.size() > 0){
		_logfile->startIndent("The following read groups will be kept as is:");
		for(auto& s : unaffectedReadGroups){
			_logfile->list(s);
		}
		_logfile->endIndent();
	}

	//open a bam file for writing
	filename = _outputName + "_mergedRG.bam";
	_logfile->list("Writing new BAM file with merged read groups to file '" + filename + "'.");
	BAM::TOutputBamFile outBam(filename, _bamFile);

	//now parse through bam file and write alignments
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignmentThatPassesFilters()){
		_bamFile.curSetNewReadGroup(readGroupMap[_bamFile.curReadGroupID()]);
		_bamFile.writeCurAlignment(outBam);

		//report
		_bamFile.printProgress();
	}
	_bamFile.printEndWithSummary();

	//close bam writer
	outBam.close(_logfile);
};

//---------------------------------------------------------------
// TGenome_filtered
// A base class without genotype likelihoods but BAM filters enabled
//---------------------------------------------------------------
TGenome_filtered::TGenome_filtered(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_basic(Params, Logfile, RandomGenerator){
	_bamFile.setFilters(Params, _logfile);
};

//---------------------------------------------------------------
// TGenome_recalibrated
// A base class with BAM filters and recalibration
//---------------------------------------------------------------
TGenome_recalibrated::TGenome_recalibrated(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_filtered(Params, Logfile, RandomGenerator){
	//initialize genotype likelihoods
	_genotypeLikelihoodCalculator.init(Params, &_bamFile.readGroups, _logfile);
};

void TGenome_recalibrated::_traverseBAMPassedQC(){
	//parse through bam file
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignmentThatPassesFilters(_alignment)){
		//parse
		_alignment.parse(_genoMap, _qualMap, _genotypeLikelihoodCalculator.getSequencingErrorModels());

		//handle alignment by derived classes
		_handleAlignment();

		//report
		_bamFile.printProgress();
	}

	//report
	_bamFile.printEndWithSummary();
};

	virtual void _handleAlignment();

void TGenome_recalibrated::printQualityDistribution(TParameters & params){


	//prepare counters
	TCountDistributionVector qualDist;


	//CHECK: what can qualityTable do better?
	/*
	//Assemble quality distribution
	int maxQinPrintQualityDistribution = alignmentParser.getMaxPhredInt() + 33;
	_logfile->list("Will assemble quality distribution up to a quality of " + toString(maxQinPrintQualityDistribution-33) + " (" + (char) maxQinPrintQualityDistribution + ").");

	//initialize tables: one overall, one per read group
	std::vector<TQualityTable> qualDist;
	for(int i=0; i<alignmentParser.numReadGroups(); ++i){
		qualDist.emplace_back(maxQinPrintQualityDistribution);
	}
	*/

	//tmp variables
	BAM::TAlignment alignment;


	//parse through bam file
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignmentThatPassesFilters(alignment)){
		//parse


		//add to quality table

		//update and write (only if alignment qualities could be calculated)
		alignment.addToQualityTable(qualDist[alignment.readGroupID], qualMap);

		//report
		_bamFile.printProgress();
	}

	//report
	_bamFile.printEndWithSummary();

	//print per read group table
	_logfile->startIndent("Writing distributions:");
	std::string outFileName;
	for(size_t i=0; i<alignmentParser.readGroups.size(); ++i){
		//open output file
		outFileName = _outputName + "_" + alignmentParser.readGroups.getName(i) + "_qualityDistribution.txt";
		_logfile->listFlush("Writing distribution for read group '" + alignmentParser.readGroups.getName(i) + "' to '" + outFileName + "' ...");
		qualDist[i].write(outFileName);
		_logfile->done();
	}

	//print overall table
	outFileName = _outputName + "_total_qualityDistribution.txt";
	_logfile->listFlush("Writing total distribution to '" + outFileName + " ...");
	TQualityTable allQualDist(maxQinPrintQualityDistribution);
	for(int i=0; i<alignmentParser.numReadGroups(); ++i)
		allQualDist.add(qualDist[i]);
	allQualDist.write(outFileName);
	_logfile->done();
	_logfile->endIndent();
}

void TGenome_recalibrated::printQualityTransformation(TParameters & params){
	//prepare alignment
	TAlignment alignment(bamFile.maxReadLength());
	alignmentParser.setParsingToTrue();
	int maxPhredInt = params.getParameterIntWithDefault("maxQ", 100);

	//create table to store counts
	TQualityTransformTables QTtables(alignmentParser.readGroups, maxPhredInt);

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, _logfile);

	//check what we compare
	if(params.parameterExists("recal2")){
		std::string recalstring = params.getParameterString("recal2");
		_logfile->startIndent("Comparing to alternative recalibration parameters:");

		TSequencingErrorModels otherSeqErrorModels;
		otherSeqErrorModels.createModels(recalstring, &bamFile.readGroups, _logfile);

		while(alignmentParser.readNextAlignment(alignment)){
			alignment.addSitesToQualityTransformTable(otherSeqErrorModels, QTtables);

			//report
			reporter.printProgress();
		}
	} else {
		_logfile->startIndent("Comparing to original qualities:");
		while(alignmentParser.readNextAlignment(alignment)){
			alignment.addSitesToQualityTransformTable(QTtables);
			//report
			reporter.printProgress();
		}
	}

	_logfile->endIndent();

	//print tables
	QTtables.writeTables(_outputName, _logfile);

	//report end
	reporter.printEnd();
};

//-------------------------------------------------------
//TGenome
//-------------------------------------------------------
TGenomeWindows::TGenomeWindows(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_basic(Params, Logfile, RandomGenerator){


	//filters for bam file


	//initialize genotype likelihoods
	genotypeLikelihoodCalculator.init(params, &bamFile.readGroups, _logfile);

	//initialize alignment parser
	alignmentParser.init(params, _logfile);

	//open FASTA reference
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		_logfile->list("Reading reference sequence from '" + fastaFile + "'.");
		reference.initialize(fastaFile);
		alignmentParser.addReference(&reference);
	}

	//outputname
	alignmentParser.setOutName(_outputName);
};

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
	if(params.parameterExists("thetaGenomeWide") || alignmentParser.considerRegions){
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
				out.write(window.chrName, window.start, window.end);
			}
			_logfile->endIndent();

			//finish
			_logfile->list("Total computation time for this window was ", timer.seconds(), "s");
		} else _logfile->list("No relevant positions -> skipping this window.");
	}
};

void TGenomeWindows::estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool onlyReadData, int numBootstraps){
	if(alignmentParser.considerRegions)
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

	if(alignmentParser.considerRegions)
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
			filename = _outputName + window.chrName + "_" + toString(window.start) + "_LLsurface.txt";
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
		region1.setChr(window.chrName);
		region2.setChr(window.chrName);
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
			timer.start();

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
			thetaOut.write(window_full.chrName, window_full.start, window_full.end);

			//finish
			_logfile->list("Total computation time for this window was ", timer.seconds(), "s");
			_logfile->endIndent();
		} else _logfile->list("No relevant positions -> skipping this window.");
	}

	//clean up
	thetaOut.close();
};

//------------------------------------------
//Callers
//------------------------------------------
GenotypeLikelihoods::TGenotypePrior* TGenomeWindows::initializeGenotypePrior(TParameters & params){
	TGenotypePrior* prior;
	_logfile->startIndent("Initializing genotype prior:");
	//read prior from parameters
	std::string priorMethod = params.getParameterStringWithDefault("prior", "theta");
	if(priorMethod == "unif"){
		prior = new TGenotypePriorUniform();
		_logfile->list("Will use a uniform prior with equal weights for all genotypes.");
	} else if(priorMethod == "theta"){
		if(params.parameterExists("fixedTheta")){
			double theta = params.getParameterDouble("fixedTheta");
			_logfile->list("Will use a fixed theta = " + toString(theta));
			bool equalBaseFreq = params.parameterExists("equalBaseFreq");
			if(equalBaseFreq)
				_logfile->list("Will use equal base frequencies.");
			else
				_logfile->list("Will estimate base frequencies individually for each window.");
			prior = new TGenotypePriorFixedTheta(theta, equalBaseFreq, _logfile, _randomGenerator);
		} else {
			_logfile->list("Will use a prior based on theta and base frequencies estimated individually for each window.");
			std::string thetaOuputName = _outputName + "_theta_estimates.txt.gz";
			if(params.parameterExists("defaultTheta")){
				double defaultTheta = params.getParameterDouble("defaultTheta");
				_logfile->list("Will use a default theta of ", defaultTheta, " for windows with limited data.");
				prior = new TGenotypePriorTheta(params, thetaOuputName, defaultTheta, _logfile, _randomGenerator);
			} else
				prior = new TGenotypePriorTheta(params, thetaOuputName, _logfile, _randomGenerator);
		}
	} else throw "Unknown prior type '" + priorMethod + "'!";
	_logfile->endIndent();

	return prior;
}

void TGenomeWindows::callGenotypes(TParameters & params){
	//make sure FASTA is open
	if(!reference.hasReference()) throw "A FASTA reference must be provided to call!";

	//--------------------------
	//initialize caller
	//--------------------------
	_logfile->startIndent("Initializing caller:");
	TCaller* caller;
	std::string method = params.getParameterStringWithDefault("method", "MLE");
	if(method == "randomBase"){
		caller = new TCallerRandomBase(_randomGenerator);
	} else if(method == "majorityBase"){
		caller = new TCallerMajorityBase(_randomGenerator);
	} else if(method == "allelePresence"){
		caller = new TCallerAllelePresence(_randomGenerator);
	} else if(method == "MLE"){
		caller = new TCallerMLE(_randomGenerator);
	} else if(method == "Bayesian"){
		caller = new TCallerBayes(_randomGenerator);
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
	caller->reportSettings(_logfile);

	//prior setting
	TGenotypePrior* prior;
	if(caller->usesPrior()){
		prior = initializeGenotypePrior(params);
	} else {
		prior = new TGenotypePrior();
	}
	caller->setPrior(prior->getPointerToPrior());

	//open output file
	std::string sampleName = params.getParameterStringWithDefault("indName", _outputName);
	caller->openVCF(_outputName, sampleName);
	_logfile->endIndent();

	//prepare windows
	//Allow for haploid windows for some callers?
	TWindow window;

	//---------------------------------------------------------------------
	// Now call, either all sites or limiting to sites with known alleles.
	//---------------------------------------------------------------------
	if(params.parameterExists("alleles")){
		//Limit to sites with known alleles
		_logfile->startIndent("Will limit calls to sites with known alleles:");
		unsigned int windowSize = alignmentParser.getWindowSize();
		TSiteSubset subset(params.getParameterString("alleles"), windowSize, _logfile, false, reference, bamFile.chromosomes);
		_logfile->endIndent();

		while(alignmentParser.readDataInNextWindow(window)){
			subset.setChr(window.chrName);
			//read data for current window
			if(window.passedFilters || caller->printSitesWithNoData()){
				//update genotype prior
				prior->update(&window, _logfile);

				//now call using known alleles
				_logfile->listFlushTime("Calling genotypes ...");
				window.callKnwonAlleles(*caller, genotypeLikelihoodCalculator, subset);
				_logfile->doneTime();
			}
		}
	} else {
		//not limiting to sites with known alleles
		//Use all sites and identify alleles
		while(alignmentParser.readDataInNextWindow(window)){
			//read data for current window
			if(window.passedFilters || caller->printSitesWithNoData()){
				//update genotype prior
				prior->update(&window, _logfile);

				//now call
				_logfile->listFlushTime("Calling genotypes ...");
				window.call(*caller, genotypeLikelihoodCalculator, reference);
				_logfile->doneTime();
			}
		}
	}

	//clean up
	delete caller;
	delete prior;
};

//---------------------------------------------------
// I/O
//---------------------------------------------------
void TGenomeWindows::printPileup(TParameters & params){
	//initialize recalibration
//	initializeRecalibration(params);

	//prepare windows
	TWindow window;

	bool printOnlySitesWithData = false;
	if(params.parameterExists("printOnlySitesWithData")){
		printOnlySitesWithData = true;
		_logfile->list("Will print only sites with data. (parameter 'printOnlySitesWithData')");
	}

	//open output
	std::string filename = _outputName + "_pileup.txt.gz";
	_logfile->list("Writing pileup to file '" + filename + "'.");
	TOutputFileZipped out(filename);

	//write header
	TGenotypeMap genoMap;
	std::vector<std::string> header = {"Chr", "position", "ref", "depth", "refDepth", "bases"};
	for(int i=0; i<10; ++i)
		header.push_back("P(D|" + genoMap.getGenotypeString(i) + ")");
	out.writeHeader(header);

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		_logfile->listFlushTime("Writing pileup ...");
		window.printPileup(genotypeLikelihoodCalculator, out, printOnlySitesWithData);
		_logfile->doneTime();
	}

	//clean up
	out.close();
};

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
		if(alignmentParser.chrChangedWindow){
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
			window.createDepthMask(minDepthForMask, maxDepthForMask, output, window.chrName);
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
			subset.setChr(window.chrName);

			//read data for current window
			if(window.passedFilters && subset.hasPositionsInWindow(window.start))
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

void TGenomeWindows::recalibrateBamFile(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(bamFile.maxReadLength());
	alignmentParser.setParsingToTrue();

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = _outputName + "_recalibrated.bam";
	BamTools::RefVector references = bamFile.bamReader.GetReferenceData();
	_logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//do we also account for PMD?
	bool withPMD = params.parameterExists("withPMD");
	if(!withPMD && alignmentParser.hasPMD) _logfile->list("Note: PMD will not be reflected in the quality scores (preferred option when using ATLAS). If you want the quality scores to reflect pmd, use \"withPMD\"!");
	else if(withPMD && alignmentParser.hasPMD) _logfile->list("Probability of PMD will be reflected in new quality scores");
	else if(withPMD && !alignmentParser.hasPMD) throw "Probability of PMD is unknown. Provide PMD patterns or remove \"withPMD\"";
	if(withPMD && !alignmentParser.hasReference) throw "Cannot run recalBAM withPMD without reference!";

//	//should we include reads that don't pass filter?
//	bool allReads = false;
//	if(params.parameterExists("allReads")) allReads = true;

	//other tmp variables
	long counter = 0;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, _logfile);

    //now parse through bam file and write alignments
	if(withPMD){
		while(alignmentParser.readNextAlignment(alignment)){
			++counter;
			alignment.recalibrateWithPMD();
			alignment.save(bamWriter, alignmentParser.genoMap, alignmentParser.qualMap);
			reporter.printProgress();
        }
	} else {
		while(alignmentParser.readNextAlignment(alignment)){
			++counter;
			alignment.save(bamWriter, alignmentParser.genoMap, alignmentParser.qualMap);
			reporter.printProgress();
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	reporter.printEnd();

	//create index of new bam file
	indexBamFile(filename);
}

void TGenomeWindows::binQualityScores(TParameters & params){
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = _outputName + "_binnedQualityScores.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	_logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//initialize alignment reading
	TAlignment alignment(bamFile.maxReadLength());
	alignmentParser.setParsingToTrue();

	//other temp variables
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	long counter = 0;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, _logfile);

    //now parse through bam file and write alignments
	while(alignmentParser.readNextAlignment(alignment)){
		++counter;

		//update and write (only if alignment qualities could be calculated)
		alignment.binQualityScores(qualMap);
		alignment.save(bamWriter, genoMap, qualMap);

		//report
		reporter.printProgress();
	}

	//close bam writer
	bamWriter.Close();

	//report
	reporter.printEnd();
}


void TGenomeWindows::allelicDepth(TParameters & params){
	//allocate table
	//std::cout << "maxDepth " << alignmentParser.getMaxDepth() << std::endl;
	params.getParameterInt("maxDepth");
	if(alignmentParser.getMaxDepth() > 100){
		_logfile->warning("Allocating count table for a max depth of " + toString(alignmentParser.getMaxDepth()) + " uses a lot of memory! Use argument maxDepth to limit.");
	}
	TAllelicDepthCounts counts(alignmentParser.getMaxDepth());

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			_logfile->listFlush("Adding imbalance values to table ...");
			window.countAlleles(counts);
			_logfile->write(" done!");
		}
	}

	//write to file
	std::string outputFileName = _outputName + "_allelicDepth.txt.gz";
	_logfile->list("Writing allelic imbalance table to '" + outputFileName + "'");
	bool writeEmpty = params.parameterExists("printAll");
	counts.write(outputFileName, writeEmpty);
};

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
			if(window.depth == -1.0) output << alignmentParser.getCurChrName() << "\t" << window.start << "\t" << window.end << "\t" << "0" << "\n";
			else output << alignmentParser.getCurChrName() << "\t" << window.start << "\t" << window.end << "\t" << window.depth << "\n";
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

void TGenomeWindows::estimateDuplicationCounts(TParameters & params){
	//assembles distribution of how often a read is duplicated
	//now: just how many reads start at the same positions

	//initialize alignment reading
	TAlignment alignment(bamFile.maxReadLength());

	//create storage
	int maxCounts = params.getParameterIntWithDefault("maxCount", 20);
	TDistributionOfCounts counts(maxCounts, "readStarts");

	//iterate through windows
	int curChrLength = 0;
	unsigned int curPos = 0;
	int countsAtPos = 0;

	while (alignmentParser.readNextAlignment(alignment)){
		if(alignmentParser.chrChangedAlignment){
			//add last pos with data
			if(curChrLength > 0){
				counts.add(countsAtPos);
				countsAtPos = 0;

				//add all positions until chromosome end to structure
				counts.add(0, curChrLength - curPos);
			}

			curChrLength = alignmentParser.getCurChrLength();
			curPos = 0;
		}

		if(alignment.getPosition() > curPos){
			//add last pos with data
			counts.add(countsAtPos);

			//add zero for all positions until here
			counts.add(0, alignment.getPosition() - curPos);

			//set counts at current position
			curPos = alignment.getPosition();
			countsAtPos = 1;
		} else if(alignment.getPosition() == curPos){
			++countsAtPos;
		} else
			throw "Bam file is not sorted!";
	}

	//write output
	std::string filename = _outputName + "_readStartsPerSite.txt";
	_logfile->listFlush("Writing distribution of read starts per site to '" + filename + "' ...");
	counts.writeCounts(filename);
	_logfile->done();
};

//---------------------------------------------------
//PMD
//---------------------------------------------------
void TGenomeWindows::estimatePMD(TParameters & params){
	//make sure FASTA is open
	if(!alignmentParser.hasReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//prepare maps
	TReadGroupMap readGroupMap(&alignmentParser.readGroups, params.getParameterString("poolReadGroups", false), _logfile);
	TGenotypeMap genoMap;

	//prepare PMD table
	int maxLengthForInference = params.getParameterIntWithDefault("length", 50);
	_logfile->list("Estimating PMD at the first " + toString(maxLengthForInference) + " positions.");
	GenotypeLikelihoods::TPMDTables pmdTables(alignmentParser.readGroups, maxLengthForInference, maxReadLength, readGroupMap);

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, _logfile);

	//iterate through BAM file
	while(alignmentParser.readNextAlignment(alignment)){
		//alignment is only filled if filters are passed
		alignment.addToPMDTables(pmdTables, genoMap);

		//report
		reporter.printProgress();
	}
	//report
	reporter.printEnd();

	//print tables and data
	std::string filename = _outputName + "_PMD_Table.txt";
	_logfile->listFlush("Writing PMD table to '" + filename + "' ...");
	pmdTables.writeTable(filename);
	_logfile->done();
	filename = _outputName + "_PMD_Table_counts.txt";
	_logfile->listFlush("Writing PMD table of counts to '" + filename + "' ...");
	pmdTables.writeTableWithCounts(filename);
	_logfile->done();
	filename = _outputName + "_PMD_input_Empiric.txt";
	_logfile->listFlush("Writing PMD input file to '" + filename + "' ...");
	pmdTables.writePMDFile(filename);
	_logfile->done();

	//estimate exponential model
	if(!params.parameterExists("onlyEmpiric")){
		filename = _outputName + "_PMD_input_Exponential.txt";
		_logfile->listFlush("Estimating PMD exponential models and writing them to '" + filename + "' ...");
		int numNRIterations = params.getParameterIntWithDefault("numNRIterations", 100);
		double eps = params.getParameterDoubleWithDefault("eps", 0.001);
		pmdTables.fitExponentialModel(numNRIterations, eps, filename, _logfile);
		_logfile->done();
	} else {
		_logfile->list("Not fitting exponential model due to user specification (parameter 'onlyEmpiric')");
	}
};

void TGenomeWindows::runPMDS(TParameters & params){
	//parse bam file and calculate PMDS for each read (seeSkoglund et al. 2014)
	//write new bam file with PMDS score added
	//parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	if(!alignmentParser.hasReference) throw "Cannot run PMDS without reference!";

	//get parameters
	double pi = params.getParameterDoubleWithDefault("pi", 0.001);
	_logfile->list("Running PMDS with rate of polymorphism (pi) = " + toString(pi));
	double minPMDS = params.getParameterDoubleWithDefault("minPMDS", -10000);
	double maxPMDS = params.getParameterDoubleWithDefault("maxPMDS", 10000);
	_logfile->list("Filtering out reads with " + toString(minPMDS) + " > PMDS > " + toString(maxPMDS));

	//other tmp
	TQualityMap qualMap;
	TGenotypeMap genoMap;

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = _outputName + "_PMDS.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	_logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, _logfile);
	long numKept = 0;

	//now parse through bam file and write alignments
	while(alignmentParser.readNextAlignment(alignment)){
		//calc PMD
		double PMDS = alignment.calculatePMDS(pi, alignmentParser.pmdObjects);

		//update and write
		if(PMDS > minPMDS && PMDS < maxPMDS){
			alignment.updateOptionalSamField("DS", PMDS);
			alignment.save(bamWriter, alignmentParser.genoMap, alignmentParser.qualMap);
			++numKept;
		}

		//report progress
		reporter.printProgress();
	}

	//close bam writer
	bamWriter.Close();

	//report
	reporter.printEndNoEndIndent();
	_logfile->conclude("Kept " + toString(numKept) + " reads");
	_logfile->endIndent();
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

void TGenomeWindows::contextStats(TParameters & params){
	//prepare table
	TContextStats table(alignmentParser.maxQualityAsPhredInt);

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//now parse through bam file and write alignments
	while(alignmentParser.readNextAlignment(alignment)){


	}

	std::string outputFileName = _outputName + "_contextInformation.txt.gz";
		_logfile->list("Writing context information to file '" + outputFileName + "'.");


    }



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
			std::cout << window.chrName << "\t" << window.start + i + 1 << "\t" << window.sites[i].getBases() << std::endl;

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

