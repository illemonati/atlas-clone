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
TGenome::TGenome(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	randomGenerator = RandomGenerator;

	//initialize alignment parser
	maxReadLength = params.getParameterIntWithDefault("maxReadLength", 1000);
	logfile->list("Will only consider reads up to " + toString(maxReadLength) + " bp. (parameter 'maxReadLength')");

	alignmentParser.init(maxReadLength, params, logfile);

	//outputname
	outputName = params.getParameterStringWithDefault("out", "");
	if(outputName == ""){
		//guess from filename. note that the genome has outputName in case we want to have multiple parsers in the future
		outputName = alignmentParser.filename;
		outputName = extractBeforeLast(outputName, ".");
	}
	alignmentParser.setOutName(outputName);
	logfile->list("Writing output files with prefix '" + outputName + "'. (parameter 'out')");

	//open FASTA reference
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		std::string fastaIndex = fastaFile + ".fai";
		logfile->list("Reading reference sequence from '" + fastaFile + "'");
		if(!reference.Open(fastaFile, fastaIndex)) throw "Failed to open FASTA file '" + fastaFile + "'! Is index file present?";
		alignmentParser.addReference(&reference);
	}

	//trimming ends
	if(params.parameterExists("trim3") || params.parameterExists("trim5")){
		int trim3 = params.getParameterIntWithDefault("trim3", 0);
		if(trim3 < 0) throw "trimming distance trim3 must be >= 0!";
		int trim5 = params.getParameterIntWithDefault("trim5", 0);
		if(trim5 < 0) throw "trimming distance trim5 must be >= 0!";
		if(trim3>0 || trim5>0){
			alignmentParser.setReadTrimming(trim3, trim5);
			logfile->list("Will trim first " + toString(trim3) + " and " + toString(trim5) + " bases from the 3' and 5' end, respectively. (parameters 'trim3', 'trim5')");
		}
	}
};

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
};

//-----------------------------------------------------
//Functions for theta estimation
//-----------------------------------------------------
bool TGenome::estimateTheta(TThetaEstimator & thetaEstimator, TWindow_base & window){
	//adding sites to estimator
	logfile->listFlush("Calculating emission probabilities ...");
	thetaEstimator.clear();
	window.addSitesToThetaEstimator(thetaEstimator.pointerToDataContainer());
	logfile->done();

	//estimate Theta
	return thetaEstimator.estimateTheta();
};

void TGenome::estimateTheta(TParameters & params){
	//Theta estimator
	TThetaEstimator thetaEstimator(params, logfile, randomGenerator);

	//open output
	std::string filename = outputName + "_theta_estimates.txt.gz";
	TThetaOutputFile thetaOut(&thetaEstimator, filename, logfile);

	//check for which segements theta is to be estimated
	if(params.parameterExists("thetaGenomeWide") || alignmentParser.considerRegions){
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
	} else {
		estimateThetaWindows(thetaEstimator, thetaOut, params.parameterExists("printAll"));
	}

	//clean up
	thetaOut.close();
};

void TGenome::estimateThetaWindows(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool printAll){
	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters || printAll){
			//measure runtime
			struct timeval startTime, endTime;
			gettimeofday(&startTime, NULL);

			logfile->startIndent("Estimating Theta:");
			if(estimateTheta(thetaEstimator, window) || printAll){
				out.write(alignmentParser.getCurChrName(), window.start, window.end);
			}
			logfile->endIndent();

			//finish
			gettimeofday(&endTime, NULL);
			logfile->list("Total computation time for this window was ", endTime.tv_sec  - startTime.tv_sec, "s");
		} else logfile->list("No relevant positions -> skipping this window.");
	}
};

void TGenome::estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool onlyReadData, int numBootstraps){
	if(alignmentParser.considerRegions)
		logfile->startIndent("Estimating theta at specific sites:");

	//prepare windows
	TWindow window;
	thetaEstimator.clear();

	//add sites to estimator
	logfile->startIndent("Adding sites to data structure:");
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters){
			//adding sites to estimator
			logfile->listFlush("Calculating emission probabilities ...");
			try{
				window.addSitesToThetaEstimator(thetaEstimator.pointerToDataContainer());
			} catch(...){
				throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size, selecting fewer regions or limiting to sites with a minimal depth (>=2 recommended).";
			}
			logfile->done();
		}
	}
	logfile->endIndent();

	//estimate Theta
	//HACK!!!!
	if(!onlyReadData){
		logfile->startIndent("Estimate theta based on a total of " + toString(thetaEstimator.sizeWithData()) + " sites:");
		thetaEstimator.estimateTheta();
	}

	if(alignmentParser.considerRegions)
		out.write("regions", "-", "-");
	else
		out.write("genome-wide", "-", "-");
	if(numBootstraps == 0)
		thetaEstimator.clear();
};

void TGenome::bootstrapTetaEstimation(int numBootstraps, TThetaEstimator & thetaEstimator){
	if(numBootstraps < 1) throw "Number of bootstraps must be > 1!";
	logfile->startIndent("Generating " + toString(numBootstraps) + " bootstrap estimates of theta:");

	//measure runtime
	struct timeval startTime, endTime;
	gettimeofday(&startTime, NULL);

	//open output file
	TOutputFileZipped bootstrapOut;
	std::string bootstrapFilename = outputName + "_theta_bootstraps.txt.gz";
	logfile->list("Writing theta bootstraps to '" + bootstrapFilename + "'");
	bootstrapOut.open(bootstrapFilename.c_str());

	//write header
	bootstrapOut.setPrecision(9);
	std::vector<std::string> header = {"Bootstrap"};
	thetaEstimator.addToHeader(header);
	bootstrapOut.writeHeader(header);

	//loop over bootstraps
	for(int s=0; s<numBootstraps; ++s){
		logfile->startIndent("Bootstrap " + toString(s+1) + " of " + toString(numBootstraps) + ":");

		//run bootstrap
		bootstrapOut << s+1;
		thetaEstimator.bootstrapTheta(&bootstrapOut);

		logfile->endIndent();
	}

	//finish
	gettimeofday(&endTime, NULL);
	logfile->list("Total computation time for theta bootstrapping was ", round((endTime.tv_sec  - startTime.tv_sec) / 6.0)/10.0, " min");
	logfile->endIndent();
};

void TGenome::calcThetaLikelihoodSurfaces(TParameters & params){
	//read params
	int steps = params.getParameterIntWithDefault("steps", 100);

	//prepare windows
	TWindow window;

	//Theta estimator
	TThetaEstimator estimator(logfile, randomGenerator);

	//iterate through windows
	std::string filename;
	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//read data for current window
		if(window.passedFilters){
			//check if we have data -> can be extended to ensure
			logfile->startIndent("Calculating likelihood surface for Theta:");

			//measure runtime
			struct timeval startTime, endTime;
			gettimeofday(&startTime, NULL);

			//adding sites to estimator
			logfile->listFlush("Calculating emission probabilities ...");
			window.addSitesToThetaEstimator(estimator.pointerToDataContainer());
			logfile->done();

			//open file
			gz::ogzstream out;
			filename = outputName + alignmentParser.getCurChrName() + "_" + toString(window.start) + "_LLsurface.txt";
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
};

void TGenome::estimateThetaRatio(TParameters & params){
	//Theta estimator
	TThetaEstimatorRatio thetaEstimatorRatio(params, logfile, randomGenerator);

	//read the two regions to be used
	logfile->startIndent("Reading regions:");

	int windowSize = alignmentParser.getWindowSize();

	//region 1
	std::string regionsFile1 = params.getParameterString("regions1");
	int siteLimit = -1;
	if(params.parameterExists("siteLimit1")){
		siteLimit = params.getParameterInt("siteLimit1");
		if(siteLimit < 0)
			throw "site limit cannot be smaller than 0!";
		logfile->startIndent("Reading first " + toString(siteLimit) + " sites from regions 1 BED file '" + regionsFile1 + "':");
	} else {
		logfile->listFlush("Reading regions 1 from BED file '" + regionsFile1 + "' ...");
	}
	TBedReader region1(regionsFile1, windowSize, alignmentParser.bamHeader.Sequences, siteLimit,logfile);
	logfile->done();

	//region 2
	siteLimit = -1;
	std::string regionsFile2 = params.getParameterString("regions2");
	if(params.parameterExists("siteLimit2")){
		siteLimit = params.getParameterInt("siteLimit2");
		if(siteLimit < 0)
			throw "site limit cannot be smaller than 0!";
		logfile->startIndent("Reading first " + toString(siteLimit) + " sites from regions 1 BED file '" + regionsFile1 + "':");
	} else {
		logfile->listFlush("Reading regions 2 from BED file '" + regionsFile2 + "' ...");
	}
	TBedReader region2(regionsFile2, windowSize, alignmentParser.bamHeader.Sequences, siteLimit, logfile);
	logfile->done();
	logfile->endIndent();

	//prepare windows
	TWindow window;

	//add sites to estimator
	logfile->startIndent("Adding sites to data structures:");
	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		region1.setChr(alignmentParser.getCurChrName());
		region2.setChr(alignmentParser.getCurChrName());
		if(window.passedFilters){
			//adding sites to estimator
			logfile->listFlush("Calculating emission probabilities for sites within defined regions ...");
			try{
				window.addSitesToThetaEstimator(thetaEstimatorRatio.pointerToDataContainer(), region1);
				window.addSitesToThetaEstimator(thetaEstimatorRatio.pointerToDataContainer2(), region2);
			} catch(...){
				throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size, selecting fewer regions or limiting to sites with a minimal depth (>=2 recommended).";
			}
			logfile->done();
		}
	}
	logfile->endIndent();

	//estimate Theta ratio
	thetaEstimatorRatio.estimateRatio(outputName);
};

void TGenome::performDownsamplingThetaQC(TParameters & params){
	//parse downsampling parameters
	std::vector<double> downSampleProbVector;
	fillVectorOfDownsamplingProbabilities(params.getParameterStringWithDefault("prob", "1.0,0.5,0.2,0.1,0.05,0.02,0.01"), downSampleProbVector);

	//report probabilities
	logfile->list("Will estimate theta after downsampling reads with probabilities (parameter 'prob'): " + concatenateString(downSampleProbVector, ", "));

	//check if full data is to be used (i.e. if prob = 1.0 is specified)
	bool printFull = false;
	if(*downSampleProbVector.begin() == 1.0){
		printFull = true;
		downSampleProbVector.erase(downSampleProbVector.begin());
	}

	//create window and estimator for full data
	TThetaEstimator estimator_full(params, logfile, randomGenerator);
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
	std::string filename = outputName + "_theta_estimates.txt.gz";
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
	thetaOut.open(filename, logfile);

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window_full)){
		if(window_full.passedFilters){
			//measure runtime
			struct timeval startTime, endTime;
			gettimeofday(&startTime, NULL);

			//estimate on full data
			if(printFull){
				logfile->startIndent("Estimating Theta on full data:");
				estimateTheta(estimator_full, window_full);
				logfile->endIndent();
			}

			if(windows.size() > 0){
				for(size_t i = 0; i<downSampleProbVector.size(); ++i){
					logfile->startIndent("Estimating Theta on downsampled data (p = " + toString(downSampleProbVector[i]) + "):");

					//downsample
					alignmentParser.downsampleWindow(*windows[i], window_full, downSampleProbVector[i], randomGenerator);

					//and estimate
					estimateTheta(estimators[i], *windows[i]);
					logfile->endIndent();
				}
			}

			//write output
			thetaOut.write(alignmentParser.getCurChrName(), window_full.start, window_full.end);

			//finish
			gettimeofday(&endTime, NULL);
			logfile->list("Total computation time for this window was ", endTime.tv_sec  - startTime.tv_sec, "s");
			logfile->endIndent();
		} else logfile->list("No relevant positions -> skipping this window.");
	}

	//clean up
	thetaOut.close();
};

//------------------------------------------
//Callers
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
			prior = new TGenotypePriorFixedTheta(theta, equalBaseFreq, logfile, randomGenerator);
		} else {
			logfile->list("Will use a prior based on theta and base frequencies estimated individually for each window.");
			std::string thetaOuputName = outputName + "_theta_estimates.txt.gz";
			if(params.parameterExists("defaultTheta")){
				double defaultTheta = params.getParameterDouble("defaultTheta");
				logfile->list("Will use a default theta of ", defaultTheta, " for windows with limited data.");
				prior = new TGenotypePriorTheta(params, thetaOuputName, defaultTheta, logfile, randomGenerator);
			} else
				prior = new TGenotypePriorTheta(params, thetaOuputName, logfile, randomGenerator);
		}
	} else throw "Unknown prior type '" + priorMethod + "'!";
	logfile->endIndent();

	return prior;
}

void TGenome::callGenotypes(TParameters & params){
	//make sure FASTA is open
	if(!alignmentParser.fastaReference) throw "A FASTA reference must be provided to call!";

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
	TWindow window;

	//---------------------------------------------------------------------
	// Now call, either all sites or limiting to sites with known alleles.
	//---------------------------------------------------------------------
	if(params.parameterExists("alleles")){
		//Limit to sites with known alleles
		logfile->startIndent("Will limit calls to sites with known alleles:");
		int windowSize = alignmentParser.getWindowSize();
		TSiteSubset subset(params.getParameterString("alleles"), reference, alignmentParser.bamHeader, windowSize, logfile, false);
		logfile->endIndent();

		while(alignmentParser.readDataInNextWindow(window)){
			subset.setChr(alignmentParser.getCurChrName());
			//read data for current window
			if(window.passedFilters || caller->printSitesWithNoData()){
				//update genotype prior
				prior->update(&window, alignmentParser.getCurChrName(), logfile);

				//now call using known alleles
				logfile->listFlush("Calling genotypes ...");
				window.callKnwonAlleles(*caller, *alignmentParser.recalObject, subset);
				logfile->done();
			}
		}
	} else {
		//not limiting to sites with known alleles
		//Use all sites and identify alleles
		while(alignmentParser.readDataInNextWindow(window)){
			//read data for current window
			if(window.passedFilters || caller->printSitesWithNoData()){
				//update genotype prior
				prior->update(&window, alignmentParser.getCurChrName(), logfile);

				//now call
				logfile->listFlush("Calling genotypes ...");
				window.call(*caller, *alignmentParser.recalObject, reference);
				logfile->done();
			}
		}
	}

	//clean up
	delete caller;
	delete prior;
}

//---------------------------------------------------
// I/O
//---------------------------------------------------
void TGenome::writeGLF(TParameters & params){
	//print all?
	bool printIfNoData = params.parameterExists("printAll");
	if(printIfNoData)
		logfile->list("Will print all sites, even those without data. (parameter 'printAll')");

	//open GLF file
	std::string outputFileName = outputName + ".glf.gz";
	logfile->list("Will write data to GLF file '" + outputFileName + "'");
	TGlfWriter writer(outputFileName);

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		if(alignmentParser.chrChangedWindow){
			writer.newChromosome(alignmentParser.getCurChrName(), (uint32_t) alignmentParser.getCurChrLength(), (uint8_t) alignmentParser.getCurChrPloidy());
		} if(window.passedFilters){
			//write to GLF
			logfile->listFlush("Adding window to GLF file ...");
			window.calculateEmissionProbabilities();
			window.addToGLF(writer, alignmentParser.getCurChrPloidy(), printIfNoData);
			logfile->done();
		}
	}
	//clean up
	writer.close();
};

void TGenome::printPileup(TParameters & params){
	//initialize recalibration
//	initializeRecalibration(params);

	//prepare windows
	TWindow window;

	bool printOnlySitesWithData = false;
	if(params.parameterExists("printOnlySitesWithData")){
		printOnlySitesWithData = true;
		logfile->list("Will print only sites with data. (parameter 'printOnlySitesWithData')");
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
	while(alignmentParser.readDataInNextWindow(window)){
		window.printPileup(alignmentParser.recalObject, out, printOnlySitesWithData);
	}

	//clean up
	out.close();
};

void TGenome::generatePSMCInput(TParameters & params){
	//read in parameters required
	double theta = params.getParameterDoubleWithDefault("theta", 0.001);
	logfile->list("Using theta = " + toString(theta));
	TThetaEstimator thetaEstimator(logfile, randomGenerator);
	thetaEstimator.setTheta(theta);

	double confidence = params.getParameterDoubleWithDefault("confidence", 0.99);
	logfile->list("Calling heterozygosity state with confidence > " + toString(confidence) + ". (parameter 'confidence')");
	int blockSize = params.getParameterIntWithDefault("block", 100);
	//make sure window size is a multiple of block length!
	if(alignmentParser.getWindowSize() % blockSize != 0) throw "Window size is not a multiple of block size!";

	//open output file
	std::ofstream output;
	std::string outputFileName = outputName + ".psmcfa";
	logfile->list("Writing PSMC input file to '" + outputFileName + "'");
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
				output << '>' << alignmentParser.getCurChrName() << '\n';
		}
		if(window.passedFilters){
			//set base frequencies
			logfile->listFlush("Calculating emission probabilities ...");
			window.calculateEmissionProbabilities();
			window.estimateBaseFrequencies();
			TBaseFrequencies baseFreq =  window.getBaseFreq();
			thetaEstimator.setBaseFreq(baseFreq);
			logfile->done();

			//create PSMC input
			logfile->listFlush("Estimating heterozygosity status ...");
			window.generatePSMCInput(thetaEstimator, blockSize, confidence, output, nCharOnLine);
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
	if(params.parameterExists("maxDepth") || params.parameterExists("minDepth"))
		throw "Cannot mask sites for sequencing depth while creating the mask!";

	std::ofstream output;
	std::string outputFileName = outputName + "_minDepth"+ toString(minDepthForMask) + "_maxDepth" + toString(maxDepthForMask) + "_depthMask.bed";
	logfile->list("Writing sites with depth < " + toString(minDepthForMask) + " and with depth > " + toString(maxDepthForMask) + " to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	//prepare windows
	TWindow window;

	//iterate through windows
	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//read data for current window
		if(window.passedFilters){
			logfile->listFlush("Writing sites to mask to output file ...");
			window.createDepthMask(minDepthForMask, maxDepthForMask, output, alignmentParser.getCurChrName());
			logfile->done();
		}
	}
	output.close();
}

//---------------------------------------------------
// Recalibration
//---------------------------------------------------
void TGenome::estimateErrorCalibrationEM(TParameters & params){
	if(alignmentParser.qualitiesScoresAreRecalibrated() && !params.parameterExists("rerecalibrate"))
		throw "Can not estimate recalibration: quality scores are already recalibrated while reading! (Use argument 'rerecalibrate' to overwrite this error)";

	//initialize maps
	TReadGroupMap readGroupMap(alignmentParser.readGroups, params.getParameterString("poolReadGroups", false), logfile);
	TQualityMap qualityMap;

	//create recalibration object
	bool writeTmpTables = false;
	if(params.parameterExists("writeTmpTables")){
		writeTmpTables = true;
		logfile->list("Will write intermediate estimates of EM and Newton-Raphson to file.");
	}
	TRecalibrationEMEstimator recalObjectEM(params, &alignmentParser.readGroups, logfile, &readGroupMap);

	//prepare windows
	TWindow window;

	//add sites to EM object
	if(params.parameterExists("alleles")){
		//Limit to sites with known alleles
		logfile->startIndent("Will limit analysis to homozygous sites with known genotypes:");
		int windowSize = alignmentParser.getWindowSize();
		TSiteSubset subset(params.getParameterString("alleles"), windowSize, logfile, true);

		//now parse through windows
		while(alignmentParser.readDataInNextWindow(window)){
			//make sure subset is on the correct chromosome
			subset.setChr(alignmentParser.getCurChrName());

			//read data for current window
			if(window.passedFilters)
				window.addToRecalibrationEM(recalObjectEM, &subset, qualityMap);
			else logfile->list("No positions in this window.");
		}
		logfile->endIndent();

		//clean up memory
		window.clear();

		//now estimate
		recalObjectEM.performEstimationKnownGenotypes(outputName, writeTmpTables);

	} else {
		logfile->startIndent("Reading data from windows:");
		while(alignmentParser.readDataInNextWindow(window)){
			//read data for current window
			if(window.passedFilters)
				window.addToRecalibrationEM(recalObjectEM, qualityMap);
			else logfile->list("No positions in this window.");
		}
		logfile->endIndent();
		logfile->endIndent();

		//clean up memory
		window.clear();

		//now estimate
		recalObjectEM.performEstimation(outputName, writeTmpTables);
	}
};

//TODO: remove? Does not currently work.
void TGenome::calculateLikelihoodErrorCalibrationEM(TParameters & params){
	//create recalibration object
	TReadGroupMap readGroupMap(alignmentParser.readGroups, params.getParameterString("poolReadGroups", false), logfile);
	TRecalibrationEMEstimator recalObjectEM(params, &alignmentParser.readGroups, logfile, &readGroupMap);
	recalObjectEM.initializeFromString(params.getParameterString("recalForLL"));

	//prepare windows
	TWindow window;

	//other tmp variables
	TQualityMap qualMap;

	//add sites to EM object
	logfile->startIndent("Reading data from windows:");
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters)
			window.addToRecalibrationEM(recalObjectEM, qualMap);
	}
	//clean up memory
	window.clear();
	logfile->endIndent();

	//calc likelihood surface
	//int numMarginalGridPoint = params.getParameterIntWithDefault("numGridPoints", 51);
	//recalObjectEM.calcLikelihoodSurface(outputName + "_LLsurface.txt", numMarginalGridPoint);

	logfile->list("LL = " + toString(recalObjectEM.calcLL()));
}

void TGenome::BQSR(TParameters & params){
	if(alignmentParser.qualitiesScoresAreRecalibrated())
		throw "Can not estimate recalibration: quality scores are already recalibrated while reading!";

	//make sure FASTA is open
	if(!alignmentParser.hasReference) throw "Can not run BQSR recalibration without a provided FASTA reference!";

	//prepare windows
	TWindow window;

	//create BQSR object
	TReadGroupMap readGroupMap(alignmentParser.readGroups, params.getParameterString("poolReadGroups", false), logfile);
	TRecalibrationBQSREstimator bqsr(params, logfile, &alignmentParser.readGroups, &readGroupMap);

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

	//loop over bam until BQSR converges
	while(!hasConverged){
		++loopNumber;
		logfile->startIndent("Running recalibration loop " + toString(loopNumber) + ":");
		//loop over all windows

		if(!bqsr.dataHasBeenStored()){
			logfile->startIndent("Reading data from BAM file:");
			while(alignmentParser.readDataInNextWindow(window)){
				if(window.passedFilters){
					//add the base to BQSR
					if(invariantSites) window.addSitesToBQSR(bqsr, alignmentParser.subset, logfile, alignmentParser.qualMap);
					else window.addSitesToBQSR(bqsr, logfile, alignmentParser.qualMap);
				} else logfile->list("No positions in this window.");
			}
			//clean up memory
			window.clear();
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
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//Assemble quality distribution
	int maxQinPrintQualityDistribution = alignmentParser.getMaxPhredInt() + 33;
	logfile->list("Will assemble quality distribution up to a quality of " + toString(maxQinPrintQualityDistribution-33) + " (" + (char) maxQinPrintQualityDistribution + ").");

	//initialize tables: one overall, one per read group
	std::vector<TQualityTable> qualDist;
	for(int i=0; i<alignmentParser.numReadGroups(); ++i){
		qualDist.emplace_back(maxQinPrintQualityDistribution);
	}

	//other tmp variables
	TQualityMap qualMap;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	while(alignmentParser.readNextAlignment(alignment)){
		//update and write (only if alignment qualities could be calculated)
		alignment.addToQualityTable(qualDist[alignment.readGroupId], qualMap);

		//report
		reporter.printProgress();
	}

	//report
	reporter.printEnd();

	//print per read group table
	logfile->startIndent("Writing distributions:");
	std::string outFileName;
	for(size_t i=0; i<alignmentParser.readGroups.size(); ++i){
		//open output file
		outFileName = outputName + "_" + alignmentParser.readGroups.getName(i) + "_qualityDistribution.txt";
		logfile->listFlush("Writing distribution for read group '" + alignmentParser.readGroups.getName(i) + "' to '" + outFileName + "' ...");
		qualDist[i].write(outFileName);
		logfile->done();
	}

	//print overall table
	outFileName = outputName + "_total_qualityDistribution.txt";
	logfile->listFlush("Writing total distribution to '" + outFileName + " ...");
	TQualityTable allQualDist(maxQinPrintQualityDistribution);
	for(int i=0; i<alignmentParser.numReadGroups(); ++i)
		allQualDist.add(qualDist[i]);
	allQualDist.write(outFileName);
	logfile->done();
	logfile->endIndent();
}

void TGenome::printQualityTransformation(TParameters & params){
	//prepare alignment
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();
	int maxPhredInt = params.getParameterIntWithDefault("maxQ", 100);

	//create table to store counts
	TQualityTransformTables QTtables(alignmentParser.readGroups, maxPhredInt);

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

	//check what we compare
	bool compareToOtherRecalibration = false;
	TRecalibration* otherRecalObject;
	if(alignmentParser.recalibrationType() == "BQSR"){
		//do we compare to recal or to raw quality scores?
		if(params.parameterExists("recal")){
			otherRecalObject = new TRecalibrationEM(params.getParameterString("recal"), &alignmentParser.readGroups, logfile);
			compareToOtherRecalibration = true;
		}
	} else if(alignmentParser.recalibrationType() == "recal"){
		if(params.parameterExists("recal2")){
			otherRecalObject = new TRecalibrationEM(params.getParameterString("recal2"), &alignmentParser.readGroups, logfile);
			compareToOtherRecalibration = true;
		}
	}
	//add alignments to tables
//	logfile->startIndent();
	logfile->listFlush("Adding sites to quality transformation tables ...");
	logfile->done();

	if(compareToOtherRecalibration){
		while(alignmentParser.readNextAlignment(alignment)){
			alignmentParser.addSitesToQualityTransformTable(alignment, otherRecalObject, QTtables);
			//report
			reporter.printProgress();
		}
	} else {
		while(alignmentParser.readNextAlignment(alignment)){
			alignmentParser.addSitesToQualityTransformTable(alignment, QTtables);
			//report
			reporter.printProgress();
		}
	}

	logfile->endIndent();

	//print tables
	QTtables.writeTables(outputName, logfile);

	//report end
	reporter.printEnd();
};

void TGenome::reportProgressParsingBamFile(const long & counter, const struct timeval & start){
	if(counter % 1000000 == 0){
		reportProgressParsingBamFileNoCheck(counter, start);
	}
}
void TGenome::reportProgressParsingBamFileNoCheck(const long & counter, const struct timeval & start){
	static struct timeval end;
	gettimeofday(&end, NULL);
	float runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
}

//---------------------------------------------------
//BAM manipulation / statistics
//---------------------------------------------------

void TGenome::recalibrateBamFile(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_recalibrated.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//do we also account for PMD?
	bool withPMD = params.parameterExists("withPMD");
	if(!withPMD && alignmentParser.hasPMD) logfile->list("Note: PMD will not be reflected in the quality scores (preferred option when using ATLAS). If you want the quality scores to reflect pmd, use \"withPMD\"!");
	else if(withPMD && alignmentParser.hasPMD) logfile->list("Probability of PMD will be reflected in new quality scores");
	else if(withPMD && !alignmentParser.hasPMD) throw "Probability of PMD is unknown. Provide PMD patterns or remove \"withPMD\"";
	if(withPMD && !alignmentParser.hasReference) throw "Cannot run recalBAM withPMD without reference!";

//	//should we include reads that don't pass filter?
//	bool allReads = false;
//	if(params.parameterExists("allReads")) allReads = true;

	//other tmp variables
	long counter = 0;
	TGenotypeMap genoMap;
	TQualityMap qualMap;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	if(withPMD){
		while(alignmentParser.readNextAlignment(alignment)){
			++counter;
			alignment.recalibrateWithPMD(alignmentParser.recalObject, qualMap);
			alignment.save(bamWriter, genoMap, alignmentParser.minQual, alignmentParser.maxQual, qualMap);
			reporter.printProgress();
        }
	} else {
		while(alignmentParser.readNextAlignment(alignment)){
			++counter;
			alignment.save(bamWriter, genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, qualMap);
			reporter.printProgress();
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	reporter.printEnd();

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
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//other temp variables
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	long counter = 0;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	while(alignmentParser.readNextAlignment(alignment)){
		++counter;

		//update and write (only if alignment qualities could be calculated)
		alignment.binQualityScores(qualMap);
		alignment.save(bamWriter, genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, qualMap);

		//report
		reporter.printProgress();
	}

	//close bam writer
	bamWriter.Close();

	//report
	reporter.printEnd();
}


void TGenome::assessSoftClipping(TParameters & params){
	//build table ??

	//initialize alignment reading
	TAlignment alignment(maxReadLength);

	//limit input / output
	bool printAll = false;
	if(params.parameterExists("printAll")){
		printAll = true;
		logfile->list("Writing soft clipping stats for all alignments to file.");
	} else
		logfile->list("Writing only stats for reads with soft clipping to file.");

	bool printSequences = false;
	if(params.parameterExists("printSequences")){
		printSequences = true;
		logfile->list("Writing soft clipped bases to file.");
	} else
		logfile->list("Writing only counts of soft clipped bases to file.");

	//open output file
	TSoftClippingStatsFile statFile(outputName, printSequences);

	//create soft clipping matrix
	TSoftClippingMatrix matrix(maxReadLength);

	//other temp variables
	TSoftClippingData softClippingData;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

	//now parse through bam file and write alignments
	while(alignmentParser.readNextAlignment(alignment)){
		softClippingData.assessSoftClipping(alignment.bamAlignment);

		//report
		if(softClippingData.hasSoftClipping || printAll){
			statFile.write(alignment.name(), alignment.getPosition(), softClippingData);
		}

		//add count to matrix
		matrix.add(softClippingData);

		//report
		reporter.printProgress();
	}

	//write matrix
	matrix.write(outputName + "_clippingStats_matrix.txt");

	//report end
	reporter.printEnd();
};

void TGenome::removeSoftClippedBasesFromReads(TParameters & params){
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_softClippedBasesRemoved.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	TAlignment alignment;
	TSoftClippingData softClippingData;
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	while (alignmentParser.readNextAlignment(alignment)){
		softClippingData.assessSoftClipping(alignment.bamAlignment);
		alignment.removeSoftClippedBases(softClippingData);

		//write
		bamWriter.SaveAlignment(alignment.bamAlignment);
//		alignment.bamAlignment.save(bamWriter, alignmentParser.genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, alignmentParser.qualMap);

		//report
		reporter.printProgress();
	}

	//close bam writer
	bamWriter.Close();

	//now generate bam index
	indexBamFile(filename);

	//report
	reporter.printEnd();
}

void TGenome::assessOverlap(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//initialize table
	int* counts = new int[maxReadLength]();

	//other variables
	float numProperPairs = 0.0;

	//open output file
	std::string filename = outputName + "_overlapStats.txt";
	std::ofstream out(filename.c_str());
	if(!out)
		throw "Failed to open file '" + filename + "' for writing!";
	out << "overlap\tcount\tproportion\n";

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

	//now parse through bam file and write alignments
    while(alignmentParser.readNextAlignment(alignment)){
		int overlap = alignment.measureOverlap();
		if(overlap >= 0){
			++counts[overlap];
			++numProperPairs;
		}

		//report
		reporter.printProgress();
	}
    reporter.printEnd();

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
			readGroupId = alignmentParser.readGroups.find(vec[0]);
			len = stringToInt(vec[1]);
			if(len < 1) throw "Max length of read group '" + vec[0] + "' is < 1!";

			//add a new readgroup for the truncated reads to the header
			readGroup = vec[0] + "_truncated";
			alignmentParser.bamHeader.ReadGroups.Add(readGroup);
			alignmentParser.readGroups.fill(alignmentParser.bamHeader);
			truncatedReadGroupId = alignmentParser.readGroups.find(readGroup);

			//copy original tags to truncated read groups
			trunc = alignmentParser.bamHeader.ReadGroups.Begin()+truncatedReadGroupId;
			orig = alignmentParser.bamHeader.ReadGroups.Begin()+readGroupId;
			trunc->Library = orig->Library;
			trunc->PlatformUnit = orig->PlatformUnit;
			trunc->PredictedInsertSize = orig->PredictedInsertSize;
			trunc->ProductionDate = orig->ProductionDate;
			trunc->Program = orig->Program;
			trunc->Sample = orig->Sample;
			trunc->SequencingCenter = orig->SequencingCenter;
			trunc->SequencingTechnology = orig->SequencingTechnology;

			//add to map
			singleEndRG.insert(std::pair<int, TReadGroupMaxLength>(readGroupId, TReadGroupMaxLength(len, truncatedReadGroupId, readGroup, 0)));
		}
	}
	logfile->done();
	logfile->conclude("read " + toString(singleEndRG.size()) + " read groups to be split.");

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_splitRG.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;
	TAlignment alignment;

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	while (alignmentParser.readNextAlignment(alignment)){
		//check if this RG needs to be parsed
		readGroupId = alignment.readGroupId;

		singleEndRGIT = singleEndRG.find(readGroupId);
		if(singleEndRGIT != singleEndRG.end()){
			//check length
			if(alignment.getBamAlignmentLength() == singleEndRGIT->second.maxLen){
				alignment.updateOptionalSamField("RG", singleEndRGIT->second.truncatedOrMergedReadGroup);
//				bamAlignment.EditTag("RG", "Z", singleEndRGIT->second.truncatedReadGroup);
			} else if(alignment.getBamAlignmentLength() > singleEndRGIT->second.maxLen){
				if(allowForLarger)
					alignment.updateOptionalSamField("RG", singleEndRGIT->second.truncatedOrMergedReadGroup);
				else {
					logfile->warning("Length of read " + alignment.name() + " in read group '" + readGroup + "' is > max length provided! Ignoring read.");
				}
			}
		}

		//write
		alignment.save(bamWriter, alignmentParser.genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, alignmentParser.qualMap);

		//report
		reporter.printProgress();
	}

	//close bam writer
	bamWriter.Close();

	//now generate bam index
	indexBamFile(filename);

	//report
	reporter.printEnd();
}


void TGenome::mergeReadGroups(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);

	//read read groups and their expected lengths
	std::string filename = params.getParameterString("readGroups");
	logfile->listFlush("Reading read groups to be merged from file '" + filename + "' ...");
	std::vector< std::vector<std::string> > readGroupsToMerge;
	std::vector< std::vector<std::string> >::reverse_iterator rIt;
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "!";

	//construct new read groups in new header object
	BamTools::SamHeader newHeader(alignmentParser.bamHeader);
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
	int* readGroupMap = new int[alignmentParser.readGroups.size()];
	for(size_t i=0; i<alignmentParser.readGroups.size(); ++i) readGroupMap[i] = -1;
	logfile->startIndent("The following read groups will be merged:");
	std::vector< std::vector<std::string> >::iterator mergeIt = readGroupsToMerge.begin();
	int oldId;
	for(size_t rg = 0; rg < newReadGroupObject.size(); ++rg, ++mergeIt){
		logfile->startIndent("New read group '" + newReadGroupObject.getName(rg) + "' will contain read groups:");
		for(std::vector<std::string>::iterator it = mergeIt->begin(); it != mergeIt->end(); ++it){
			logfile->list(*it);
			oldId = alignmentParser.readGroups.find(*it);
			if(readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
			readGroupMap[oldId] = rg;
		}
		logfile->endIndent();
	}
	logfile->endIndent();

	//now add read groups that will not be merged
	bool printed = false;
	std::string name;
	for(size_t i = 0; i < alignmentParser.readGroups.size(); ++i){
		//check if it is mapped, otherwise add
		if(readGroupMap[i] < 0){
			if(!printed){
				logfile->startIndent("The following read groups will be kept as is:");
				printed = true;
			}
			name = alignmentParser.readGroups.getName(i);
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
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, newHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//prepare reporting
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	while (alignmentParser.readNextAlignment(alignment)){
		//get read group info
		alignment.bamAlignment.GetTag("RG", readGroup);
		oldId = alignmentParser.readGroups.find(readGroup);

		//save as new RG
		alignment.updateOptionalSamField("RG", newReadGroupObject.getName(readGroupMap[oldId]));
//		alignment.bamAlignment.EditTag("RG", "Z", newReadGroupObject.getName(readGroupMap[oldId]));
		alignment.save(bamWriter, alignmentParser.genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, alignmentParser.qualMap);

		//report
		reporter.printProgress();
	}

	//close bam writer
	bamWriter.Close();

	//report
	reporter.printEnd();
};

void TGenome::parseSplitMergeReadGroupSettings(TParameters & params, std::map<int, TReadGroupMaxLength> & RGSettings){
	//do we have to ignore read groups?
	std::map<std::string,int> readGroupsToIgnore;
	if(params.parameterExists("ignoreReadGroups")){
		std::string ignoredReadGroupsFile = params.getParameterString("ignoreReadGroups");
		logfile->listFlush("Reading read groups to ignore from file '" + ignoredReadGroupsFile + "' ...");
		std::ifstream file(ignoredReadGroupsFile.c_str());
		if(!file) throw "Failed to open file '" + ignoredReadGroupsFile + "!";
		int lineNum = 0;

		while(file.good() && !file.eof()){
			std::vector<std::string> vec;
			fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
			++lineNum;
			if(!vec.empty() && vec.size() != 1)
				throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + ignoredReadGroupsFile + "'!";
			else if(vec.size() == 1){
				readGroupsToIgnore.emplace(vec[0],1);
			}
		}
	}

	std::string readGroupSettingsFile = params.getParameterString("readGroupSettings");
	logfile->listFlush("Reading single end read groups from file '" + readGroupSettingsFile + "' ...");
	std::ifstream file(readGroupSettingsFile.c_str());
	if(!file) throw "Failed to open file '" + readGroupSettingsFile + "!";
	int lineNum = 0;

	std::vector<std::string> paired;
	std::vector<std::string> mixed;
	std::vector<std::string> single;

	//read settings file
	while(file.good() && !file.eof()){
		std::vector<std::string> vec;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		++lineNum;

		//parse line
		if(!vec.empty()){
			if(readGroupsToIgnore.size() == 0 || readGroupsToIgnore.find(vec[0]) == readGroupsToIgnore.end()){
				//get RG info
				std::string sequencingType = vec[1];
				int readGroupId = alignmentParser.readGroups.find(vec[0]);

				//add needed new RG
				if(sequencingType == "paired"){
					if(vec.size() != 2)
						throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + readGroupSettingsFile + "'!";
					paired.push_back(vec[0]);
					RGSettings.emplace(readGroupId, TReadGroupMaxLength(-1, readGroupId, vec[0], 2));

				} else if(sequencingType == "mixed"){
					if(vec.size() != 3)
						throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + readGroupSettingsFile + "'!";
					std::string readGroupTruncated = vec[0] + "_truncated";
					mixed.push_back(vec[0]);
					mixed.push_back(readGroupTruncated);
					alignmentParser.bamHeader.ReadGroups.Add(readGroupTruncated);

					int truncatedReadGroupId = alignmentParser.readGroups.find(readGroupTruncated);
					int len = stringToInt(vec[2]);
					if(len < 1) throw "Max length of read group '" + vec[0] + "' is < 1!";
					RGSettings.emplace(readGroupId, TReadGroupMaxLength(len, truncatedReadGroupId, readGroupTruncated, 1));

				} else if(sequencingType == "single"){
					if(vec.size() != 3)
						throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + readGroupSettingsFile + "'!";
					std::string readGroupTruncated = vec[0] + "_truncated";
					single.push_back(vec[0]);
					single.push_back(readGroupTruncated);
					alignmentParser.bamHeader.ReadGroups.Add(readGroupTruncated);

					int truncatedReadGroupId = alignmentParser.readGroups.addTruncatedOrMergedRG(alignmentParser.bamHeader, vec[0], readGroupTruncated);
					int len = stringToInt(vec[2]);
					if(len < 1) throw "Max length of read group '" + vec[0] + "' is < 1!";
					RGSettings.emplace(readGroupId, TReadGroupMaxLength(len, truncatedReadGroupId, readGroupTruncated, 0));

				} else {
					throw "Unknown sequencing type '" + sequencingType + "' on line " + toString(lineNum) + " in file '" + readGroupSettingsFile + "'! Expected 'single', 'mixed' or 'paired'.";
				}
			}
		}
	}
	logfile->done();
	logfile->conclude("read " + toString(single.size() / 2) + " single-end read groups to be split.");
	logfile->conclude("read " + toString(mixed.size() / 2) + " mixed read groups to be split and merged.");
	logfile->conclude("read " + toString(paired.size()) + " paired read groups to be merged.");

	unsigned int predictedSize = single.size() + mixed.size() + paired.size();
	if(predictedSize > alignmentParser.readGroups.size())
		throw "Number of read groups in header incorrect!";
	else if(predictedSize < alignmentParser.readGroups.size()){
		logfile->conclude("Will write " + toString(alignmentParser.readGroups.size() - (single.size() + mixed.size() + paired.size())) + " read group(s) to BAM unchanged (due their presence in ignoreReadGroups or their absence from readGroupSettingsFile).");
	}
};

void TGenome::findPairedReadGroupsToMergeReads(TParameters & params, std::vector<bool> & pairedReadGroups){
	std::string pairedRG = params.getParameterStringWithDefault("pairedReadGroups", "all");
	if(pairedRG == "all"){
		//all are used, initialize to true
		logfile->list("Will merge pairs in all read groups");

		for(size_t i=0; i<pairedReadGroups.size(); ++i)
			pairedReadGroups[i] = true;
	} else {
		//change the paired to true
		std::vector<std::string> vec;
		fillVectorFromString(pairedRG, vec, ',');
		logfile->startIndent("Will only merge pairs in the following read groups:");
		for(unsigned int i=0; i<vec.size(); ++i){
			pairedReadGroups[alignmentParser.readGroups.find(vec.at(i))] = true;
			logfile->list(vec.at(i));
		}
		logfile->endIndent();
	}
};

void TGenome::filterBAM(TParameters & params){

	//soft clipped bases will translated to 'N'

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();
	alignmentParser.setUpdateBlacklistToTrue();
	alignmentParser.setWriteBlacklistToFileToTrue();

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_filtered.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	if(alignmentParser.hasPMD)
		logfile->warning("PMD is given but not relevant for filtering!");

	//create alignment storage
	int maxDist = params.getParameterIntWithDefault("acceptedDistance", 2000);
	TAlignmentMerger merger(&bamWriter, &alignmentParser, maxDist);
	logfile->list("Mates that are farther than " + toString(maxDist) + " apart will be considered orphans. (parameter 'acceptedDistance')");
	if(params.parameterExists("keepOrphans")){
		logfile->list("Will keep keep orphaned reads.");
		merger.keepOrphans();
	} else {
		logfile->list("Will ignore orphaned reads and not write them to BAM (use 'keepOrphans' to keep them).");
	}

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	int curChr = 0;

	while(alignmentParser.readNextAlignment(alignment) && alignmentParser.getNumAlignmentsRead()){
		//if on new chromosome, empty storage
		if(curChr != alignment.chrNumber){
			//write all ready currently in storage
			merger.clear();
			curChr = alignment.chrNumber;
		}

		//check if first alignment in storage is too far away from current read (after checking for chr change)
		//if yes, first alignment in storage is considered an orphan
		merger.writeUpTo(alignment.getPosition());

		//attempt merging of paired reads
		if(alignment.isPaired){
			//Ignore reads in black list
			if(alignmentParser.isInBlacklist(alignment.name()) || !alignment.isProperPair){
				merger.addAsImproperPair(alignment);
			} else {
				//is a proper pair: attempt merging
				merger.checkForMateAndWriteUnmerged(alignment);
			}
		}

		//read is in single-end read group
		else {
			//Ignore reads in black list
			if(alignmentParser.isInBlacklist(alignment.name())){
				alignmentParser.removeFromBlacklist(alignment, "read was in the blacklist");
			} else {
				merger.addReadyToBeWritten(alignment);
			}
		}

		//report
		reporter.printProgress();
	}

	//write unwritten reads
	merger.clear();

	//report end
	reporter.printEnd();

	//close bam writer
	bamWriter.Close();

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
};

void TGenome::setMergerSettings(TParameters & params, TAlignmentMerger & merger){
	//which data to keep
	if(params.parameterExists("keepOrphans")){
		logfile->list("Will keep keep orphaned reads.");
		merger.keepOrphans();
	} else {
		logfile->list("Will ignore orphaned reads and not write them to BAM (use 'keepOrphans' to keep them).");
	}

	//set merging method
	std::string method = params.getParameterStringWithDefault("mergingMethod", "keepRandomRead");
	if(method == "keepRandomRead"){
		logfile->list("Will keep random read for all of overlapping positions");
		merger.keepRandomRead();
	} else if(method == "keepRandomBase"){
		logfile->list("Will keep random base at overlapping positions.");
		merger.keepRandomBase();
	} else if(method == "keepHighestQualBase"){
		logfile->list("Will keep base with higher quality score at overlapping positions.");
	} else
		throw "Unknown merging method " + method + "! Use 'keepRandomRead', 'keepRandomBase' or 'keepHighestQualBase'.";


	//decide if we update quality score
	if(params.parameterExists("updateQuality")){
		logfile->list("Will update quality scores of prefered bases to reflect information from overlapping bases.");
	} else {
		merger.keepOriginalQuality();
		logfile->list("Will keep original quality scores of the prefered bases (use updateQuality to update quality scores).");
	}
}

void TGenome::splitMerge(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();
	alignmentParser.setUpdateBlacklistToTrue();
	alignmentParser.setWriteBlacklistToFileToTrue();

	if(alignmentParser.hasPMD)
		logfile->warning("PMD is given but not relevant for read merging.");

	//which read groups are paired-end?
	std::map<int, TReadGroupMaxLength> RGSettings;
	parseSplitMergeReadGroupSettings(params, RGSettings);

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_mergedReads.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";


	//create alignment storage
	int maxDist = params.getParameterIntWithDefault("acceptedDistance", 2000);
	logfile->list("Mates that are farther than " + toString(maxDist) + " apart will be considered orphans. (parameter 'acceptedDistance')");

	TAlignmentMerger merger(&bamWriter, &alignmentParser, maxDist);
	setMergerSettings(params, merger);

	//for splitter: allow for reads longer than max provided?
	if(params.parameterExists("allowForLarger")){
		merger.allowForLarger();
		logfile->list("Adding single end reads that are longer than maxLength to 'truncated' read group without throwing an error. (parameter 'allowForLarger')");
	}

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	int curChr = 0;
	while (alignmentParser.readNextAlignment(alignment) && alignmentParser.getNumAlignmentsRead()){
		//if on new chromosome, empty storage
		if(curChr != alignment.chrNumber){
			//write all ready currently in storage
			merger.clear();
			curChr = alignment.chrNumber;
		}

		//check if first alignment in storage is too far away from current read (after checking for chr change)
		//if yes, first alignment in storage is considered an orphan
		merger.writeUpTo(alignment.getPosition());

		//attempt merging of paired reads
		std::map<int, TReadGroupMaxLength>::iterator RGSettingsIt;
		RGSettingsIt = RGSettings.find(alignment.readGroupId);

		if(RGSettingsIt != RGSettings.end()){

			if(alignment.isPaired && (RGSettingsIt->second.sequencingType == 1 || RGSettingsIt->second.sequencingType == 2)){
				//Ignore reads in black list
				if(alignmentParser.isInBlacklist(alignment.name()) || !alignment.isProperPair){
					merger.addAsImproperPair(alignment);
				} else {
					//is a proper pair: attempt merging
					merger.addToBeMerged(alignment, randomGenerator);
				}

			//attempt splitting of single end reads
			} else if(!alignment.isPaired && (RGSettingsIt->second.sequencingType == 1 || RGSettingsIt->second.sequencingType == 0)){
				//Ignore reads in black list
				if(alignmentParser.isInBlacklist(alignment.name())){
					alignmentParser.removeFromBlacklist(alignment, "read was in the blacklist");
				} else {
					merger.addToBeSplit(alignment, RGSettingsIt);
				}
			} else throw "nonsensical settings found for read " + alignment.name();

		//read is in read group without a specified type -> should not be split or merged
		} else {
			//Ignore reads in black list
			if(alignmentParser.isInBlacklist(alignment.name())){
				alignmentParser.removeFromBlacklist(alignment, "read was in the blacklist");
			} else {
				merger.addReadyToBeWritten(alignment);
			}
		}

		//report
		reporter.printProgress();
	}

	//write unwritten reads
	merger.clear();

	//report end
	reporter.printEnd();

	//close bam writer
	bamWriter.Close();

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
};

void TGenome::mergePairedEndReads(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();
	alignmentParser.setUpdateBlacklistToTrue();
	alignmentParser.setWriteBlacklistToFileToTrue();

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_mergedReads.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	if(alignmentParser.hasPMD)
		logfile->warning("PMD is given but not relevant for read merging.");

	//which read groups are paired-end?
	std::vector<bool> mergeThisReadGroup(alignmentParser.numReadGroups(), false);
	findPairedReadGroupsToMergeReads(params, mergeThisReadGroup);

	//create alignment storage

	int maxDist = params.getParameterIntWithDefault("acceptedDistance", 2000);
	logfile->list("Mates that are farther than " + toString(maxDist) + " apart will be considered orphans. (parameter 'acceptedDistance')");
	TAlignmentMerger merger(&bamWriter, &alignmentParser, maxDist);
	setMergerSettings(params, merger);

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	int curChr = 0;

	while (alignmentParser.readNextAlignment(alignment) && alignmentParser.getNumAlignmentsRead()){
		//if on new chromosome, empty storage
		if(curChr != alignment.chrNumber){
			//write all ready currently in storage
			merger.clear();
			curChr = alignment.chrNumber;
		}

		//check if first alignment in storage is too far away from current read (after checking for chr change)
		//if yes, first alignment in storage is considered an orphan
		merger.writeUpTo(alignment.getPosition());

		//attempt merging of paired reads
		if(alignment.isPaired && mergeThisReadGroup[alignment.readGroupId]){
			//Ignore reads in black list
			if(alignmentParser.isInBlacklist(alignment.name()) || !alignment.isProperPair){
				merger.addAsImproperPair(alignment);
			} else {
				//is a proper pair: attempt merging
				merger.addToBeMerged(alignment, randomGenerator);
			}
		}

		//read is in single-end read group
		else {
			//Ignore reads in black list
			if(alignmentParser.isInBlacklist(alignment.name())){
				alignmentParser.removeFromBlacklist(alignment, "read was in the blacklist");
			} else {
				merger.addReadyToBeWritten(alignment);
			}
		}

		//report
		reporter.printProgress();
	}

	//write unwritten reads
	merger.clear();

	//report end
	reporter.printEnd();

	//close bam writer
	bamWriter.Close();

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
};

void TGenome::fillVectorOfDownsamplingProbabilities(std::string prob, std::vector<double> & downSampleProbVector){
	if(!stringContainsOnly(prob, "-0123456789.,")){
		throw "Wrong format on probability list: use floating point numbers delimited by commas (e.g. 0.1,0.2,0.5).";
	}

	fillVectorFromString(prob, downSampleProbVector, ',');

	//make sure it is sorted
	std::sort(downSampleProbVector.begin(), downSampleProbVector.end());
	std::reverse(downSampleProbVector.begin(), downSampleProbVector.end());

	//check if all numbers are within (0,1]
	for(double it: downSampleProbVector){
		if(it <= 0.0 || it > 1.0){
			throw "All downsampling probabilities must be within (0,1.0], which is not the case for " + toString(it) + "!";
		}
	}
};

void TGenome::downSampleBamFile(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setUpdateBlacklistToTrue();
	alignmentParser.setWriteBlacklistToFileToTrue();
	TReadList keep;

	//read downsampling rate
	std::vector<double> downSampleProbVector;
	fillVectorOfDownsamplingProbabilities(params.getParameterString("prob"), downSampleProbVector);
	if(downSampleProbVector.size() < 1)
		throw "No downsampling probabilities provided!";

	//read how many replicates
	int numProbs = downSampleProbVector.size();

	//check if probs are between 0 and 1, save in array and print them
	logfile->list("Will accept reads with probabilities (parameter 'prob'): " + concatenateString(downSampleProbVector, ", "));
	if(*downSampleProbVector.begin() == 1.0) logfile->warning("Probability of 1 will result in identical file!");

	//open bam files for writing
	BamTools::BamWriter* bamWriter = new BamTools::BamWriter[numProbs];
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->startIndent("Writing results to the following files:");

	for(int i=0; i<numProbs; ++i){
		//construct and print filename
		std::string filename = outputName + "_downsampled_" + toString(downSampleProbVector[i]) + ".bam";
		logfile->list(filename);
		//open file
		if(!bamWriter[i].Open(filename, alignmentParser.bamHeader, references))	throw "Failed to open BAM file '" + filename + "'!";
	}
	logfile->endIndent();

	//other temp variables
	long counter = 0;
	double r;
	TGenotypeMap genoMap;
	TQualityMap qualMap;

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	while (alignmentParser.readNextAlignment(alignment)){
		++counter;

		//accept read or not?
		for(int i=0; i<numProbs; ++i){
			if(alignmentParser.isInBlacklist(alignment.name())){
				alignmentParser.removeFromBlacklist(alignment, "was in blacklist");
				continue;
			} if(keep.isInReadList(alignment.name())){
				alignment.save(bamWriter[i], genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, qualMap);
			} else {
				r = randomGenerator->getRand(); //inside loop to avoid correlation when multiple probs
				if(r < downSampleProbVector[i]){
					alignment.save(bamWriter[i], genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, qualMap);
					if(alignment.isProperPair)
						keep.addToReadList(alignment, "passed downsampling");
				} else {
					if(alignment.isProperPair){
						alignmentParser.addToBlacklist(alignment, "did not pass downsampling");
					}
				}
			}
		}

		//report
		reporter.printProgress();
	}

	//close bam writer and clean up memory
	for(int i=0; i<numProbs; ++i){
		bamWriter[i].Close();

		std::string filename = outputName + "_downsampled_" + toString(downSampleProbVector[i]) + ".bam";

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

	delete[] bamWriter;

	//report end
	reporter.printEnd();
};

void TGenome::downSampleReads(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//read parameters
	double fraction = params.getParameterDoubleWithDefault("prob", 0.1);
	logfile->list("Each base has a probability of " + toString(fraction)+ " of being masked. (parameter 'prob')");

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_downsampledReads_" + toString(fraction) + ".bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	TGenotypeMap genoMap;
	TQualityMap qualMap;

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and write alignments
	while (alignmentParser.readNextAlignment(alignment)){
		alignment.downsampleAlignment(fraction, *randomGenerator, qualMap);
		alignment.save(bamWriter, genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting, qualMap);

		//report
		reporter.printProgress();
	}

	//close bam writer
	bamWriter.Close();

	//report end
	reporter.printEnd();
};

void TGenome::diagnoseBamFile(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);

	//get max params
	int maxMQ = params.getParameterIntWithDefault("maxMQ", 100);
	int maxReadLength = params.getParameterIntWithDefault("maxReadLength", 1000);

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

    //calculate length of genome
    double totLength = (double) alignmentParser.calcReferenceLength();

    //other temp variables
    std::vector<double> depth;
    double totalDepth = 0.0;
    int numProperPairs = 0;
    long sumFragLen = 0;
    long sumSquaredFragLen = 0;
    int numReadGroups = alignmentParser.readGroups.size();

    long** mappingQuality = new long*[numReadGroups];
    long** readLength = new long*[numReadGroups];
    for(int i = 0; i < numReadGroups; ++i){
    	depth.push_back(0);
    	mappingQuality[i] = new long[maxMQ + 1]; //+1 for zero bin
    	readLength[i] = new long[maxReadLength + 1];
    	for(int j=0; j<100; ++j) mappingQuality[i][j]=0;
    	for(int j=0; j<500; ++j) readLength[i][j]=0;
    }

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

    //now parse through bam file and sum number of aligned bases
	while (alignmentParser.readNextAlignment(alignment)){
        //fragment length
        if(alignment.isProperPair){
        	if(!alignment.isReverseStrand){
        		++numProperPairs;
        		int32_t insSize = alignment.getInsertSize();
        		sumFragLen += abs(insSize);
        		sumSquaredFragLen += (insSize * insSize);
        	}
        }

        //depth
        int length = alignment.getUsableLength(alignmentParser.minQual, alignmentParser.maxQual);
        totalDepth += length;
        depth[alignment.readGroupId] += length;

        //mapping quality
        if(alignment.mappingQuality > maxMQ)
        	throw "Mapping quality of alignment " + alignment.name() + " is larger than maxMQ (" + toString(alignment.mappingQuality) + ">" + toString(maxMQ) +")";
        ++mappingQuality[alignment.readGroupId][alignment.mappingQuality];

        //read length
        if(alignment.getBamAlignmentLength() > maxReadLength)
    	   throw "Read length of alignment " + alignment.name() + " is larger than maxReadLength (" + toString(alignment.getParsedLength()) + ">" + toString(maxReadLength) +")";

        ++readLength[alignment.readGroupId][alignment.getBamAlignmentLength()];

        //report
        reporter.printProgress();
    }

	//report end
	reporter.printEnd();
	logfile->list("Approximate sequencing depth was estimated at " + toString(totalDepth/totLength));

	//writing output files
	logfile->listFlush("Writing to output files ...");

    //depth
    outputDepth << "readGroup\tApproximate_depth";
    outputDepth << "\nallReadGroups\t" << totalDepth/totLength;
    for(int r=0; r<numReadGroups; ++r){
        outputDepth << "\n" << alignmentParser.readGroups.getName(r) << "\t" << depth[r]/totLength;
    }
    outputDepth << "\n";

    //MQ
    long tot;
    outputMQ << "readGroup\tMapping_quality\tCount";
    for(int i=0; i<100; ++i){
    	tot = 0;
    	for(int r=0; r<numReadGroups; ++r) tot += mappingQuality[r][i];
		outputMQ << "\nallReadGroups\t" << i << "\t" << tot;

    }
    for(int r=0; r<numReadGroups; ++r){
        for(int i=0; i<100; ++i){
            outputMQ << "\n" << alignmentParser.readGroups.getName(r) << "\t" << i << "\t" << mappingQuality[r][i];
        }
    }
    outputMQ << "\n";

    //RL
    outputReadLen << "readGroup\tRead_length\tCount";
    for(int i=0; i<500; ++i){
    	tot = 0;
    	for(int r=0; r<numReadGroups; ++r) tot += readLength[r][i];
		outputReadLen << "\nallReadGroups\t" << i << "\t" << tot;
    }
    for(int r=0; r<numReadGroups; ++r){
        for(int i=0; i<500; ++i){
            outputReadLen << "\n" << alignmentParser.readGroups.getName(r)<< "\t" << i << "\t" << readLength[r][i];
        }
    }
    outputReadLen << "\n";

    //FL
    float mean = float(sumFragLen)/float(numProperPairs);
    float var = float(sumSquaredFragLen) / float(numProperPairs) - (mean*mean);
    fragmentStats << "mean: " << mean << "\n" << "variance: " << var << "\n";

    logfile->done();

    //clena up
    outputDepth.close();
    outputMQ.close();
    outputReadLen.close();
    fragmentStats.close();

    for(int i = 0; i < numReadGroups; ++i){
    	delete[] mappingQuality[i];
    	delete[] readLength[i];
    }
    delete [] mappingQuality;
    delete [] readLength;
}

void TGenome::allelicDepth(TParameters & params){
	//allocate table
	if(alignmentParser.getMaxDepth() > 100){
		logfile->warning("Allocating count table for a max depth of " + toString(alignmentParser.getMaxDepth()) + " uses a lot of memory! Use argument maxDepth to limit.");
	}
	TAllelicDepthCounts counts(alignmentParser.getMaxDepth());

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			window.countAlleles(counts);
			logfile->listFlush("Adding imbalance values to table ...");
			logfile->write(" done!");
		}
	}

	//write to file
	std::string outputFileName = outputName + "_allelicDepth.txt";
	logfile->list("Writing allelic imbalance table to '" + outputFileName + "'");
	bool writeEmpty = params.parameterExists("printAll");
	counts.write(outputFileName, writeEmpty);
};

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
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			//write to file
			logfile->listFlush("Writing sequencing depth estimates to file ...");
			if(window.depth == -1.0) output << alignmentParser.getCurChrName() << "\t" << window.start << "\t" << window.end << "\t" << "0" << "\n";
			else output << alignmentParser.getCurChrName() << "\t" << window.start << "\t" << window.end << "\t" << window.depth << "\n";
			logfile->done();
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
};

void TGenome::estimateDepthPerSite(TParameters & params){
	//initialize count object
	int maxDepth = params.getParameterIntWithDefault("maxDepth", 20);
	TDistributionOfCounts counts(maxDepth, "depth");

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			logfile->listFlush("Adding depth to table ...");
			window.countDepthPerSite(counts);
			logfile->done();
		}
	}

	//write
	std::string filename = outputName + "_distDepthPerSite.txt";
	logfile->listFlush("Writing depth distribution to '" + filename + "' ...");
	counts.writeCounts(filename);
	logfile->done();

	filename = outputName + "_cumulativeDepthPerSite.txt";
	logfile->listFlush("Writing normalized cumulative depth distribution to '" + filename + "' ...");
	counts.writeNormalizedCumulativeCounts(filename);
	logfile->done();

	filename = outputName + "_quantilesDepthPerSite.txt";
	logfile->listFlush("Writing quantiles of depth distribution '" + filename + "' ...");
	counts.writeQuantiles(filename);
	logfile->done();
};

void TGenome::writeDepthPerSite(TParameters & params){
	gz::ogzstream out;

	std::string outputFileName = outputName + "_depthPerSite.txt.gz";
	logfile->list("Writing per site depth to '" + outputFileName + "'");

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
			logfile->listFlush("Writing depth per site ...");
			window.printDepthPerSite(out, alignmentParser.getCurChrName());
			logfile->done();
		}
	}

	//clean up
	out.close();
};

void TGenome::estimateDuplicationCounts(TParameters & params){
	//assembles distribution of how often a read is duplicated
	//now: just how many reads start at the same positions

	//initialize alignment reading
	TAlignment alignment(maxReadLength);

	//create storage
	int maxCounts = params.getParameterIntWithDefault("maxCount", 20);
	TDistributionOfCounts counts(maxCounts, "readStarts");

	//iterate through windows
	int curChr = 0;
	int curChrLength = alignmentParser.chrNumberToLength(curChr);
	int curPos = 0;
	int countsAtPos = 0;
	while (alignmentParser.readNextAlignment(alignment)){
		if(alignment.chrNumber != curChr){
			//add last pos with data
			counts.add(countsAtPos);
			countsAtPos = 0;

			//add all positions until chromosome end to structure
			counts.add(0, curChrLength - curPos);
			curChr = alignment.chrNumber;
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
			countsAtPos = countsAtPos + 1;
		} else
			throw "Bam file is not sorted!";
	}

	//write output
	std::string filename = outputName + "_readStartsPerSite.txt";
	logfile->listFlush("Writing distribution of read starts per site to '" + filename + "' ...");
	counts.writeCounts(filename);
	logfile->done();
};

//---------------------------------------------------
//PMD
//---------------------------------------------------
void TGenome::estimatePMD(TParameters & params){
	//make sure FASTA is open
	if(!alignmentParser.hasReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//prepare maps
	TReadGroupMap readGroupMap(alignmentParser.readGroups, params.getParameterString("poolReadGroups", false), logfile);
	TGenotypeMap genoMap;

	//prepare PMD table
	int maxLengthForInference = params.getParameterIntWithDefault("length", 50);
	logfile->list("Estimating PMD at the first " + toString(maxLengthForInference) + " positions.");
	TPMDTables pmdTables(alignmentParser.readGroups, maxLengthForInference, maxReadLength, readGroupMap);

	//measure progress and runtime
	TBamProgressReporter reporter(&alignmentParser, logfile);

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
	if(!params.parameterExists("onlyEmpiric")){
		filename = outputName + "_PMD_input_Exponential.txt";
		logfile->listFlush("Estimating PMD exponential models and writing them to '" + filename + "' ...");
		int numNRIterations = params.getParameterIntWithDefault("numNRIterations", 100);
		double eps = params.getParameterDoubleWithDefault("eps", 0.001);
		pmdTables.fitExponentialModel(numNRIterations, eps, filename, logfile);
		logfile->done();
	} else {
		logfile->list("Not fitting exponential model due to user specification (parameter 'onlyEmpiric')");
	}
}


void TGenome::runPMDS(TParameters & params){
	//parse bam file and calculate PMDS for each read (seeSkoglund et al. 2014)
	//write new bam file with PMDS score added
	//parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	if(!alignmentParser.hasReference) throw "Cannot run PMDS without reference!";

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

	//other tmp
	TQualityMap qualMap;
	TGenotypeMap genoMap;

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_PMDS.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//now parse through bam file and write alignments
	while(alignmentParser.readNextAlignment(alignment)){
		++counter;

		//calc PMD
		double PMDS = alignment.calculatePMDS(pi, alignmentParser.pmdObjects);

		//update and write
		if(PMDS > minPMDS && PMDS < maxPMDS){
			alignment.updateOptionalSamField("DS", PMDS);
			alignment.save(bamWriter, genoMap, alignmentParser.minQual, alignmentParser.maxQual, qualMap);
		} else ++counterF;

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

};

void TGenome::printMateInformationPerSite(TParameters & params){
	//open output file
	std::string outputFileName = outputName + "_mateInformation.txt.gz";
	logfile->list("Writing mate information to file '" + outputFileName + "'.");
	TOutputFileZipped out(outputFileName);
	out.writeHeader({"chr", "pos", "depth", "numA", "numC", "numG", "numT", "numAlleles", "numFirstMate", "numSecondMate", "numFwd", "numRev"});

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			logfile->listFlush("Writing mate info per site ...");
			window.printMateInformationPerSite(out, alignmentParser.getCurChrName());
			logfile->done();
		}
	}

	//clean up
	out.close();
};

void TGenome::contextStats(TParameters & params){
	std::string outputFileName = outputName + "_contextInformation.txt.gz";
	logfile->list("Writing context information to file '" + outputFileName + "'.");
	TOutputFileZipped out(outputFileName);
	out.writeHeader({"quality","cAA", "cAC", "cAG", "cAT", "cCA", "cCC", "cCG", "cCT", "cGA", "cGC", "cGG", "cGT", "cTA", "cTC", "cTG", "cTT", "cNA", "cNC", "cNG", "cNT", "cAN", "cCN", "cGN", "cTN", "cNN"}); //N means unknwon base or "nothing", i.e. end of read or del
	int numContext = 25;

	//prepare windows
	TWindow window;

	//prepare table
    int** contextCounts = new int*[alignmentParser.maxQual];
    for(int i = 0; i < alignmentParser.maxQual; ++i){
    	contextCounts[i] =  new int[numContext];
    	for(int j=0; j<numContext; ++j)
    		contextCounts[i][j]=0;
    }

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters){
			window.contextStats(contextCounts, alignmentParser.qualMap);
		}
	}

    //write to file
    for(int i=0; i<alignmentParser.maxQual; ++i){
    	out << i;
    	for(int j=0; j<numContext; ++j){
    		out << contextCounts[i][j];
    	}
		out << std::endl;
    }

    for(int i = 0; i < alignmentParser.maxQual; ++i){
    	delete[] contextCounts[i];
    }
    delete [] contextCounts;
}


