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
TGenome::TGenome(TLog* Logfile, TParameters & params){
	logfile = Logfile;

	//read parameters
	filename = params.getParameterString("bam");
	windowSize = params.getParameterDouble("window");
	//if(windowSize < 1000) throw "Window size should be at least 1Kb!";
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 0.95);
	initializePostMortemDamage(params);

	//outputname
	outputName = params.getParameterStringWithDefault("out", "");
	if(outputName == ""){
		//guess from filename
		outputName = filename;
		outputName = extractBeforeLast(outputName, ".");
	}

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;
	curStart = -1;
	curEnd = -1;
	oldPos = -1;

	//open BAM file
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";

	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + filename + "'!";

	//read header
	bamHeader = bamReader.GetHeader();
	readGroups.fill(bamHeader);
	chrIterator = bamHeader.Sequences.End();

	//initialize recalibration
	initializeRecalibration(params);

	//check if we mask sites
	if(params.parameterExists("mask")){
		doMasking = true;
		std::string maskFile = params.getParameterString("mask");
		logfile->list("Will mask all sites listed in BED file '" + maskFile + "'");
		mask = new TBedReader(maskFile, windowSize);
		//mask->print();
	} else doMasking = false;
};

bool TGenome::iterateChromosome(TWindowPair & windowPair){
	if(chrIterator == bamHeader.Sequences.End()){
		chrIterator = bamHeader.Sequences.Begin();
		chrNumber = 0;
	} else {
		logfile->endNumbering();
		//move to next
		++chrIterator;
		++chrNumber;
	}

	//did we reach end?
	if(chrIterator == bamHeader.Sequences.End()){
		curEnd = 0;
		return false;
	}

	//jump reader
	bamReader.Jump(chrNumber, 0);

	//restart windows
	chrLength = stringToLong(chrIterator->Length);
	curStart = 0;
	curEnd = 0;
	oldPos = -1;
	windowPair.nextPointer->move(0, windowSize);

	//advance mask
	if(doMasking) mask->setChr(chrIterator->Name);

	//write progress
	logfile->startNumbering("Parsing chromosome '" + chrIterator->Name + "':");
	return true;
}

bool TGenome::iterateWindow(TWindowPair & windowPair){
	if(curEnd > 0) logfile->endIndent();

	//swap window pairs
	windowPair.swap();

	//move to next region
	curStart = curEnd;
	curEnd += windowSize;
	if(curEnd > chrLength){
		curEnd = chrLength + 1;
		windowPair.curPointer->end = curEnd;
	}
	if(curStart >= chrLength) return false;

	//move next
	windowPair.nextPointer->move(curEnd, curEnd + windowSize);

	//report
	logfile->number("Window [" + toString(curStart) + ", " + toString(curEnd) + ") on '" + chrIterator->Name + "':");
	logfile->addIndent();

	return true;
};

bool TGenome::readData(TWindowPair & windowPair){
	logfile->listFlush("Reading data ...");

	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//parse through reads
	int numReads = 0;
	while(bamReader.GetNextAlignment(bamAlignement) && bamAlignement.RefID==chrNumber){
		//std::cout << "REF ID = " << bamAlignement.RefID << "\tpos = " << bamAlignement.Position << std::endl;
		//only take those reads that pass QC
		if(!bamAlignement.IsFailedQC() && !bamAlignement.IsDuplicate() && bamAlignement.Position >= curStart){
			//check if bam file is sorted
			if(bamAlignement.Position < oldPos){
				throw "BAM file must be sorted by position!";
			}
			oldPos = bamAlignement.Position;

			//check if still within current window and add to window
			if(bamAlignement.Position >= curEnd){
				//nextWindow->addFromRead(bamAlignement);
				break;
			} else {
				++numReads;
				if(windowPair.curPointer->addFromRead(bamAlignement, &readGroups)){
					//add also to next window in case reads overhangs current window -> function returns true
					windowPair.nextPointer->addFromRead(bamAlignement, &readGroups);
				}
			}
		}
	}

	//apply mask
	if(doMasking) windowPair.curPointer->applyMask(mask);

	//calc coverage
	windowPair.curPointer->calcCoverage();

	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");
	logfile->conclude("read data from " + toString(numReads) + " reads.");
	logfile->conclude("coverage is " + toString(windowPair.curPointer->coverage));
	logfile->conclude(toString(windowPair.curPointer->fractionsitesCoverageAtLeastTwo * 100) + "% of all sites are covered at least twice");
	logfile->conclude(toString(windowPair.curPointer->fractionSitesNoData * 100) + "% of all sites have no data");

	return true;
};

void TGenome::initializePostMortemDamage(TParameters & params){
	std::string pmdString;
	if(params.parameterExists("pmd")){
		pmdString = params.getParameterString("pmd");
		logfile->list("Initializing post mortem damage for both C->T and G->A with function '" + pmdString +"'.");
		pmdObject.initializeFunction(pmdString, pmdCT);
		pmdObject.initializeFunction(pmdString, pmdGA);
		logfile->conclude(pmdObject.getFunctionString(pmdCT));
	} else {
		if(!params.parameterExists("pmdCT") && !params.parameterExists("pmdGA")){
			//no post mortem damage
			pmdString = "none";
			pmdObject.initializeFunction(pmdString, pmdGA);
			pmdObject.initializeFunction(pmdString, pmdCT);
		} else {
			//separate model for C->T and G->A
			//first C->T
			if(!params.parameterExists("pmdCT")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdCT' has to be provided!";
			pmdString = params.getParameterString("pmdCT");
			logfile->list("Initializing post mortem C->T damage with function '" + pmdString +"'.");
			pmdObject.initializeFunction(pmdString, pmdCT);
			logfile->conclude(pmdObject.getFunctionString(pmdCT));

			//second G->A
			if(!params.parameterExists("pmdGA")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdGA' has to be provided!";
			pmdString = params.getParameterString("pmdGA");
			logfile->list("Initializing post mortem G->A damage with function '" + pmdString +"'.");
			pmdObject.initializeFunction(pmdString, pmdGA);
			logfile->conclude(pmdObject.getFunctionString(pmdGA));
		}
	}
}

void TGenome::initializeRecalibration(TParameters & params){
	if(params.parameterExists("BQSRQuality")){
		recalObject = new TRecalibrationBQSR(&bamHeader, params, &pmdObject, logfile);
	} else {
		recalObject = new TRecalibration();
	}
}

void TGenome::estimateTheta(TParameters & params){
	//EM params
 	EMParameters EMParams;
	EMParams.numThetaOnlyUpdates = params.getParameterIntWithDefault("iterationsThetaOnly", 10);
	EMParams.numIterations = params.getParameterIntWithDefault("iterations", 100);

	EMParams.maxEpsilon = params.getParameterDoubleWithDefault("maxEps", 0.000001);
	EMParams.NewtonRalphsonNumIterations = params.getParameterIntWithDefault("NRiterations", 10);
	EMParams.NewtonRalphsonMaxF = params.getParameterDoubleWithDefault("maxF", 0.00001);

	//params regarding initial search
	EMParams.initalTheta = params.getParameterDoubleWithDefault("initTheta", 0.01);
	EMParams.initThetaSearchFactor = params.getParameterDoubleWithDefault("initThetaSearchFactor", 100);
	EMParams.initThetaNumSearchIterations = params.getParameterDoubleWithDefault("initThetaNumSearchIterations", 10);

	//open output
	std::ofstream out;
	out.open((outputName + "_theta_estimates.txt").c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	out << std::setprecision(9) << "Chr\t";
	out << "start\tend\tcoverage\tmissing\ttwoOrMore\tpi(A)\tpi(C)\tpi(G)\tpi(T)\ttheta_MLE\ttheta_C95_l\ttheta_C95_u\tLL";
	out << "\n";

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//check if we have data -> can be extended to ensure
			if(windows.cur->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				//estimate Theta
				out << chrIterator->Name << "\t";
				windows.cur->estimateTheta(EMParams, pmdObject, recalObject, out, logfile);
			}
		}
	}

	//clean up
	out.close();
}

void TGenome::calcLikelihoodSurfaces(TParameters & params){
	//read params
	int steps = params.getParameterIntWithDefault("steps", 100);
	int numWindows = params.getParameterIntWithDefault("numWindows", 1);

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	int windowsCalculated = 0;
	std::string filename;
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//check if we have data -> can be extended to ensure
			if(windows.cur->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				//open file
				std::ofstream out;
				filename = outputName + chrIterator->Name + "_" + toString(windows.cur->start) + "_LLsurface.txt";
				out.open(filename.c_str());
				if(!out) throw "Failed to open output file '" + outputName + "'!";

				//calc surface
				logfile->listFlush("Calculating likelihood surface ...");
				windows.cur->calcLikelihoodSurface(pmdObject, recalObject, out, steps);
				logfile->write(" done!");

				//close output
				out.close();

				//check if we break
				++windowsCalculated;
				if(windowsCalculated >= numWindows) break;
			}
		}
		if(windowsCalculated >= numWindows) break;
	}
}

void TGenome::callMLEGenotypes(TParameters & params){
	//open output
	std::string filename = outputName + "_MLEGenotypes.txt.gz";
	logfile->list("Writing MLE genotypes to '" + filename + "'");
	gz::ogzstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";

	//do we print sites with no data?
	bool printIfNoData = params.parameterExists("printAll");

	//write header
	out << std::setprecision(3);
	out << "chr\tpos\tcoverage\tL(AA)\tL(AC)\tL(AG)\tL(AT)\tL(CC)\tL(CG)\tL(CT)\tL(GG)\tL(GT)\tL(TT)\tMLE\tQ\n";

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//call genotypes
			logfile->listFlush("Calling MLE genotypes ...");
			windows.cur->callMLEGenotype(pmdObject, recalObject, out, chrIterator->Name, printIfNoData);
			logfile->write(" done!");
		}
	}

	//clean up
	out.close();
}

void TGenome::printPileup(){
	//prepare windows
	TWindowPairDiploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_pileup.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	TGenotypeMap genoMap;
	out << "Chr\tposition\tcoverage\tbases";
	for(int i=0; i<10; ++i)
		out << "\tPem(" << genoMap.getGenotypeString(i) << ")";
	out << "\n";

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//print pileup
			windows.cur->printPileup(pmdObject, recalObject, out, chrIterator->Name);
		}
	}

	//clean up
	out.close();
}


void TGenome::estimateErrorCalibrationEM(TParameters & params){
	//Initialize calibration parameters.
	//We assume a linear model log(error) = eta = beta_0 + beta_1 * quality + beta_2 * position in reads
	//                                            + gamma_A * Ind(d = A) + gamma_C * Ind(d = C) + gamma_G * Ind(d = G) + gamma_T * Ind(d = T)

	/*

	//read EM parameters
	int numEMIterations = params.getParameterIntWithDefault("iterations", 10);
	logfile->list("Will perform at max " + toString(numEMIterations) + " EM iterations.");
	double maxEpsilon = params.getParameterDoubleWithDefault("maxEps", 0.000001);
	logfile->list("Will stop EM when deltaLL < " + toString(maxEpsilon));
	double NewtonRalphsonNumIterations = params.getParameterIntWithDefault("NRiterations", 10);
	logfile->list("Will conduct at max " + toString(NewtonRalphsonNumIterations) + " Newton-Ralphson iterations");
	double NewtonRalphsonMaxF = params.getParameterDoubleWithDefault("maxF", 0.00001);
	logfile->list("Will stop Newton-Ralphson when F < " + toString(NewtonRalphsonMaxF));

	//create recalibration object
	TRecalibrationEM recalObject(&params, logfile);

	//prepare windows
	TWindowPairHaploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_calibration.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";
	out << "iteration";
	recalObject.writeHeader(out);

	//run EM iterations
	for(int i=0; i < numEMIterations; ++i){
		logfile->startIndent("EM Iteration " + toString(i+1) + ":");

		//prepare recal object for next iteration
		recalObject.initEMStep();

		//run Newton-Ralphosn iteration
		for(int j=0; j<NewtonRalphsonNumIterations; ++j){
			logfile->startIndent("Calculating Jacobian for Newton-Ralphson iteration " + toString(j+1) + ":");
			recalObject.initNetwonRalphsonStep();

			while(iterateChromosome(windows)){
				while(iterateWindow(windows)){
					//read data for current window
					readData(windows);
					windows.cur->estimateBaseFrequencies();
					windows.cur->addToJacobian(&recalObject);
				}
			}

			//clean up memory
			windows.clear();

			//perform parameter update
			logfile->endIndent();
			logfile->listFlush("Updating parameters ...");
			recalObject.runNewtonRalphson();
			logfile->write(" done!");

			//check if we break
			if(recalObject.maxF < NewtonRalphsonMaxF){
				logfile->list("Stopping Newton-Ralphson since maxF < " + toString(NewtonRalphsonMaxF));
				break;
			}
		}

		//save current parameters
		recalObject.saveParams();

		//Report to file
		recalObject.writeParams(out);
		logfile->endIndent();

		//check if we break EM

	}
	logfile->endIndent();

	//clean up
	out.close();

	*/
}

void TGenome::fillSequence(std::vector<double> & vec, std::string & str){
	//it is either a number, or a sequence min-max:num steps
	std::string::size_type posDash = str.find_first_of('-', 1);
	if(posDash != std::string::npos){
		std::string::size_type posColon = str.find_first_of(':');
		if(posColon != std::string::npos){
			//fill sequence
			double min = stringToDoubleCheck(str.substr(0, posDash));
			double max = stringToDoubleCheck(str.substr(posDash+1, posColon-posDash-1));
			int numSteps = stringToIntCheck(str.substr(posColon+1));
			double step = (max - min) / (double) (numSteps - 1.0);
			for(int i=0; i<numSteps; ++i) vec.push_back(min + (double) i * step);
		} else throw "Unable to understand sequence '" + str + "'!";
	} else {
		//should be a number
		vec.push_back(stringToDoubleCheck(str));
	}
}

void TGenome::calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params){
	/*
	//read vectors of betas to test
	logfile->startIndent("Will calculate likelihood on a grid with these marginal ranges:");
	logfile->startIndent("Will calculate LL for these parameter values:");

	//beta0
	std::vector<double> beta0;
	std::string tmp = params.getParameterString("beta0");
	fillSequence(beta0, tmp);
	logfile->list("beta0 = " + concatenateString(beta0, ", "));

	//beta1
	std::vector<double> beta1;
	tmp = params.getParameterString("beta1");
	fillSequence(beta1, tmp);
	logfile->list("beta1 = " + concatenateString(beta1, ", "));

	//beta2
	std::vector<double> beta2;
	tmp = params.getParameterString("beta2");
	fillSequence(beta2, tmp);
	logfile->list("beta2 = " + concatenateString(beta2, ", "));
	logfile->endIndent();

	//create recalibration object
	TRecalibrationEM recalObject(logfile);

	//prepare windows
	TWindowPairHaploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_calibration_LL_surface.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";
	recalObject.writeHeader(out);
	out << "\tLL\n";

	//run loop over all combinations
	double* recalParams = new double[3];
	for(std::vector<double>::iterator itBeta0=beta0.begin(); itBeta0!=beta0.end(); ++itBeta0){
		recalParams[0] = *itBeta0;
		for(std::vector<double>::iterator itBeta1=beta1.begin(); itBeta1!=beta1.end(); ++itBeta1){
			recalParams[1] = *itBeta1;
			for(std::vector<double>::iterator itBeta2=beta2.begin(); itBeta2!=beta2.end(); ++itBeta2){
				recalParams[2] = *itBeta2;
				logfile->startIndent("Calculating LL for beta0 = " + toString(*itBeta0) + ", beta1 = " + toString(*itBeta1) + ", beta2 = " + toString(*itBeta2) + ":");

				//reset
				recalObject.resetLikelihood();

				//set parameters
				recalObject.setParams(recalParams);

				//loop over all windows
				while(iterateChromosome(windows)){
					while(iterateWindow(windows)){
						//read data for current window
						readData(windows);
						//calc LL
						windows.cur->estimateBaseFrequencies();
						windows.cur->addToLikelihoodRecalibration(&recalObject);
					}
				}

				//clean up memory
				windows.clear();

				//output
				logfile->conclude("LL = " + toString(recalObject.logLikelihood));
				recalObject.writeParams(out);
				out << "\t" << recalObject.logLikelihood << std::endl;
				logfile->endIndent();
			}
		}
	}

	//clean up
	out.close();
	*/
}

void TGenome::BQSR(TParameters & params){
	//read vectors of betas to test
	logfile->startIndent("Estimating recalibration parameters:");

	//prepare windows
	TWindowPairHaploid windows;

	//open FASTA reference
	std::string fastaFile = params.getParameterString("fasta");
	std::string fastaIndex = fastaFile + ".fai";
	BamTools::Fasta reference;
	reference.Open(fastaFile, fastaIndex);

	//create BQSR object
	TRecalibrationBQSR bqsr(&bamHeader, params, &pmdObject, logfile);
	if(bqsr.allConverged()){
		logfile->list("No need to estimate any BQSR cells. Aborting Program.");
		return;
	}

	//loop over bam until BQSR converges
	bool hasConverged = false;
	int loopNumber = 0;
	while(!hasConverged){
		++loopNumber;
		logfile->startIndent("Running recalibration loop " + toString(loopNumber) + ":");
		//loop over all windows
		while(iterateChromosome(windows)){
			while(iterateWindow(windows)){
				//read data for current window
				readData(windows);

				//add reference data
				windows.cur->addReferenceBaseToSites(reference, chrNumber);

				//add the base to BQSR
				windows.cur->addSitesToBQSR(bqsr, logfile);
			}
		}

		//clean up memory
		windows.clear();

		//estimate epsilon
		hasConverged = bqsr.estimateEpsilon();
		logfile->endIndent();
	}

	//write results to file
	bqsr.writeToFile(outputName);

}
