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
	if(windowSize < 1000) throw "Window size should be at least 1Kb!";
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 0.95);
	initializePostMortemDamage(params, logfile);

	//outputname
	outname = params.getParameterStringWithDefault("out", "");
	if(outname == ""){
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
	chrIterator = bamHeader.Sequences.End();
};

bool TGenome::iterateChromosome(TWindowPair & windowPair){
	if(chrIterator == bamHeader.Sequences.End()){
		chrIterator = bamHeader.Sequences.Begin();
		chrNumber = 0;
	} else {
		logfile->endIndent();
		logfile->endNumbering();
		//move to next
		++chrIterator;
		++chrNumber;
	}

	//did we reach end?
	if(chrIterator == bamHeader.Sequences.End()){
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

	//write progress
	logfile->startNumbering("Parsing chromosome '" + chrIterator->Name + "':");
	return true;
}

bool TGenome::iterateWindow(TWindowPair & windowPair){
	if(curEnd > 0) logfile->endIndent();

	//move to next region
	curStart = curEnd;
	curEnd += windowSize;
	if(curEnd > chrLength) curEnd = chrLength + 1;
	if(curStart > chrLength) return false;

	//swap windows and move next
	windowPair.swap();
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
				if(windowPair.curPointer->addFromRead(bamAlignement))
					//add also to next window in case reads overhangs current window -> function returns true
					windowPair.nextPointer->addFromRead(bamAlignement);

			}
		}
	}

	windowPair.curPointer->calcCoverage();

	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");
	logfile->conclude("read data from " + toString(numReads) + " reads.");
	logfile->conclude("coverage is " + toString(windowPair.curPointer->coverage));
	logfile->conclude(toString(windowPair.curPointer->fractionsitesCoverageAtLeastTwo) + " of all sites are covered at least twice");
	logfile->conclude(toString(windowPair.curPointer->fractionSitesNoData) + " of all sites have no data");

	return true;
};

void TGenome::initializePostMortemDamage(TParameters & params, TLog* logfile){
	std::string pmdString;
	if(params.parameterExists("pmd")){
		pmdString = params.getParameterString("pmd");
		logfile->list("Initializing post mortem damage for both C->T and G->A with function '" + pmdString +"'.");
		pmdObject.initializeFunction(pmdString, pmdCT);
		pmdObject.initializeFunction(pmdString, pmdGA);
		logfile->conclude(pmdObject.getFunctionString(pmdCT));
	} else {
		//separate model for C->T and G->A
		//first C->T
		if(!params.parameterExists("pmdCT")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdCT' has to be provided!";
		pmdString = params.getParameterString("pmdCT");
		logfile->list("Initializing post mortem C->T damage with function '" + pmdString +"'.");
		pmdObject.initializeFunction(pmdString, pmdCT);
		logfile->conclude(pmdObject.getFunctionString(pmdCT));

		//second G->A
		if(!params.parameterExists("pmdGA")) throw "Problem initializing post mortem damaga: argument 'pmd' or 'pmdGA' has to be provided!";
		pmdString = params.getParameterString("pmdGA");
		logfile->list("Initializing post mortem G->A damage with function '" + pmdString +"'.");
		pmdObject.initializeFunction(pmdString, pmdGA);
		logfile->conclude(pmdObject.getFunctionString(pmdGA));
	}
}

void TGenome::estimateTheta(TParameters & params){
	//Error rate recalibration
	Recalibration recal(params.getParameterString("recal", false));

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
				windows.cur->estimateTheta(EMParams, pmdObject, recal, out, logfile);
			}
		}
	}

	//clean up
	out.close();
}

void TGenome::calcLikelihoodSurfaces(TParameters & params){
	//Error rate recalibration
	Recalibration recal(params.getParameterString("recal", false));

	//read params
	int steps = params.getParameterIntWithDefault("steps", 1000);
	int numWindows = params.getParameterIntWithDefault("numWindows", 1);

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	int windowsCalculated = 0;
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
				out.open((outputName + chrIterator->Name + "_" + toString(windows.cur->start) + "_LLsurface.txt").c_str());
				if(!out) throw "Failed to open output file '" + outputName + "'!";

				//calc surface
				logfile->listFlush("Calculating likelihood surface ...");
				windows.cur->calcLikelihoodSurface(pmdObject, recal, out, steps);
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

void TGenome::printPileup(){
	//prepare windows
	TWindowPairDiploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_pileup.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	GenotypeMap genoMap;
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
			windows.cur->printPileup(pmdObject, out, chrIterator->Name);
		}
	}

	//clean up
	out.close();
}

void TGenome::estimateErrorCalibration(TParameters & params){
	//read params
	double minA = params.getParameterDoubleWithDefault("minA", 0.0);
	if(minA < 0.0) throw "minA has to be > 0.0!";
	double maxA = params.getParameterDoubleWithDefault("maxA", 0.2);
	if(maxA < minA) throw "minA has to be < maxA!";
	double minB = params.getParameterDoubleWithDefault("minB", 0.0);
	if(minB <= 0.0) throw "minB has to be > 0.0!";
	double maxB = params.getParameterDoubleWithDefault("maxB", 1.0);
	if(maxB < minB) throw "minB has to be < maxB!";
	if(maxB > 1.0) throw "maxB has to be < 1.0!";

	int steps = params.getParameterIntWithDefault("steps", 10);
	int numIterations = params.getParameterIntWithDefault("iterations", 2);

	//prepare windows
	TWindowPairHaploid windows;

	//open output
	std::ofstream out;
	std::string filename = outputName + "_calibration.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";
	out << "a\tb\tLL\n";

	//prepare storage for grid search
	std::vector<double> a, b, LL;
	std::vector<double>::iterator aIt, bIt, LLIt;

	//add a=0, b=0
	a.push_back(0.0);
	b.push_back(0.0);
	LL.push_back(0.0);

	//prepare grid for b != 0
	double sizeA = (maxA - minA) / (double)(steps - 1.0);
	double sizeB = (maxB - minB) / (double)(steps - 1.0); //without b=0
	double tmpA, tmpB;
	for(int i=0; i<steps; ++i){
		tmpA = minA + i*sizeB;
		if(tmpA > 0.0){
			for(int j=0; j<steps; ++j){
				//only one with b=0
				tmpB = minB + j*sizeB;
				if(tmpB > 0.0){
					a.push_back(minA + i*sizeA);
					b.push_back(minB + j*sizeB);
					LL.push_back(0.0);
				}
			}
		}
	}

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//calc LL for a=0 b=0
			aIt = a.begin();
			bIt = b.begin();
			LLIt = LL.begin();
			windows.cur->estimateBaseFrequencies();
			windows.cur->calculateEmissionProbabilities(pmdObject);
			*LLIt += windows.cur->calcLogLikelihood();
			++aIt; ++bIt; ++LLIt;

			//calc LL for all combinations of a and b
			for(; aIt != a.end(); ++aIt, ++bIt, ++LLIt){
				windows.cur->calculateEissionProbabilitiesWithScaledError(pmdObject, *aIt, *bIt);
				*LLIt += windows.cur->calcLogLikelihood();
			}
		}
	}

	//Report
	aIt = a.begin();
	bIt = b.begin();
	LLIt = LL.begin();
	for(; aIt != a.end(); ++aIt, ++bIt, ++LLIt){
		out << *aIt << "\t" << *bIt << "\t" << *LLIt << "\n";
	}

	out.close();

}

