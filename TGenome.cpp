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
	initializeRandomGenerator(params);

	//initialize alignment parser
	maxReadLength = params.getParameterIntWithDefault("maxReadLength", 1000);
	alignmentParser.init(maxReadLength, params, logfile);

	//outputname
	outputName = params.getParameterStringWithDefault("out", "");
	if(outputName == ""){
		//guess from filename
		outputName = alignmentParser.filename;
		outputName = extractBeforeLast(outputName, ".");
	}
	logfile->list("Writing output files with prefix '" + outputName + "'.");

	//open FASTA reference
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		std::string fastaIndex = fastaFile + ".fai";
		logfile->list("Reading reference sequence from '" + fastaFile + "'");
		if(!reference.Open(fastaFile, fastaIndex)) throw "Failed to open FASTA file '" + fastaFile + "'! Is index file present?";
		fastaReference = true;
		alignmentParser.addReference(&reference);
	} else fastaReference = false;

/*	//initialize post mortem damage
	hasPMD = false;
	initializePostMortemDamage(params);
	doRecalibration = false;
	recalObjectInitialized = false;*/



	//trimming ends
	if(params.parameterExists("trim3") || params.parameterExists("trim5")){
		int trim3 = params.getParameterIntWithDefault("trim3", 0);
		if(trim3 < 0) throw "trimming distance trim3 must be >= 0!";
		int trim5 = params.getParameterIntWithDefault("trim5", 0);
		if(trim5 < 0) throw "trimming distance trim5 must be >= 0!";
		if(trim3>0 || trim5>0){
			alignmentParser.setReadTrimming(trim3, trim5);
			logfile->list("Will trim first " + toString(trim3) + " and " + toString(trim5) + " bases from the 3' and 5' end, respectively.");
		}
	}

};

void TGenome::initializeRandomGenerator(TParameters & params){
	logfile->listFlush("Initializing random generator ...");

	if(params.parameterExists("fixedSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else {
		randomGenerator=new TRandomGenerator();
	}
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
	randomGeneratorInitialized = true;
}

void TGenome::indexBamFile(std::string & filename){
	logfile->listFlush("Creating index of BAM file '" + filename + "' ...");
	BamTools::BamReader reader;
	if(!reader.Open(alignmentParser.filename))
		throw "Failed to open BAM file '" + filename + "' for indexing!";

	// create index for BAM file
	reader.CreateIndex(BamTools::BamIndex::STANDARD);

	//close BAM file
	reader.Close();
	logfile->done();
}

//-----------------------------------------------------
//Functions for theta estimation
//-----------------------------------------------------
void TGenome::openThetaOutputFile(std::ofstream & out, TThetaEstimator & estimator){
	std::string filename = outputName + "_theta_estimates.txt";
	logfile->list("Writing theta estimates to '" + filename + "'");
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";

	//write header
	out << std::setprecision(9) << "Chr\t";
	out << "start\tend";
	estimator.writeHeader(out);
	out << "\n";
}


void TGenome::estimateTheta(TParameters & params){
//	//initialize recalibration
//	alignmentParser.initializeRecalibration(params);

	//Theta estimator
	TThetaEstimator thetaEstimator(params, logfile);

	//open output
	std::ofstream out; openThetaOutputFile(out, thetaEstimator);

	//check for which segements theta is to be estimated
	if(params.parameterExists("thetaGenomeWide") || alignmentParser.considerRegions){
		if(params.parameterExists("thetaGenomeWide"))
			logfile->startIndent("Estimating theta genome-wide:");
		else logfile->startIndent("Estimating theta at specific sites:");

		//HACK!!
		bool onlyBootstrap = params.parameterExists("onlyBootstrap");
//		removed this for building
//		estimateThetaGenomeWide(thetaEstimator, out, onlyBootstrap);
		logfile->endIndent();
		if(params.parameterExists("bootstraps")){
			int numBootstraps = params.getParameterInt("bootstraps");
			//		removed this for building

	//		bootstrapTetaEstimation(numBootstraps, thetaEstimator);
		}
	} else
		estimateThetaWindows(thetaEstimator, out);

	//clean up
	out.close();
}

void TGenome::estimateThetaWindows(TThetaEstimator & thetaEstimator, std::ofstream & out){
	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters){

			logfile->startIndent("Estimating Theta:");

			//measure runtime
			struct timeval startTime, endTime;
			gettimeofday(&startTime, NULL);

			//adding sites to estimator
			logfile->listFlush("Calculating emission probabilities ...");
			thetaEstimator.clear();
			window.addSitesToThetaEstimator(alignmentParser.recalObject, thetaEstimator);
			logfile->done();

			//estimate Theta
			if(thetaEstimator.estimateTheta()){
				out << alignmentParser.chrIterator->Name << "\t" << window.start << "\t" << window.end;
				thetaEstimator.writeResultsToFile(out);
			}

			//clear theta estimator
			thetaEstimator.clear();

			//finish
			gettimeofday(&endTime, NULL);
			logfile->list("Total computation time for this window was ", endTime.tv_sec  - startTime.tv_sec, "s");
			logfile->endIndent();
		}
	}
};

void TGenome::estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, std::ofstream & out, bool onlyReadData){
	if(alignmentParser.considerRegions)
		logfile->startIndent("Estimating theta at specific sites:");

	//prepare windows
	TWindow window;

	//add sites to estimator
	logfile->startIndent("Adding sites to data structure:");
	while(alignmentParser.readDataInNextWindow(window)){
		if(window.passedFilters){
			//adding sites to estimator
			logfile->listFlush("Calculating emission probabilities ...");
			try{
				thetaEstimator.clear();
				window.addSitesToThetaEstimator(alignmentParser.recalObject, thetaEstimator);
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
		out  << "\t-\t-"; //chromosome, start, end
	else
		out  << "genome-wide\t-\t-"; //chromosome, start, end
	thetaEstimator.writeResultsToFile(out);
}

void TGenome::bootstrapTetaEstimation(int numBootstraps, TThetaEstimator & thetaEstimator){
	if(numBootstraps < 1) throw "Number of bootstraps must be > 1!";
	logfile->startIndent("Generating " + toString(numBootstraps) + " bootstrap estimates of theta:");

	//measure runtime
	struct timeval startTime, endTime;
	gettimeofday(&startTime, NULL);

	//open output file
	std::ofstream bootstrapOut;
	std::string bootstrapFilename = outputName + "_theta_bootstraps.txt";
	logfile->list("Writing theta bootstraps to '" + bootstrapFilename + "'");
	bootstrapOut.open(bootstrapFilename.c_str());
	if(!bootstrapOut) throw "Failed to open output file '" + bootstrapFilename + "'!";

	//write header
	bootstrapOut << std::setprecision(9) << "Bootstrap";
	thetaEstimator.writeHeader(bootstrapOut);
	bootstrapOut << "\n";

	//loop over bootstraps
	for(int s=0; s<numBootstraps; ++s){
		logfile->startIndent("Bootstrap " + toString(s+1) + " of " + toString(numBootstraps) + ":");

		//run bootstrap
		bootstrapOut << s+1;
		thetaEstimator.bootstrapTheta(*randomGenerator, bootstrapOut);

		logfile->endIndent();
	}

	//finish
	gettimeofday(&endTime, NULL);
	logfile->list("Total computation time for theta bootstrapping was ", round((endTime.tv_sec  - startTime.tv_sec) / 6.0)/10.0, " min");
	logfile->endIndent();
}

void TGenome::calcLikelihoodSurfaces(TParameters & params){
	//read params
	int steps = params.getParameterIntWithDefault("steps", 100);

	//prepare windows
	TWindow window;

	//Theta estimator
	TThetaEstimator estimator(logfile);

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
			window.addSitesToThetaEstimator(alignmentParser.recalObject, estimator);
			logfile->done();

			//open file
			std::ofstream out;
			filename = outputName + alignmentParser.chrIterator->Name + "_" + toString(window.start) + "_LLsurface.txt";
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
}
/*
//------------------------------------------
//Callers
//------------------------------------------
void TGenome::callMLEGenotypes(TParameters & params){
	//set all booleans in the booleans class
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = false;
	bool noAltIfHomoRef = false;

	bool gVCF = false;
	bool writeVCF = false;
	bool beagle = false, printOnlyGL = false;


	std::string indName;
	TSiteSubset* subset = NULL;

	//only call at specific sites?
	if(params.parameterExists("sites")){
		bool invariantSites = false;

		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, alignmentParser.bamHeader, alignmentParser.windowSize, logfile, invariantSites);
		else subset = new TSiteSubset(params.getParameterString("sites"), alignmentParser.windowSize, logfile, invariantSites);
		limitToSitesWithKnownAlleles = true;

	//if not, how much information should be printed?
	} else {
		if(params.parameterExists("printAll")){
			printIfNoData = true;
			logfile->list("Will print all sites, even those without data");
		}
		if(params.parameterExists("noAltIfHomoRef")){
			noAltIfHomoRef = true;
			logfile->list("Will not print alternative alleles when genotype is 0/0");
		}
		if(params.parameterExists("gVCF")){
			gVCF = true;
			if(!printIfNoData) throw "gVCF format includes calls for all sites. Use parameter \"printAll\".";
			if(noAltIfHomoRef) throw "gVCF format includes printing alternative alleles even if genotype is 0/0. Remove \"printIfNoData\".";
			if(!fastaReference) throw "Can not print VCF file without reference!";
			logfile->list("Will print output in gVCF format");
		}
	}

	if(params.parameterExists("beagle")){
		if(limitToSitesWithKnownAlleles == false) throw "Need sites file specifying major and minor alleles for beagle format!";
		beagle=true;
		logfile->list("Will print output in beagle format");
		printOnlyGL = params.parameterExists("printOnlyGL");
		if(printOnlyGL) logfile->list("Will print only genotype likelihoods");
		indName = params.getParameterStringWithDefault("indName", outputName);
	}

	if(params.parameterExists("vcf")){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;
	}

	if((writeVCF + gVCF + beagle) > 1) throw "More than one output format specified!";

	//open output file
	gz::ogzstream out;
	std::string outputFileName;

//	boolsForVCF boolaeans(writeVCF, noAltIfHomoRef, gVCF);



	if(params.parameterExists("vcf") || gVCF){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;
		std::string sName = params.getParameterStringWithDefault("sampleName", outputName);


		//open file
		if(gVCF) outputFileName = outputName + "_MLEGenotypes.gvcf.gz";
		else outputFileName = outputName + "_MLEGenotypes.vcf.gz";
		if(gVCF) logfile->list("Writing MLE genotypes in gVCF format to '" + outputFileName + "'");
		else logfile->list("Writing MLE genotypes in VCF format to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		out << "##fileformat=VCFv4.2\n";
		out << "##source=ATLAS\n";
		out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		out << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
//		if(!limitToSitesWithKnownAlleles) out << "##INFO=<ID=GG,Number=10,Type=Integer,Description=\"Phred-scaled relative likelihoods of all genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT and TT\">\n";
		out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		out << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n";
		out << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
		out << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n";
		out << "##FORMAT=<ID=GG,Number=10,Type=Integer,Description=\"Phred-scaled likelihoods for all genotypes in alphabetical order\">\n";
		out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sName << "\n";
	} else if(beagle){
		//open file
		outputFileName = outputName + "_MLEGenotypes.beagle.gz";
		logfile->list("Writing MLE genotypes to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		if(printOnlyGL) out << indName << "\t" << indName << "\t" << indName << "\n";
		else out << "#marker\talleleA\talleleB\t" << indName << "\t" << indName << "\t" << indName << "\n";

	} else {
		//open file
		outputFileName = outputName + "_MLEGenotypes.txt.gz";
		logfile->list("Writing MLE genotypes to '" + outputFileName + "'");
		out.open(outputFileName.c_str());
		if(!out) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		out << "chr\tpos\tRef";
		if(limitToSitesWithKnownAlleles) out << "\tAlt\tdepth\tL(RR)\tL(RA)\tL(AA)\tMLE\tQ\n";
		else out << "\tdepth\tL(AA)\tL(AC)\tL(AG)\tL(AT)\tL(CC)\tL(CG)\tL(CT)\tL(GG)\tL(GT)\tL(TT)\tMLE\tQ\n";
	}

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles || beagle) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			if((!limitToSitesWithKnownAlleles && !beagle) || subset->hasPositionsInWindow(windows.cur->start)){
				//read data for current window
				if(readData(windows) || printIfNoData){
					//call genotypes
					logfile->listFlush("Calling MLE genotypes ...");
					if(limitToSitesWithKnownAlleles || beagle){
						windows.cur->addReferenceBaseToSites(subset);
						windows.cur->callMLEGenotypeKnownAlleles(recalObject, subset, *randomGenerator, out, chrIterator->Name, writeVCF, noAltIfHomoRef, beagle, printOnlyGL);
					} else {
						if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
						windows.cur->callMLEGenotype(recalObject, *randomGenerator, out, chrIterator->Name, printIfNoData, fastaReference, writeVCF, gVCF, noAltIfHomoRef);
					}
					logfile->done();
				}
			}
		}
	}
	//clean up
	out.close();
}

bool TGenome::initThetaEstimatorForCallers(TParameters & params, TThetaEstimator* & thetaEstimator){
	if(params.parameterExists("theta")){
		double theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
		thetaEstimator = new TThetaEstimator(logfile);
		thetaEstimator->setTheta(theta);
		return false;
	} else {
		//prepare theta estimator
		thetaEstimator = new TThetaEstimator(params, logfile);
		return true;
	}
}

void TGenome::callBayesianGenotypes(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	TThetaEstimator* thetaEstimator;
	bool estimateTheta = initThetaEstimatorForCallers(params, thetaEstimator);
	std::ofstream outTheta;
	if(estimateTheta) openThetaOutputFile(outTheta, *thetaEstimator);

	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = true;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(alignmentParser.windowsPredefined) throw "Using site subsets is currently not implemented if windows are predefined from a BED file.";
		bool invariantSites = false;
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile, invariantSites);
		limitToSitesWithKnownAlleles = true;
	} else {
		printIfNoData = params.parameterExists("printAll");
		if(printIfNoData) logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream output;
	std::string outputFileName;
	if(params.parameterExists("vcf")){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;
		//open file
		outputFileName = outputName + "_BayesianGenotypes.vcf.gz";
		logfile->list("Writing Bayesian genotypes in VCF format to '" + outputFileName + "'");
		output.open(outputFileName.c_str());
		if(!output) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		output << "##fileformat=VCFv4.2\n";
		output << "##source=estimHet\n";
		output << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		if(!limitToSitesWithKnownAlleles) output << "##INFO=<ID=PP,Number=10,Type=Integer,Description=\"Phred-scaled posterior probabilities of all genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT and TT\">\n";
		output << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		output << "##FORMAT=<ID=GP,Number=1,Type=Integer,Description=\"Genotype posterior probabilities (phred-scaled)\">\n";
		output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		outputFileName = outputName + "_BayesianGenotypes.txt.gz";
		logfile->list("Writing Bayesian genotypes to '" + outputFileName + "'");
		output.open(outputFileName.c_str());
		if(!output) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		output << "chr\tpos";
		if(fastaReference) output << "\tRef";
		if(limitToSitesWithKnownAlleles) output << "\talt\tdepth\tP(RR|D)\tP(RA|D)\tP(AA|D)\tMAP\tQ\n";
		else output << "\tdepth\tP(AA|D)\tP(AC|D)\tP(AG|D)\tP(AT|D)\tP(CC|D)\tP(CG|D)\tP(CT|D)\tP(GG|D)\tP(GT|D)\tP(TT|D)\tMAP\tQ\n";
	}

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			//read data for current window
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				if(readData(windows) || printIfNoData){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing && estimateTheta){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN && estimateTheta){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//estimate Theta?
						if(estimateTheta){
							//adding sites to estimator
							logfile->listFlush("Calculating emission probabilities ...");
							windows.cur->addSitesToThetaEstimator(recalObject, *thetaEstimator);
							logfile->done();

							//estimate Theta
							thetaEstimator->estimateTheta();

							//write results to file
							outTheta << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t";
							thetaEstimator->writeResultsToFile(outTheta);

							//clear theta estimator
							(*thetaEstimator).clear();

						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							thetaEstimator->setBaseFreq(windows.cur->baseFreq);
						}

						//call Bayesian genotypes
						logfile->listFlush("Calling Bayesian Genotypes ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callBayesianGenotypeKnownAlleles(subset, *thetaEstimator, *randomGenerator, output, chrIterator->Name, writeVCF);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callBayesianGenotype(*thetaEstimator, *randomGenerator, output, chrIterator->Name, printIfNoData, fastaReference, writeVCF);
						}
						logfile->done();
					}
				}
			}
		}
	}

	//clean up
	if(estimateTheta){
		outTheta.close();
		delete thetaEstimator;
	}
	if(limitToSitesWithKnownAlleles) delete subset;
}

void TGenome::callAllelePresence(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	TThetaEstimator* thetaEstimator;
	bool estimateTheta = initThetaEstimatorForCallers(params, thetaEstimator);
	std::ofstream outTheta;
	if(estimateTheta) openThetaOutputFile(outTheta, *thetaEstimator);

	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = false;
	bool noAltIfHomoRef = false;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(windowsPredefined) throw "Using site subsets is currently not implemented if windows are predefined from a BED file.";
		bool invariantSites = false;
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile, invariantSites);
		limitToSitesWithKnownAlleles = true;
	} if(params.parameterExists("noAltIfHomoRef")){
		noAltIfHomoRef = true;
		logfile->list("Will not print alternative alleles when genotype is 0/0");
	} else {
		if(params.parameterExists("printAll")){
			printIfNoData = true;
			logfile->list("Will print all sites, even those without data");
		}
		if(params.parameterExists("noAltIfHomoRef")){
			noAltIfHomoRef = true;
			logfile->list("Will not print alternative alleles when genotype is 0/0");
		}
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
	gz::ogzstream outAllelePresence;
	std::string outputFileName;
	if(params.parameterExists("vcf")){
		if(!fastaReference) throw "Can not print VCF file without reference!";
		writeVCF = true;

		//open file
		outputFileName = outputName + "_AllelePresence.vcf.gz";
		logfile->list("Writing estimates of allele presence in VCF format to '" + outputFileName + "'");
		outAllelePresence.open(outputFileName.c_str());
		if(!outAllelePresence) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		outAllelePresence << "##fileformat=VCFv4.2\n";
		outAllelePresence << "##source=ATLAS\n";
		outAllelePresence << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		outAllelePresence << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
		outAllelePresence << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		outAllelePresence << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
		outAllelePresence << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
		if (limitToSitesWithKnownAlleles) outAllelePresence << "##FORMAT=<ID=PP,Number=2,Type=Integer,Description=\"Phred-scaled posterior probabilities of allele presence in the order A, C, G and T\">\n";
		else outAllelePresence << "##FORMAT=<ID=PP,Number=4,Type=Integer,Description=\"Phred-scaled posterior probabilities of allele presence in the order A, C, G and T\">\n";
		outAllelePresence << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << outputName << "\n";
	} else {
		//open file
		outputFileName = outputName + "_AllelePresence.txt.gz";
		logfile->list("Writing estimates of allele presence to '" + outputFileName + "'");
		outAllelePresence.open(outputFileName.c_str());
		if(!outAllelePresence) throw "Failed to open output file '" + outputFileName + "'!";

		//write header
		outAllelePresence << "chr\tpos";
		if(limitToSitesWithKnownAlleles) outAllelePresence << "\tRef\tAlt\tdepth\tP(Ref|D)\tP(Alt|D)\tMAP\tQ\n";
		else {
			if(fastaReference) outAllelePresence << "\tRef";
			outAllelePresence << "\tdepth\tP(A|D)\tP(C|D)\tP(G|D)\tP(T|D)\tMAP\tQ\n";
		}
	}

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				//read data for current window
				//TODO: ask Dan if following line should actually be if(readData || printIfNoData)
				if(readData(windows) || printIfNoData){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing && estimateTheta){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN && estimateTheta){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//estimate Theta?
						if(estimateTheta){
							//adding sites to estimator
							logfile->listFlush("Calculating emission probabilities ...");
							(*thetaEstimator).clear();
							windows.cur->addSitesToThetaEstimator(recalObject, *thetaEstimator);
							logfile->done();

							//estimate Theta
							thetaEstimator->estimateTheta();

							//write results to file
							outTheta << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t";
							thetaEstimator->writeResultsToFile(outTheta);

							//clear theta estimator
							(*thetaEstimator).clear();

						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							thetaEstimator->setBaseFreq(windows.cur->baseFreq);
						}

						//call allele presence
						logfile->listFlush("Calling allele presence ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callAllelePresenceKnwonAlleles(subset, *thetaEstimator, *randomGenerator, outAllelePresence, chrIterator->Name, writeVCF, noAltIfHomoRef);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callAllelePresence(*thetaEstimator, *randomGenerator, outAllelePresence, chrIterator->Name, printIfNoData, fastaReference, writeVCF, noAltIfHomoRef);
						}
						logfile->done();
					}
				}
			} else logfile->list("No positions in this window.");
		}
	}

	//clean up
	if(estimateTheta){
		outTheta.close();
		delete thetaEstimator;
	}
	if(limitToSitesWithKnownAlleles) delete subset;
}

void TGenome::randomBaseCaller(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	bool printIfNoData = false;
	if(params.parameterExists("printAll")){
		printIfNoData = true;
		logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	gz::ogzstream randomBases;
	std::string outputFileName;

	//open file
	outputFileName = outputName + "_randomBase.txt.gz";
	logfile->list("Writing random base calls to '" + outputFileName + "'");
	randomBases.open(outputFileName.c_str());
	if(!randomBases) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	randomBases << "chr\tpos\tref\tdepth\tpileup\trandom_base\n";

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows) || printIfNoData){
				//call random allele
				logfile->listFlush("Calling random base ...");
				if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
				windows.cur->callRandomBase(*randomGenerator, randomBases, chrIterator->Name, printIfNoData);
				logfile->done();

			} else logfile->list("No positions in this window.");
		}
	}
}

void TGenome::majorityBaseCaller(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	bool printIfNoData = false;
	if(params.parameterExists("printAll")){
		printIfNoData = true;
		logfile->list("Will print all sites, even those without data");
	}

	//open output: vcf or flat file?
	gz::ogzstream randomBases;
	std::string outputFileName;

	//open file
	outputFileName = outputName + "_majorityCalls.txt.gz";
	logfile->list("Writing random base calls to '" + outputFileName + "'");
	randomBases.open(outputFileName.c_str());
	if(!randomBases) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	randomBases << "chr\tpos\tref\tdepth\tpileup\trandom_base\n";

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows) || printIfNoData){
				//call random allele
				logfile->listFlush("Calling random base ...");
				if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
				windows.cur->majorityCall(*randomGenerator, randomBases, chrIterator->Name, printIfNoData);
				logfile->done();

			} else logfile->list("No positions in this window.");
		}
	}
}


//---------------------------------------------------
// I/O
//---------------------------------------------------
void TGenome::writeGLF(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//print all?
	bool printIfNoData = params.parameterExists("printAll");
	if(printIfNoData)
		logfile->list("Will print all sites, even those without data");

	//open GLF file
	std::string outputFileName = outputName + ".glf.gz";
	logfile->list("Will write data to GLF file '" + outputFileName + "'");
	TGlfWriter writer(outputFileName);

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		writer.newChromosome(chrIterator->Name, stringToLong(chrIterator->Length));
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				//write to GLF
				logfile->listFlush("Adding window to GLF file ...");
				windows.cur->calculateEmissionProbabilities(recalObject);
				windows.cur->addToGLF(writer, printIfNoData);
				logfile->done();
			}
		}
	}
	//clean up
	writer.close();
}

void TGenome::combineBeagleFiles(TParameters & params){
	std::string list = params.getParameterString("beagleList");
	std::string sites = params.getParameterString("sites");

	std::string bamName;
	std::string beagleLine;
	std::string sitesLine;
	std::vector<std::string> sitesLineVec;
	//std::vector<std::string> beagleLineVec;
	std::string marker;

	std::ifstream listFile(list.c_str());
	if(!listFile) throw "Failed to open file '" + list + "!";


//	std::string *beagleLines = new std::string[numBeagle];
	std::vector<gz::igzstream*> input;
	std::vector<gz::igzstream*>::iterator inputIt;


	//Open all files and store them in ifs
	while(listFile.good() && !listFile.eof()){
		std::getline(listFile, bamName);
		gz::igzstream* f = new gz::igzstream(bamName.c_str(), std::ios::in); // create in free store
		input.push_back(f);
		delete f;
	}


	std::ifstream sitesFile(sites.c_str());

	int lineNum = 0;
	while(sitesFile.good() && !sitesFile.eof()){
		++lineNum;
		std::getline(sitesFile, sitesLine);
		fillVectorFromStringWhiteSpaceSkipEmpty(sitesLine, sitesLineVec);
		marker = ""; marker += sitesLineVec[0] + "_" + sitesLineVec[1];
		for(inputIt=input.begin(); inputIt!=input.end(); ++inputIt){
//			std::getline(*(input.at(inputIt)), beagleLine);
//			fillVectorFromStringWhiteSpaceSkipEmpty(beagleLine, sitesLineVec);
//			if(input.at(inputIt) == marker) std::cout << marker << std::endl;
		}
	}
}

*/
void TGenome::printPileup(TParameters & params){
	//initialize recalibration
//	initializeRecalibration(params);

	//prepare windows
	TWindow window;

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
		window.printPileup(alignmentParser.recalObject, out, alignmentParser.chrNumberToName(window.chrNumber));
	}

	//clean up
	out.close();
}
/*
void TGenome::generatePSMCInput(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//read in parameters required
	double theta = params.getParameterDoubleWithDefault("theta", 0.001);
	logfile->list("Using theta = " + toString(theta));
	TThetaEstimator thetaEstimator(logfile);
	thetaEstimator.setTheta(theta);

	double confidence = params.getParameterDoubleWithDefault("confidence", 0.99);
	logfile->list("Calling heterozygosity state with confidence > " + toString(confidence));
	int blockSize = params.getParameterIntWithDefault("block", 100);
	//make sure window size is a multiple of block length!
	if(windowSize % blockSize != 0) throw "Window size is not a multiple of block size!";

	//open output file
	std::ofstream output;
	std::string outputFileName = outputName + ".psmcfa";
	logfile->list("Writing PSMC input file to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int nCharOnLine = 0;

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		if(nCharOnLine > 0) output << '\n';
		output << '>' << chrIterator->Name << '\n';
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//set base frequencies
			logfile->listFlush("Calculating emission probabilities ...");
			windows.cur->calculateEmissionProbabilities(recalObject);
			windows.cur->estimateBaseFrequencies();
			thetaEstimator.setBaseFreq(windows.cur->baseFreq);
			logfile->done();

			//create PSMC input
			logfile->listFlush("Estimating heterozygosity status ...");
			if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
			windows.cur->generatePSMCInput(thetaEstimator, blockSize, confidence, output, nCharOnLine);
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
	if(params.parameterExists("maxDepth") || params.parameterExists("minDepth")) throw "Cannot mask sites for sequencing depth while creating the mask!";

	std::ofstream output;
	std::string outputFileName = outputName + "_minDepth"+ toString(minDepthForMask) + "_maxDepth" + toString(maxDepthForMask) + "_depthMask.bed";
	logfile->list("Writing sites with depth < " + toString(minDepthForMask) + " and with depth > " + toString(maxDepthForMask) + " to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			logfile->listFlush("Writing sites to mask to output file ...");
			windows.cur->createDepthMask(minDepthForMask, maxDepthForMask, output, chrIterator->Name);
			logfile->done();
		}
	}
	output.close();
}
*/
//---------------------------------------------------
//recalibration
//---------------------------------------------------
void TGenome::estimateErrorCalibrationEM(TParameters & params){
	//create recalibration object
	std::string filename;
	if(params.parameterExists("recal")) filename = params.getParameterString("recal");
	else filename = "empty";

	bool writeTmpTables = false;
	if(params.parameterExists("writeTmpTables")){
		writeTmpTables = true;
		logfile->list("Will write intermediate estimates of EM and Newton-Raphson to file.");
	}
	TReadGroupMap readGroupMap(&alignmentParser.bamHeader, params, logfile);
	TQualityMap qualityMap;

	TRecalibrationEM recalObjectEM(&alignmentParser.bamHeader, filename, params, logfile, readGroupMap);
	if(!recalObjectEM.estimatetionRequired){
		logfile->list("No need to estimate anything. Aborting Program.");
		return;
	}

	//prepare windows
	TWindow window;

	//add sites to EM object
	logfile->startIndent("Reading data from windows:");
	while(alignmentParser.readDataInNextWindow(window)){
		//read data for current window
		if(window.passedFilters)
			window.addToRecalibrationEM(recalObjectEM, qualityMap);
		else logfile->list("No positions in this window.");
	}
	//clean up memory
	window.clear();
	logfile->endIndent();

	//run EM iterations
	recalObjectEM.runEM(outputName, writeTmpTables);
}
/*
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

void TGenome::calculateLikelihoodErrorCalibrationEM(TParameters & params){
	//create recalibration object
	std::string filename = params.getParameterString("recal");
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);
	TRecalibrationEM recalObjectEM(&bamHeader, filename, params, logfile, readGroupMap);

	//prepare windows
	TWindowPair windows;

	//add sites to EM object
	logfile->startIndent("Reading data from windows:");
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				windows.cur->addToRecalibrationEM(recalObjectEM);
			}
		}
	}
	//clean up memory
	windows.clear();
	logfile->endIndent();

	//calc likelihood surface
	//int numMarginalGridPoint = params.getParameterIntWithDefault("numGridPoints", 51);
	//recalObjectEM.calcLikelihoodSurface(outputName + "_LLsurface.txt", numMarginalGridPoint);

	logfile->list("LL = " + toString(recalObjectEM.calcLL()));
}

void TGenome::BQSR(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not run BQSR recalibration without a provided FASTA reference!";

	//prepare windows
	TWindowPair windows;

	//create BQSR object
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);
	TRecalibrationBQSR bqsr(&bamHeader, params, logfile, readGroupMap);
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
	TSiteSubset* subset = NULL;

	if(params.parameterExists("sites")){
		subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		invariantSites = true;
	}

	//loop over bam until BQSR converges
	while(!hasConverged){
		++loopNumber;
		logfile->startIndent("Running recalibration loop " + toString(loopNumber) + ":");
		//loop over all windows

		if(!bqsr.dataHasBeenStored()){
			logfile->startIndent("Reading data from BAM file:");
			while(iterateChromosome(windows)){
				if(invariantSites) subset->setChr(chrIterator->Name);
				while(iterateWindow(windows)){
					if((invariantSites && subset->hasPositionsInWindow(windows.cur->start)) || !invariantSites){
						//read data for current window
						if(readData(windows)){
							//add reference data
							windows.cur->addReferenceBaseToSites(reference, chrNumber);

							//add the base to BQSR
							if(invariantSites) windows.cur->addSitesToBQSR(bqsr, subset, logfile);
							else windows.cur->addSitesToBQSR(bqsr, logfile);
						}
					} else logfile->list("No positions in this window.");
				}
			}
			//clean up memory
			windows.clear();
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
	int maxQinPrintQualityDistribution = maxPhredInt;
//			params.getParameterIntWithDefault("maxQ", 100);

	maxQinPrintQualityDistribution += 33; //internally, quality is phred(error) + 33!
	logfile->list("Will assemble quality distribution up to a quality of " + toString(maxQinPrintQualityDistribution-33) + " (" + (char) maxQinPrintQualityDistribution + ").");

	//initialize tables: one overall, one per read group
	std::vector<TQualityTable> qualDist;
	for(int i=0; i<readGroups.size(); ++i)
		qualDist.emplace_back(maxQinPrintQualityDistribution);

	//vars
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
    long counter = 0;

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader, alignment)){
		++counter;
        if(useChromosome[alignment.chrNumber] && alignment.passedFilters && readGroups.readGroupInUse(alignment.readGroupId)){

			//update and write (only if alignment qualities could be calculated)
			alignment.addToQualityTable(qualDist[alignment.readGroupId]);

			//report
			reportProgressParsingBamFile(counter, start);
        }
	}

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

	//print per read group table
	logfile->startIndent("Writing distributions:");
	std::string outFileName;
	for(int i=0; i<readGroups.size(); ++i){
		//open output file
		outFileName = outputName + "_" + readGroups.getName(i) + "_qualityDistribution.txt";
		logfile->listFlush("Writing distribution for read group '" + readGroups.getName(i) + "' to '" + outFileName + "' ...");
		qualDist[i].write(outFileName);
		logfile->done();
	}

	//print overall table
	outFileName = outputName + "_total_qualityDistribution.txt";
	logfile->listFlush("Writing total distribution to '" + outFileName + " ...");
	TQualityTable allQualDist(maxQinPrintQualityDistribution);
	for(int i=0; i<readGroups.size(); ++i)
		allQualDist.add(qualDist[i]);
	allQualDist.write(outFileName);
	logfile->done();
	logfile->endIndent();
}
*/
void TGenome::printQualityTransformation(TParameters & params){
	//prepare alignment
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();
	int maxQ = params.getParameterIntWithDefault("maxQ", 100);

	//create table to store counts
	std::vector<TQualityTransformTable*> QTtables;
	for(int i=0; i<alignmentParser.readGroups.size() + 1; ++i){
		TQualityTransformTable* QT = new TQualityTransformTable(maxQ);
		QTtables.push_back(QT);
	}

	//prepare output
	std::ofstream out;
	std::string filename;

	//add alignments to tables
	logfile->listFlush("Adding sites to quality transformation tables ...");
	while(alignmentParser.readNextAligment(alignment)){
		if(alignmentParser.recalObjectInitialized2) alignmentParser.addSitesToQualityTransformTable(alignment, alignmentParser.recalObject, alignmentParser.recalObject2, QTtables, logfile);
		else alignmentParser.addSitesToQualityTransformTable(alignment, alignmentParser.recalObject, QTtables, logfile);
	}
	logfile->done();

	//print final tables for read groups
	for(int i=0; i<alignmentParser.readGroups.size(); ++i){
		filename = outputName + "_" + alignmentParser.readGroups.getName(i) + "_qualityTransformation.txt";
		out.open(filename.c_str());
		if(!out) throw "Failed to open output file '" + filename + "'!";
		QTtables[i]->printTable(out);
		out.close();

		//clean up vector
		delete QTtables.at(i);
	}
	//print table for total data
	filename = outputName + "_total_qualityTransformation.txt";
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	QTtables.at(QTtables.size()-1)->printTable(out);
	out.close();
	delete QTtables.at(QTtables.size()-1);
}

void TGenome::reportProgressParsingBamFile(const long & counter, const struct timeval & start){
	if(counter % 1000000 == 0){
		static struct timeval end;
		gettimeofday(&end, NULL);
		float runtime = (end.tv_sec  - start.tv_sec)/60.0;
		logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	}
}

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
	if(withPMD && !fastaReference) throw "Cannot run recalBAM withPMD without reference!";

	//other tmp variables
	long counter = 0;
	TGenotypeMap genoMap;
	TQualityMap qualMap;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	if(withPMD){
		while(alignmentParser.readNextAligment(alignment)){
			++counter;
			if(!fastaReference) throw "Cannot take PMD into account in BAM quality scores without reference!";
			alignmentParser.addReference(&reference);

			alignmentParser.recalibrateWithPMD(alignment);
			alignment.save(bamWriter, genoMap, alignmentParser.minQual, alignmentParser.maxQual, qualMap);

			reportProgressParsingBamFile(counter, start);
        }
	} else {
		while(alignmentParser.readNextAligment(alignment)){
			++counter;

			alignmentParser.recalibrate(alignment);
			alignment.save(bamWriter, genoMap, alignmentParser.maxQualForPrinting, alignmentParser.maxQualForPrinting, qualMap);

			reportProgressParsingBamFile(counter, start);
		}
	}

	//close bam writer
	bamWriter.Close();
	logfile->done();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

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

/*

void TGenome::binQualityScores(TParameters & params){
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_binnedQualityScores.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//other temp variables
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader, alignment)){
		++counter;

		//update and write (only if alignment qualities could be calculated)
		alignment.binQualityScores(qualMap);
		alignment.save(bamWriter, genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting);

		//report
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

//---------------------------------------------------
//BAM manipulation / statistics
//---------------------------------------------------
void TGenome::assessSoftClipping(TParameters & params){
	//build table ??

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//open output file
	std::string filename = outputName + "_clippingStats.txt.gz";
	gz::ogzstream out(filename.c_str());
	if(!out)
		throw "Failed to open file '" + filename + "' for writing!";
	out << "Read\tposition\tClippedLeft\tnNotClipped\tnClippedRight\n";

	//other temp variables
	std::vector<BamTools::CigarOp>::iterator it;
	int S_left, S_right, middle;
	long counter = 0;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
	gettimeofday(&start, NULL);

	//now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader, alignment)){
		alignment.assessSoftClipping(S_left, middle, S_right);

		//report
		out << alignmentParser.bamAlignment.Name << "\t" << alignmentParser.bamAlignment.Position << "\t" << S_left << "\t" << middle << "\t" << S_right << "\n";

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close output file
	out.close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
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
			readGroupId = readGroups.find(vec[0]);
			len = stringToInt(vec[1]);
			if(len < 1) throw "Max length of read group '" + vec[0] + "' is < 1!";

			//add a new readgroup for the truncated reads to the header
			readGroup = vec[0] + "_truncated";
			bamHeader.ReadGroups.Add(readGroup);
			readGroups.fill(bamHeader);
			truncatedReadGroupId = readGroups.find(readGroup);

			//copy original tags to truncated read groups
			trunc = bamHeader.ReadGroups.Begin()+truncatedReadGroupId;
			orig = bamHeader.ReadGroups.Begin()+readGroupId;
			trunc->Library = orig->Library;
			trunc->PlatformUnit = orig->PlatformUnit;
			trunc->PredictedInsertSize = orig->PredictedInsertSize;
			trunc->ProductionDate = orig->ProductionDate;
			trunc->Program = orig->Program;
			trunc->Sample = orig->Sample;
			trunc->SequencingCenter = orig->SequencingCenter;
			trunc->SequencingTechnology = orig->SequencingTechnology;

			//add to map
			singleEndRG.insert(std::pair<int, TReadGroupMaxLength>(readGroupId, TReadGroupMaxLength(len, truncatedReadGroupId, readGroup)));
		}
	}
	logfile->done();
	logfile->conclude("read " + toString(singleEndRG.size()) + " read groups to be split.");

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_splitRG.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;
 	BamTools::BamAlignment bamAlignment;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;


    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		//check if this RG needs to be parsed
		singleEndRGIT = singleEndRG.find(readGroupId);
		if(singleEndRGIT != singleEndRG.end()){
			//check length
			if(bamAlignment.Length == singleEndRGIT->second.maxLen)
				bamAlignment.EditTag("RG", "Z", singleEndRGIT->second.truncatedReadGroup);
			else if(bamAlignment.Length > singleEndRGIT->second.maxLen){
				if(allowForLarger) bamAlignment.EditTag("RG", "Z", singleEndRGIT->second.truncatedReadGroup);
				else logfile->warning("Length of read in read group '" + readGroup + "' is > max length provided! Ignoring read.");
			}
		}

		//write
		bamWriter.SaveAlignment(bamAlignment);

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//now generate bam index
	indexBamFile(filename);

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}


void TGenome::mergeReadGroups(TParameters & params){
	//read read groups and their expected lengths
	std::string filename = params.getParameterString("readGroups");
	logfile->listFlush("Reading read groups to be merged from file '" + filename + "' ...");
	std::vector< std::vector<std::string> > readGroupsToMerge;
	std::vector< std::vector<std::string> >::reverse_iterator rIt;
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "!";

	//construct new read groups in new header object
	BamTools::SamHeader newHeader(bamHeader);
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
	int* readGroupMap = new int[bamHeader.ReadGroups.Size()];
	for(int i=0; i<bamHeader.ReadGroups.Size(); ++i) readGroupMap[i] = -1;
	logfile->startIndent("The following read groups will be merged:");
	std::vector< std::vector<std::string> >::iterator mergeIt = readGroupsToMerge.begin();
	int oldId;
	for(int rg = 0; rg < newReadGroupObject.size(); ++rg, ++mergeIt){
		logfile->startIndent("New read group '" + newReadGroupObject.getName(rg) + "' will contain read groups:");
		for(std::vector<std::string>::iterator it = mergeIt->begin(); it != mergeIt->end(); ++it){
			logfile->list(*it);
			oldId = readGroups.find(*it);
			if(readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
			readGroupMap[oldId] = rg;
		}
		logfile->endIndent();
	}
	logfile->endIndent();

	//now add read groups that will not be merged
	bool printed = false;
	std::string name;
	for(int i = 0; i < readGroups.size(); ++i){
		//check if it is mapped, otherwise add
		if(readGroupMap[i] < 0){
			if(!printed){
				logfile->startIndent("The following read groups will be kept as is:");
				printed = true;
			}
			name = readGroups.getName(i);
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
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, newHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;
 	BamTools::BamAlignment bamAlignment;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		//get read group info
		bamAlignment.GetTag("RG", readGroup);
		oldId = readGroups.find(readGroup);

		//save as new RG
		bamAlignment.EditTag("RG", "Z", newReadGroupObject.getName(readGroupMap[oldId]));
		bamWriter.SaveAlignment(bamAlignment);

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

*/
void TGenome::mergePairedEndReads(TParameters & params){
	//TODO: make merge function that accounts for indels!
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_mergedReads.bam";
	BamTools::RefVector references = alignmentParser.bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, alignmentParser.bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	if(alignmentParser.hasPMD) logfile->warning("PMD is given but relevant for read merging.");

	//should we omit some reads that did not pass validateSamFile?
	bool blacklistGiven = false;
	std::map <std::string, int> readsToOmit;

	if(params.parameterExists("blacklist")){
		blacklistGiven = true;
		//open blacklist file
		std::string blacklist = params.getParameterString("blacklist");
		logfile->listFlush("Reading reads to be omitted from '" + blacklist + "...");
		std::ifstream file(blacklist.c_str());
		if(!file) throw "Failed to open file '" + blacklist + "!";

		int lineNum = 0;
		std::vector<std::string> vec;

		//fill list of reads to omit
		while(file.good() && !file.eof()){
			++lineNum;
			fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
			if(!vec.empty()) readsToOmit.insert(std::pair<std::string,int>(vec[0], 1));
		}
		logfile->write("done! Read " + toString(lineNum) + " read names");
	}


	//other temp variables
	long counter = 0;
	TGenotypeMap genoMap;
	TQualityMap qualMap;

 	//BamTools::BamAlignment bamAlignment;

	//create storage for reads until their mate was found
	std::vector< std::pair<TAlignment*, bool> > alignmentStorage;
	std::vector< std::pair<TAlignment*, bool> >::iterator it;
	std::vector< std::pair<TAlignment*, bool> >::iterator IT;

	//BamTools::BamAlignment* alignmentPointer;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);
	int curChr = -1;

	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

    //now parse through bam file and write alignments
	while(alignmentParser.readNextAligment(alignment)){
		++counter;
		std::cout << "parsed " << alignment.alignmentName << std::endl;
		if((readsToOmit.count(alignment.alignmentName) > 0)){
			continue;

//		} else if(abs(alignment.bamAlignment.InsertSize) < alignment.bamAlignment.AlignedBases.size()){
//			logfile->warning("read " + bamAlignment.Name + " was filtered out because it was longer than the insert size");
//			readsToOmit.insert(std::pair<std::string,int>(bamAlignment.Name, 1));
//			continue;
//
		} else if(!alignment.bamAlignment.IsProperPair() || (alignment.bamAlignment.IsReverseStrand() && alignment.bamAlignment.IsMateReverseStrand()) || (!alignment.bamAlignment.IsReverseStrand() && !alignment.bamAlignment.IsMateReverseStrand())){
			continue;
		}

		else {
			//if on new chromosome, empty storage
			if(curChr != alignment.chrNumber){
				for(it = alignmentStorage.begin(); it != alignmentStorage.end(); ++it){
					bamWriter.SaveAlignment((it->first)->bamAlignment);
					delete it->first;
				}
				alignmentStorage.clear();
				curChr = alignment.chrNumber;
			}
			//add alignment to storage
				if(!alignment.bamAlignment.IsReverseStrand()) {
					//mate on forward strand is always first in bam file
					alignmentStorage.push_back(std::pair<TAlignment*, bool>(new TAlignment(alignment), false));
					std::cout << "saved alignment in storage" << std::endl;
				}
				else if(alignment.bamAlignment.IsReverseStrand()){
					std::cout << "found that it is reverse strand" << std::endl;
					//find first mate -> should be in storage
					for(it=alignmentStorage.begin(); it!=alignmentStorage.end(); ++it){
						std::cout << "going through all alignments in storage:" << std::endl;
						std::cout << "it->first->alignmentName " << it->first->alignmentName << std::endl;
						if(it->first->alignmentName == alignment.alignmentName){
							//check if this read accepts mate
							if(it->second) throw "First read of '" + alignment.alignmentName + "' is not paired or has already been merged!";
							//merge reads
							std::cout << "accepted mate" << std::endl;
							TAlignment* alignmentPointer = it->first;
							if(alignment.bamAlignment.Position > alignmentPointer->position + alignmentPointer->bamAlignment.AlignedBases.size()){
								std::cout << "reads dont overlap" << std::endl;
								//reads do not overlap -> add Ns in between
								//int numN = bamAlignment.Position - (alignmentPointer->Position + alignmentPointer)->AlignedBases.size() - 1;
								int numN = abs(alignment.bamAlignment.InsertSize) - alignment.bamAlignment.AlignedBases.size() - alignmentPointer->bamAlignment.AlignedBases.size();
								std::cout << "numN " << numN << std::endl;
								for(int i=0; i<numN; ++i){
									alignmentPointer->bamAlignment.AlignedBases += 'N';
									alignmentPointer->bamAlignment.AlignedQualities += '!';
								}
								alignmentPointer->bamAlignment.AlignedBases += alignment.bamAlignment.AlignedBases;
								alignmentPointer->bamAlignment.AlignedQualities += alignment.bamAlignment.AlignedQualities;
							} else if(alignment.position < alignmentPointer->position) std::cout << "Second mate of read '" << alignment.alignmentName << "' has pos < pos of first mate!";
							else {
								std::cout << "reads overlap" << std::endl;
								//reads do overlap -> merge them
								std::string alignmentString;
								std::string quality;
								int firstOverlap = alignment.position - alignmentPointer->position;
								int lastOverlapPlusOne = -1;
								if((alignmentPointer->position + alignmentPointer->bamAlignment.AlignedBases.size()) <= (alignment.position + alignment.bamAlignment.AlignedBases.size())) lastOverlapPlusOne = alignmentPointer->bamAlignment.AlignedBases.size();
								else lastOverlapPlusOne = alignment.bamAlignment.AlignedBases.size();

								//first copy from first mate
								alignmentString = alignmentPointer->bamAlignment.AlignedBases.substr(0, firstOverlap);
								quality = alignmentPointer->bamAlignment.AlignedQualities.substr(0, firstOverlap);
								//decide which alignment has higher quality in overlap
								for(int i=firstOverlap; i<lastOverlapPlusOne; ++i){
									//decide which quality is higher
									if(alignmentPointer->bamAlignment.AlignedQualities.at(i) > alignment.bamAlignment.AlignedQualities.at(i - firstOverlap)){
										alignmentString += alignmentPointer->bamAlignment.AlignedBases.at(i);
										quality += alignmentPointer->bamAlignment.AlignedQualities.at(i);
									} else {
										alignmentString += alignment.bamAlignment.AlignedBases.at(i - firstOverlap);
										quality += alignment.bamAlignment.AlignedQualities.at(i - firstOverlap);
									}
								}

								//add rest from second
								if(alignmentPointer->position + alignmentPointer->bamAlignment.AlignedBases.size() < alignment.position + alignment.bamAlignment.AlignedBases.size()){
									alignmentString += alignment.bamAlignment.AlignedBases.substr(lastOverlapPlusOne - firstOverlap);
									quality += alignment.bamAlignment.AlignedQualities.substr(lastOverlapPlusOne - firstOverlap);
								}

								if(alignmentString.length() != abs(alignment.bamAlignment.InsertSize) + 1  && alignmentString.length() != abs(alignment.bamAlignment.InsertSize)) logfile->warning("merged alignment length of reads " + alignment.alignmentName + " is not equal to original insert size (+1)!");

								//set
								alignmentPointer->bamAlignment.AlignedBases = alignmentString;
								alignmentPointer->bamAlignment.AlignedQualities = quality;
							}

							std::cout << "got to updating part!" << std::endl;
							//update

//							alignmentPointer->QueryBases = alignmentPointer->AlignedBases;
//							alignmentPointer->Qualities = alignmentPointer->AlignedQualities;
//							alignmentPointer->Length = alignmentPointer->AlignedBases.size();
//							alignmentPointer->CigarData.clear();
//							alignmentPointer->CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MATCH_CHAR, alignmentPointer->Length));
//							alignmentPointer->SetIsFirstMate(false);
//							alignmentPointer->SetIsSecondMate(false);
//							alignmentPointer->SetIsPaired(false);
//							alignmentPointer->SetIsProperPair(false);
//							alignmentPointer->SetIsMateReverseStrand(false);
//							alignmentPointer->SetIsReverseStrand(false); //the read that comes first in BAM is always fwd strand
//							alignmentPointer->MateRefID = -1;
//							alignmentPointer->MatePosition = -1;
//							alignmentPointer->InsertSize = 0;

							alignmentPointer->setToSingleEnd();
							std::cout << "done setting to single end" << std::endl;

							it->second = true;
//								if(bamAlignment.Name == "HWI-ST558:352:C7T9BACXX:1:1201:7676:94794") std::cout << "aligment " << bamAlignment.Name << "was updated" << std::endl;

							//write if is first in vector
							if(it == alignmentStorage.begin()){
								std::cout << "we are at alignmentstorage begin" << std::endl;

								//write all that are OK
								int counter = 0;
								for(; it != alignmentStorage.end(); ++it){
									std::cout << "counter " << counter << std::endl;
									++counter;
									if(it->second){
										std::cout << "in if" << std::endl;
										alignmentPointer->save(bamWriter, genoMap, alignmentParser.minQual, alignmentParser.maxQual, qualMap);
										std::cout << "outside save function" << std::endl;
//										bamWriter.SaveAlignment(*(it->first)); //saves the alignment to the bam file
										delete it->first;
										std::cout << "was able to delete" << std::endl;
									} else {
										//first that can not be written -> erase all previous ones!
										it = alignmentStorage.erase(alignmentStorage.begin(), it);
										break;
									}
								}
								if(it == alignmentStorage.end()){
									alignmentStorage.clear();
								}
							} break;
						}
					}

					if(!alignmentStorage.empty() && it == alignmentStorage.end()) throw "One read of '" + alignment.alignmentName + "' is reverse mate, but forward one has not been read!";
				} else{
					for(it=alignmentStorage.begin(); it != alignmentStorage.end(); ++it){
						delete it->first;
					}
					alignmentStorage.clear();
					throw "One read of '" + alignment.alignmentName + "' is paired, but neither first nor second mate!";
//				} else {
//					//read is not paired: add to storage or write
//					if(alignmentStorage.empty()) bamWriter.SaveAlignment(alignment.bamAlignment);
//					else alignmentStorage.push_back(std::pair<TAlignment*, bool>(new alignment, true));
				}
		}

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

/*
void TGenome::downSampleBamFile(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);

	//read downsampling rate
	std::string prob = params.getParameterString("prob");
	//check if prob is a vector of multiple probabilities
	std::vector<double> downSampleProbVector;
	if(!stringContainsOnly(prob, "-0123456789.,")) throw "Wrong format on probability list: use floating point numbers delimited by commas (e.g. 0.1,0.2,0.5).";
	fillVectorFromString(prob, downSampleProbVector, ',');

	//read how many replicates
	int numProbs = downSampleProbVector.size();
	int times = params.getParameterIntWithDefault("times", 1);
	if(times > 1 && numProbs > 1) throw "Replicated downsampling is not implemented for more than one prob";
	if(times > 1) numProbs = times;

	//check if probs are between 0 and 1, save in array and print them
	double* downSampleProb = new double[numProbs];
	logfile->listFlush("Will accept reads with probabilities");
	bool first = true;
	int i=0;
	if(times == 1){
		for(std::vector<double>::iterator it=downSampleProbVector.begin(); it!=downSampleProbVector.end(); ++it, ++i){
			if(first) first = false;
			else logfile->flush(",");
			logfile->flush(" " + toString(*it));
			if(*it <= 0.0 || *it > 1.0) throw "All probabilities have to be between  0 and 1!";
			if(*it == 1.0) logfile->warning("Probability of 1 will result in identical file");
			if(*it == 0.0) logfile->warning("Probability of 0 will result in empty file");
			downSampleProb[i] = *it;
		}
	} else {
		for(int i=0; i<times; ++i){
			downSampleProb[i] = *(downSampleProbVector.begin());
		}
	}
	logfile->newLine();

	//open bam files for writing
	BamTools::BamWriter* bamWriter = new BamTools::BamWriter[numProbs];
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->startIndent("Writing results to the following files:");
	if(times > 1){
		for(i=0; i<numProbs; ++i){
			//construct and print filename
			filename = outputName + "_downsampled_" + toString(downSampleProb[i]) + "_" + toString(i) + ".bam";
			logfile->list(filename);
			//open file
			if(!bamWriter[i].Open(filename, bamHeader, references))	throw "Failed to open BAM file '" + filename + "'!";
		}
	} else {
		for(i=0; i<numProbs; ++i){
			//construct and print filename
			filename = outputName + "_downsampled_" + toString(downSampleProb[i]) + ".bam";
			logfile->list(filename);
			//open file
			if(!bamWriter[i].Open(filename, bamHeader, references))	throw "Failed to open BAM file '" + filename + "'!";
		}
	}
	logfile->endIndent();

	//other temp variables
	long counter = 0;
	double r;
	TGenotypeMap genoMap;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader, alignment)){
		//filters
        if(!readGroups.readGroupInUse(alignment.readGroupId)) continue;
        if(!useChromosome[alignment.chrNumber]) continue;
		++counter;

		//accept read or not?
		for(i=0; i<numProbs; ++i){
			r = randomGenerator->getRand(); //inside loop to avoid correlation when multiple probs
			if(r < downSampleProb[i])
				alignment.save(bamWriter[i], genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting);
		}

		//report
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer and clean up memory
	for(i=0; i<numProbs; ++i){
		bamWriter[i].Close();
	}
	delete[] downSampleProb;
	delete[] bamWriter;

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::downSampleReads(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//read parameters
	double fraction = params.getParameterDoubleWithDefault("fraction", 0.1);
	logfile->list("Each base has a probability of " + toString(fraction)+ " of being masked.");

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_downsampledReads" + toString(fraction) + ".bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	long counter = 0;
	TGenotypeMap genoMap;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start;
    gettimeofday(&start, NULL);

    //now parse through bam file and write alignments
	while(alignmentParser.readAlignment(bamReader, alignment)){

		alignment.downsampleAlignment(fraction, *randomGenerator);
		alignment.save(bamWriter, genoMap, alignmentParser.minQualForPrinting, alignmentParser.maxQualForPrinting);

		//report
		++counter;
		reportProgressParsingBamFile(counter, start);
	}

	//close bam writer
	bamWriter.Close();

	//report
	reportProgressParsingBamFile(counter, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::diagnoseBamFile(TParameters & params){
	//initialize alignment reading
	TAlignment alignment(maxReadLength);

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
    double totLength = 0.0;
    int chrNum = 0;
    for(chrIterator = bamHeader.Sequences.Begin(); chrIterator!=bamHeader.Sequences.End(); ++chrIterator, ++chrNum)
        if(useChromosome[chrNum]) totLength += stringToLong(chrIterator->Length);

    //prepare reporting
    logfile->startIndent("Parsing through BAM file:");
    struct timeval start;
    gettimeofday(&start, NULL);

    //other temp variables
    long counter = 0;
    int RGInd = -1;
    std::vector<double> cov;
    double totCov = 0.0;
    int numProperPairs = 0;
    long sumFragLen = 0;
    long sumSquaredFragLen = 0;

    long** MQ = new long*[readGroups.size()];
    long** RL = new long*[readGroups.size()];
    for(int i = 0; i < readGroups.size(); ++i){
    	cov.push_back(0);
    	MQ[i] = new long[100];
    	RL[i] = new long[500];
    	for(int j=0; j<100; ++j) MQ[i][j]=0;
    	for(int j=0; j<500; ++j) RL[i][j]=0;
    }

    //now parse through bam file and sum number of aligned bases
    //TODO: avoid getting properties from bamAlignment, use alignmentParser
    while(alignmentParser.readAlignment(bamReader, alignment)){
    	//filters
        if(!readGroups.readGroupInUse(alignment.readGroupId)) continue;
        if(!useChromosome[alignment.chrNumber]) continue;
        if(!alignment.passedFilters) continue;
        if(alignment.isProperPair){
        	if(!alignment.isReverseStrand){
        		++numProperPairs;
        		sumFragLen += abs(alignmentParser.bamAlignment.InsertSize);
        		sumSquaredFragLen += (alignmentParser.bamAlignment.InsertSize * alignmentParser.bamAlignment.InsertSize);
        	}
        }

        RGInd = alignment.readGroupId;
        totCov += alignmentParser.bamAlignment.AlignedBases.length();
        cov[RGInd] += alignmentParser.bamAlignment.AlignedBases.length();
        ++MQ[RGInd][alignment.mappingQuality];
        ++RL[RGInd][alignmentParser.bamAlignment.Length];

        //report
        ++counter;
        reportProgressParsingBamFile(counter, start);
    }

    //report
    reportProgressParsingBamFile(counter, start);logfile->list("Reached end of BAM file!");
    logfile->removeIndent();
    logfile->list("Approximate sequencing depth was estimated at " + toString(totCov/totLength));

    logfile->listFlush("Writing to output files ...");

    //cov
    outputDepth << "RG\tApproximate_depth";
    outputDepth << "\nallReadGroups\t" << totCov/totLength;
    for(int r=0; r<readGroups.size(); ++r){
        outputDepth << "\n" << readGroups.getName(r) << "\t" << cov[r]/totLength;
    }
    outputDepth << "\n";


    //MQ
    long tot;
    outputMQ << "RG\tMapping_quality\tCount";
    for(int i=0; i<100; ++i){
    	tot = 0;
    	for(int r=0; r<readGroups.size(); ++r) tot += MQ[r][i];
		outputMQ << "\nallReadGroups\t" << i << "\t" << tot;

    }
    for(int r=0; r<readGroups.size(); ++r){
        for(int i=0; i<100; ++i){
            outputMQ << "\n" << readGroups.getName(r) << "\t" << i << "\t" << MQ[r][i];
        }
    }
    outputMQ << "\n";


    //RL
    outputReadLen << "RG\tRead_length\tCount";
    for(int i=0; i<500; ++i){
    	tot = 0;
    	for(int r=0; r<readGroups.size(); ++r) tot += RL[r][i];
		outputReadLen << "\nallReadGroups\t" << i << "\t" << tot;
    }
    for(int r=0; r<readGroups.size(); ++r){
        for(int i=0; i<500; ++i){
            outputReadLen << "\n" << readGroups.getName(r)<< "\t" << i << "\t" << RL[r][i];
        }
    }
    outputReadLen << "\n";


    //FL
    float mean = float(sumFragLen)/float(numProperPairs);
    float var = float(sumSquaredFragLen) / float(numProperPairs) - (mean*mean);
    fragmentStats << "mean: " << mean << "\n" << "variance: " << var << "\n";

    logfile->done();

    outputDepth.close();
    outputMQ.close();
    outputReadLen.close();
    fragmentStats.close();

    for(int i = 0; i < readGroups.size(); ++i){
    	delete MQ[i];
    	delete RL[i];
    }
    delete [] MQ;
    delete [] RL;
}

void TGenome::allelicDepth(TParameters & params){
	std::ofstream output;
	std::string outputFileName = outputName + "_allelicDepth.txt";
	logfile->list("Writing allelic imbalance table to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	//int maxCov = params.getParameterIntWithDefault("maxCov", 20);
	int maxAllelicDepth = params.getParameterInt("maxAllelicDepth");
	int size = maxAllelicDepth+1; // need 0 bin
	int nCharOnLine = 0;

	//prepare array
	long**** siteCounts = new long***[size];
	for(int i=0; i<size; ++i){
		siteCounts[i] = new long**[size];
		for(int j=0; j<size; ++j){
			siteCounts[i][j] = new long*[size];
			for(int k=0; k<size; ++k){
				siteCounts[i][j][k] = new long[size];
				for(int l=0; l<size; ++l){
					siteCounts[i][j][k][l] = 0;
				}
			}
		}
	}

	//write header
	output << "A\tC\tG\tT\tCounts\tDepth" << std::endl;

	//prepare windows
	TWindowPair windows;
	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			windows.cur->countAlleles(siteCounts, maxAllelicDepth);
			logfile->listFlush("Adding imbalance values to table ...");
			logfile->write(" done!");
		}
	}

	//write to file
	for(int i=0; i<(size); ++i){
		for(int j=0; j<(size); ++j){
			for(int k=0; k<(size); ++k){
				for(int l=0; l<(size); ++l){
					output << i << "\t" << j << "\t" << k << "\t" << l << "\t" << siteCounts[i][j][k][l] << "\t" << i + j + k + l;
					output << std::endl;
				}
			}
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
	delete[] siteCounts;
}

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
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//write to file
			logfile->listFlush("Writing sequencing depth estimates to file ...");
			if(windows.cur->depth == -1.0) output << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t" << "0" << "\n";
			else output << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t" << windows.cur->depth << "\n";
			logfile->done();
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
}

void TGenome::estimateDepthPerSite(TParameters & params){
	std::ofstream output;
	std::ofstream outputNormalized;
	std::ofstream outputQuantiles;

	std::string outputFileName = outputName + "_depthPerSite.txt";
	std::string outputFileNameNormalized = outputName + "_normalizedCumulativeDepthPerSite.txt";
	std::string outputFileNameQuantiles = outputName + "_quantilesDepthPerSite.txt";


	logfile->list("Writing cumulative depth distribution to '" + outputFileName + "'");
	logfile->list("Writing normalized cumulative depth distribution to '" + outputFileNameNormalized + "'");
	logfile->list("Writing quantiles of cumulative depth distribution to '" + outputFileNameQuantiles + "'");


	output.open(outputFileName.c_str());
	outputNormalized.open(outputFileNameNormalized.c_str());
	outputQuantiles.open(outputFileNameQuantiles.c_str());


	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	if(!outputNormalized) throw "Failed to open output file '" + outputFileNameNormalized + "'!";
	if(!outputQuantiles) throw "Failed to open output file '" + outputFileNameQuantiles + "'!";


	int maxCov = params.getParameterIntWithDefault("maxDepth", 20);
	if(!maxCov) throw "No maximum depth specified!";
	int size = maxCov + 2; // need 0 bin and >maxCov bin
	int nCharOnLine = 0;
	double totalSites = 0.0;

	//prepare arrays
	long * siteDepth = new long[size];
	for(int i=0; i<size; ++i){
		siteDepth[i] = 0;
	}

	long * siteDepthNormalized = new long[size];
	for(int i=0; i<size; ++i){
		siteDepthNormalized[i] = 0;
	}

	//write headers
	output << "depth\tcounts" << std::endl;
	outputNormalized << "depth\tcounts" << std::endl;
	outputQuantiles << "percent\tquantile" << std::endl;

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			logfile->listFlush("Adding depth to table ...");
			windows.cur->calcDepthPerSite(siteDepth, maxCov);
			logfile->done();
		}
	}

	//write to files
	//read counts
	for(int i=0; i<(size-1); ++i){
		output << i << "\t" << siteDepth[i] << "\n";
		totalSites += (double) siteDepth[i];
	}
	totalSites += siteDepth[size - 1];
	output << ">" << maxCov << "\t" << siteDepth[size - 1] << std::endl;

	//normalized cumulative distribution and quantiles
	double cumul = 0.0;
	double norm = 0.0;
	float percentages[14] = {0.001, 0.005, 0.01, 0.025, 0.05, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999};
	int numberOfPercentages = 14;
	int quantiles[numberOfPercentages];
	std::fill_n(quantiles, numberOfPercentages, -1);
	for(int i=0; i<(size-1); ++i){
		cumul += (double) (siteDepth[i]);
		if(i > 0 && norm == 1.0) outputNormalized << i << "\t" << 0 << "\n";
		else {
			norm = cumul / totalSites;
			outputNormalized << i << "\t" << norm << "\n";
		}
		for(int p=0; p<numberOfPercentages; ++p){
			//if smallest quantiles = 0
			if(i == 0 && norm > percentages[p]) quantiles[p] = i;
			else if(quantiles[p] == -1 && norm >= percentages[p]){
				quantiles[p] = i;
			}
		}
		if(quantiles[numberOfPercentages] == -1 && i >= percentages[numberOfPercentages]) quantiles[numberOfPercentages] = i;
	}
	if(norm == 1.0) outputNormalized << ">" << maxCov << "\t" << 0 << "\n";
	else {
		cumul += (double) siteDepth[size - 1];
		norm = cumul / totalSites;
		outputNormalized << ">" << maxCov << "\t" << norm << std::endl;
	}

	for(int p=0; p<numberOfPercentages; ++p){
		if(quantiles[p] == -1) quantiles[p] = maxCov;
	}

	for(int i=0; i < numberOfPercentages; ++i){
		outputQuantiles << percentages[i] << "\t" << quantiles[i] << std::endl;
	}

	//clean up
	if(nCharOnLine > 0){
		output << '\n';
		outputNormalized << '\n';
	}
	output.close();
	outputNormalized.close();
	outputQuantiles.close();

	delete[] siteDepth;
	delete[] siteDepthNormalized;

}

void TGenome::writeDepthPerSite(TParameters & params){
	gz::ogzstream out;

	std::string outputFileName = outputName + "_depthPerSite.txt.gz";
	logfile->list("Writing per site depth to '" + outputFileName + "'");

	out.open(outputFileName.c_str());
	if(!out) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	out << "chr\tpos\tdepth" << std::endl;

	//prepare windows
	TWindowPair windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			logfile->listFlush("Writing depth per site ...");
			windows.cur->printDepthPerSite(out, chrIterator->Name);
			logfile->done();
		}
	}

	//clean up
	out.close();
}

//---------------------------------------------------
//PMD
//---------------------------------------------------
void TGenome::estimatePMD(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	//prepare maps
	TReadGroupMap readGroupMap(&bamHeader, params, logfile);
	TGenotypeMap genoMap;

	//prepare PMD table
	int maxLength = params.getParameterIntWithDefault("length", 50);
	logfile->list("Estimating PMD at the first " + toString(maxLength) + " positions.");
	TPMDTables pmdTables(readGroups, maxLength, readGroupMap);

	//measure progress and runtime
	struct timeval start;
	long numreadsAdded = 0;
	gettimeofday(&start, NULL);

	//iterate through BAM file
	while(alignmentParser.readAlignment(bamReader, alignment)){
        if(useChromosome[alignment.chrNumber] && alignment.passedFilters && readGroups.readGroupInUse(alignment.readGroupId))
			alignment.addToPMDTables(pmdTables, genoMap);

		//report
		++numreadsAdded;
		reportProgressParsingBamFile(numreadsAdded, start);
	}

	//report
	reportProgressParsingBamFile(numreadsAdded, start);
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();

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
	filename = outputName + "_PMD_input_Exponential.txt";
	logfile->listFlush("Estimating PMD exponential models and writing them to '" + filename + "' ...");
	int numNRIterations = params.getParameterIntWithDefault("numNRIterations", 100);
	double eps = params.getParameterDoubleWithDefault("eps", 0.001);
	pmdTables.fitExponentialModel(numNRIterations, eps, filename, logfile);
	logfile->done();
}


void TGenome::runPMDS(TParameters & params){
	//parse bam file and calculate PMDS for each read (seeSkoglund et al. 2014)
	//write new bam file with PMDS score added
	//parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)

	//initialize alignment reading
	TAlignment alignment(maxReadLength);
	alignmentParser.setParsingToTrue();

	initializeRecalibration(params);
	if(!fastaReference) throw "Cannot run PMDS without reference!";

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
	filename = outputName + "_PMDS.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//now parse through bam file and write alignments
	double PMDS;
	while(alignmentParser.readAlignment(bamReader, alignment)){
		++counter;

        if(useChromosome[alignment.chrNumber] && alignment.passedFilters && readGroups.readGroupInUse(alignment.readGroupId)){
			//recalibrate quality scores
			alignment.recalibrate(*recalObject, qualMap);

			//calc PMD
			PMDS = alignment.calculatePMDS(pi, pmdObjects);

			//update and write
			if(PMDS > minPMDS && PMDS < maxPMDS){
				alignment.updateOptionalSamField("DS", PMDS);
				alignment.save(bamWriter, genoMap, alignmentParser.minQual, alignmentParser.maxQual);
			} else ++counterF;
		} else alignment.save(bamWriter, genoMap, alignmentParser.minQual, alignmentParser.maxQual);

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

}

*/
