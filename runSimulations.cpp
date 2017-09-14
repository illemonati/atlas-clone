/*
 * runSimulations.cpp
 *
 *  Created on: May 19, 2017
 *      Author: wegmannd
 */

#include "runSimulations.h"

//TODO: this function is horribly long. Some parts shoudl either be functions, or pushed to TSimulator

void runSimulations(TParameters & params, TLog* logfile){
	//initialize random generator
	TRandomGenerator* randomGenerator;
	logfile->listFlush("Initializing random generator ...");
	if(params.parameterExists("fixedSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator=new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");

	//obtain base frequencies
	std::vector<float> baseFreq;
	std::string tmp = params.getParameterStringWithDefault("baseFreq", "0.25,0.25,0.25,0.25");
	fillVectorFromString(tmp, baseFreq, ',');
	if(baseFreq.size() != 4) throw "baseFreq vector must have size = 4!";
	logfile->list("simulating chromosomes with base frequencies A:" + toString(baseFreq[0]) + " C:" + toString(baseFreq[1])+ " G:" + toString(baseFreq[1])+ " T:" + toString(baseFreq[1]));

	//initialize simulator
	TSimulator simulator(logfile, randomGenerator, baseFreq);

	//read other parameters
	logfile->startIndent("Reading simulation parameters:");
	std::string outname = params.getParameterStringWithDefault("out", "ATLAS_simulations");
	logfile->list("Will write output files with tag '" + outname + "'.");
	int sampleSize = params.getParameterIntWithDefault("sampleSize", 1);
	logfile->list("Will simulate data from " + toString(sampleSize) + " individuals.");
	if(sampleSize < 1)
		throw "Sample size needs to be at least 1!";
	double referenceDivergence = params.getParameterDoubleWithDefault("refDiv", 0.01);
	logfile->list("Will simulate data with reference divergence = " + toString(referenceDivergence) + ".");

	//read length & depth
	float depth = params.getParameterDoubleWithDefault("depth", 10.0);
	simulator.setDepth(depth);
	int readLength = params.getParameterIntWithDefault("readLength", 100);
	simulator.setReadLength(readLength);
	logfile->list("Will simulate reads of length " + toString(readLength) + " to a depth of " + toString(depth) + ".");

	//chromosomes
	std::vector<std::string> string_vec;
	std::vector<long> chrLength;
	params.fillParameterIntoVectorWithDefault("chrLength", string_vec, ',', "1000000");
	repeatIndexes(string_vec, chrLength);
	std::vector<int> ploidy;
	if(params.parameterExists("ploidy")){
		params.fillParameterIntoVector("ploidy", string_vec, ',');
		repeatIndexes(string_vec, ploidy);
	} else {
		for(size_t i=0; i<chrLength.size(); ++i)
			ploidy.push_back(2);
	}
	if(ploidy.size() != chrLength.size())
		throw "List of chromosome lengths and ploidies differ in length!";
	std::vector<bool> haploid;
	for(std::vector<int>::iterator it=ploidy.begin(); it!=ploidy.end(); ++it){
		if(*it == 1) haploid.push_back(true);
		else if(*it == 2) haploid.push_back(false);
		else throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!";
	}

	if(chrLength.size() < 1)
		throw "Issue understanding length of chromosomes!";
	if(chrLength.size() == 1){
		int numChr = params.getParameterIntWithDefault("numChr", 1);
		std::string text = "Will simulate " + toString(numChr) ;
		if(haploid[0]) text += " haploid";
		else text += " diploid";
		text += " chromosome(s) of length " + toString(chrLength[0]) + " each.";
		logfile->list(text);
		simulator.initializeChromosomes(numChr, chrLength[0], haploid[0]);
	} else {
		logfile->startIndent("Will simulate " + toString(chrLength.size()) + " chromosome(s) of the following length:");
		std::vector<bool>::iterator hIt=haploid.begin();
		std::string text;
		for(std::vector<long>::iterator it=chrLength.begin(); it!=chrLength.end(); ++it, ++hIt){
			text = toString(*it) + " (";
			if(*hIt) text += "haploid)";
			else text += "diploid)";
			logfile->list(text);
		}
		simulator.initializeChromosomes(chrLength, haploid);
		logfile->endIndent();
	}

	//quality distribution
	double meanQual = params.getParameterDoubleWithDefault("meanQual", 30);
	double sdQual = params.getParameterDoubleWithDefault("sdQual", 10);
	logfile->list("Will simulate normal distributed quality scores with mean = " + toString(meanQual) + " and sd = " + toString(sdQual));
	simulator.setQualityDistribution(meanQual, sdQual);
	int maxQual = params.getParameterDoubleWithDefault("maxQual", 500);
	simulator.setMaxQual(maxQual);
	logfile->list("Will cap qualities at " + toString(maxQual));

	//quality transformation
	if(params.parameterExists("qualTransform")){
		params.fillParameterIntoVector("qualTransform", string_vec, ',');
		std::vector<double> beta;
		repeatIndexes(string_vec, beta);
		if(beta.size() != 24)
			throw "Wrong number of beta values for quality transformation (" + toString(beta.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";
		std::string s = concatenateString(beta, ",");
		logfile->list("Will transform qualities with beta = {" + s + "}");
		simulator.setQualityTransformation(beta);
	} else if(params.parameterExists("recal")){
		std::string filename = params.getParameterString("recal");
		logfile->listFlush("Reading recalibration parameters from '" + filename + "' ...");
		std::ifstream file(filename.c_str());
		if(!file) throw "Failed to open file '" + filename + "' for reading!";

		//tmp variables for reading
		std::string tmp;
		std::vector<std::string> vec;
		std::vector<double> beta;

		//skip header
		std::getline(file, tmp);

		//parse file to read details for each read group
		std::getline(file, tmp);
		fillVectorFromString(tmp, vec, "\t");
		logfile->done();

		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 25) throw "Found " + toString(vec.size()) + " instead of 25 columns in '" + filename + "' on line 2!";
			for(int i=1; i < 25; ++i) beta.push_back(stringToFloat(vec[i]));
			logfile->done();
			if(beta.size() != 24)
				throw "Wrong number of beta values for quality transformation (" + toString(beta.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";
			std::string s = concatenateString(beta, ",");
			logfile->list("Will transform qualities with beta = {" + s + "}");
			simulator.setQualityTransformation(beta);
		}
	}

	//initialize PMD
	TPMD pmdObject;
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		if(params.parameterExists("pmd")){
			std::string pmdString = params.getParameterString("pmd");
			logfile->list("Initializing PMD for both C->T and G->A with function '" + pmdString +"'.");
			pmdObject.initializeFunction(pmdString, pmdGA, logfile);
			pmdObject.initializeFunction(pmdString, pmdCT, logfile);
			logfile->conclude(pmdObject.getFunctionString(pmdCT));
			if(params.parameterExists("pmdCT")) logfile->warning("Ignoring argument 'pmdCT'!");
			if(params.parameterExists("pmdGA")) logfile->warning("Ignoring argument 'pmdGA'!");
		} else {
			//first C->T
			if(!params.parameterExists("pmdCT")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdCT' has to be provided!";
			std::string pmdStringCT = params.getParameterString("pmdCT");
			logfile->list("Initializing post mortem C->T damage with function '" + pmdStringCT +"'.");
			pmdObject.initializeFunction(pmdStringCT, pmdCT, logfile);
			logfile->conclude(pmdObject.getFunctionString(pmdCT));

			//second G->A
			if(!params.parameterExists("pmdGA")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdGA' has to be provided!";
			std::string pmdStringGA = params.getParameterString("pmdGA");
			logfile->list("Initializing post mortem G->A damage with function '" + pmdStringGA +"'.");
			pmdObject.initializeFunction(pmdStringGA, pmdGA, logfile);
			logfile->conclude(pmdObject.getFunctionString(pmdGA));
		}
		simulator.setPMD(&pmdObject);
	}
	logfile->endIndent();

	//simulate differently depending on number of individuals
	if(sampleSize == 1){
		logfile->startIndent("Simulating a single individual:");

		std::vector<std::string> tmp;
		params.fillParameterIntoVectorWithDefault("theta", tmp, ',', "0.001");
		std::vector<double> thetas;
		repeatIndexes(tmp, thetas);
		if(thetas.size() == 1){
			logfile->list("Will simulate data with theta = " + toString(thetas[0]) + ".");
			for(int i=1; i<chrLength.size(); ++i)
				thetas.push_back(thetas[0]);
		} else {
			if(thetas.size() != chrLength.size())
				throw "Number of theta values does not match number of chromosomes!";
			logfile->list("Will simulate data with chromosome specific thetas " + concatenateString(thetas, ", "));
		}
		simulator.simulateSingleIndividual(thetas, referenceDivergence, outname);
	} else if(sampleSize == 2 && params.parameterExists("phi")){
		logfile->startIndent("Simulating two individuals with predefined genetic distance:");
		//simulate according to genetic distance
		//parse distance
		std::vector<std::string> tmp;
		std::vector<double> phis;
		params.fillParameterIntoVector("phi", tmp, ',');
		repeatIndexes(tmp, phis);
		//params.fillParameterIntoVector("phi", phis, ',');
		if(phis.size() != 9)
			throw "Wrong number of phi! Required are nine values for 00/00, 00/01, 01/00, 00/11, 01/01, 01/02, 00/12, 01/22, 01/23";
		//normalize phis
		double sum = 0.0;
		for(std::vector<double>::iterator it=phis.begin(); it!=phis.end(); ++it)
			sum += *it;
		if(sum != 1.0){
			logfile->list("Normalizing phi to sum to one (currently summing to " + toString(sum) + ").");
			for(std::vector<double>::iterator it=phis.begin(); it!=phis.end(); ++it)
				*it /= sum;
		}

		logfile->list("Used phi are: " + concatenateString(phis, ", "));

		//now run simulations
		simulator.simulateIndividualPair(phis, referenceDivergence, outname);
	} else {
		logfile->startIndent("Simulating multiple individuals based on an SFS:");
		//prepare SFS
		//TODO: think about ploidy!
		logfile->startIndent("Preparing SFS:");
		std::vector<SFS*> sfs;
		if(params.parameterExists("sfs")){
			std::vector<std::string> tmp;
			std::vector<std::string> sfsFileNames;
			params.fillParameterIntoVector("sfs", tmp, ',');
			repeatIndexes(tmp, sfsFileNames);
			if(sfsFileNames.size() != chrLength.size())
				throw "List of chromosome lengths and SFS files differ in length!";

			//now read those files
			std::vector<int>::iterator ploidyIt=ploidy.begin();
			for(std::vector<std::string>::iterator it=sfsFileNames.begin(); it!=sfsFileNames.end(); ++it, ++ploidyIt){
				logfile->listFlush("Reading the sfs from file '" + *it + "' ...");
				sfs.push_back(new SFS(*it));
				logfile->done();
				if((*sfs.rbegin())->numChromosomes != *ploidyIt * sampleSize)
					throw "SFS does not match sample size! It contains data for " + toString((*sfs.rbegin())->numChromosomes) + " instead of " + toString(2*sampleSize) + " chromosomes.";
			}
		} else if(params.parameterExists("theta")){
			//parse theta from command line
			std::vector<std::string> tmp;
			params.fillParameterIntoVector("theta", tmp, ',');
			std::vector<double> thetas;
			repeatIndexes(tmp, thetas);

			if(thetas.size() != chrLength.size())
				throw "Number of theta values does not match number of chromosomes!";
			logfile->list("Will simulate data with chromosome specific thetas " + concatenateString(thetas, ", "));

			//generate SFS for each chromosome
			logfile->listFlush("Generating SFS from provided theta values ...");
			int chr = 1;
			std::string filename;
			logfile->listFlush("Writing true SFS to '" + filename + "' ...");
			std::vector<double>::iterator thetaIt = thetas.begin();
			for(std::vector<int>::iterator it=ploidy.begin(); it!=ploidy.end(); ++it, ++thetaIt, ++chr){
				sfs.push_back(new SFS(*it * sampleSize, (float) *thetaIt));

				//save true SFS
				filename = outname + "_trueSFS_chr" + toString(chr) + ".txt";
				(*sfs.rbegin())->writeToFile(filename);
			}
			logfile->done();
		} else throw "Specifying either sfs or theta is required to initialize the SFS!";

		//run simulations
		simulator.simulatePopulationFromSFS(sfs, sampleSize, referenceDivergence, outname);

		//deleting SFS
		for(std::vector<SFS*>::iterator it=sfs.begin(); it!=sfs.end(); ++it)
			delete *it;

		logfile->endIndent();
	}
	logfile->endIndent();

	//clean up
	delete randomGenerator;
}
