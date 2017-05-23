/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */


#include "TGenome.h"
#include "TBase.h"
#include <typeinfo>
#include <map>


//-------------------------------------------------------
//TGenome
//-------------------------------------------------------
TGenome::TGenome(TLog* Logfile, TParameters & params){
	logfile = Logfile;
	initializeRandomGenerator(params);

	//open BAM file
	filename = params.getParameterString("bam");
	logfile->list("Reading data from BAM file '" + filename + "'.");
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";
	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + filename + "'!";

	maxReadLength = params.getParameterIntWithDefault("maxReadLength", 1000);
	logfile->list("Will only consider reads up to " + toString(maxReadLength) + " bp.");

	//read window parameter
	if(!params.parameterExists("window") && params.parameterExists("windows")) logfile->warning("Argument 'windows' specified, but unknown. Did you mean 'window'?");
	std::string tmp = params.getParameterStringWithDefault("window", "1000000");
	//check if it is a number
	if(stringContainsOnly(tmp, "1234567890.Ee-+")){
		windowsPredefined = false;
		windowSize = stringToInt(tmp);
		logfile->list("Setting window size to " + toString(windowSize));
		if(windowSize < maxReadLength) throw "Window size " + tmp + " out of range! Windows must be at least as large as the max read length (" + toString(maxReadLength) + " bp)!";
	} else {
		windowsPredefined = true;
		logfile->listFlush("Limiting analysis to windows defined in '" + tmp + "'...");
		windows = new TBed(tmp);
		logfile->write(" done!");
		logfile->conclude("read " + toString(windows->size()) + " on " + toString(windows->getNumChromosomes()) + " chromosomes");
	}
	numWindowsOnChr = 0;

	//outputname
	outputName = params.getParameterStringWithDefault("out", "");
	if(outputName == ""){
		//guess from filename
		outputName = filename;
		outputName = extractBeforeLast(outputName, ".");
	}
	logfile->list("Writing output files with prefix '" + outputName + "'.");

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;
	curStart = -1;
	curEnd = -1;
	oldPos = -1;
	oldAlignementMustBeConsidered = false;

	//read header
	bamHeader = bamReader.GetHeader();
	readGroups.fill(bamHeader);
	chrIterator = bamHeader.Sequences.End();

	//open FASTA reference
	if(params.parameterExists("fasta")){
		std::string fastaFile = params.getParameterString("fasta");
		std::string fastaIndex = fastaFile + ".fai";
		logfile->list("Reading reference sequence from '" + fastaFile + "'");
		if(!reference.Open(fastaFile, fastaIndex)) throw "Failed to open FASTA file '" + fastaFile + "'! Is index file present?";
		fastaReference = true;
	} else fastaReference = false;

	//initialize post mortem damage
	hasPMD = false;
	initializePostMortemDamage(params);
	doRecalibration = false;
	recalObjectInitialized = false;

	//check if we mask sites
	if(params.parameterExists("mask")){
		if(windowsPredefined) throw "Masking is currently not implemented if windows are predefined from a BED file.";
		if(params.parameterExists("sites")) throw "Masking is currently not implemented if variant positions are also specified with \"sites\"";
		if(params.parameterExists("regions")) throw "Cannot use mask and regions at the same time";
		doMasking = true;
		std::string maskFile = params.getParameterString("mask");
		logfile->startIndent("Will mask all sites listed in BED file '" + maskFile + "':");
		logfile->listFlush("Reading file ...");
		mask = new TBedReader(maskFile, windowSize, bamHeader.Sequences, logfile);
		logfile->write(" done!");
		logfile->endIndent();
		//mask->print();
	} else doMasking = false;

	if(params.parameterExists("maskCpG")){
		if(!fastaReference) throw "Cannot mask CpG sites without reference!";
		doCpGMasking = true;
		std::string maskFile = params.getParameterString("maskCpG");
		logfile->list("Will mask all CpG sites");
	} else doCpGMasking = false;

	if(params.parameterExists("regions")){
		if(windowsPredefined) throw "Regions is currently not implemented if windows are predefined from a BED file.";
		if(params.parameterExists("sites")) throw "Regions is currently not implemented if variant positions are also specified with \"sites\"";
		considerRegions = true;
		std::string regionsFile = params.getParameterString("regions");
		logfile->startIndent("Will limit analysis to all regions listed in BED file '" + regionsFile + "':");
		logfile->listFlush("Reading file ...");
		mask = new TBedReader(regionsFile, windowSize, bamHeader.Sequences, logfile);
		logfile->write(" done!");
		logfile->endIndent();
	} else considerRegions = false;

	//filters
	if(params.parameterExists("minCoverage") || params.parameterExists("maxCoverage")){
		applyCoverageFilter = true;
		minCoverage = params.getParameterIntWithDefault("minCoverage", 0);
		maxCoverage = params.getParameterIntWithDefault("maxCoverage", 1000000);
		logfile->list("Will filter out sites with coverage < " + toString(minCoverage) + " or > " + toString(maxCoverage));
	} else {
		applyCoverageFilter = false;
		minCoverage = 0;
		maxCoverage = 1000000;
	}
	if(params.parameterExists("minQual") || params.parameterExists("maxQual")){
		applyQualityFilter = true;
		minQuality = params.getParameterInt("minQual") + 33; //bamAlignment.Qualities is in ascii
		maxQuality = params.getParameterInt("maxQual") + 33;
		logfile->list("Will filter out bases with quality <= " + toString(minQuality-33) + " or >= " + toString(maxQuality-33));
	} else {
		applyQualityFilter = false;
		minQuality = 34; //filter out quality 0 by default
		maxQuality = 1000000;
	}

	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 1.0);
	if(maxMissing > 1.0) throw "maxMissing must be smaller or equal to 1.0!";

	maxRefN = params.getParameterDoubleWithDefault("maxRefN", 1.0);
	if(maxRefN > 1.0) throw "maxRefN must be smaller or equal to 1.0!";
	if(maxRefN < 1.0 && fastaReference == false) throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided.";


	//limit chrs and / or windows
	useChromosome = new bool[bamHeader.Sequences.Size()];
	if(params.parameterExists("chr")){
		logfile->startIndent("Will limit analysis to the following chromosomes:");

		//set all chromosomes to false
		for(int i=0; i<bamHeader.Sequences.Size(); ++i)
				useChromosome[i] = false;

		//parse chromosome names
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("chr"), vec, ',');
		int num;
		for(std::vector<std::string>::iterator it=vec.begin(); it!=vec.end(); ++it){
			//find chromosome
			num = 0;
			for(chrIterator = bamHeader.Sequences.Begin(); chrIterator != bamHeader.Sequences.End(); ++chrIterator, ++num){
				if(chrIterator->Name == *it){
					useChromosome[num] = true;
					logfile->list(*it);
					break;
				}
			}
			if(chrIterator == bamHeader.Sequences.End()) throw "Chromosome '" + *it + "' is not present in the bam header!";
		}
		chrIterator = bamHeader.Sequences.End();
		limitChr = 1000000;
		logfile->endIndent();
	} else {
		limitChr = params.getParameterIntWithDefault("limitChr", 1000000);
		if(params.parameterExists("limitChr")) logfile->list("Will limit analysis to the first " + toString(limitChr) + " chromosomes.");
		for(int i=0; i<bamHeader.Sequences.Size(); ++i)
			useChromosome[i] = true;
	}
	limitWindows = params.getParameterLongWithDefault("limitWindows", 1000000000);
	if(params.parameterExists("limitWindows")) logfile->list("Will limit analysis to the first " + toString(limitWindows) + " windows per chromosome.");

	//limit readGroups
	if(params.parameterExists("readGroup")){
		limitReadGroups = true;
		fillVectorFromString(params.getParameterString("readGroup"), readGroupsInUse, ',');
		logfile->startIndent("Will limit analysis to the following read groups:");
		for(int i=0; i < readGroups.numGroups; i++){
			if(std::find(readGroupsInUse.begin(), readGroupsInUse.end(), readGroups.getName(i)) != readGroupsInUse.end()){
				readGroups.inUse[i] = true;
				logfile->list(readGroups.getName(i));
			} else readGroups.inUse[i] = false;
		}
		logfile->endIndent();
	} else limitReadGroups = false;

};

void TGenome::jumpToEnd(){
	chrIterator = bamHeader.Sequences.End();
	chrNumber = -1;
}

void TGenome::restartChromosome(TWindowPair & windowPair){
	chrIterator = bamHeader.Sequences.Begin();
	chrNumber = 0;

	moveChromosome(windowPair);
}

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

	//do we use this chromosome? if not, move on!
	while(chrIterator != bamHeader.Sequences.End() && !useChromosome[chrNumber]){
		++chrIterator;
		++chrNumber;
	}

	//did we reach end?
	if(chrIterator == bamHeader.Sequences.End() || chrNumber >= limitChr){
		curEnd = 0;
		chrIterator = bamHeader.Sequences.End();
		return false;
	}

	moveChromosome(windowPair);
	return true;
}

void TGenome::moveChromosome(TWindowPair & windowPair){
	//jump reader
	bamReader.Jump(chrNumber, 0);
	oldAlignementMustBeConsidered = false;
	chrLength = stringToLong(chrIterator->Length);

	//restart windows
	curStart = 0;
	curEnd = 0;
	oldPos = -1;
	windowNumber = 0;
	if(windowsPredefined){
		windows->setChr(chrIterator->Name);
		numWindowsOnChr = windows->getNumWindowsOnCurChr();
		int nextEnd = windows->curWindowEnd();
		if(nextEnd > chrLength) nextEnd = chrLength + 1;
		else windowPair.nextPointer->move(windows->curWindowStart(), nextEnd);
	} else {
		numWindowsOnChr = ceil(chrLength / (double) windowSize);
		int nextEnd = windowSize;
		if(nextEnd > chrLength) nextEnd = chrLength + 1;
		windowPair.nextPointer->move(0, nextEnd);
	}

	//advance mask
	if(doMasking || considerRegions) mask->setChr(chrIterator->Name);

	//write progress
	logfile->startNumbering("Parsing chromosome '" + chrIterator->Name + "':");
}

bool TGenome::iterateWindow(TWindowPair & windowPair){
	if(curEnd > 0) logfile->endIndent();

	//swap window pairs
	windowPair.swap();

	//move to next region
	curStart = windowPair.curPointer->start;
	curEnd = windowPair.curPointer->end;
	if(curStart >= chrLength || windowNumber >= limitWindows) return false;

	//move next
	if(windowsPredefined){
		if(numWindowsOnChr < 1){
			logfile->conclude("No windows on this chromosome.");
			return false;
		}
		//jump reader if large gap to previous window
		if(windowPair.curPointer->start - windowPair.nextPointer->end > maxReadLength)
			bamReader.Jump(chrNumber, curStart);

		//now move coordinates of next window
		if(windows->nextWindow()){
			int nextEnd = windows->curWindowEnd();
			if(nextEnd > chrLength) nextEnd = chrLength + 1;
			windowPair.nextPointer->move(windows->curWindowStart(), nextEnd);
		} else {
			windowPair.nextPointer->move(chrLength + 1, chrLength + 2);
		}
	} else {
		long nextEnd = curEnd + windowSize;
		if(nextEnd > chrLength) nextEnd = chrLength + 1;
		windowPair.nextPointer->move(curEnd, nextEnd);
	}

	++windowNumber;

	//report
	logfile->number("Window [" + toString(curStart) + ", " + toString(curEnd) + "] of " + toString(numWindowsOnChr) + " on '" + chrIterator->Name + "':");
	logfile->addIndent();

	return true;
};

bool TGenome::addAlignementToWindows(BamTools::BamAlignment & alignement, TWindowPair & windowPair){
	//std::cout << "REF ID = " << bamAlignement.RefID << "\tpos = " << bamAlignement.Position << std::endl;
	//only take those reads that pass QC
	if(!alignement.IsFailedQC() && !alignement.IsDuplicate() && alignement.Position >= curStart){
		//check if bam file is sorted
		if(alignement.Position < oldPos){
			throw "BAM file must be sorted by position!";
		}
		oldPos = alignement.Position;

		//check if still within current window and add to window
		if(alignement.Position >= curEnd) return false;
		else {
			if(windowPair.curPointer->addFromRead(alignement, pmdObjects, &readGroups, minQuality, maxQuality)){
				//add also to next window in case reads overhangs current window -> function returns true
				windowPair.nextPointer->addFromRead(alignement, pmdObjects, &readGroups, minQuality, maxQuality);
			}
		}
	}
	return true; //continue
}

bool TGenome::readData(TWindowPair & windowPair){
	logfile->listFlush("Reading data ...");

	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//parse through reads
	if(oldAlignementMustBeConsidered){
		oldAlignementMustBeConsidered = false;
		if(!addAlignementToWindows(bamAlignment, windowPair)){
			//next read is for a later window
			oldAlignementMustBeConsidered = true;
			gettimeofday(&end, NULL);
			logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");
			logfile->conclude("No data in this window.");
			return false; //still only in next window
		}
	}

	while(bamReader.GetNextAlignment(bamAlignment) && bamAlignment.RefID==chrNumber){

		if(readGroups.readGroupInUse(bamAlignment)){

			//filter out unmapped reads and those that did not pass QC
			if(bamAlignment.IsMapped() && !bamAlignment.IsFailedQC()){

				if(bamAlignment.IsPrimaryAlignment()){

					//check if read is paired and reject reads with pairs on different chromosomes (maybe too harsh?)
					//if(!bamAlignment.IsPaired() || bamAlignment.MateRefID == bamAlignment.RefID){

						//check if insert size is shorter than read, this means we are reading the adaptor sequence
						if(!bamAlignment.IsPaired() || abs(bamAlignment.InsertSize) >= bamAlignment.AlignedBases.length()){

							if(bamAlignment.AlignedBases.size() > maxReadLength)
								throw "Alignment '" +  bamAlignment.Name + "' is long than the max read length! Please change max read length to parse this data.";

							if(!addAlignementToWindows(bamAlignment, windowPair)){

								//read is beyond window and should be reconsidered
								oldAlignementMustBeConsidered = true;
								break;
							}
						}
						else logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
					//}
				}
			}
		}
	}
	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");

	if(windowPair.curPointer->numReadsInWindow > 0){
		//apply masks and filters
		if(doMasking){
			logfile->listFlush("Masking sites ...");
			windowPair.curPointer->applyMask(mask, considerRegions);
			logfile->write(" done!");
		} else if(considerRegions){
			logfile->listFlush("Masking sites outside regions ...");
			windowPair.curPointer->applyMask(mask, considerRegions);
			logfile->write(" done!");
		} else if(doCpGMasking){
			logfile->listFlush("Masking CpG sites ...");
			windowPair.curPointer->maskCpG(reference, chrNumber);
			logfile->write(" done!");
		} if(applyCoverageFilter){
			windowPair.curPointer->applyCoverageFilter(minCoverage, maxCoverage);
		} if(maxRefN < 1.0 && fastaReference == true){
			windowPair.curPointer->addReferenceBaseToSites(reference, chrNumber);
			windowPair.curPointer->calcFracN();
		}

		//calc coverage
		windowPair.curPointer->calcCoverage();

		//report
		logfile->conclude("read data from " + toString(windowPair.curPointer->numReadsInWindow) + " reads.");
		logfile->conclude("coverage is " + toString(windowPair.curPointer->coverage));
		logfile->conclude(toString(windowPair.curPointer->fractionsitesCoverageAtLeastTwo * 100) + "% of all sites are covered at least twice");
		logfile->conclude(toString(windowPair.curPointer->fractionSitesNoData * 100) + "% of all sites have no data");
		if(maxRefN < 1.0 && fastaReference == true) logfile->conclude(toString(windowPair.curPointer->fractionRefIsN * 100) + "% of all reference bases are 'N'");
		return true;
	} else {
		logfile->conclude("No data in this window.");
		return false;
	}
};

void TGenome::initializePostMortemDamage(TParameters & params){
	logfile->startIndent("Initializing Post Mortem Damage (PMD):");
	//create an array of TPMD objects for each read group
	pmdObjects = new TPMD[readGroups.numGroups];

	//now fill them!
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		//all read groups have the same pmd
		logfile->list("Initializing one PMD function for all read groups.");
		if(params.parameterExists("pmd")){
			std::string pmdString = params.getParameterString("pmd");
			logfile->list("Initializing PMD for both C->T and G->A with function '" + pmdString +"'.");
			for(int i=0; i<readGroups.numGroups; ++i){
				pmdObjects[i].initializeFunction(pmdString, pmdGA, logfile);
				pmdObjects[i].initializeFunction(pmdString, pmdCT, logfile);
			}
			logfile->conclude(pmdObjects[0].getFunctionString(pmdCT));
			if(params.parameterExists("pmdCT")) logfile->warning("Ignoring argument 'pmdCT'!");
			if(params.parameterExists("pmdGA")) logfile->warning("Ignoring argument 'pmdGA'!");
		} else {
			if(!params.parameterExists("pmdCT")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdCT' has to be provided!";
			std::string pmdStringCT = params.getParameterString("pmdCT");
			logfile->list("Initializing post mortem C->T damage with function '" + pmdStringCT +"'.");
			for(int i=0; i<readGroups.numGroups; ++i){
				pmdObjects[i].initializeFunction(pmdStringCT, pmdCT, logfile);
			}
			logfile->conclude(pmdObjects[0].getFunctionString(pmdCT));

			//second G->A
			if(!params.parameterExists("pmdGA")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdGA' has to be provided!";
			std::string pmdStringGA = params.getParameterString("pmdGA");
			logfile->list("Initializing post mortem G->A damage with function '" + pmdStringGA +"'.");
			for(int i=0; i<readGroups.numGroups; ++i){
				pmdObjects[i].initializeFunction(pmdStringGA, pmdGA, logfile);
			}
			logfile->conclude(pmdObjects[0].getFunctionString(pmdGA));
		}
		hasPMD = true;
	} else if(params.parameterExists("pmdFile")){
		//read from file for each read group
		std::string filename = params.getParameterString("pmdFile");
		logfile->list("Reading PMD from  from file '" + filename + "'.");
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
				//get read group
				if(readGroups.readGroupExists(vec[0])){ //ignore if it does not exist
					readGroupId = readGroups.find(vec[0]);
					//initialize functions
					pmdObjects[readGroupId].initializeFunction(vec[1], pmdCT, logfile);
				//	logfile->conclude("For read group '" + vec[0] + "', C->T: " + pmdObjects[readGroupId].getFunctionString(pmdCT));
					pmdObjects[readGroupId].initializeFunction(vec[2], pmdGA, logfile);
				//	logfile->conclude("For read group '" + vec[0] + "', G->A: " + pmdObjects[readGroupId].getFunctionString(pmdGA));
				}
			}
		}

		//close file
		file.close();

		//test if we have a function for all read groups
		for(int i=0; i<readGroups.numGroups; ++i){
			if(!pmdObjects[i].functionInitialized(pmdCT)) throw "PMD C->T for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
			if(!pmdObjects[i].functionInitialized(pmdGA)) throw "PMD G->A for read group '" + readGroups.getName(i) + "' is missing in file '" + filename + "'!";
		}
		hasPMD = true;
	} else {
		//no post mortem damage
		logfile->list("Assuming there is no PMD in the data.");
		std::string pmdString = "none";
		for(int i=0; i<readGroups.numGroups; ++i){
			pmdObjects[i].initializeFunction(pmdString, pmdGA, logfile);
			pmdObjects[i].initializeFunction(pmdString, pmdCT, logfile);
		}
	}
	logfile->endIndent();
}

void TGenome::initializeRecalibration(TParameters & params){
	if(params.parameterExists("recal")){
		std::string filename = params.getParameterString("recal");
		recalObject = new TRecalibrationEM(&bamHeader, filename, params, logfile);
		doRecalibration = true;
	} else if(params.parameterExists("BQSRQuality")){
		recalObject = new TRecalibrationBQSR(&bamHeader, params, logfile);
		doRecalibration = true;
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
		doRecalibration = false;
		recalObject = new TRecalibration();
	}
	recalObjectInitialized = true;

	//check if estimation is required, in which case throw an error!
	if(recalObject->requiresEstimation()) throw "Can not use provided recalibration: estimation is required!";
}

void TGenome::openThetaOutputFile(std::ofstream & out){
	std::string filename = outputName + "_theta_estimates.txt";
	logfile->list("Writing theta estimates to '" + filename + "'");
	out.open(filename.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	out << std::setprecision(9) << "Chr\t";
	out << "start\tend\tcoverage\tmissing\ttwoOrMore\tpi(A)\tpi(C)\tpi(G)\tpi(T)\ttheta_MLE\ttheta_C95_l\ttheta_C95_u\tLL";
	out << "\n";
}

void TGenome::initializeRandomGenerator(TParameters & params){
	logfile->listFlush("Initializing random generator ...");

	if(params.parameterExists("fixedSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator=new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator=new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
	randomGeneratorInitialized = true;
}



void TGenome::estimateTheta(TParameters & params){
	//read parameters
	bool thetaGenomeWide = false;
	if(params.parameterExists("thetaGenomeWide")){
		if(considerRegions) throw "thetaGenomeWide can presently not be used in combination with regions!";
		logfile->list("estimating theta for all sites with depth >= 2");
		thetaGenomeWide = true;
	}
	//initialize recalibration
	initializeRecalibration(params);

	//EM params
	EMParameters EMParams(params, logfile);

	//open output
	std::ofstream out; openThetaOutputFile(out);

	//prepare windows
	TWindowPairDiploid windows;

	if(considerRegions){
		TWindowDiploidSiteSubset* windowSitesSubset = new TWindowDiploidSiteSubset(mask);
		while(iterateChromosome(windows)){
			mask->setChr(chrIterator->Name);
			while(iterateWindow(windows)){
				if(readData(windows)){
					if(windows.cur->fractionSitesNoData > maxMissing){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//copy sites to a fake window
						logfile->listFlush("Adding relevant sites to data structure ...");
						windowSitesSubset->copySites(windows.cur);
						logfile->done();
						//estimate Theta
						windowSitesSubset->estimateTheta(EMParams, recalObject, out, logfile);
					}
				} else logfile->list("No relevant positions -> skipping this window.");
			}
		} delete windowSitesSubset;

	} else if(thetaGenomeWide){
		std::vector<TSiteDiploid*> siteVec;
		while(iterateChromosome(windows)){
			while(iterateWindow(windows)){
				if(readData(windows)){
					if(windows.cur->fractionSitesNoData > maxMissing){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//add informative sites to siteVec
						logfile->listFlush("Adding relevant sites to data structure ...");
						try{
						windows.cur->addSitesWithDepthTwoOrMoreToVector(siteVec);
						} catch(...){
							throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size or selecting fewer sites.";
						}
						logfile->done();
					}
				} else logfile->list("No relevant positions -> skipping this window.");
			}
		}
		//estimate Theta
		logfile->list("will estimate theta based on a total of " + toString(siteVec.size()) + " sites");
		TWindowDiploidSpecificSites specificSites =  TWindowDiploidSpecificSites(siteVec);
		specificSites.estimateTheta(EMParams, recalObject, out, logfile);

	} else {
		//iterate through windows
		while(iterateChromosome(windows)){
			while(iterateWindow(windows)){
				if(readData(windows)){
					if(windows.cur->fractionSitesNoData > maxMissing){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//estimate Theta
						out << chrIterator->Name << "\t";
						windows.cur->estimateTheta(EMParams, recalObject, out, logfile);
					}
				} else logfile->list("No relevant positions -> skipping this window.");
			}
		}
	}


	//clean up
	out.close();
}



void TGenome::calcLikelihoodSurfaces(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//read params
	int steps = params.getParameterIntWithDefault("steps", 100);
	int limitWindows = params.getParameterIntWithDefault("limitWindows", 1);

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	int windowsCalculated = 0;
	std::string filename;
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				//check if we have data -> can be extended to ensure
				if(windows.cur->fractionSitesNoData > maxMissing){
					logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
				} if(windows.cur->fractionRefIsN > maxRefN){
					logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
				} else {
					//open file
					std::ofstream out;
					filename = outputName + chrIterator->Name + "_" + toString(windows.cur->start) + "_LLsurface.txt";
					out.open(filename.c_str());
					if(!out) throw "Failed to open output file '" + outputName + "'!";

					//calc surface
					logfile->listFlush("Calculating likelihood surface ...");
					windows.cur->calcLikelihoodSurface(recalObject, out, steps);
					logfile->write(" done!");

					//close output
					out.close();

					//check if we break
					++windowsCalculated;
					if(windowsCalculated >= limitWindows) break;
				}
			}
		}
	}
}

void TGenome::callMLEGenotypes(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//set all booleans in the booleans class
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = false;
	bool gVCF = false;
	bool noAltIfHomoRef = false;
	bool beagle = false, printOnlyGL = false;
	std::string indName;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("beagle")){
			bool invariantSites = false;
			subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile, invariantSites);
			beagle=true;
			logfile->list("Will print output in beagle format");
			printOnlyGL = params.parameterExists("printOnlyGL");
			if(printOnlyGL) logfile->list("Will print only genotype likelihoods");
			indName = params.getParameterStringWithDefault("indName", outputName);
		}
	else if(params.parameterExists("sites")){
		bool invariantSites = false;
		if(fastaReference) subset = new TSiteSubset(params.getParameterString("sites"), reference, bamHeader, windowSize, logfile, invariantSites);
		else subset = new TSiteSubset(params.getParameterString("sites"), windowSize, logfile, invariantSites);
		limitToSitesWithKnownAlleles = true;
		if(params.parameterExists("noAltIfHomoRef")){
			noAltIfHomoRef = true;
			logfile->list("Will not print alternative alleles when genotype is 0/0");
		}
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
			logfile->list("Will print output in gVCF format");
		}
	}

	//open output: vcf or flat file?
	bool writeVCF = false;
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
		if(limitToSitesWithKnownAlleles) out << "\tAlt\tcoverage\tL(RR)\tL(RA)\tL(AA)\tMLE\tQ\n";
		else out << "\tcoverage\tL(AA)\tL(AC)\tL(AG)\tL(AT)\tL(CC)\tL(CG)\tL(CT)\tL(GG)\tL(GT)\tL(TT)\tMLE\tQ\n";
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles || beagle) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			if((!limitToSitesWithKnownAlleles && !beagle) || subset->hasPositionsInWindow(windows.cur->start)){
				//read data for current window
				if(readData(windows)){
					//call genotypes
					logfile->listFlush("Calling MLE genotypes ...");
					if(limitToSitesWithKnownAlleles || beagle){
						windows.cur->addReferenceBaseToSites(subset);
						windows.cur->callMLEGenotypeKnownAlleles(recalObject, subset, *randomGenerator, out, chrIterator->Name, writeVCF, noAltIfHomoRef, beagle, printOnlyGL);
					} else {
						if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
						windows.cur->callMLEGenotype(recalObject, *randomGenerator, out, chrIterator->Name, printIfNoData, fastaReference, writeVCF, gVCF, noAltIfHomoRef);
					}
					logfile->write(" done!");
				}
			}
		}
	}
	//clean up
	out.close();
}

void TGenome::callBayesianGenotypes(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	double theta;
	bool estimateTheta;
	EMParameters* EMParams = NULL;
	std::ofstream outTheta;
	if(params.parameterExists("theta")){
		estimateTheta = false;
		theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
	} else {
		estimateTheta = true;
		//read EM params
		EMParams = new EMParameters(params, logfile);
		openThetaOutputFile(outTheta);
	}

	//limit to a set of sites? Print all sites, even those without data?
	bool limitToSitesWithKnownAlleles = false;
	bool printIfNoData = true;
	TSiteSubset* subset = NULL;
	if(params.parameterExists("sites")){
		if(windowsPredefined) throw "Using site subsets is currently not implemented if windows are predefined from a BED file.";
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
		if(limitToSitesWithKnownAlleles) output << "\talt\tcoverage\tP(RR|D)\tP(RA|D)\tP(AA|D)\tMAP\tQ\n";
		else output << "\tcoverage\tP(AA|D)\tP(AC|D)\tP(AG|D)\tP(AT|D)\tP(CC|D)\tP(CG|D)\tP(CT|D)\tP(GG|D)\tP(GT|D)\tP(TT|D)\tMAP\tQ\n";
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			//read data for current window
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				if(readData(windows)){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing && estimateTheta){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN && estimateTheta){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//set Theta
						if(estimateTheta){
							outTheta << chrIterator->Name << "\t";
							windows.cur->estimateTheta((*EMParams), recalObject, outTheta, logfile);
						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							windows.cur->setTheta(theta);
						}

						//call Bayesian genotypes
						logfile->listFlush("Calling Bayesian Genotypes ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callBayesianGenotypeKnownAlleles(subset, *randomGenerator, output, chrIterator->Name, writeVCF);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callBayesianGenotype(*randomGenerator, output, chrIterator->Name, printIfNoData, fastaReference, writeVCF);
						}
						logfile->write(" done!");
					}
				}
			}
		}
	}

	//clean up
	if(estimateTheta){
		outTheta.close();
		delete EMParams;
	}
	if(limitToSitesWithKnownAlleles) delete subset;
}

void TGenome::callAllelePresence(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//do we estimate theta or is it given?
	double theta;
	bool estimateTheta;
	EMParameters* EMParams = NULL;
	std::ofstream outTheta;
	if(params.parameterExists("theta")){
		estimateTheta = false;
		theta = params.getParameterDouble("theta");
		logfile->list("Using theta = " + toString(theta));
	} else {
		estimateTheta = true;
		openThetaOutputFile(outTheta);
		//read EM params
		EMParams = new EMParameters(params, logfile);
	}

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
		if(limitToSitesWithKnownAlleles) outAllelePresence << "\tRef\tAlt\tcoverage\tP(Ref|D)\tP(Alt|D)\tMAP\tQ\n";
		else {
			if(fastaReference) outAllelePresence << "\tRef";
			outAllelePresence << "\tcoverage\tP(A|D)\tP(C|D)\tP(G|D)\tP(T|D)\tMAP\tQ\n";
		}
	}

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		if(limitToSitesWithKnownAlleles) subset->setChr(chrIterator->Name);
		while(iterateWindow(windows)){
			if(!limitToSitesWithKnownAlleles || subset->hasPositionsInWindow(windows.cur->start)){
				//read data for current window
				//TODO: ask Dan if following line should actually be if(readData || printIfNoData)
				if(readData(windows)){
					//check if we have data -> can be extended to ensure
					if(windows.cur->fractionSitesNoData > maxMissing && estimateTheta){
						logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
					} if(windows.cur->fractionRefIsN > maxRefN && estimateTheta){
						logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
					} else {
						//set Theta
						if(estimateTheta){
							outTheta << chrIterator->Name << "\t";
							windows.cur->estimateTheta((*EMParams), recalObject, outTheta, logfile);
						} else {
							windows.cur->calculateEmissionProbabilities(recalObject);
							windows.cur->estimateBaseFrequencies();
							windows.cur->setTheta(theta);
						}

						//call allele presence
						logfile->listFlush("Calling allele presence ...");
						if(limitToSitesWithKnownAlleles){
							windows.cur->addReferenceBaseToSites(subset);
							windows.cur->callAllelePresenceKnwonAlleles(subset, *randomGenerator, outAllelePresence, chrIterator->Name, writeVCF, noAltIfHomoRef);
						} else {
							if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
							windows.cur->callAllelePresence(*randomGenerator, outAllelePresence, chrIterator->Name, printIfNoData, fastaReference, writeVCF, noAltIfHomoRef);
						}
						logfile->write(" done!");
					}
				}
			} else logfile->list("No positions in this window.");
		}
	}

	//clean up
	if(estimateTheta){
		outTheta.close();
		delete EMParams;
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
	randomBases << "chr\tpos\tref\tcoverage\tpileup\trandom_base\n";

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows) || params.parameterExists("printAll")){
				//call random allele
				logfile->listFlush("Calling random base ...");
				if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
				windows.cur->callRandomBase(*randomGenerator, randomBases, chrIterator->Name, printIfNoData);
				logfile->write(" done!");

			} else logfile->list("No positions in this window.");
		}
	}
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


/*
	if(!bamName.empty()){

	}

	fillVectorFromStringWhiteSpaceSkipEmpty(trueLine, trueVec);
*/

}

void TGenome::printPileup(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

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
			if(readData(windows)){
				//print pileup
				windows.cur->printPileup(recalObject, out, chrIterator->Name);
			}
		}
	}

	//clean up
	out.close();
}

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
	TRecalibrationEM recalObjectEM(&bamHeader, filename, params, logfile);
	if(!recalObjectEM.estimatetionRequired){
		logfile->list("No need to estimate anything. Aborting Program.");
		return;
	}

	//prepare windows
	TWindowPairHaploid windows;

	//add sites to EM object
	logfile->startIndent("Reading data from windows:");
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)) windows.cur->addToRecalibrationEM(recalObjectEM);
			else logfile->list("No positions in this window.");
		}
	}
	//clean up memory
	windows.clear();
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
}*/

void TGenome::calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params){
	//create recalibration object
	std::string filename = "empty";
	TRecalibrationEM recalObjectEM(&bamHeader, filename, params, logfile);

	//prepare windows
	TWindowPairHaploid windows;

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
	int numMarginalGridPoint = params.getParameterIntWithDefault("numGridPoints", 51);
	recalObjectEM.calcLikelihoodSurface(outputName + "_LLsurface.txt", numMarginalGridPoint);
}

void TGenome::BQSR(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not run BQSR recalibration without a provided FASTA reference!";

	//prepare windows
	TWindowPairHaploid windows;

	//create BQSR object
	TRecalibrationBQSR bqsr(&bamHeader, params, logfile);
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

	/*
	//check if we use first chromosome for initial convergence
	if(params.parameterExists("preConverge")){
		int numPreConvLoops = params.getParameterInt("preConverge");
		int numChrPreConv = params.getParameterIntWithDefault("limitChrPreConv", 1) - 1;
		int numWindowsPreConv = params.getParameterIntWithDefault("limitWindowsPreConv", 100);
		logfile->startIndent("Running initial BQSR (pre-convergence):");
		logfile->list("Limiting pre-convergence to the first " + toString(numChrPreConv) + " chromosomes.");
		logfile->list("Limiting pre-convergence to the first " + toString(numWindowsPreConv) + " windows on each chromosome.");

		//run until it converges
		while(!hasConverged && loopNumber < numPreConvLoops){
			++loopNumber;
			logfile->startIndent("Running pre-recalibration loop " + toString(loopNumber) + ":");

			//jump to first chromosome
			while(chrNumber < numChrPreConv && iterateChromosome(windows)){
				//iterate over all windows
				while(windowNumber < numWindowsPreConv && iterateWindow(windows)){
					//read data for current window
					if(readData(windows)){
						//add reference data
						windows.cur->addReferenceBaseToSites(reference, chrNumber);

						//add the base to BQSR
						windows.cur->addSitesToBQSR(bqsr, logfile);
					}
				}
				logfile->endIndent();
			}
			logfile->endIndent();

			//clean up memory and restart from chr 1
			windows.clear();
			jumpToEnd();

			//estimate epsilon
			hasConverged = bqsr.estimateEpsilon();

			//write results to file
			bqsr.writeToFile(outputName + "_preConvergence_Loop" + toString(loopNumber));
			logfile->endIndent();
		}

		//reset counters and such
		bqsr.reopenEstimation();
		hasConverged = false;
		loopNumber = 0;
		logfile->endIndent();
		jumpToEnd();
	} */


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
}

void TGenome::printQualityTransformation(TParameters & params){
	//initialize recalibration
	//compare to a second recalibration definition?
	if(params.parameterExists("recal2") && !params.parameterExists("recal")) throw "use recal instead of recal2 for comparison of recalibrated qualities to original qualities!";
	recalObjectInitialized2 = false;
	if(params.parameterExists("recal")){
		std::string nameRecal = params.getParameterString("recal");
		recalObject = new TRecalibrationEM(&bamHeader, nameRecal, params, logfile);
		if(params.parameterExists("recal2")){
			std::string nameRecal2 = params.getParameterString("recal2");
			recalObject2 = new TRecalibrationEM(&bamHeader, nameRecal2, params, logfile);
			doRecalibration2 = true;
			recalObjectInitialized2 = true;
		} else if(params.parameterExists("BQSRQuality")){
			recalObject2 = new TRecalibrationBQSR(&bamHeader, params, logfile);
			doRecalibration2 = true;
			recalObjectInitialized2 = true;
		}
	} else if(params.parameterExists("BQSRQuality")){
		recalObject = new TRecalibrationBQSR(&bamHeader, params, logfile);
		doRecalibration = true;
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
		doRecalibration = false;
		recalObject = new TRecalibration();
	}
	//recalObjectInitialized = true;

	//check if estimation is required, in which case throw an error!
	if(recalObject->requiresEstimation()) throw "Can not use provided recalibration: estimation is required!";



	//prepare windows
	TWindowPairHaploid windows;

	int maxQ = params.getParameterIntWithDefault("maxQ", 100);



	//create table to store counts
	std::vector<TQualityTransformTable*> QTtables;
	for(int i=0; i<readGroups.numGroups; ++i){
		TQualityTransformTable* QT = new TQualityTransformTable(maxQ);
		QTtables.push_back(QT);
	}

	//prepare output
	std::ofstream out;
	std::string filename;

	//loop over windows and add to table
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			//read data for current window
			if(readData(windows)){
				if(recalObjectInitialized2) windows.cur->addSitesToQualityTransformTable(recalObject, recalObject2, QTtables, logfile);
				else windows.cur->addSitesToQualityTransformTable(recalObject, QTtables, logfile);
			}
		}
	}

	//print final tables
	for(int i=0; i<readGroups.numGroups; ++i){
		filename = outputName + "_" + readGroups.getName(i) + "_qualityTransformation.txt";
		out.open(filename.c_str());
		if(!out) throw "Failed to open output file '" + filename + "'!";
		QTtables[i]->printTable(out);
		out.close();

		//clean up vector
		delete QTtables.at(i);
	}

	//clean up memory
	windows.clear();

}

void TGenome::createBase(TBase** basePointer, char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId){
	if(base == 'A') *basePointer = new TBaseDiploidA(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);
	else if(base == 'C') *basePointer = new TBaseDiploidC(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);
	else if(base == 'G') *basePointer = new TBaseDiploidG(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);
	else *basePointer = new TBaseDiploidT(quality, posInRead, revPosInRead,  pmdCT, pmdGA, context, readGroupId);
}

char TGenome::returnBaseQualityAsChar(char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId){
	TBase* basePointer;
	createBase(&basePointer, base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);

	char qual = recalObject->getQualityAsChar(basePointer); //is in ascii, already has filter

	delete basePointer;
	return qual;
}

double TGenome::returnBaseQualityWithPMDAsCharRevMapping(char & base, char & refBase, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId){
	double qual = returnBaseQuality(base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId); //error rate
	if(base == 'T' && refBase == 'C') qual = 1.0 - ((1.0 - qual)*(1.0 - pmdGA)); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
	else if(base == 'A' && refBase == 'G') qual = 1.0 - ((1.0 - qual)*(1.0 - pmdCT)); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
	qual = round(-10.0 * log10(qual)) + 33.0; //cutoffs are in ascii
	if(qual < 33.0) qual = 33.0;
	else if(qual > 126.0) qual = 126.0;
	return qual;
}

double TGenome::returnBaseQualityWithPMDAsCharFwdMapping(char & base, char & refBase, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId){
	double qual = returnBaseQuality(base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);
	if(base == 'T' && refBase == 'C') qual = 1.0 - ((1.0 - qual)*(1.0 - pmdCT)); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
	else if(base == 'A' && refBase == 'G')	qual = 1.0 - ((1.0 - qual)*(1.0 - pmdGA)); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
	qual = round(-10.0 * log10(qual)) + 33.0; //cutoffs are in ascii
	if(qual < 33.0) qual = 33.0;
	else if(qual > 126.0) qual = 126.0;
	return qual;
}

double TGenome::returnBaseQuality(char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId){
	TBase* basePointer;
	createBase(&basePointer, base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);
	double qual = recalObject->getErrorRate(basePointer);
	delete basePointer;
	return qual;
}

bool TGenome::recalibrateAlignment(BamTools::BamAlignment & alignment, std::string & qual, TGenotypeMap & genoMap, bool withPMD, int & begin, std::string & ref, std::map <std::string, int> & mateTooLong){
	//variables
	char base, refBase, quality, newQual;
	BaseContext context;
	int posInRead, revPosInRead;
	double pmdCT, pmdGA;
	qual.clear();

	std::string bases = alignment.AlignedBases;
	std::string qualities = alignment.AlignedQualities;

	int len = bases.size();

	/*if(withPMD){
		//withPMD requires knowledge about distance from read ends. If there are indels, these distances can only be calculated by
		//taking cigar string into account. Since
		if(alignment.AlignedBases.size() != alignment.QueryBases.size()) withPMD = false;
	}*/
	int readGroupId = readGroups.find(alignment);

	//parse into bases
	if(alignment.IsProperPair()){
		if(abs(alignment.InsertSize) >= alignment.AlignedBases.size() && mateTooLong.count(alignment.Name) == 0 ){
			//now recalibrate
			if(alignment.IsReverseStrand()){
				// hence it is second in bam file and maps on reverse strand -> FLIP BASES
				//hence P(C->T) is given by  f(insert size - len + pos) (add this to the reverse table)
				//and P(G->A) is given as f(read len - pos - 1) (add this to forward table)
				for(int pos = 0; pos < len; ++pos){
					base = bases.at(pos);
					quality = qualities.at(pos);
					if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
						if((int)quality >= minQuality && (int)quality <= maxQuality){
							if(pos == (len - 1)) context = genoMap.getContextReverseRead('N', base);
							else context = genoMap.getContextReverseRead(alignment.AlignedBases.at(pos + 1), base);

							posInRead = abs(alignment.InsertSize)-len+pos;
							revPosInRead = len - pos - 1;
							pmdCT = pmdObjects[readGroupId].getProbGA(posInRead);
							pmdGA = pmdObjects[readGroupId].getProbCT(revPosInRead);

							if(withPMD){
								refBase = ref[bamAlignment.Position-begin+pos];
								newQual = returnBaseQualityWithPMDAsCharRevMapping(base, refBase, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);
							} else newQual = returnBaseQualityAsChar(base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);

							qual += newQual;

						} else qual += quality;
					} else if(base == '-'){
						//this means there was a deletion, don't want to add quality '!' to new quality string!
						continue;
					} else qual += quality;
				}
			} else {
				//Hence it is first in the bam file and maps on forward strand
				//Hence P(C->T) is given as a function of pos (add this to the in the forward table)
				//And P(G->A) is given by (insert size) - pos -1 (add this to the reverse table)
				for(int pos = 0; pos < len; ++pos){
					base = bases.at(pos);
					quality = qualities.at(pos);
					if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
						if((int)quality >= minQuality && (int)quality <= maxQuality){
							if(pos == 0) context = genoMap.getContext('N', base);
							else context = genoMap.getContext(alignment.AlignedBases.at(pos - 1), base);
							posInRead = pos;
							revPosInRead = abs(alignment.InsertSize) - pos - 1;

							pmdCT = pmdObjects[readGroupId].getProbCT(posInRead);
							pmdGA = pmdObjects[readGroupId].getProbGA(revPosInRead);

							//get new quality
							if(withPMD){
								refBase = ref[bamAlignment.Position-begin+pos];
								newQual = returnBaseQualityWithPMDAsCharFwdMapping(base, refBase, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);
							} else newQual = returnBaseQualityAsChar(base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);

							qual += newQual;
						} else 	qual += quality;
					} else if(base == '-'){
						//this means there was a deletion, don't want to add quality '!' to new quality string!
						continue;
					} else qual += quality;
				}
			}
		} else {
			logfile->warning("One mate of the following alignment pair is longer than its insert size: " + alignment.Name);
			mateTooLong.insert(make_pair(alignment.Name, 1));
			return false;
		}
	} else { //single end
		if(alignment.IsReverseStrand()){
			for(int pos = 0; pos < len; ++pos){
				base = bases.at(pos);
				quality = qualities.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
					if((int)quality >= minQuality && (int)quality <= maxQuality){
						if(pos == (len - 1)) context = genoMap.getContextReverseRead('N', base);
						else context = genoMap.getContextReverseRead(alignment.AlignedBases.at(pos + 1), base);
						posInRead = len - pos - 1;
						revPosInRead = pos;
						pmdCT = pmdObjects[readGroupId].getProbGA(posInRead);
						pmdGA = pmdObjects[readGroupId].getProbCT(revPosInRead);

						//get new quality
						if(withPMD){
							refBase = ref[bamAlignment.Position-begin+pos];
							newQual = returnBaseQualityWithPMDAsCharRevMapping(base, refBase, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);
						}
						else newQual = returnBaseQualityAsChar(base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);

						qual += newQual;

					} else qual += quality;
				} else if(base == '-'){
					//this means there was a deletion, don't want to add quality '!' to new quality string!
					continue;
				} else qual += quality;
			}
		} else {
			for(int pos = 0; pos < len; ++pos){
				base = bases.at(pos);
				quality = qualities.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
					if((int)quality >= minQuality && (int)quality <= maxQuality){
						if(pos == 0) context = genoMap.getContext('N', base);
						else context = genoMap.getContext(alignment.AlignedBases.at(pos - 1), base);
						posInRead = pos;
						revPosInRead = len - pos - 1;
						pmdCT = pmdObjects[readGroupId].getProbCT(posInRead);
						pmdGA = pmdObjects[readGroupId].getProbGA(revPosInRead);

						//get new quality
						if(withPMD){
							refBase = ref[bamAlignment.Position-bamAlignment.Position+pos];
							newQual = returnBaseQualityWithPMDAsCharFwdMapping(base, refBase, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);
						}
						else newQual = returnBaseQualityAsChar(base, quality, posInRead, revPosInRead, pmdCT, pmdGA, context, readGroupId);

						qual += newQual;

					} else qual += quality;
				} else if(base == '-'){
					//this means there was a deletion, don't want to add quality '!' to new quality string!
					continue;
				} else qual += quality;
			}
		}
	}
	return true;
}

void TGenome::recalibrateBamFile(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_recalibrated.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	TGenotypeMap genoMap;
	std::string qual;
	long counter = 0;
	std::map <std::string, int> mateTooLong;
	bool withPMD = params.parameterExists("withPMD");

	if(!withPMD && hasPMD) logfile->warning("The pmd pattern will not have any affect on the quality scores. If you want the quality scores to reflect pmd, use \"withPMD\"!");
	else if(withPMD && hasPMD) logfile->list("Probability of PMD will be reflected in new quality scores");
	else if(withPMD && !hasPMD) throw "Probability of PMD is unknown. Provide PMD patterns or remove \"withPMD\"";
	if(withPMD && !fastaReference) throw "Cannot run recalBAM withPMD without reference!";

	//declare reference (this is necessary for withPMD option)
	int curChr = -1;
	int len, begin = 0, windowSize = 1000000;
	int stop = begin + windowSize; //note that end is last position + 1
	std::string ref; //fasta object fills string

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		len = bamAlignment.AlignedBases.size();
		if(withPMD && (bamAlignment.Position + len >= stop || curChr!=bamAlignment.RefID)){
			curChr = bamAlignment.RefID;
			begin = bamAlignment.Position;
			stop = begin + windowSize;
			reference.GetSequence(curChr, begin, stop, ref);
		}
		++counter;
		//update and write (only if alignment qualities could be calculated)
		if(recalibrateAlignment(bamAlignment, qual, genoMap, withPMD, begin, ref, mateTooLong)){
			bamAlignment.Qualities = qual;
			bamAlignment.QueryBases = bamAlignment.AlignedBases;
			bamAlignment.CigarData.clear();
			bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MATCH_CHAR, bamAlignment.Qualities.size()));

			bamWriter.SaveAlignment(bamAlignment);
		}
		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::binQualityScores(TParameters & params){
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	std::string filename = outputName + "_binnedQualityScores.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//other temp variables
	TGenotypeMap genoMap;
	long counter = 0;
	int len;
	int q;
	char qChar;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		std::string qual;

		++counter;
		len = bamAlignment.Qualities.size();
		for(int i=0; i<len; ++i){
			q = (int) bamAlignment.Qualities[i] - 33;
			if(q >= 2 && q <= 9) qChar = 6 + 33;
			else if(q >= 10 && q <= 19) qChar = 15 + 33;
			else if(q >= 20 && q <= 24) qChar = 22 + 33;
			else if(q >= 25 && q <= 29) qChar = 27 + 33;
			else if(q >= 30 && q <= 34) qChar = 33 + 33;
			else if(q >= 35 && q <= 39) qChar = 37 + 33;
			else if(q >= 40) qChar = 40 + 33;
			else {
				logfile->warning("Quality in alignment " + bamAlignment.Name + " has quality " + toString(q) + ". Will keep it as is.");
				qChar = q + 33;
			}
			qual += qChar;
		}

		//update and write (only if alignment qualities could be calculated)
		bamAlignment.Qualities = qual;
		bamWriter.SaveAlignment(bamAlignment);

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
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

			//ad to map
			singleEndRG.insert(std::pair<int, TReadGroupMaxLength>(readGroupId, TReadGroupMaxLength(len, truncatedReadGroupId, readGroup)));
		}
	}
	logfile->write(" done!");
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

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		++counter;
		//get read group info
		bamAlignment.GetTag("RG", readGroup);
		readGroupId = readGroups.find(readGroup);

		//check if this RG needs to be parse
		singleEndRGIT = singleEndRG.find(readGroupId);
		if(singleEndRGIT != singleEndRG.end()){
			//check length
			if(bamAlignment.Length == singleEndRGIT->second.maxLen)
				bamAlignment.EditTag("RG", "Z", singleEndRGIT->second.truncatedReadGroup);
			else if(bamAlignment.Length > singleEndRGIT->second.maxLen && !allowForLarger) logfile->warning("Length of read in read group '" + readGroup + "' is > max length provided!");
		}

		//write
		bamWriter.SaveAlignment(bamAlignment);

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
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
	logfile->write(" done!");

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

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;
	std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		++counter;

		//get read group info
		bamAlignment.GetTag("RG", readGroup);
		oldId = readGroups.find(readGroup);

		//save as new RG
		bamAlignment.EditTag("RG", "Z", newReadGroupObject.getName(readGroupMap[oldId]));
		bamWriter.SaveAlignment(bamAlignment);

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::addReadToPMD(TWindowDiploid* window, TGenotypeMap & genoMap, std::string & ref, TPMDTables & pmdTables){
	//Extract Read Group Info
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
	int readGroupId = readGroups.find(readGroup);
	int length = bamAlignment.AlignedBases.size();
	char base;
	int quality;
	Base readBase, refBase;

	//add to PMD
	//distinguish between cases
	int internalPos = bamAlignment.Position - window->start;
	//paired end
	if(readGroups.inUse[readGroupId] == true){
		if(!bamAlignment.IsDuplicate()){
			if(bamAlignment.IsPrimaryAlignment()){
				if(bamAlignment.IsProperPair()){
					if(abs(bamAlignment.InsertSize) >= bamAlignment.AlignedBases.length()){
						if(bamAlignment.IsReverseStrand()){
							// hence it is second in bam file and maps on reverse strand -> FLIP BASES
							//hence P(C->T) is given by  f(insert size - len + pos) (add this to the reverse table)
							//and P(G->A) is given as f(read len - pos - 1) (add this to forward table)
							for(int pos = 0; pos < length; ++pos, ++internalPos){
								base = bamAlignment.AlignedBases[pos];
								if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip ann other
									quality = bamAlignment.AlignedQualities[pos];
									if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality d0es not make sense
										readBase = genoMap.flipBase(base);
										//std::cout << " " << internalPos << "," << ref[internalPos] << std::flush;
										refBase = genoMap.flipBase(ref[internalPos]);

										pmdTables.addForward(readGroupId, length - pos - 1, refBase, readBase);
										pmdTables.addReverse(readGroupId, abs(bamAlignment.InsertSize)-length+pos, refBase, readBase);
									}
								}
							}
						} else {
							//Hence it is first in the bam file and maps on forward strand
							//Hence P(C->T) is given as a function of pos (add this to the in the forward table)
							//And P(G->A) is given by (insert size) - pos -1 (add this to the reverse table)
							for(int pos = 0; pos < length; ++pos, ++internalPos){
								base = bamAlignment.AlignedBases[pos];
								if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
									quality = bamAlignment.AlignedQualities[pos];
									if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
										readBase = genoMap.getBase(base);
										refBase = genoMap.getBase(ref[internalPos]);

										pmdTables.addForward(readGroupId, pos, refBase, readBase);
										pmdTables.addReverse(readGroupId, bamAlignment.InsertSize - pos - 1, refBase, readBase);
									}
								}
							}
						}
					} else logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);

				//single end
				} else {
					if(bamAlignment.IsReverseStrand()){
						//single end & reverse
						//forward position = len - pos - 1
						//reverse position = pos
						//FLIP BASES!
						for(int pos = 0; pos < length; ++pos, ++internalPos){
							base = bamAlignment.AlignedBases[pos];
							if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip ann other
								quality = bamAlignment.AlignedQualities[pos];
								if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
									readBase = genoMap.flipBase(base);
									//std::cout << " " << internalPos << "," << ref[internalPos] << std::flush;
									refBase = genoMap.flipBase(ref[internalPos]);

									pmdTables.addForward(readGroupId, length - pos - 1, refBase, readBase);
									pmdTables.addReverse(readGroupId, pos, refBase, readBase);
								}
							}
						}
					} else {
						//single end & forward
						//forward position = pos
						//reverse position = len - pos -1
						for(int pos = 0; pos < length; ++pos, ++internalPos){
							base = bamAlignment.AlignedBases[pos];
							if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
								quality = bamAlignment.AlignedQualities[pos];
								if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
									readBase = genoMap.getBase(base);
									refBase = genoMap.getBase(ref[internalPos]);

									pmdTables.addForward(readGroupId, pos, refBase, readBase);
									pmdTables.addReverse(readGroupId, length - pos - 1, refBase, readBase);
								}
							}
						}
					}
				}
			}
		}
	}
}

void TGenome::estimatePMD(TParameters & params){
	//make sure FASTA is open
	if(!fastaReference) throw "Can not estimate PMD without a provided FASTA reference!";

	//prepare windows
	TWindowPairDiploid windows;

	//prepare PMD table
	int maxLength = params.getParameterIntWithDefault("length", 50);
	logfile->list("Estimating PMD at the first " + toString(maxLength) + " positions.");
	TPMDTables pmdTables(&readGroups, maxLength);

	//measure runtime
	struct timeval start, end;

	//tmp variables
	int fastaEnd;
	std::string ref;
	TGenotypeMap genoMap;
	long numreadsAdded;

	//iterate through windows
	while(iterateChromosome(windows)){
		while(iterateWindow(windows)){
			gettimeofday(&start, NULL);
			logfile->listFlush("Adding reads to PMD tables ...");
			numreadsAdded = 0;

			//get fasta reference
			fastaEnd = windows.cur->end + 500; //note that end is last position + 1
			reference.GetSequence(chrNumber, windows.cur->start, fastaEnd, ref);

			//still have old alignment to condider?
			if(oldAlignementMustBeConsidered && bamAlignment.Position >= windows.cur->end){
					logfile->write(" done!");
					logfile->conclude("No data in this window.");
			} else {
				if(oldAlignementMustBeConsidered){
					//add this read to window
					oldAlignementMustBeConsidered = false;
					++numreadsAdded;
					addReadToPMD(windows.cur, genoMap, ref, pmdTables);
				}

				//parse additional reads
				while(bamReader.GetNextAlignment(bamAlignment) && bamAlignment.RefID==chrNumber){
					//check if this read is beyond window
					if(bamAlignment.Position >= windows.cur->end){
						oldAlignementMustBeConsidered = true;
						break;
					}
					++numreadsAdded;
					addReadToPMD(windows.cur, genoMap, ref, pmdTables);
				}

				//report
				gettimeofday(&end, NULL);
				logfile->write(" done (added " + toString(numreadsAdded) + " reads in " + toString(end.tv_sec  - start.tv_sec) + "s)!");
			}
		}
	}

	//print tables and data
	std::string filename = outputName + "_PMD_Table.txt";
	logfile->listFlush("Writing PMD table to '" + filename + "' ...");
	pmdTables.writeTable(filename);
	logfile->write(" done!");
	filename = outputName + "_PMD_Table_counts.txt";
	logfile->listFlush("Writing PMD table of counts to '" + filename + "' ...");
	pmdTables.writeTableWithCounts(filename);
	logfile->write(" done!");
	filename = outputName + "_PMD_input_Empiric.txt";
	logfile->listFlush("Writing PMD input file to '" + filename + "' ...");
	pmdTables.writePMDFile(filename);
	logfile->write(" done!");

	//estimate exponential model
	filename = outputName + "_PMD_input_Exponential.txt";
	logfile->listFlush("Estimating PMD exponential models and writing them to '" + filename + "' ...");
	int numNRIterations = params.getParameterIntWithDefault("numNRIterations", 100);
	double eps = params.getParameterDoubleWithDefault("eps", 0.001);
	pmdTables.fitExponentialModel(numNRIterations, eps, filename, logfile);
	logfile->write(" done!");
}

float TGenome::calculatePMDS(int readGroup, char & ref, char & read, double & pmdCT, double & pmdGA, double & errorRate, double & pi, float & probPMD, float & probNoPMD){
	double epsThird = errorRate / 3.0;
	double fourEpsThird = 4.0*epsThird;
	//std::cout << "error: " << errorRate << " pmdGA " << pmdGA << " pmdCT " << pmdCT << std::endl;


	if(ref == 'A'){
		if(read == 'A'){
			probPMD = 1.0 - errorRate - pi + fourEpsThird*pi + pmdGA*pi/3.0*(1.0-fourEpsThird);
			probNoPMD = 1.0 - errorRate - pi + fourEpsThird*pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'C'){ //ok
			probPMD = errorRate - fourEpsThird*pi + pi - pi*pmdCT*(fourEpsThird-1.0);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'G'){
			probPMD = errorRate - fourEpsThird*pi + pi + pi*pmdGA*(fourEpsThird-1.0);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'T'){
			probPMD = errorRate - fourEpsThird*pi + pi + pi*pmdCT*(1.0-fourEpsThird);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
	}
	else if (ref == 'C'){
		if(read == 'A'){
			probPMD = errorRate + pi - 2.0*errorRate*pi + pi*pmdGA*(1.0-fourEpsThird);
			probNoPMD = errorRate + pi - 2.0*errorRate*pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'C'){
			probPMD = 1.0 - pi - errorRate + fourEpsThird*pi + (1.0-pi)*pmdCT*(fourEpsThird-1.0);
			probNoPMD = 1.0 - pi - errorRate + fourEpsThird*pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'G'){
			probPMD = errorRate - fourEpsThird*pi + pi + pmdGA*(fourEpsThird*pi - pi);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'T'){
			probPMD = epsThird + (1.0-pi)*pmdCT*(1.0-fourEpsThird);
			probNoPMD = epsThird;
			return probPMD/probNoPMD;
		}
	}
	else if (ref == 'G'){
		if(read == 'A'){
			probPMD = pmdGA*(3.0-3.0*pi+4.0*errorRate+4.0*errorRate*pi) + errorRate - fourEpsThird*pi + pi;
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'C'){
			probPMD = errorRate - fourEpsThird*pi + pi + pi*pmdCT*(fourEpsThird - 1.0);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'G'){
			probPMD = 1.0 - errorRate - pi + fourEpsThird*pi + (1.0-pi)*pmdGA*(fourEpsThird-1.0);
			probNoPMD = 1.0 - errorRate - pi + fourEpsThird*pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'T'){
			probPMD = errorRate - fourEpsThird*pi + pi + pi*pmdCT*(1.0-fourEpsThird);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
	}
	else if(ref == 'T'){
		if(read == 'A'){
			probPMD = errorRate - fourEpsThird*pi + pi - epsThird*pi*pmdCT + pi*pmdGA*(1.0-errorRate);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'C'){
			probPMD = errorRate - fourEpsThird*pi + pi + pi*pmdCT*(fourEpsThird - 1.0);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'G'){
			probPMD = errorRate - fourEpsThird*pi + pi + pi*pmdGA*(fourEpsThird - 1.0);
			probNoPMD = errorRate - fourEpsThird*pi + pi;
			return probPMD/probNoPMD;
		}
		else if(read == 'T'){
			probPMD = 1.0 - errorRate - pi + fourEpsThird*pi + pmdCT*(pi/3.0-fourEpsThird*pi/3.0);
			probNoPMD = 1.0 - errorRate - pi + fourEpsThird*pi;
			return probPMD/probNoPMD;
		}
	}
	return -1.0;
}

void TGenome::runPMDS(TParameters & params){
	//parse bam file and calculate PMDS for each read (seeSkoglund et al. 2014)
	//write new bam file with PMDS score added
	//parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)

	initializeRecalibration(params);

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
	int curChr = -1;
	long counter = 0, counterF = 0;

	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_PMDS.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	//tmp variables
	int len;
	char base, refBase, quality;
	BaseContext context;
	int distFrom5Prime, distFrom3Prime;
	double pmdCT, pmdGA, qual=-1.0;
	int readGroupId;
	TGenotypeMap genoMap;
	std::string readGroup;
	float PMDS = 0.0, probNoPMD = -1.0, probPMD = -1.0;


	if(!fastaReference) throw "Cannot run PMDS without reference!";

	//get reference
	int begin = 0;
	int windowSize = 1000000;
	int stop = begin + windowSize; //note that end is last position + 1
	std::string ref; //fasta object fills string

	//now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		++counter;
		PMDS = 0;
		len = bamAlignment.AlignedBases.size();

		//get readgroup info
		bamAlignment.GetTag("RG", readGroup);
		readGroupId = readGroups.find(readGroup);

		//update reference
		if(bamAlignment.Position + len >= stop || curChr!=bamAlignment.RefID){
			curChr = bamAlignment.RefID;
			begin = bamAlignment.Position;
			stop = begin + windowSize;
			reference.GetSequence(curChr, begin, stop, ref);
		}

		//paired-end
		if(bamAlignment.IsProperPair()){
			if(abs(bamAlignment.InsertSize) >= bamAlignment.AlignedBases.size()){
				if(bamAlignment.IsReverseStrand()){
					for(int pos = 0; pos < len; ++pos){
						base = bamAlignment.AlignedBases.at(pos);
						refBase = ref[bamAlignment.Position-begin+pos];
						if((base == 'A' || base == 'C' || base == 'G' || base == 'T') && (refBase == 'A' || refBase == 'C' || refBase == 'G' || refBase == 'T')){ //skip any other
							quality = bamAlignment.AlignedQualities[pos];
							if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
								if(pos == (len - 1)) context = genoMap.getContextReverseRead('N', base);
								else context = genoMap.getContextReverseRead(bamAlignment.AlignedBases.at(pos+1), base);
								distFrom5Prime = len - pos - 1; //is this 5' of the read?
								distFrom3Prime = abs(bamAlignment.InsertSize) - len + pos;
								//bases not flipped -> flip pmd pattern
								pmdCT = pmdObjects[readGroupId].getProbCT(distFrom5Prime); //the main pattern we want to correct for is pmdCT. why not distance from 3'?
								pmdGA = pmdObjects[readGroupId].getProbGA(distFrom3Prime);
								qual = returnBaseQuality(base, quality, distFrom5Prime, distFrom3Prime, pmdCT, pmdGA, context, readGroupId);

								PMDS += log(calculatePMDS(readGroupId, refBase, base, pmdCT, pmdGA, qual, pi, probPMD, probNoPMD));

							}
						}

					}
				} else {
					for(int pos = 0; pos < len; ++pos){
						base = bamAlignment.AlignedBases[pos];
						refBase = ref[bamAlignment.Position-begin+pos];
						if((base == 'A' || base == 'C' || base == 'G' || base == 'T') && (refBase == 'A' || refBase == 'C' || refBase == 'G' || refBase == 'T')){ //skip any other
							quality = bamAlignment.AlignedQualities[pos];
							if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
								if(pos == 0) context = genoMap.getContext('N', base);
								else context = genoMap.getContext(bamAlignment.AlignedBases.at(pos-1), base);
								distFrom5Prime = pos;
								distFrom3Prime = abs(bamAlignment.InsertSize) - pos - 1;
								pmdCT = pmdObjects[readGroupId].getProbCT(distFrom5Prime);
								pmdGA = pmdObjects[readGroupId].getProbGA(distFrom3Prime);
								qual = returnBaseQuality(base, quality, distFrom5Prime, distFrom3Prime, pmdCT, pmdGA, context, readGroupId);

								PMDS += log(calculatePMDS(readGroupId, refBase, base, pmdCT, pmdGA, qual, pi, probPMD, probNoPMD));
							}
						}
					}
				}


			} else logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
		}
		//single-end
		else{
			if(bamAlignment.IsReverseStrand()){
				for(int pos = 0; pos < len; ++pos){
					base = bamAlignment.AlignedBases.at(pos);
					refBase = ref.at(bamAlignment.Position-begin+pos);
					if((base == 'A' || base == 'C' || base == 'G' || base == 'T') && (refBase == 'A' || refBase == 'C' || refBase == 'G' || refBase == 'T')){ //skip any other
				//	if((refBase == 'C' && (base=='C'||base=='T')) || (refBase == 'G' && (base == 'G' || base == 'A'))){
						quality = bamAlignment.AlignedQualities.at(pos);
						if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
							if(pos == (len - 1)) context = genoMap.getContextReverseRead('N', base);
							else context = genoMap.getContextReverseRead(bamAlignment.AlignedBases.at(pos+1), base);
							distFrom5Prime = len - pos - 1;
							distFrom3Prime = pos;

							//bases not flipped -> flip pmd pattern
							pmdCT = pmdObjects[readGroupId].getProbGA(distFrom3Prime);
							pmdGA = pmdObjects[readGroupId].getProbCT(distFrom5Prime);

							qual = returnBaseQuality(base, quality, distFrom5Prime, distFrom3Prime, pmdCT, pmdGA, context, readGroupId);

							PMDS += log(calculatePMDS(readGroupId, refBase, base, pmdCT, pmdGA, qual, pi, probPMD, probNoPMD));
						}
					}
				}
			} else {
				for(int pos = 0; pos < len; ++pos){
					base = bamAlignment.AlignedBases.at(pos);
					refBase = ref.at(bamAlignment.Position-begin+pos);
					if((base == 'A' || base == 'C' || base == 'G' || base == 'T') && (refBase == 'A' || refBase == 'C' || refBase == 'G' || refBase == 'T')){ //skip any other
				//	if((refBase == 'C' && (base=='C'||base=='T')) || (refBase == 'G' && (base == 'G' || base == 'A'))){
						quality = bamAlignment.AlignedQualities[pos];
						if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense

							if(pos == 0) context = genoMap.getContext('N', base);
							else context = genoMap.getContext(bamAlignment.AlignedBases.at(pos-1), base);
							distFrom5Prime = pos;
							distFrom3Prime = len - pos - 1;

							pmdCT = pmdObjects[readGroupId].getProbCT(distFrom5Prime);
							pmdGA = pmdObjects[readGroupId].getProbGA(distFrom3Prime);
							qual = returnBaseQuality(base, quality, distFrom5Prime, distFrom3Prime, pmdCT, pmdGA, context, readGroupId);

							PMDS += log(calculatePMDS(readGroupId, refBase, base, pmdCT, pmdGA, qual, pi, probPMD, probNoPMD));
						}
					}
				}
			}
		}
		if(bamAlignment.HasTag("DS") == false) bamAlignment.AddTag("DS", "f", PMDS);
		else bamAlignment.EditTag("DS", "f", PMDS);

		//update and write
		if(PMDS > minPMDS && PMDS < maxPMDS) bamWriter.SaveAlignment(bamAlignment);
		else ++counterF;

		//report progress
		if(counter % 1000000 == 0){
		gettimeofday(&end, NULL);
		logfile->list("Analyzed " + toString(counter) + " reads in " + toString(end.tv_sec  - start.tv_sec) + "s and filtered out " + toString(counterF) + " of them!");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report

	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Analyzed " + toString(counter) + " reads in " + toString(runtime) + " min. and filtered out " + toString(counterF) + " of them!");

}

void TGenome::mergePairedEndReads(TParameters & params){
	//open a bam file for writing
	BamTools::BamWriter bamWriter;
	filename = outputName + "_mergedReads.bam";
	BamTools::RefVector references = bamReader.GetReferenceData();
	logfile->list("Writing results to '" + filename + "'.");
	if (!bamWriter.Open(filename, bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";

	if(hasPMD) logfile->warning("PMD is given but relevant for read merging.");

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

	//create storage for reads until their mate was found
	std::vector< std::pair<BamTools::BamAlignment*, bool> > alignmentStorage;
	std::vector< std::pair<BamTools::BamAlignment*, bool> >::iterator it;
	std::vector< std::pair<BamTools::BamAlignment*, bool> >::iterator IT;

	BamTools::BamAlignment* alignmentPointer;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;
	int curChr = -1;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		++counter;
		if((readsToOmit.count(bamAlignment.Name) > 0)){
			continue;
			//alignmentStorage.push_back(std::pair<BamTools::BamAlignment*, bool>(new BamTools::BamAlignment(bamAlignment), true));
		} else if(abs(bamAlignment.InsertSize) < bamAlignment.AlignedBases.size()){
			logfile->warning("filtered out because of adapter: " + bamAlignment.Name);
			readsToOmit.insert(std::pair<std::string,int>(bamAlignment.Name, 1));
			continue;
		} else if(!bamAlignment.IsProperPair() || !bamAlignment.IsPrimaryAlignment()|| bamAlignment.IsDuplicate() || (bamAlignment.IsReverseStrand() && bamAlignment.IsMateReverseStrand()) || (!bamAlignment.IsReverseStrand() && !bamAlignment.IsMateReverseStrand())){
			readsToOmit.insert(std::pair<std::string,int>(bamAlignment.Name, 1));
			continue;
		}
		else {
			//if on new chromosome, empty storage
			if(curChr != bamAlignment.RefID){
				for(it = alignmentStorage.begin(); it != alignmentStorage.end(); ++it){
					bamWriter.SaveAlignment(*(it->first));
					delete it->first;
				}
				alignmentStorage.clear();
				curChr = bamAlignment.RefID;
			}
			//add alignment to storage
					if(bamAlignment.IsProperPair()){
						if(!bamAlignment.IsReverseStrand()) {
							//mate on forward strand is always first in bam file
							alignmentStorage.push_back(std::pair<BamTools::BamAlignment*, bool>(new BamTools::BamAlignment(bamAlignment), false));
						}
						else if(bamAlignment.IsReverseStrand()){
							//find first mate -> should be in storage
							for(it=alignmentStorage.begin(); it!=alignmentStorage.end(); ++it){
								if(it->first->Name == bamAlignment.Name){
									//check if this read accepts mate
									if(it->second) throw "First read of '" + bamAlignment.Name + "' is not paired or has already been merged!";
									//merge reads
									alignmentPointer = it->first;
									if(bamAlignment.Position > alignmentPointer->Position + alignmentPointer->AlignedBases.size()){
										//reads do not overlap -> add Ns in between
										//int numN = bamAlignment.Position - (alignmentPointer->Position + alignmentPointer)->AlignedBases.size() - 1;
										int numN = abs(bamAlignment.InsertSize) - bamAlignment.AlignedBases.size() - alignmentPointer->AlignedBases.size();
										for(int i=0; i<numN; ++i){
											alignmentPointer->AlignedBases += 'N';
											alignmentPointer->AlignedQualities += '!';
										}
										alignmentPointer->AlignedBases += bamAlignment.AlignedBases;
										alignmentPointer->AlignedQualities += bamAlignment.AlignedQualities;
									} else if(bamAlignment.Position < alignmentPointer->Position) std::cout << "Second mate of read '" << bamAlignment.Name << "' has pos < pos of first mate!";
									else {
										//reads do overlap -> merge them
										std::string alignment;
										std::string quality;
										int firstOverlap = bamAlignment.Position - alignmentPointer->Position;
										int lastOverlapPlusOne = -1;
										if((alignmentPointer->Position + alignmentPointer->AlignedBases.size()) <= (bamAlignment.Position + bamAlignment.AlignedBases.size())) lastOverlapPlusOne = alignmentPointer->AlignedBases.size();
										else lastOverlapPlusOne = bamAlignment.AlignedBases.size();

										//first copy from first mate
										alignment = alignmentPointer->AlignedBases.substr(0, firstOverlap);
										quality = alignmentPointer->AlignedQualities.substr(0, firstOverlap);
										//decide which alignment has higher quality in overlap
										for(int i=firstOverlap; i<lastOverlapPlusOne; ++i){
											//decide which quality is higher
											if(alignmentPointer->AlignedQualities.at(i) > bamAlignment.AlignedQualities.at(i - firstOverlap)){
												alignment += alignmentPointer->AlignedBases.at(i);
												quality += alignmentPointer->AlignedQualities.at(i);
											} else {
												alignment += bamAlignment.AlignedBases.at(i - firstOverlap);
												quality += bamAlignment.AlignedQualities.at(i - firstOverlap);
											}
										}

										//add rest from second
										if(alignmentPointer->Position + alignmentPointer->AlignedBases.size() < bamAlignment.Position + bamAlignment.AlignedBases.size()){
											alignment += bamAlignment.AlignedBases.substr(lastOverlapPlusOne - firstOverlap);
											quality += bamAlignment.AlignedQualities.substr(lastOverlapPlusOne - firstOverlap);
										}

										if(alignment.length() != abs(bamAlignment.InsertSize) + 1  && alignment.length() != abs(bamAlignment.InsertSize)) throw "merged alignment length of reads " + bamAlignment.Name + " is not equal to original insert size (+1)!";

										//set
										alignmentPointer->AlignedBases = alignment;
										alignmentPointer->AlignedQualities = quality;
									}
									//update

									alignmentPointer->QueryBases = alignmentPointer->AlignedBases;
									alignmentPointer->Qualities = alignmentPointer->AlignedQualities;
									alignmentPointer->Length = alignmentPointer->AlignedBases.size();
									alignmentPointer->CigarData.clear();
									alignmentPointer->CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MATCH_CHAR, alignmentPointer->Length));
									alignmentPointer->SetIsFirstMate(false);
									alignmentPointer->SetIsSecondMate(false);
									alignmentPointer->SetIsPaired(false);
									alignmentPointer->SetIsProperPair(false);
									alignmentPointer->SetIsMateReverseStrand(false);
									alignmentPointer->SetIsReverseStrand(false); //the read that comes first in BAM is always fwd strand
									alignmentPointer->MateRefID = -1;
									alignmentPointer->MatePosition = -1;
									alignmentPointer->InsertSize = 0;

									it->second = true;
	//								if(bamAlignment.Name == "HWI-ST558:352:C7T9BACXX:1:1201:7676:94794") std::cout << "aligment " << bamAlignment.Name << "was updated" << std::endl;

									//write if is first in vector
									if(it == alignmentStorage.begin()){

										//write all that are OK
										for(; it != alignmentStorage.end(); ++it){
											if(it->second){
												bamWriter.SaveAlignment(*(it->first)); //saves the alignment to the bam file
	//											if(bamAlignment.Name == "HWI-ST558:352:C7T9BACXX:1:2101:1301:32927") std::cout << "aligment " << bamAlignment.Name << " was output" << std::endl;
												delete it->first;
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

							if(!alignmentStorage.empty() && it == alignmentStorage.end()) throw "One read of '" + bamAlignment.Name + "' is reverse mate, but forward one has not been read!";
						} else{
							for(it=alignmentStorage.begin(); it != alignmentStorage.end(); ++it){
								delete it->first;
							}
							alignmentStorage.clear();
							throw "One read of '" + bamAlignment.Name + "' is paired, but neither first nor second mate!";
						}
					} else {
						//read is not paired: add to storage or write
						if(alignmentStorage.empty()) bamWriter.SaveAlignment(bamAlignment);
						else alignmentStorage.push_back(std::pair<BamTools::BamAlignment*, bool>(new BamTools::BamAlignment(bamAlignment), true));
					}


		}

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}


void TGenome::generatePSMCInput(TParameters & params){
	//initialize recalibration
	initializeRecalibration(params);

	//read in parameters required
	double theta = params.getParameterDoubleWithDefault("theta", 0.001);
	logfile->list("Using theta = " + toString(theta));
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
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		if(nCharOnLine > 0) output << '\n';
		output << '>' << chrIterator->Name << '\n';
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			//set Theta
			logfile->listFlush("Calculating emission probabilities ...");
			windows.cur->calculateEmissionProbabilities(recalObject);
			windows.cur->estimateBaseFrequencies();
			windows.cur->setTheta(theta);
			logfile->write(" done!");

			//create PSMC input
			logfile->listFlush("Estimating heterozygosity status ...");
			if(fastaReference) windows.cur->addReferenceBaseToSites(reference, chrNumber);
			windows.cur->generatePSMCInput(blockSize, confidence, output, nCharOnLine);
			logfile->write(" done!");
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
}


void TGenome::downSampleBamFile(TParameters & params){
	//read downsampling rate
	std::string prob = params.getParameterString("prob");
	int times = params.getParameterIntWithDefault("times", 1);
	//check if prob is a vector of multiple probabilities
	std::vector<double> downSampleProbVector;
	if(!stringContainsOnly(prob, "-0123456789.,")) throw "Wrong format on probability list: use floating point numbers delimited by commas (e.g. 0.1,0.2,0.5).";
	fillVectorFromString(prob, downSampleProbVector, ',');

	std::vector<double>::iterator it;
	int numProbs = downSampleProbVector.size();

	if(times > 1 && numProbs > 1) throw "Replicated downsampling is not implemented for more than one prob";
	if(times > 1) numProbs = times;

	//check if probs are between 0 and 1, save in array and print them
	double* downSampleProb = new double[numProbs];
	logfile->listFlush("Will accept reads with probabilities");
	bool first = true;
	int i=0;
	if(times == 1){
		for(it=downSampleProbVector.begin(); it!=downSampleProbVector.end(); ++it, ++i){
			if(first) first = false;
			else logfile->flush(",");
			logfile->flush(" " + toString(*it));
			if(*it <= 0.0 || *it >= 1.0) throw "All probabilities have to be between >0 and <1!";
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
			filename = outputName + "_downsampled" + toString(downSampleProb[i]) + "_" + toString(i) + ".bam";
			logfile->list(filename);
			//open file
			if(!bamWriter[i].Open(filename, bamHeader, references))	throw "Failed to open BAM file '" + filename + "'!";
		}
	} else {
		for(i=0; i<numProbs; ++i){
			//construct and print filename
			filename = outputName + "_downsampled" + toString(downSampleProb[i]) + ".bam";
			logfile->list(filename);
			//open file
			if(!bamWriter[i].Open(filename, bamHeader, references))	throw "Failed to open BAM file '" + filename + "'!";
		}
	}
	logfile->endIndent();

	//other temp variables
	long counter = 0;
	double r;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		++counter;

		//accept read or not?
		r = randomGenerator->getRand();
		for(i=0; i<numProbs; ++i){
			if(r < downSampleProb[i]){
				bamWriter[i].SaveAlignment(bamAlignment);
			}
		}

		//report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer and clean up memory
	for(i=0; i<numProbs; ++i){
		bamWriter[i].Close();
	}
	delete[] downSampleProb;
	delete[] bamWriter;

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::downSampleReads(TParameters & params){

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
	double r;

	//prepare reporting
	logfile->startIndent("Parsing through BAM file:");
	struct timeval start, end;
    gettimeofday(&start, NULL);
	float runtime;

    //now parse through bam file and write alignments
	while (bamReader.GetNextAlignment(bamAlignment)){
		++counter;
		for(int i=0; i<bamAlignment.Length; ++i){
			r = randomGenerator->getRand();
			if(r < fraction){
				bamAlignment.QueryBases.at(i) = 'N';
				bamAlignment.Qualities.at(i) = '!';
			}
		}
		bamWriter.SaveAlignment(bamAlignment);

		// intermetiate report report
		if(counter % 1000000 == 0){
			gettimeofday(&end, NULL);
			runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
		}
	}

	//close bam writer
	bamWriter.Close();

	//report
	gettimeofday(&end, NULL);
	runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
	logfile->list("Reached end of BAM file!");
	logfile->removeIndent();
}

void TGenome::estimateApproximateCoverage(TParameters & params){    //get genome length
    //open output files
    std::ofstream outputCoverage;
    std::string outputFileNameCov = outputName + "_approximateCoverage.txt";
    logfile->list("Writing coverage estimates to '" + outputFileNameCov + "'");
    outputCoverage.open(outputFileNameCov.c_str());
    if(!outputCoverage) throw "Failed to open output file '" + outputFileNameCov + "'!";

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

    double totLength = 0.0;
    for(chrIterator = bamHeader.Sequences.Begin(); chrIterator!=bamHeader.Sequences.End(); ++chrIterator)
        totLength += stringToLong(chrIterator->Length);

    //prepare reporting
    logfile->startIndent("Parsing through BAM file:");
    struct timeval start, end;
    gettimeofday(&start, NULL);
    float runtime;

    //other temp variables
    long counter = 0;
    //double toNumAlignedBases = 0.0;
    //std::vector<int*> MQ;
    //std::vector<int*> RL;
    int RGInd = -1;

/*
    long MQ[readGroups.numGroups][100];
    long RL[readGroups.numGroups][500];
    double cov[readGroups.numGroups];
*/

    std::vector<double> cov;
    double totCov =0.0;

    long** MQ = new long*[readGroups.numGroups];
    long** RL = new long*[readGroups.numGroups];
    for(int i = 0; i < readGroups.numGroups; ++i){
    	cov.push_back(0);
    	MQ[i] = new long[100];
    	RL[i] = new long[500];
//    	memset(MQ[i], 0, 100);
 //   	memset(RL[i], 0, 100);
    	for(int j=0; j<100; ++j) MQ[i][j]=0;
    	for(int j=0; j<500; ++j) RL[i][j]=0;
    }


    //now parse through bam file and sum number of aligned bases
    while (bamReader.GetNextAlignment(bamAlignment)){
    	//filters
        if(!readGroups.readGroupInUse(bamAlignment)) continue;
        if(!useChromosome[bamAlignment.RefID]) continue;
        if(bamAlignment.IsDuplicate()) continue;
        if(!bamAlignment.IsPrimaryAlignment()) continue;
        if(bamAlignment.IsPaired() && !bamAlignment.IsProperPair()) continue;

        ++counter;
        RGInd = readGroups.find(bamAlignment);
        totCov += bamAlignment.AlignedBases.length();
        cov[RGInd] += bamAlignment.AlignedBases.length();
        ++MQ[RGInd][bamAlignment.MapQuality];
        ++RL[RGInd][bamAlignment.Length];

        //report
        if(counter % 1000000 == 0){
            gettimeofday(&end, NULL);
            runtime = (end.tv_sec  - start.tv_sec)/60.0;
            logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
        }
    }

    //report
    gettimeofday(&end, NULL);
    runtime = (end.tv_sec  - start.tv_sec)/60.0;
    logfile->list("Parsed " + toString(counter) + " reads in " + toString(runtime) + " min.");
    logfile->list("Reached end of BAM file!");
    logfile->removeIndent();
    logfile->list("Approximate coverage was estimated at " + toString(totCov/totLength));

    logfile->listFlush("Writing to output files ...");

    //cov
    outputCoverage << "RG\tApproximate_coverage";
    outputCoverage << "\nallReadGroups\t" << totCov/totLength;
    for(int r=0; r<readGroups.numGroups; ++r){
        outputCoverage << "\n" << readGroups.getName(r) << "\t" << cov[r]/totLength;
    }

    //MQ
    long tot;
    outputMQ << "RG\tMapping_quality\tCount";
    for(int i=0; i<100; ++i){
    	tot = 0;
    	for(int r=0; r<readGroups.numGroups; ++r) tot += MQ[r][i];
		outputMQ << "\nallReadGroups\t" << i << "\t" << tot;

    }
    for(int r=0; r<readGroups.numGroups; ++r){
        for(int i=0; i<100; ++i){
            outputMQ << "\n" << readGroups.getName(r) << "\t" << i << "\t" << MQ[r][i];
        }
    }

    //RL
    outputReadLen << "RG\tRead_length\tCount";
    for(int i=0; i<500; ++i){
    	tot = 0;
    	for(int r=0; r<readGroups.numGroups; ++r) tot += RL[r][i];
		outputReadLen << "\nallReadGroups\t" << i << "\t" << tot;
    }
    for(int r=0; r<readGroups.numGroups; ++r){
        for(int i=0; i<500; ++i){
            outputReadLen << "\n" << readGroups.getName(r)<< "\t" << i << "\t" << RL[r][i];
        }
    }

    logfile->write(" done!");

    outputCoverage.close();
    outputMQ.close();
    outputReadLen.close();

    for(int i = 0; i < readGroups.numGroups; ++i){
    	delete MQ[i];
    	delete RL[i];
    }
    delete MQ;
    delete RL;
}



void TGenome::estimateApproximateCoveragePerWindow(TParameters & params){
	//open output file
	std::ofstream output;
	std::string outputFileName = outputName + "_coverage.txt";
	logfile->list("Writing coverage estimates to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int nCharOnLine = 0;

	//write header
	output << "chr\tstart\tend\tcoverage" << std::endl;

	//prepare windows
	TWindowPairDiploid windows;

	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);

			//write to file
			logfile->listFlush("Writing coverage to file ...");
			if(windows.cur->coverage == -1.0) output << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t" << "0" << "\n";
			else output << chrIterator->Name << "\t" << windows.cur->start << "\t" << windows.cur->end << "\t" << windows.cur->coverage << "\n";
			logfile->write(" done!");
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
}

void TGenome::estimateCoveragePerSite(TParameters & params){
	std::ofstream output;
	std::string outputFileName = outputName + "_coveragePerSite.txt";
	logfile->list("Writing coverage estimates to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int maxCov = params.getParameterIntWithDefault("maxCov", 20);
	if(!maxCov) throw "No maximum coverage specified!";
	int size = maxCov + 2; // need 0 bin and >maxCov bin
	int nCharOnLine = 0;

	//prepare array
	long * siteCoverage = new long[size];
	for(int i=0; i<size; ++i){
		siteCoverage[i] = 0;
	}

	//write header
	output << "coverage\tcounts" << std::endl;

	//prepare windows
	TWindowPairDiploid windows;


	//iterate through windows
	while(iterateChromosome(windows)){
		//write chromosome to file
		while(iterateWindow(windows)){
			//read data for current window
			readData(windows);
			windows.cur->calcCoveragePerSite(siteCoverage, maxCov);

			logfile->listFlush("Adding coverages to table ...");
			logfile->write(" done!");
		}
	}

	//write to file
	for(int i=0; i<(size-1); ++i){
		output << i << "\t" << siteCoverage[i] << "\n";
	}
	output << ">" << maxCov << "\t" << siteCoverage[size - 1] << std::endl;

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
	delete[] siteCoverage;
}
