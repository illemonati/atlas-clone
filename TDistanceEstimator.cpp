/*
 * TDistanceCalculator.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#include "TDistanceEstimator.h"

TDistanceEstimator::TDistanceEstimator(TLog* Logfile){
	logfile = Logfile;
	phredToLikelihood = NULL;
}

void TDistanceEstimator::printGLF(TParameters & params){

	//test first to parse GLF files
	std::string glf = params.getParameterString("glf");
	TGlfReader reader(glf);

	//print file
	reader.printToEnd();
}

void TDistanceEstimator::initializePhredToLikelihoodTable(){
	phredToLikelihood = new double[256];
	phredToLikelihood[0] = 1.0;
	for(int i=0; i<256; ++i){
		phredToLikelihood[i] = pow(10.0, (double) -i/10.0);
	}
	phredToLikelihoodInitialized = true;
}

void TDistanceEstimator::deletePhredToLikelihoodTable(){
	if(phredToLikelihoodInitialized)
		delete[] phredToLikelihood;
}


//------------------------------------------------------------------
void TDistanceEstimator::estimateDistances(TParameters & params){
	//prepare stuff
	initializePhredToLikelihoodTable();

	//open all GLF files specified
	std::vector<std::string> glfNames;
	params.fillParameterIntoVector("glf", glfNames, ',');
	int numGLFs = glfNames.size();
	if(numGLFs < 2)
		throw "At least two GLF files have to be provided to estimate distances!";

	//open files
	TGlfReader* glfs = new TGlfReader[numGLFs];
	logfile->startIndent("Opening GLF files:");
	int g1 = 0;
	for(std::vector<std::string>::iterator it=glfNames.begin(); it != glfNames.end(); ++it, ++g1){
		logfile->list("Opening GLF '" + *it + "' ...");
		glfs[g1].open(*it);
		logfile->done();
	}
	logfile->endIndent();

	//in windows or whole genome?
	long windowLen = params.getParameterLongWithDefault("window", 1000000);
	logfile->list("Will estimate genetic distance in windows of length " + toString(windowLen) + ".");

	//now calculate all pairwise distances
	int g2;
	std::string filename;
	for(g1=0; g1<(numGLFs-1); ++g1){
		for(g2 = g1+1; g2 < numGLFs; ++g2){
			//estimate distance
			std::cout << "RUN " << g1 << " and " << g2 << ".... " << std::endl;

			filename = glfNames[g1] + "_" + glfNames[g2] + "_distanceEstimates.txt.gz";

			estimateDistanceInWindows(filename, glfs[g1], glfs[g2], windowLen);

			//report

			//rewind GLFs
			glfs[g1].rewind();
			glfs[g2].rewind();

		}
	}

	//close all glf handlers
	for(g1=0; g1<numGLFs; ++g1)
		glfs[g1].close();
	delete[] glfs;
}

/*
void TDistanceEstimator::estimateDistanceInWindows(TGlfReader & g1, TGlfReader & g2, long windowLen){
	//initialize variables
	std::string curChr = "";
	long lastPosOfCurWindow = windowLen;
	bool keepReading = true;
	bool isGood1 = true;
	bool isGood2 = true;

	//initialize storage for two windows


	//read first positions
	if(!g1.readNext()) throw "Failed to read first position from GLF '" + g1.name() + "'!";
	if(!g2.readNext()) throw "Failed to read first position from GLF '" + g2.name() + "'!";

	//parse GLFs in windows
	while(keepReading){

		//if on new chromosome, window is done
		if(g1.chr != curChr){
			//something happens
		}

		if(isGood1 && isGood2){
			//keep adding to window
			if(g2.position == g1.position){
				//add data

				//advance both
				isGood1 = g1.readNext();
				isGood2 = g2.readNext();
			} else if(g2.position < g1.position){
				//add data

				//advance g2
				isGood2 = g2.readNext();
			} else {
				//add data

				//advance g1
				isGood1 = g1.readNext();
			}
		} else if(isGood1 && !isGood2){
			//think about case in which one  file has reached end
		} else if(!isGood1 && isGood2){

		} else {
			//reached end!

		}

		//is window done?
		//TODO: needs to also be called after a chromosom jump!!
		if(g1.position > lastPosOfCurWindow && g2.position > lastPosOfCurWindow){
			//do calculations for this window
		}



	}

	//do calculations for last window
}
*/

void TDistanceEstimator::writeDistanceEstimates(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd, int & numsitesWithData){
	out << chr << "\t" << windowStart << "\t" << windowEnd << "\t" << numsitesWithData << "\n";
}

void TDistanceEstimator::writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd){
	out << chr << "\t" << windowStart << "\t" << windowEnd << "\t0\n";
}

void TDistanceEstimator::estimateDistanceInWindows(std::string filename, TGlfReader & g1, TGlfReader & g2, long windowLen){
	//initialize variables
	bool keepReading = true;
	bool isGood1 = true;
	bool isGood2 = true;

	//initialize storage for two windows
	//TODO: share across pairs? Do all pairs at once?
	int** genoQual1 = new int*[windowLen];
	int** genoQual2 = new int*[windowLen];
	for(int i=0; i<windowLen; ++i){
		genoQual1[i] = new int[10];
		genoQual2[i] = new int[10];
	}

	//open output file
	gz::ogzstream out(filename.c_str());
	if(!out) throw "Failed to ipen file '" + filename + "' for writing!";

	//prepare variables
	std::string curChr;
	long curChrLen;
	long windowStart;
	long windowEnd;

	int numSitesWithData = 100;

	//parse GLFs in windows
	while(keepReading){
		//move to new chromosome
		curChr = g1.chr();
		curChrLen = g1.chrLength();
		long windowStart = 1;
		long windowEnd = windowLen+1;

		//parse all windows of chromosome
		while(windowStart < curChrLen){

			std::cout << "Chr '" << curChr << "' [" << windowStart  << ", " << windowEnd << "):" << std::endl;

			//read data
			isGood1 = g1.readNextWindow(genoQual1, curChr, windowStart, windowEnd);
			if(isGood1 || g1.eof()){
				isGood2 = g2.readNextWindow(genoQual2, curChr, windowStart, windowEnd);
				if(isGood2){
					//estimate distance

					//write to file
					writeDistanceEstimates(out, curChr, windowStart, windowEnd, numSitesWithData);


				} else writeDistanceEstimatesNoData(out, curChr, windowStart, windowEnd);
			} writeDistanceEstimatesNoData(out, curChr, windowStart, windowEnd);

/*
			//print data
			for(int i=0; i<windowLen; ++i){
				std::cout << "Sample 1: " << genoQual1[i][0];
				for(int g=1; g<10; ++g){
					std::cout << "," << genoQual1[i][g];
				}
				std::cout << "\tSample 2:" << genoQual2[i][0];
				for(int g=1; g<10; ++g){
					std::cout << "," << genoQual2[i][g];
				}
				std::cout << std::endl;
			}
*/

			//move window
			windowStart = windowEnd;
			windowEnd = windowStart + windowLen;

		}

		std::cout << "END OF WINDOW LOOP ...." << std::endl;

		std::cout << "G1.EOF = " << g1.eof() << std::endl;
		std::cout << "G2.EOF = " << g2.eof() << std::endl;

		if(g1.eof() && g2.eof())
			keepReading = false;
	}

}
























