/*
 * TDistanceCalculator.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#include "TDistanceEstimator.h"

//-------------------------------------------------
//TPhiToGenoMap
//-------------------------------------------------
TGenoToPhiMap::TGenoToPhiMap(){
	genoToPhiMap = new int*[10];
	for(int i=0; i<10; ++i)
		genoToPhiMap[i] = new int[10];

	//case aa/aa
	genoToPhiMap[AA][AA] = 0;
	genoToPhiMap[CC][CC] = 0;
	genoToPhiMap[GG][GG] = 0;
	genoToPhiMap[TT][TT] = 0;

	//case aa/ab
	genoToPhiMap[AA][AC] = 1;
	genoToPhiMap[AA][AG] = 1;
	genoToPhiMap[AA][AT] = 1;
	genoToPhiMap[CC][AC] = 1;
	genoToPhiMap[CC][CG] = 1;
	genoToPhiMap[CC][CT] = 1;
	genoToPhiMap[GG][AG] = 1;
	genoToPhiMap[GG][CG] = 1;
	genoToPhiMap[GG][GT] = 1;
	genoToPhiMap[TT][AT] = 1;
	genoToPhiMap[TT][CT] = 1;
	genoToPhiMap[TT][GT] = 1;

	//case ab/aa
	genoToPhiMap[AC][AA] = 2;
	genoToPhiMap[AG][AA] = 2;
	genoToPhiMap[AT][AA] = 2;
	genoToPhiMap[AC][CC] = 2;
	genoToPhiMap[CG][CC] = 2;
	genoToPhiMap[CT][CC] = 2;
	genoToPhiMap[AG][GG] = 2;
	genoToPhiMap[CG][GG] = 2;
	genoToPhiMap[GT][GG] = 2;
	genoToPhiMap[AT][TT] = 2;
	genoToPhiMap[CT][TT] = 2;
	genoToPhiMap[GT][TT] = 2;

	//case aa/bb
	genoToPhiMap[AA][CC] = 3;
	genoToPhiMap[AA][GG] = 3;
	genoToPhiMap[AA][TT] = 3;
	genoToPhiMap[CC][AA] = 3;
	genoToPhiMap[CC][GG] = 3;
	genoToPhiMap[CC][TT] = 3;
	genoToPhiMap[GG][AA] = 3;
	genoToPhiMap[GG][CC] = 3;
	genoToPhiMap[GG][TT] = 3;
	genoToPhiMap[TT][AA] = 3;
	genoToPhiMap[TT][CC] = 3;
	genoToPhiMap[TT][GG] = 3;

	//case ab/ab
	genoToPhiMap[AC][AC] = 4;
	genoToPhiMap[AG][AG] = 4;
	genoToPhiMap[AT][AT] = 4;
	genoToPhiMap[CG][CG] = 4;
	genoToPhiMap[CT][CT] = 4;
	genoToPhiMap[GT][GT] = 4;

	//case ab/ac
	genoToPhiMap[AC][AG] = 5;
	genoToPhiMap[AC][AT] = 5;
	genoToPhiMap[AC][CG] = 5;
	genoToPhiMap[AC][CT] = 5;
	genoToPhiMap[AG][AC] = 5;
	genoToPhiMap[AG][AT] = 5;
	genoToPhiMap[AG][CG] = 5;
	genoToPhiMap[AG][GT] = 5;
	genoToPhiMap[AT][AC] = 5;
	genoToPhiMap[AT][AG] = 5;
	genoToPhiMap[AT][CT] = 5;
	genoToPhiMap[AT][GT] = 5;
	genoToPhiMap[CG][AC] = 5;
	genoToPhiMap[CG][CT] = 5;
	genoToPhiMap[CG][AG] = 5;
	genoToPhiMap[CG][GT] = 5;
	genoToPhiMap[CT][AC] = 5;
	genoToPhiMap[CT][CG] = 5;
	genoToPhiMap[CT][AT] = 5;
	genoToPhiMap[CT][GT] = 5;
	genoToPhiMap[GT][AG] = 5;
	genoToPhiMap[GT][CG] = 5;
	genoToPhiMap[GT][AT] = 5;
	genoToPhiMap[GT][CT] = 5;

	//case aa/bc
	genoToPhiMap[GG][AC] = 6;
	genoToPhiMap[TT][AC] = 6;
	genoToPhiMap[CC][AG] = 6;
	genoToPhiMap[TT][AG] = 6;
	genoToPhiMap[CC][AT] = 6;
	genoToPhiMap[GG][AT] = 6;
	genoToPhiMap[AA][CG] = 6;
	genoToPhiMap[TT][CG] = 6;
	genoToPhiMap[AA][CT] = 6;
	genoToPhiMap[GG][CT] = 6;
	genoToPhiMap[AA][GT] = 6;
	genoToPhiMap[CC][GT] = 6;

	//case ab/cc
	genoToPhiMap[AC][GG] = 7;
	genoToPhiMap[AC][TT] = 7;
	genoToPhiMap[AG][CC] = 7;
	genoToPhiMap[AG][TT] = 7;
	genoToPhiMap[AT][CC] = 7;
	genoToPhiMap[AT][GG] = 7;
	genoToPhiMap[CG][AA] = 7;
	genoToPhiMap[CG][TT] = 7;
	genoToPhiMap[CT][AA] = 7;
	genoToPhiMap[CT][GG] = 7;
	genoToPhiMap[GT][AA] = 7;
	genoToPhiMap[GT][CC] = 7;

	//case ab/cd
	genoToPhiMap[AC][GT] = 8;
	genoToPhiMap[AG][CT] = 8;
	genoToPhiMap[AT][CG] = 8;
	genoToPhiMap[CG][AT] = 8;
	genoToPhiMap[CT][AG] = 8;
	genoToPhiMap[GT][AC] = 8;

	/*
	phiToGenoMap = new Genotype**[9];

	//Case aa/aa
	phiToGenoMap[0] = new Genotype*[];
	phiToGenoMap[0][0] = new Genotype[2];
	phiToGenoMap[0][0][0] = AA;
	*/

};
//----------------------------------------------------
//TDistanceEstimate
//----------------------------------------------------
TDistanceEstimate::TDistanceEstimate(){
	phi = new double[9];
	old_phi = new double[9];
	for(int i=0; i<9; ++i){
		phi[i] = 0.0;
		old_phi[i] = 0.0;
	}
	LL = 0.0;
	old_LL = 0.0;
	probGeno = new double*[10];
	for(int i=0; i<10; ++i)
		probGeno[i] = new double[10];
};

void TDistanceEstimate::guessPi(int** genoQual1, int** genoQual2, long numSites){
	//just estimate pi as average posterior probability
	pi1.clear();
	pi2.clear();
	double sum1, sum2;
	int i;

	//now loop over sites
	for(long s=0; s<numSites; ++s){
		sum1 = 0.0; sum2 = 0.0;
		for(i=0; i<10; ++i){
			sum1 += phredToLik[genoQual1[s][i]];
			sum2 += phredToLik[genoQual2[s][i]];
		}
		for(i=0; i<10; ++i){
			pi1.addNoRef(genoMap.genotypeToBase[i][0], phredToLik[genoQual1[s][i]] / sum1);
			pi1.addNoRef(genoMap.genotypeToBase[i][0], phredToLik[genoQual1[s][i]] / sum1);
			pi1.addNoRef(genoMap.genotypeToBase[i][1], phredToLik[genoQual2[s][i]] / sum2);
			pi1.addNoRef(genoMap.genotypeToBase[i][1], phredToLik[genoQual2[s][i]] / sum2);
		}
	}

	//normalize
	pi1.normalize();
	pi = &pi1;
	old_pi = &pi2;
}

void TDistanceEstimate::guessPhi(int** genoQual1, int** genoQual2, long numSites){
	//set to zero
	for(int i=0; i<9; ++i)
		phi[i] = 0.0;

	//now loop over sites and add posterior probs
	double sum1, sum2;
	int i;
	for(long s=0; s<numSites; ++s){
		sum1 = 0.0; sum2 = 0.0;
		for(i=0; i<10; ++i){
			sum1 += phredToLik[genoQual1[s][i]];
			sum2 += phredToLik[genoQual2[s][i]];
		}
		for(g1 = 0; g1<10; ++g1){
			for(g2 = 0; g2<10; ++g2){
				phi[genoToPhiMap(g1,g2)] += (phredToLik[genoQual1[s][g1]] / sum1) * (phredToLik[genoQual2[s][g2]] / sum2);
			}
		}
	}

	//normalize
	sum1 = 0.0;
	for(i=0; i<9; ++i)
		sum1 += phi[i];
	for(i=0; i<9; ++i)
		phi[i] /= sum1;

}

void TDistanceEstimate::fill_K(TBaseFrequencies  & thesePi){
	//normalizing constant for each phi class
	//case of one base
	K[0] = 1.0;

	//cases of two bases
	K[1] = 0.0;
	for(int i=0; i<3; ++i){
		for(int j=i; j<4; ++j)
			K[1] += thesePi[i] * thesePi[j];
	}
	K[2] = 2.0 * K[1]; //account for AC vs CA
	K[3] = K[1];
	K[4] = K[1];

	//cases of three bases
	K[5] = thesePi[1] * thesePi[2] * thesePi[3]
		 + thesePi[0] * thesePi[2] * thesePi[3]
		 + thesePi[0] * thesePi[1] * thesePi[3]
		 + thesePi[0] * thesePi[1] * thesePi[2];
	K[5] = 3.0 * K[5]; //account for ways to distribute

	K[6] = K[5];
	K[7] = 2.0 * K[5]; //twice as many cases than other cases with three bases!

	//case of four bases: each of the 6 cases is equally likely
	K[8] = 1.0 / 6.0;
}

void TDistanceEstimate::fill_P_g_given_phi_pi(double* thesePhi, TBaseFrequencies & thesePi){
	//case aa/aa
	probGeno[AA][AA] = thesePhi[0] * thesePi[A];
	probGeno[CC][CC] = thesePhi[0] * thesePi[C];
	probGeno[GG][GG] = thesePhi[0] * thesePi[G];
	probGeno[TT][TT] = thesePhi[0] * thesePi[T];

	//cases aa/ab
	double tmp = thesePhi[1] * K[1];
	double tmp2 = tmp * thesePi[A];
	probGeno[AA][AC] = tmp2 * thesePi[C];
	probGeno[AA][AG] = tmp2 * thesePi[G];
	probGeno[AA][AT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[CC][AC] = tmp2 * thesePi[A];
	probGeno[CC][CG] = tmp2 * thesePi[G];
	probGeno[CC][CT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[G];
	probGeno[GG][AG] = tmp2 * thesePi[A];
	probGeno[GG][CG] = tmp2 * thesePi[C];
	probGeno[GG][GT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[T];
	probGeno[TT][AT] = tmp2 * thesePi[A];
	probGeno[TT][CT] = tmp2 * thesePi[C];
	probGeno[TT][GT] = tmp2 * thesePi[G];

	//case ab/aa
	tmp = thesePhi[2] * K[2];
	tmp2 = tmp * thesePi[A];
	probGeno[AC][AA] = tmp2 * thesePi[C];
	probGeno[AG][AA] = tmp2 * thesePi[G];
	probGeno[AT][AA] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[AC][CC] = tmp2 * thesePi[A];
	probGeno[CG][CC] = tmp2 * thesePi[G];
	probGeno[CT][CC] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[G];
	probGeno[AG][GG] = tmp2 * thesePi[A];
	probGeno[CG][GG] = tmp2 * thesePi[C];
	probGeno[GT][GG] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[T];
	probGeno[AT][TT] = tmp2 * thesePi[A];
	probGeno[CT][TT] = tmp2 * thesePi[C];
	probGeno[GT][TT] = tmp2 * thesePi[G];

	//case aa/bb
	tmp = thesePhi[3] * K[3];
	tmp2 = tmp * thesePi[A];
	probGeno[AA][CC] = tmp2 * thesePi[C];
	probGeno[AA][GG] = tmp2 * thesePi[G];
	probGeno[AA][TT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[CC][AA] = tmp2 * thesePi[A];
	probGeno[CC][GG] = tmp2 * thesePi[G];
	probGeno[CC][TT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[G];
	probGeno[GG][AA] = tmp2 * thesePi[A];
	probGeno[GG][CC] = tmp2 * thesePi[C];
	probGeno[GG][TT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[T];
	probGeno[TT][AA] = tmp2 * thesePi[A];
	probGeno[TT][CC] = tmp2 * thesePi[C];
	probGeno[TT][GG] = tmp2 * thesePi[G];

	//case ab/ab
	tmp = thesePhi[4] * K[4];
	tmp2 = tmp * thesePi[A];
	probGeno[AC][AC] = tmp2 * thesePi[C];
	probGeno[AG][AG] = tmp2 * thesePi[G];
	probGeno[AT][AT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[CG][CG] = tmp2 * thesePi[G];
	probGeno[CT][CT] = tmp2 * thesePi[T];
	probGeno[GT][GT] = tmp * thesePi[G] * thesePi[T];

	//case ab/ac
	tmp = thesePhi[5] * K[5];
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[G];
	probGeno[AC][AG] = tmp2;
	probGeno[AC][CG] = tmp2;
	probGeno[AG][AC] = tmp2;
	probGeno[AG][CG] = tmp2;
	probGeno[CG][AC] = tmp2;
	probGeno[CG][AG] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[T];
	probGeno[AC][AT] = tmp2;
	probGeno[AC][CT] = tmp2;
	probGeno[AT][AC] = tmp2;
	probGeno[AT][CT] = tmp2;
	probGeno[CT][AC] = tmp2;
	probGeno[CT][AT] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[G] * thesePi[T];
	probGeno[AG][GT] = tmp2;
	probGeno[AG][AT] = tmp2;
	probGeno[AT][AG] = tmp2;
	probGeno[AT][GT] = tmp2;
	probGeno[GT][AG] = tmp2;
	probGeno[GT][AT] = tmp2;
	tmp2 = tmp * thesePi[C] * thesePi[G] * thesePi[T];
	probGeno[CG][CT] = tmp2;
	probGeno[CG][GT] = tmp2;
	probGeno[CT][CG] = tmp2;
	probGeno[CT][GT] = tmp2;
	probGeno[GT][CG] = tmp2;
	probGeno[GT][CT] = tmp2;


	//case aa/bc
	tmp = thesePhi[6] * K[6];
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[G];
	probGeno[AA][CG] = tmp2;
	probGeno[CC][AG] = tmp2;
	probGeno[GG][AC] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[T];
	probGeno[AA][CT] = tmp2;
	probGeno[CC][AT] = tmp2;
	probGeno[TT][AC] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[G] * thesePi[T];
	probGeno[AA][GT] = tmp2;
	probGeno[GG][AT] = tmp2;
	probGeno[TT][AG] = tmp2;
	tmp2 = tmp * thesePi[C] * thesePi[G] * thesePi[T];
	probGeno[CC][GT] = tmp2;
	probGeno[GG][CT] = tmp2;
	probGeno[TT][CG] = tmp2;

	//case ab/cc
	tmp = thesePhi[7] * K[7];
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[G];
	probGeno[AC][GG] = tmp2;
	probGeno[AG][CC] = tmp2;
	probGeno[CG][AA] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[T];
	probGeno[AC][TT] = tmp2;
	probGeno[AT][CC] = tmp2;
	probGeno[CT][AA] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[G] * thesePi[T];
	probGeno[AG][TT] = tmp2;
	probGeno[AT][GG] = tmp2;
	probGeno[GT][AA] = tmp2;
	tmp2 = tmp * thesePi[C] * thesePi[G] * thesePi[T];
	probGeno[CG][TT] = tmp2;
	probGeno[CT][GG] = tmp2;
	probGeno[GT][CC] = tmp2;

	//case ab/cd
	tmp = thesePhi[8] * K[8];
	probGeno[AC][GT] = tmp;
	probGeno[AG][CT] = tmp;
	probGeno[AT][CG] = tmp;
	probGeno[CG][AT] = tmp;
	probGeno[CT][AG] = tmp;
	probGeno[GT][AC] = tmp;
}

void TDistanceEstimate::estimatePhiWithEM(int** genoQual1, int** genoQual2, long numSites, int maxNumIterations){
	//prepare estimates
	guessPi(genoQual1, genoQual2, numSites);
	guessPhi(genoQual1, genoQual2, numSites);

	//variables


	//now run EM
	for(int iter=0; iter<maxNumIterations; ++iter){
		//swap parameters
		tmp_pi = pi;
		pi = old_pi;
		old_pi = tmp_pi;
		tmp_phi = phi;
		phi = old_phi;
		old_phi = tmp_phi;

		//calculate P(g|phi, pi)
		fill_P_g_given_phi_pi(old_phi, *old_pi);

		//loop across loci


	}

	//clean up
	delete[] probGeno;
};

//----------------------------------------------------
//TDistanceEstimator
//----------------------------------------------------
TDistanceEstimator::TDistanceEstimator(TLog* Logfile){
	logfile = Logfile;
}

void TDistanceEstimator::printGLF(TParameters & params){

	//test first to parse GLF files
	std::string glf = params.getParameterString("glf");
	TGlfReader reader(glf);

	//print file
	reader.printToEnd();
}

//------------------------------------------------------------------
void TDistanceEstimator::estimateDistances(TParameters & params){
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

//void TDistanceEstimator::estimatePhis(int** genoQual1, int** genoQual2, long numSites){
	//variables


	//estimate starting values


	//run EM on these sites

//}






















