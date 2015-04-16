/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */


#include "TGenome.h"

//-------------------------------------------------------
//TBase
//-------------------------------------------------------
double TBase::probPDM(int & pos, Constants & constants){
	//probability of a Post Mortem Damag (according to Skoglund et al. 2014)
	return constants.pdmC + pow(1.0 - constants.pdmLambda, (double) pos - constants.pdmLambda) * constants.pdmLambda;
}

void TBaseA::fillEmissionProbabilities(Constants & constants){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pdm = probPDM(pos3, constants);

	emissionProbabilities.set(AA, oneMinusError);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pdm) * errorOneThird + pdm * oneMinusError);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, 0.5 - errorOneThird);
	emissionProbabilities.set(AG, ((1.0 + pdm) * oneMinusError + (1.0 - pdm) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, (pdm * oneMinusError + (2.0 - pdm) * errorOneThird) / 2.0);
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(CG));
};

void TBaseC::fillEmissionProbabilities(Constants & constants){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pdm = probPDM(pos5, constants);

	emissionProbabilities.set(AA , errorOneThird);
	emissionProbabilities.set(CC , (1.0 - pdm) * oneMinusError + pdm * errorOneThird);
	emissionProbabilities.set(GG , errorOneThird);
	emissionProbabilities.set(TT , errorOneThird);
	emissionProbabilities.set(AC , ((1.0 - pdm) * oneMinusError + (1.0 + pdm) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG , errorOneThird);
	emissionProbabilities.set(AT , errorOneThird);
	emissionProbabilities.set(CG , emissionProbabilities.get(AC));
	emissionProbabilities.set(CT , emissionProbabilities.get(AC));
	emissionProbabilities.set(GT , errorOneThird);
};

void TBaseG::fillEmissionProbabilities(Constants & constants){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pdm = probPDM(pos3, constants);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pdm) * oneMinusError + pdm * errorOneThird);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, errorOneThird);
	emissionProbabilities.set(AG, ((1.0 - pdm) * oneMinusError + (1.0 + pdm) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AG));
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(AG));

};

void TBaseT::fillEmissionProbabilities(Constants & constants){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pdm = probPDM(pos5, constants);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, (1.0 - pdm) * errorOneThird + pdm * oneMinusError);
	emissionProbabilities.set(GG, errorOneThird);
	emissionProbabilities.set(TT, oneMinusError);
	emissionProbabilities.set(AC, (pdm * oneMinusError + (2.0 - pdm) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG, errorOneThird);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AC));
	emissionProbabilities.set(CT, ((1.0 + pdm) * oneMinusError + (1.0 - pdm) * errorOneThird) / 2.0);
	emissionProbabilities.set(GT, 0.5 - errorOneThird);
};
//-------------------------------------------------------
//TSite
//-------------------------------------------------------
void TSite::add(char & base, char & quality, int pos5, int pos3, Constants & constants){
	if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
		double error = pow(10.0, - ((double)quality - 33.0)/10.0);
		if(error < 1.0){
			if(base == 'A') bases.push_back(new TBaseA(error, pos5, pos3, constants));
			else if(base == 'C') bases.push_back(new TBaseC(error, pos5, pos3, constants));
			else if(base == 'G') bases.push_back(new TBaseG(error, pos5, pos3, constants));
			else bases.push_back(new TBaseT(error, pos5, pos3, constants));
		}
	}
};

std::string TSite::getBases(){
	std::string b = "";
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		b += (*it)->getBase();
	}
	return b;
}

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	double weight = 1.0 / bases.size();
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		(*it)->addToBaseFrequencies(frequencies, weight);
	}
}

void TSite::calcEmissionProbabilities(){
	//calculate normalized genotype probabilities according to Bayes rule
	for(int i=0; i<10; ++i){
		emissionProbabilities[i] = 1.0;
		for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
			emissionProbabilities[i] *= (*it)->getEmissionProbability(i);
		}
	}

}

std::string TSite::getEmissionProbs(){
	std::string b = toString(emissionProbabilities[0]);
	for(int i=1; i<10; ++i){
		b += ", " + toString(emissionProbabilities[i]);
	}
	return b;
}

void TSite::calculateP_g(double* genotypeProbabilities){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<10; ++i){
		P_g[i] =  emissionProbabilities[i] * genotypeProbabilities[i];
		sum += P_g[i];
	}
	for(int i=0; i<10; ++i){
		P_g[i] /= sum;
	}
}
//-------------------------------------------------------
//Twindow
//-------------------------------------------------------
TWindow::TWindow(Constants* SharedConstants){
	sharedConstants = SharedConstants;
	start = -1;
	end = -1;
	length = -1;
	sites = NULL;
	sitesInitialized = false;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
	theta = 0.0;
};

TWindow::TWindow(long Start, long End, Constants* SharedConstants){
	sharedConstants = SharedConstants;
	start = Start;
	end = End;
	initSites(end - start); //end NOT in window!
};

void TWindow::initSites(long newLength){
	if(sitesInitialized)
		delete[] sites;
	length = newLength;
	sites = new TSite[length];
	sitesInitialized = true;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
}

void TWindow::clear(){
	for(int i=0; i<length; ++i) sites[i].clear();
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
};

void TWindow::move(long Start, long End){
	start = Start;
	end = End;
	if(sitesInitialized){
		if((start - end) != length)
			initSites(end - start);
		else  clear();
	} else initSites(end - start);
};

bool TWindow::addFromRead(BamTools::BamAlignment & bamAlignement){
	/* Note:
	 * Function returns true if read also maps to next window and
	 * returns false if end of read is within this (or a previous) window
	 */
	if(bamAlignement.Position >= end) return true;

	//find first position to be within window
	double len = bamAlignement.AlignedBases.length();
	if(bamAlignement.Position + len < start) return false;

	//find which position to consider first
	int firstPos = start - bamAlignement.Position;
	if(firstPos < 0) firstPos = 0;
	int lastPos = len;
	if(bamAlignement.Position + lastPos > end) lastPos = end - bamAlignement.Position;

	//add sites
	int internalPos = bamAlignement.Position + firstPos - start;

	for(int pos = firstPos; pos < lastPos; ++pos, ++internalPos){
		/* Note:
		 * Reference is 5' -> 3'
		 * Hence for any read mapping "forward":
		 * 		- distance from 5' = pos + 1
		 * 		- distance from 3' = len - pos
		 * For any read mapping "reverse":
		 * 		- distance from 5' = len - pos
		 * 		- distance from 3' = pos + 1
		 * In both cases: 1) distance is 1-based!
		 *                2) Ignoring indels when calculating distances
		 *                3) Function add need first pos5, then pos3
		 */

		if(bamAlignement.IsReverseStrand())
			sites[internalPos].add(bamAlignement.AlignedBases.at(pos), bamAlignement.AlignedQualities.at(pos), len - pos, pos + 1, *sharedConstants);
		else
			sites[internalPos].add(bamAlignement.AlignedBases.at(pos), bamAlignement.AlignedQualities.at(pos), pos + 1, len - pos, *sharedConstants);
	}

	//return if part of the read maps to next window
	if(lastPos == len) return false;
	else return true;
}

void TWindow::runEM(){
	//estimate initial base frequencies
	//calculate per site emission probabilities
	for(int i=0; i<length; ++i){
		sites[i].calcEmissionProbabilities();
		sites[i].addToBaseFrequencies(baseFreq);
	}
	baseFreq.normalize();
	//baseFreq.print();

	//set initial parameters
	theta = 0.01;
	double expTheta = exp(-theta);
	double rho = expTheta / (1.0 - expTheta);
	double mu = length;
	double oldTheta;

	//prepare storage
	double pGenotype[10];
	double P_G[10];
	double tmp[4];
	double tmpSum;
	arma::mat Jacobian(6,6);
	arma::vec F(6);
	arma::mat JxF(6,6);
	Genotype geno;
	double maxF;

	//start EM loop
	for(int iter = 0; iter < sharedConstants->numIterations; ++iter){
		//b) calculate	substitution probabilities
		for(int i=0; i<4; ++i){
			//homozygous genotypes
			pGenotype[sharedConstants->getGenotype(i,i)] = baseFreq[i] * (expTheta + baseFreq[i] * (1.0 - expTheta));
			//heterozygous genotypes
			for(int j=i+1; j<4; ++j){
				pGenotype[sharedConstants->getGenotype(i,j)] = 2.0 * baseFreq[i] * baseFreq[j] *  (1.0 - expTheta);
			}
		}

		//c) Calculate all genotype probabilities for all sites
		//calculate P_g for each site
		for(int i=0; i<length; ++i){
			sites[i].calculateP_g(pGenotype);
		}
		//calculate P_g
		for(int g=0; g<10; ++g){
			P_G[g] = 0.0;
			for(int i=0; i<length; ++i){
				P_G[g] += sites[i].P_g[g];
			}
		}

		//get sensible starting point for rho in first iteration
		if(iter == 0){
			double sumPKK = 0.0;
			for(int i=0; i<4; ++i)
				sumPKK += P_G[sharedConstants->getGenotype(i, i)];
			rho = (sumPKK) / (length - sumPKK);
		}

		//d) Find new parameter estimates using Newton-Ralphson
		for(int n=0; n<sharedConstants->NewtonRalphsonNumIterations; ++n){
			//i) calculate F (Note: index is zero based!)
			F(4) = length;
			F(5) = 0.0;
			for(int k=0; k<4; ++k){
				geno = sharedConstants->getGenotype(k, k);
				tmpSum = 0.0;
				for(int l=0; l<4; ++l){
					if(l != k){
						tmpSum += P_G[sharedConstants->getGenotype(k, l)];
					}
				}
				F(k) = P_G[geno] * (1.0 + baseFreq[k] / (rho + baseFreq[k])) + tmpSum - mu * baseFreq[k];
				F(4) -= P_G[geno] * (rho + 1.0 ) / (rho + baseFreq[k]);
				F(5) += baseFreq[k];
			}
			F(5) = F(5) - 1.0;

			//ii) fill Jacobian (Note: index is zero based!)
			Jacobian.zeros();
			tmpSum = 0.0;
			for(int k=0; k<4; ++k){
				tmp[k] = P_G[sharedConstants->getGenotype(k, k)] / ((baseFreq[k] + rho)*(baseFreq[k] + rho));
				tmpSum += tmp[k];
			}

			for(int k=0; k<4; ++k){
				Jacobian(k,k) = tmp[k] * rho - mu;
				Jacobian(k,4) = - tmp[k];
				Jacobian(5,k) = 1.0;
				Jacobian(4,k) = tmp[k] * (rho + 1.0);
				Jacobian(k,5) = - baseFreq[k];
				Jacobian(4,4) += tmp[k] * (1.0 - baseFreq[k]);
			}

			//iii) now estimate new parameters
			JxF = inv(Jacobian) * F;
			for(int k=0; k<4; ++k){
				baseFreq[k] -= JxF(k);
			}
			rho -= JxF(4);
			mu -= JxF(5);

			//check if we break
			maxF = 0.0;
			for(int i=0; i<6; ++i){
				if(F(i) > maxF) maxF = F(i);
			}
			if(maxF < sharedConstants->NewtonRalphsonMaxF) break;

		}

		//check if we stop EM
		oldTheta = theta;
		theta = -log(rho / (1.0 + rho));
		if(fabs(theta - oldTheta) < sharedConstants->maxEpsilon) break;
	}
}

void TWindow::writeEMResults(std::ofstream & out){
	//position
	out << start << "\t" << end;
	//coverage NOTE: assumes coverage has been calculated before...
	out << "\t" << coverage << "\t" << fractionSitesNoData << "\t" << fractionsitesCoverageAtLeastTwo;
	//estimated params
	for(int i=0; i<4; ++i)
		out << "\t" << baseFreq[i];
	out << "\t" << theta << std::endl;
}


void TWindow::printPileup(){
	for(int i=0; i<length; ++i){
		std::cout << start + i << "\t" << sites[i].bases.size() << "\t" << sites[i].getBases() << "\n";
	}
}

void TWindow::calcCoverage(){
	//calculate and return coverage
	coverage = 0.0;
	long noData = 0;
	long plentyData = 0;
	for(int i=0; i<length; ++i){
		coverage += sites[i].bases.size();
		if(sites[i].bases.size() == 0) ++ noData;
		else if(sites[i].bases.size() > 1) ++ plentyData;
	}
	coverage = coverage / (double) length;
	fractionSitesNoData = (double) noData / (double) length;;
	fractionsitesCoverageAtLeastTwo = (double) plentyData / (double) length;
}

//-------------------------------------------------------
//TGenome
//-------------------------------------------------------
TGenome::TGenome(TLog* Logfile, TParameters & params){
	logfile = Logfile;

	//read parameters
	filename = params.getParameterString("bam");
	windowSize = params.getParameterDouble("window");
	if(windowSize < 1000) throw "Window size should be at least 1Kb!";
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 0.9);

	//post mortem damage
	sharedConstants.pdmC = params.getParameterDoubleWithDefault("pdmC", 0.01);
	sharedConstants.pdmLambda = params.getParameterDoubleWithDefault("pdmLambda", 0.3);

	//EM params
	sharedConstants.numIterations = params.getParameterIntWithDefault("iterations", 100);
	sharedConstants.maxEpsilon = params.getParameterDoubleWithDefault("maxEps", 0.00001);
	sharedConstants.NewtonRalphsonNumIterations = params.getParameterIntWithDefault("NRiterations", 20);
	sharedConstants.NewtonRalphsonMaxF = params.getParameterDoubleWithDefault("maxF", 0.000001);

	//output
	if(params.parameterExists("out"))
		outputName = params.getParameterString("out");
	else {
		//guess from filename
		outputName = filename;
		outputName = extractBeforeLast(outputName, ".");
		outputName += "_theta_estimates.txt";
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

	//prepare windows
	curWindow = new TWindow(&sharedConstants);
	nextWindow = new TWindow(&sharedConstants);
};

bool TGenome::iterateChromosome(){
	if(chrIterator == bamHeader.Sequences.End()){
		chrIterator = bamHeader.Sequences.Begin();
		chrNumber = 0;
	} else {
		logfile->endIndent();
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
	nextWindow->move(0, windowSize);

	//write progress
	logfile->startIndent("Parsing chromosome '" + chrIterator->Name + "':");
	return true;
}

bool TGenome::iterateWindow(){
	//move to next region
	curStart = curEnd;
	curEnd += windowSize;
	if(curEnd > chrLength) curEnd = chrLength + 1;
	if(curStart > chrLength) return false;

	//swap windows and move next
	TWindow* tmp = curWindow;
	curWindow = nextWindow;
	nextWindow = tmp;
	nextWindow->move(curEnd, curEnd + windowSize);

	return true;
};

bool TGenome::readData(){
	logfile->listFlush("Reading data on '" + chrIterator->Name + "' at [" + toString(curStart) + ", " + toString(curEnd) + ") ...");

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
				if(curWindow->addFromRead(bamAlignement))
					//add also to next window in case reads overhangs current window -> function returns true
					nextWindow->addFromRead(bamAlignement);

			}
		}
	}

	curWindow->calcCoverage();

	logfile->write(" done!");
	logfile->conclude("read data from " + toString(numReads) + " reads.");
	logfile->conclude("coverage is " + toString(curWindow->coverage));
	logfile->conclude(toString(curWindow->fractionsitesCoverageAtLeastTwo) + " of all sites are covered at least twice");
	logfile->conclude(toString(curWindow->fractionSitesNoData) + " of all sites have no data");

	return true;
};

void TGenome::estimateTheta(){
	//write header
	std::ofstream out;
	out.open(outputName.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	out << "Chr\tstart\tend\tcoverage\tmissing\ttwoOrMore\tf(A)\tf(C)\tf(G)\tf(T)\ttheta" << std::endl;

	//iterate through windows
	while(iterateChromosome()){
		while(iterateWindow()){
			//read data for current window
			readData();

			//check if we have data -> can be extended to ensure
			if(curWindow->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				//run EM
				logfile->listFlush("Performing EM estimation of theta ...");
				curWindow->runEM();
				logfile->write(" done!");
				logfile->conclude("theta was estimated at " + toString(curWindow->theta));

				//save results to file
				out << chrIterator->Name << "\t";
				curWindow->writeEMResults(out);
			}
		}
	}

	//close output
	out.close();
}


