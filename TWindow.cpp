/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"

//-------------------------------------------------------
//TWindow_base
//-------------------------------------------------------
TWindow_base::TWindow_base(){
	start = 0;
	end = 0;
	length = 0;
	chrNumber = 0;
	sites = NULL;
	sitesInitialized = false;
	depth = 0;
	fractionSitesNoData = 0.0;
	fractionRefIsN = 0.0;
	fractionDepthAtLeastTwo = 0.0;
	numSitesWithData = 0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
	passedFilters = false;
};

TWindow_base::TWindow_base(TWindow_base & other){
	//initialize coordinates and sites
	sites = NULL;
	sitesInitialized = false;
	stealFromOther(other);
};

void TWindow_base::stealFromOther(TWindow_base & other){
	//initialize coordinates and sites
	setCoordinates(other.start, other.end, other.chrNumber);

	//copy data
	for(unsigned int i=0; i<length; ++i){
		sites[i].stealFromOther(&other.sites[i]);
	}

	depth = other.depth;
	fractionSitesNoData = other.fractionSitesNoData;
	fractionRefIsN = other.fractionRefIsN;
	fractionDepthAtLeastTwo = other.fractionDepthAtLeastTwo;
	numSitesWithData = other.numSitesWithData;
	numReadsInWindow = other.numReadsInWindow;
	referenceBaseAdded = other.referenceBaseAdded;
	passedFilters = other.passedFilters;
};

TWindow_base::TWindow_base(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	//initialize coordinates and sites
	sites = NULL;
	sitesInitialized = false;
	downsampleFromOther(other, readUpToDepth, downsamplingProb, randomGenerator);
};

void TWindow_base::downsampleFromOther(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	//set coordinates
	setCoordinates(other.start, other.end, other.chrNumber);

	//fill sites by downsampling
	numReadsInWindow = other.fillSitesDownsampling(sites, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	calcDepth();
};

void TWindow_base::downsampleFromOther(TWindow & other, TSiteSubset* subset, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	//set coordinates
	setCoordinates(other.start, other.end, other.chrNumber);

	//fill sites by downsampling
	numReadsInWindow = other.fillSitesSubsetDownsampling(sites, subset, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	calcDepth();
};

TWindow_base::~TWindow_base(){
	//delete sites
	clear();
	if(sitesInitialized)
		delete[] sites;
};

void TWindow_base::initSites(long newLength){
	if(sitesInitialized){
		clear();
		delete[] sites;
	}
	length = newLength;
	if(length > 0){
		try{
			sites = new TSite[length];;
		} catch(...){
			throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size or selecting fewer sites.";
		}
	} else {
		sites = nullptr;
	}

	sitesInitialized = true;
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionDepthAtLeastTwo = -1.0;
	numReadsInWindow = 0;
}

void TWindow_base::clear(){
	if(sitesInitialized){
		for(unsigned int i=0; i<length; ++i)
			sites[i].clear();
	}
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionRefIsN = -1.0;
	fractionDepthAtLeastTwo = -1.0;
	numSitesWithData = 0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
	passedFilters = false;
};

void TWindow_base::setCoordinates(long Start, long End, int ChrNumber){
	start = Start;
	end = End;
	chrNumber = ChrNumber;
	if(sitesInitialized){
		if((end - start) != length)
			initSites(end - start);
		else
			clear();
	} else initSites(end - start);
}

//TODO: what is difference???
void TWindow_base::move(long Start, long End, int ChrNumber){
	setCoordinates(Start, End, ChrNumber);
};

void TWindow_base::addReferenceBaseToSites(BamTools::Fasta & reference){
	if(!referenceBaseAdded){
		int stop = end - 1; //note that end is last position + 1
		std::string ref; //fasta object fills string
		reference.GetSequence(chrNumber, start, stop, ref);
		for(unsigned int i=0; i<length; ++i){
			sites[i].setRefBase(ref[i]);
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::addReferenceBaseToSites(TSiteSubset* subset){
	if(!referenceBaseAdded){
		if(subset->hasPositionsInWindow(start)){
			//now only run over sites listed in that window
			std::map<unsigned int,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
			int pos;
			for(std::map<unsigned int,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				pos = it->first - start;
				sites[pos].setRefBase(it->second.first);
			}
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::applyMask(TBedReader* mask, bool doInverseMasking){
	unsigned int first = 0;
	if(doInverseMasking){
		if(mask->hasPositionsInWindow(start)){
			std::vector<unsigned int> thesePos = mask->getPositionInWindow(start);
			for(std::vector<unsigned int>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				unsigned int pos = *it - start;
				//clear sites between regions (if there are none pos==first)
				for(unsigned int i=first; i<pos; ++i){
					if(pos < length){
						sites[i].clear();
					}
				}
				first = pos + 1;
			}
			//clear rest of window if necessary
			for(unsigned int i=first; i<end-start; ++i){
				sites[i].clear();
			}
		//else clear entire window
		} else clear();
	} else {
		if(mask->hasPositionsInWindow(start)){
			std::vector<unsigned int> thesePos = mask->getPositionInWindow(start);
			//skip sites listed in mask by setting their hasData = false
			for(std::vector<unsigned int>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				unsigned int pos = *it - start;
				if(pos < length)
					sites[pos].clear();
			}
		}
	}
};

void TWindow_base::maskCpG(){
	throw "maskCpG is not functional!";
	std::string ref; //fasta object fills string
	//note that end is last position + 1
	for(unsigned int i=0; i<length; ++i){
		if(ref[i+1] == 'C' && ref[i+2] == 'G')
			sites[i].clear();
		else if(ref[i] == 'C' && ref[i+1] == 'G')
			sites[i].clear();
	}
};

void TWindow_base::estimateBaseFrequencies(){
	//estimate initial base frequencies
	baseFreq.clear();
	for(unsigned int i=0; i<length; ++i){
		sites[i].addToBaseFrequencies(baseFreq);
	}
	baseFreq.normalize();
};


void TWindow_base::calculateEmissionProbabilities(){
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData)
			sites[i].calcEmissionProbabilities();
	}
};

void TWindow_base::call(TCaller & caller, TRecalibration & recalObject, BamTools::Fasta & reference){
	//add reference to sites
	addReferenceBaseToSites(reference);

	//loop over sites and call
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData)
			sites[i].calcEmissionProbabilities();
		caller.call(chrName, start + i, sites[i]);
	}
};

void TWindow_base::callKnwonAlleles(TCaller & caller, TRecalibration & recalObject, TSiteSubset & subset){
	//check if we need to process this window
	if(subset.hasPositionsInWindow(start)){
		//add reference to sites
		addReferenceBaseToSites(&subset);

		//now only run over sites listed in that window
		std::map<unsigned int,std::pair<char,char> > thesePos = subset.getPositionInWindow(start);
		for(std::map<unsigned int,std::pair<char,char> >::iterator it = thesePos.begin(); it!=thesePos.end(); ++it){
			int pos = it->first - start;
			if(sites[pos].hasData){
				sites[pos].calcEmissionProbabilities();

				//call
				caller.call(chrName, start + pos, sites[pos], it->second.first, it->second.second);
			}
		}
	}
};

void TWindow_base::printPileup(TRecalibration* recalObject, gz::ogzstream & out, bool printOnlySitesWithData){
	//print pileup
	for(unsigned int i=0; i<length; ++i){
		if((printOnlySitesWithData && sites[i].hasData) || !printOnlySitesWithData){
			sites[i].calcEmissionProbabilities();
			out << chrName << "\t" << start + i + 1;
			sites[i].printPileup(out);
			out << "\n";
		}
	}
};

void TWindow_base::printPileupToScreen(TRecalibration* recalObject){
	//print pileup
	for(unsigned int i=0; i<length; ++i){
		sites[i].calcEmissionProbabilities();
		std::cout << chrName << "\t" << start + i + 1;
		sites[i].printPileupToScreen();
		std::cout << "\n";
	}
};

void TWindow_base::calcDepth(){
	//calculate and return coverage
	depth = 0.0;
	long noData = 0;
	long plentyData = 0;
	int cov;
	for(unsigned int i=0; i<length; ++i){
		cov = sites[i].depth();
		depth += cov;
		if(cov == 0) ++noData;
		else if(cov > 1) ++ plentyData;
	}

	depth = depth / (double) length;
	numSitesWithData = length - noData;
	fractionSitesNoData = (double) noData / (double) length;
	fractionDepthAtLeastTwo = (double) plentyData / (double) length;
};

void TWindow_base::calcFracN(){
	double numN = 0.0;
	for(unsigned int i=0; i<length; ++i)	if(sites[i].referenceBase == 'N' || sites[i].referenceBase == 'n')
		++numN;
	fractionRefIsN = numN / (double) length;
};

void TWindow_base::calcDepthPerSite(long* siteDepth, size_t maxDepth){
	//calculate and return coverage
	depth = 0.0;
	long noData = 0;
	long plentyData = 0;
	size_t cov;

	for(unsigned int i=0; i<length; ++i){
		cov = sites[i].depth();
		if(cov <= maxDepth)	siteDepth[cov] += 1;
		else siteDepth[maxDepth + 1] += 1; //else it should be in the "greater than" bin

		if(cov == 0) ++ noData;
		else if(cov > 1) ++ plentyData;
	}

	depth = depth / (double) length;
	fractionSitesNoData = (double) noData / (double) length;
	fractionDepthAtLeastTwo = (double) plentyData / (double) length;
};

void TWindow_base::countDepthPerSite(TDistributionOfCounts & counts){
	for(unsigned int i=0; i<length; ++i)
		counts.add(sites[i].depth());
};

void TWindow_base::printDepthPerSite(gz::ogzstream & out, const std::string & chr){
	//print depth for each site to file
	for(unsigned int i=0; i<length; ++i){
		out << chr << "\t" << start + i + 1 << "\t" << sites[i].depth() << "\n";
	}
};

void TWindow_base::printMateInformationPerSite(TOutputFileZipped & out, const std::string & chr){
	int* alleleCounts = new int[4];
	int* mateCounts = new int[2];
	int* frCounts = new int[2];

	for(unsigned int i=0; i<length; ++i){
		//chr, pos and depth
		out << chr << start + i + 1 << sites[i].depth();

		//num alleles
		sites[i].countAlleles(alleleCounts);
		int numAlleles = 0;
		for(int b=0; b<4; ++b){
			out << alleleCounts[b];
			if(alleleCounts[b] > 0)
				++numAlleles;
		}
		out << numAlleles;

		//mates
		sites[i].countMates(mateCounts);
		out << mateCounts[0] << mateCounts[1];

		//fwd / rev
		sites[i].countFwdRev(frCounts);
		out << frCounts[0] << frCounts[1];
		out << std::endl;
	}

	delete[] alleleCounts;
	delete[] mateCounts;
	delete[] frCounts;
};

void TWindow_base::countAlleles(TAllelicDepthCounts & counts){
	//calculate and return imbalance
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].depth() > 0){
			sites[i].countAllelesForImbalance(counts);
		} else {
			counts.addSiteZeroDepth();
		}

	}
};

void TWindow_base::writeNonConservedBed(std::ofstream & output){
	//calculate and return imbalance
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].depth() > 0 && sites[i].refDepth() < sites[i].depth()){
			output << chrName << "\t" << start + i << "\t" << start + i + 1 << "\n";
		}
	}
};

void TWindow_base::contextStats(int** contextCounts, TQualityMap & qualMap){
	for(unsigned int i=0; i<length; ++i){
		sites[i].contextStats(contextCounts, qualMap);
	}
};

void TWindow_base::applyDepthFilter(const size_t minDepth, const size_t maxDepth){
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData){
			if(sites[i].bases.size() < minDepth || sites[i].bases.size() > maxDepth)
				sites[i].clear();
		}
	}
};

void TWindow_base::createDepthMask(size_t minDepthForMask, size_t maxDepthForMask, std::ofstream & outputMaskFile, const std::string & chr){
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData){
			if(sites[i].bases.size() < minDepthForMask || sites[i].bases.size() > maxDepthForMask){
				outputMaskFile << chr << "\t" << start + i << "\t" << start + i + 1 << "\n";
			}
		}
	}
};

void TWindow_base::addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TLog* logfile, TQualityMap & qualMap){
	logfile->listFlush("Adding sites to BQSR ...");
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData){
			bqsr.addSite(sites[i], qualMap);
		}
	}
	logfile->done();
};

void TWindow_base::addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TSiteSubset* subset, TLog* logfile, TQualityMap & qualMap){
	logfile->listFlush("Adding sites to BQSR ...");
	//now only run over sites listed in that window
	std::map<unsigned int,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
	int pos;
	for(std::map<unsigned int,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
		pos = it->first - start;
		if(sites[pos].hasData){
			sites[pos].setRefBase(it->second.second);
			bqsr.addSite(sites[pos], qualMap);
		}
	}
	logfile->done();

};

void TWindow_base::addSitesToPMDTable(TPMDTables & pmdTables, TLog* logfile){
	logfile->listFlush("Adding sites to PMD tables ...");
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData){
			//pmdTables.add(sites[i]);
		}
	}
	logfile->done();
};

void TWindow_base::addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer){
	//assumes that emission probabilities were calculated
	for(unsigned int i=0; i<length; ++i){
		sites[i].calcEmissionProbabilities();
		thetaDataContainer->add(sites[i]);
	}
};

void TWindow_base::addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TBedReader & region){
	//assumes that emission probabilities were calculated
	//only add sites from regions
	if(region.hasPositionsInWindow(start)){
		std::vector<unsigned int> thesePos = region.getPositionInWindow(start);
		for(std::vector<unsigned int>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
			unsigned int pos = *it - start;
			if(pos < length){
				sites[pos].calcEmissionProbabilities();
				thetaDataContainer->add(sites[pos]);
			}
		}
	}
};

void TWindow_base::addToGLF(TGlfWriter & writer, const int ploidy, bool printAll){
	//TODO: calculate root mean squared mapping qualities for sites (now just passing 0). Would be helpful in VCFs as well
	if(printAll){
		for(unsigned int i=0; i<length; ++i){
			writer.writeSite(start + i, sites[i].depth(), 0, sites[i].emissionProbabilities);
		}
	} else {
		for(unsigned int i=0; i<length; ++i){
			if(sites[i].hasData){
				writer.writeSite(start + i, sites[i].depth(), 0, sites[i].emissionProbabilities);
			}
		}
	}
};

void TWindow_base::generatePSMCInput(TThetaEstimator & estimator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine){
	//calc prior probabilities on Genotypes
	double* pGenotype = new double[10];
	estimator.fillPGenotype(pGenotype);

	//now call heterozygosity in blocks
	int nBlocks = length / blockSize;
	int start;
	double logPHomo;
	double logConfidence = log(confidence);
	double logConfidenceHet = log(1.0 - confidence);
	double tmp;

	//loop over blocks
	for(int b=0; b<nBlocks; ++b){
		start = b*blockSize;
		logPHomo = 0.0;

		for(int i=0; i<blockSize; ++i){
			if(sites[start + i].hasData){
				tmp = sites[start + i].calculatePHomozygous(pGenotype);
				logPHomo += log(tmp);
			}
		}

		//check if we are heterozygous
		if(logPHomo > logConfidence){
			out << 'T';
		} else if(logPHomo < logConfidenceHet){
			out << 'K';
		} else {
			out << 'N';
		}

		//do we add a new line?
		if(nCharOnLine == 59){
			nCharOnLine = 0;
			out << '\n';
		} else ++nCharOnLine;
	}
	delete[] pGenotype;
};

void TWindow_base::addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TQualityMap & qualMap){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject.addSite(sites[i], qualMap);
		}
	}
};

void TWindow_base::addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TSiteSubset* subset, TQualityMap & qualMap){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	//now only run over sites listed in that window
	std::map<unsigned int,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
	int pos;
	for(std::map<unsigned int,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
		pos = it->first - start;
		if(sites[pos].hasData && it->second.first == it->second.second){
			recalObject.addSite(sites[pos], qualMap, genoMap.getBase(it->second.first));
		}
	}
};

//-------------------------------------------------------
//Twindow
//-------------------------------------------------------
TWindow::TWindow():TWindow_base(){};

TWindow::~TWindow(){
	//delete alignments
	for(std::vector<TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end(); ++alignmentIt)
		delete *alignmentIt;
	usedAlignments.clear();

	for(std::vector<TAlignment*>::iterator alignmentIt=emptyAlignments.begin(); alignmentIt != emptyAlignments.end(); ++alignmentIt)
		delete *alignmentIt;
	emptyAlignments.clear();

};

TAlignment* TWindow::swapUsedForEmptyAlignment(TAlignment* usedAlignment, const unsigned int & maxReadLength){
	//save used alignment on proper stack
	usedAlignments.push_back(usedAlignment);

	//return empty alignment, either from stack or create new
	if(emptyAlignments.size() > 0){
		TAlignment* alignment = *(emptyAlignments.rbegin());
		emptyAlignments.pop_back();
		return alignment;
	} else {
		TAlignment* alignment = new TAlignment(maxReadLength);
		return alignment;
	}
};

void TWindow::review(){
	//update pointers
	/*
	firstAlignmentwithPosOutsideWindow = usedAlignments.end()-1;
	while((*firstAlignmentwithPosOutsideWindow)->position > end && firstAlignmentwithPosOutsideWindow != usedAlignments.begin())
		--firstAlignmentwithPosOutsideWindow;
		*/

	//fillSites();
	//calcDepth();
};

void TWindow::cleanUpUsedAlignments(){
	//now check and move the rest
	for(std::vector<TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end();){
		if((*alignmentIt)->position < end && (*alignmentIt)->lastAlignedPositionWithRespectToRef >= start && (*alignmentIt)->chrNumber == chrNumber){
			++alignmentIt;
		} else{
			(*alignmentIt)->clear();
			emptyAlignments.push_back(*alignmentIt);
			alignmentIt = usedAlignments.erase(alignmentIt);
		}
	}
};

void TWindow::clearAllUsedAlignments(){
	for(std::vector<TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end();){
		(*alignmentIt)->clear();
		emptyAlignments.push_back(*alignmentIt);
		alignmentIt = usedAlignments.erase(alignmentIt);
	}
};

void TWindow::move(unsigned int Start, unsigned int End, int ChrNumber){
	setCoordinates(Start, End, ChrNumber);
	cleanUpUsedAlignments();
};

void TWindow::printStacks(){
	std::cout << "USED ALIGMENTS:";
	for(TAlignment* alignmentIt : usedAlignments)
		std::cout << " " << alignmentIt << " : " << alignmentIt->alignmentName << " pos " << alignmentIt->position;
	std::cout << std::endl;

	std::cout << "EMPTY ALIGMENTS:";
	for(std::vector<TAlignment*>::iterator alignmentIt=emptyAlignments.begin(); alignmentIt != emptyAlignments.end(); ++alignmentIt)
		std::cout << " " << *alignmentIt;
	std::cout << std::endl;
};

void TWindow::checkAlignmentForFillingSites(TAlignment* alignmentIt){
	//check if alignment start is inside window
	if(alignmentIt->position >= end){
		throw "alignment should be assigned to next window!";
	}
};

void TWindow::setFirstPositionWithinWindow(TAlignment* alignmentIt, unsigned int & firstPos, unsigned int & p){
	//genomic position of alignment as seen from window perspective
	firstPos = alignmentIt->position - start;

	//set position in read
	p = 0;

	//is the beginning of the read part of previous window? increase starting p for adding bases!
	if(firstPos < 0){
		while(p < alignmentIt->length && (firstPos + alignmentIt->bases[p].alignedPos) < 0){
			++p;
		} if(p == alignmentIt->length){
			throw "alignment should be assigned to previous window! Name: " + alignmentIt->alignmentName + ". In window " + toString(start) + "-" + toString(end) + ". with position " + toString(alignmentIt->position);
		}
	}
};

//------------------------------------------------------
//fill sites
//------------------------------------------------------
void TWindow::fillSites(TAlignment* alignmentIt, TSite* theseSites, const unsigned int & readUpToDepth){
	//genomic position of alignment as seen from window perspective
	unsigned int firstPos, p;
	setFirstPositionWithinWindow(alignmentIt, firstPos, p);

	//position in window where first one = 0
	//p is at first position of read in window
	for(; p < alignmentIt->length; ++p){
		if(alignmentIt->bases[p].aligned && alignmentIt->bases[p].base != N){
			unsigned int internalPos = firstPos + alignmentIt->bases[p].alignedPos;
			//if read extends past window length
			if(internalPos >= length){
				break; //since part of the read maps to next window
			}
			if(theseSites[internalPos].depth() < readUpToDepth){
				theseSites[internalPos].add(&alignmentIt->bases[p]);
			}
		}
	}
};

int TWindow::fillSites(TSite* theseSites, const unsigned int & readUpToDepth){
	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(TAlignment* alignmentIt : usedAlignments){
		//now fill
		fillSites(alignmentIt, theseSites, readUpToDepth);
		++counter;
	}
	return counter;
};

int TWindow::fillSitesDownsampling(TSite* theseSites, const unsigned int & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator){
	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(TAlignment* alignmentIt : usedAlignments){
		//fill if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			fillSites(alignmentIt, theseSites, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

void TWindow::fillSites(const unsigned int & readUpToDepth){
	numReadsInWindow = fillSites(sites, readUpToDepth);
};

//------------------------------------------------------
//fill sites according to subset
//------------------------------------------------------
void TWindow::fillSitesSubset(TAlignment* alignmentIt, TSite* theseSite, std::map<unsigned int,std::pair<char,char> > & thesePos, const unsigned int & readUpToDepth){
	//genomic position of alignment as seen from window perspective
	unsigned int firstPos, p;
	setFirstPositionWithinWindow(alignmentIt, firstPos, p);

	//position in window where first one = 0
	//p is at first position of read in window
	for(; p < alignmentIt->length; ++p){
		if(alignmentIt->bases[p].alignedPos && alignmentIt->bases[p].base != N){
			unsigned int internalPos = firstPos + alignmentIt->bases[p].alignedPos;
			//if read extends past window length
			if(internalPos >= length)
				break; //since part of the read maps to next window
			if(thesePos.find(internalPos) != thesePos.end() && sites[internalPos].depth() < readUpToDepth)
				sites[internalPos].add(&alignmentIt->bases[p]);
		}
	}
};

int TWindow::fillSitesSubset(TSite* theseSites, TSiteSubset* subset, const unsigned int & readUpToDepth){
	//get positions that are used
	std::map<unsigned int,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);

	//add reads in usedAlignments to sites in window
	unsigned int firstPos, p;
	int counter = 0;
	for(TAlignment* alignmentIt : usedAlignments){
		//genomic position of alignment as seen from window perspective
		setFirstPositionWithinWindow(alignmentIt, firstPos, p);

		//now fill
		fillSitesSubset(alignmentIt, theseSites, thesePos, readUpToDepth);
		++counter;
	}
	return counter;
};

int TWindow::fillSitesSubsetDownsampling(TSite* theseSites, TSiteSubset* subset, const unsigned int & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator){
	//get positions that are used
	std::map<unsigned int,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(TAlignment* alignmentIt : usedAlignments){
		//check if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			//now fill
			fillSitesSubset(alignmentIt, theseSites, thesePos, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

void TWindow::fillSitesSubset(TSiteSubset* subset, const unsigned int & readUpToDepth){
	numReadsInWindow = fillSitesSubset(sites, subset, readUpToDepth);
};
