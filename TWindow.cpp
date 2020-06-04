/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"

using namespace GenotypeLikelihoods;

//-------------------------------------------------------
//TWindow_base
//-------------------------------------------------------
TWindow_base::TWindow_base(){
	startPos = 0;
	endPos = 0;
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

/*
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
*/

TWindow_base::TWindow_base(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	//initialize coordinates and sites
	sites = NULL;
	sitesInitialized = false;
	downsampleFromOther(other, readUpToDepth, downsamplingProb, randomGenerator);
};

void TWindow_base::downsampleFromOther(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	//set coordinates
	setCoordinates(other.startPos, other.endPos, other.chrNumber);

	//fill sites by downsampling
	numReadsInWindow = other.fillSitesDownsampling(sites, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	calcDepth();
};

void TWindow_base::downsampleFromOther(TWindow & other, TSiteSubset & subset, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	//set coordinates
	setCoordinates(other.startPos, other.endPos, other.chrNumber);

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
	length = newLength;
	if(length > 0){
		try{
			sites.resize(length);
		} catch(...){
			throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size or selecting fewer sites.";
		}
	} else {
		sites.clear();
	}

	sitesInitialized = true;
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionDepthAtLeastTwo = -1.0;
	numReadsInWindow = 0;
};

void TWindow_base::clear(){
	for(auto& s : sites){
		s.clear();
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
	startPos = Start;
	endPos = End;
	chrNumber = ChrNumber;
	if(sitesInitialized){
		if((endPos - startPos) != length)
			initSites(endPos - startPos);
		else
			clear();
	} else initSites(endPos - startPos);
};

void TWindow_base::move(unsigned int Start, unsigned int End, int ChrNumber, TLog* logfile){
	setCoordinates(Start, End, ChrNumber);
};

void TWindow_base::addReferenceBaseToSites(BAM::TFastaBuffer & reference){
	if(!referenceBaseAdded && reference.hasReference()){
		int stop = endPos - 1; //note that end is last position + 1
		std::string ref; //fasta object fills string
		reference.fill(chrNumber, startPos, stop, ref);
		for(unsigned int i=0; i<length; ++i){
			sites[i].setRefBase(genoMap.getBase(ref[i]));
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::addReferenceBaseToSites(TSiteSubset & subset){
	if(!referenceBaseAdded){
		if(subset.hasPositionsInWindow(startPos)){
			//now only run over sites listed in that window
			std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(startPos);
			int pos;
			for(auto& it : thesePositions){
				pos = it.position - startPos;
				sites[pos].setRefBase(it.ref);
			}
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::applyMask(BAM::TBedReader* mask, bool doInverseMasking){
	unsigned int first = 0;
	if(doInverseMasking){
		if(mask->hasPositionsInWindow(startPos)){
			std::vector<unsigned int> thesePos = mask->getPositionInWindow(startPos);
			for(std::vector<unsigned int>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				unsigned int pos = *it - startPos;
				//clear sites between regions (if there are none pos==first)
				for(unsigned int i=first; i<pos; ++i){
					if(pos < length){
						sites[i].clear();
					}
				}
				first = pos + 1;
			}
			//clear rest of window if necessary
			for(unsigned int i=first; i<endPos-startPos; ++i){
				sites[i].clear();
			}
		//else clear entire window
		} else clear();
	} else {
		if(mask->hasPositionsInWindow(startPos)){
			std::vector<unsigned int> thesePos = mask->getPositionInWindow(startPos);
			//skip sites listed in mask by setting their hasData = false
			for(std::vector<unsigned int>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				unsigned int pos = *it - startPos;
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

/*
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
*/

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
			output << _chrName << "\t" << startPos + i << "\t" << startPos + i + 1 << "\n";
		}
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
				outputMaskFile << chr << "\t" << startPos + i << "\t" << startPos + i + 1 << "\n";
			}
		}
	}
};

void TWindow_base::addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TGenotypeLikelihoodCalculator & GL_calculator){
	//assumes that emission probabilities were calculated
	GenotypeLikelihoods::TGenotypeLikelihoods genoLik;
	for(unsigned int i=0; i<length; ++i){
		GL_calculator.calculateGenotypeLikelihoods(sites[i].bases, genoLik);
		thetaDataContainer->add(sites[i], genoLik);
	}
};

void TWindow_base::addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TGenotypeLikelihoodCalculator & GL_calculator, BAM::TBedReader & region){
	//assumes that emission probabilities were calculated
	//only add sites from regions
	if(region.hasPositionsInWindow(startPos)){
		GenotypeLikelihoods::TGenotypeLikelihoods genoLik;
		std::vector<unsigned int> thesePos = region.getPositionInWindow(startPos);
		for(std::vector<unsigned int>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
			unsigned int pos = *it - startPos;
			if(pos < length){
				GL_calculator.calculateGenotypeLikelihoods(sites[pos].bases, genoLik);
				thetaDataContainer->add(sites[pos], genoLik);
			}
		}
	}
};

void TWindow_base::addToGLF(TGlfWriter & writer, TGenotypeLikelihoodCalculator & GL_calculator, bool printAll){
	//TODO: calculate root mean squared mapping qualities for sites (now just passing 0). Would be helpful in VCFs as well
	GenotypeLikelihoods::TGenotypeLikelihoods genoLik;
	if(printAll){
		for(unsigned int i=0; i<length; ++i){
			GL_calculator.calculateGenotypeLikelihoods(sites[i].bases, genoLik);
			writer.writeSite(startPos + i, sites[i].depth(), 0, genoLik);
		}
	} else {
		for(unsigned int i=0; i<length; ++i){
			if(sites[i].hasData){
				GL_calculator.calculateGenotypeLikelihoods(sites[i].bases, genoLik);
				writer.writeSite(startPos + i, sites[i].depth(), 0, genoLik);
			}
		}
	}
};

void TWindow_base::generatePSMCInput(TThetaEstimator & thetaEstimator, TGenotypeLikelihoodCalculator & GL_calculator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine){
	//calc prior probabilities on Genotypes
	TGenotypeData prior;
	TGenotypePosteriorProbabilities posterior;
	TGenotypeLikelihoods genoLik;
	thetaEstimator.fillPGenotype(prior);

	//estimating base frequencies
	estimateBaseFrequencies();
	thetaEstimator.setBaseFreq(baseFreq);

	//now call heterozygosity in blocks
	int nBlocks = length / blockSize;
	double logPHomo;
	double logConfidence = log(confidence);
	double logConfidenceHet = log(1.0 - confidence);
	double tmp;

	//loop over blocks
	for(int b=0; b<nBlocks; ++b){
		int blockStart = b*blockSize;
		logPHomo = 0.0;

		for(int i=0; i<blockSize; ++i){
			if(sites[blockStart + i].hasData){
				GL_calculator.calculateGenotypeLikelihoods(sites[blockStart + 1].bases, genoLik);
				posterior.fill(genoLik, prior);
				logPHomo += log(posterior.probHomozygous());
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
};

void TWindow_base::addToRecalibrationEM(GenotypeLikelihoods::TRecalibrationEMEstimator & recalObject, TQualityMap & qualMap){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	for(unsigned int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject.addSite(sites[i], qualMap);
		}
	}
};

void TWindow_base::addToRecalibrationEM(GenotypeLikelihoods::TRecalibrationEMEstimator & recalObject, TSiteSubset & subset, TQualityMap & qualMap){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	//now only run over sites listed in that window
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(startPos);
	int pos;
	for(auto& it : thesePositions){
		pos = it.position - startPos;
		if(sites[pos].hasData && it.ref == it.alt){
			recalObject.addSite(sites[pos], qualMap, it.ref);
		}
	}
};

//-------------------------------------------------------
//Twindow
//-------------------------------------------------------
TWindow::TWindow():TWindow_base(){};

TWindow::~TWindow(){
	//delete alignments
	for(std::vector<BAM::TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end(); ++alignmentIt)
		delete *alignmentIt;
	usedAlignments.clear();

	for(std::vector<BAM::TAlignment*>::iterator alignmentIt=emptyAlignments.begin(); alignmentIt != emptyAlignments.end(); ++alignmentIt)
		delete *alignmentIt;
	emptyAlignments.clear();

};

BAM::TAlignment* TWindow::swapUsedForEmptyAlignment(BAM::TAlignment* usedAlignment){
	//save used alignment on proper stack
	usedAlignments.push_back(usedAlignment);

	//return empty alignment, either from stack or create new
	if(emptyAlignments.size() > 0){
		BAM::TAlignment* alignment = *(emptyAlignments.rbegin());
		emptyAlignments.pop_back();
		return alignment;
	} else {
		BAM::TAlignment* alignment = new BAM::TAlignment();
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

void TWindow::_cleanUpUsedAlignments(TLog* logfile){
	if(usedAlignments.size() > 0){
		//measure time
		logfile->listFlushTime("Cleaning up data storage ...");

		//go through alignments
		for(std::vector<BAM::TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end();){
			if((*alignmentIt)->position < endPos && (*alignmentIt)->lastAlignedPositionWithRespectToRef >= startPos && (*alignmentIt)->refID() == chrNumber){
				++alignmentIt;
			} else{
				(*alignmentIt)->clear();
				emptyAlignments.push_back(*alignmentIt);
				alignmentIt = usedAlignments.erase(alignmentIt);
			}
		}

		//report
		logfile->doneTime();
	}
};

void TWindow::_clearAllUsedAlignments(){
	for(std::vector<BAM::TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end();){
		(*alignmentIt)->clear();
		emptyAlignments.push_back(*alignmentIt);
		alignmentIt = usedAlignments.erase(alignmentIt);
	}
};

void TWindow::move(unsigned int Start, unsigned int End, int ChrNumber, TLog* logfile){
	setCoordinates(Start, End, ChrNumber);
	_cleanUpUsedAlignments(logfile);
};

void TWindow::printStacks(){
	std::cout << "USED ALIGMENTS:";
	for(BAM::TAlignment* alignmentIt : usedAlignments)
		std::cout << " " << alignmentIt << " : " << alignmentIt->name << " pos " << alignmentIt->position;
	std::cout << std::endl;

	std::cout << "EMPTY ALIGMENTS:";
	for(std::vector<BAM::TAlignment*>::iterator alignmentIt=emptyAlignments.begin(); alignmentIt != emptyAlignments.end(); ++alignmentIt)
		std::cout << " " << *alignmentIt;
	std::cout << std::endl;
};

void TWindow::_checkAlignmentForFillingSites(BAM::TAlignment* alignmentIt){
	//check if alignment start is inside window
	if(alignmentIt->position >= endPos){
		throw "alignment should be assigned to next window!";
	}
};

void TWindow::_setFirstPositionWithinWindow(BAM::TAlignment* alignmentIt, uint16_t & p){
	//genomic position of alignment as seen from window perspective
	uint32_t firstPos = alignmentIt->position - startPos;

	//set position in read
	p = 0;

	//is the beginning of the read part of previous window? increase starting p for adding bases!
	if(firstPos < 0){
		while(p < alignmentIt->parsedLength() && (alignmentIt->positionInRef(p)) < startPos){
			++p;
		}
		if(p == alignmentIt->parsedLength()){
			throw "alignment should be assigned to previous window! Name: " + alignmentIt->name + ". In window " + toString(startPos) + "-" + toString(endPos) + ". with position " + toString(alignmentIt->position);
		}
	}
};

//------------------------------------------------------
//fill sites
//------------------------------------------------------
void TWindow::_fillSites(BAM::TAlignment* alignmentIt, std::vector<TSite> & sites, const uint32_t & readUpToDepth){
	//genomic position of alignment as seen from window perspective
	uint16_t p;
	_setFirstPositionWithinWindow(alignmentIt, p);

	//position in window where first one = 0
	//p is at first position of read in window
	for(; p < alignmentIt->parsedLength(); ++p){
		if(alignmentIt->isAlignedAtInternalPos(p) && alignmentIt->base(p).base != N){
			uint32_t internalPos = alignmentIt->positionInRef(p)- startPos;

			//if read extends past window length
			if(internalPos >= length) break; //since part of the read maps to next window

			if(sites[internalPos].depth() < readUpToDepth){
				sites[internalPos].add(&alignmentIt->base(p));
			}
		}
	}
};

void TWindow::fillSites(std::vector<TSite> & sites, const uint32_t & readUpToDepth){
	sites.resize(length);

	//add reads in usedAlignments to sites in window
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//now fill
		_fillSites(alignmentIt, sites, readUpToDepth);
	}
};

int TWindow::fillSitesDownsampling(std::vector<TSite> & sites, const uint32_t & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator){
	sites.resize(length);

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//fill if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			_fillSites(alignmentIt, sites, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

void TWindow::fillSites(const unsigned int & readUpToDepth){
	fillSites(sites, readUpToDepth);
	numReadsInWindow = usedAlignments.size();
};

//------------------------------------------------------
//fill sites according to subset
//------------------------------------------------------
void TWindow::_fillSitesSubset(BAM::TAlignment* alignmentIt, std::vector<TSite> & sites, std::set<TSiteSubsetSite> & thesePos, const uint32_t & readUpToDepth){
	//genomic position of alignment as seen from window perspective
	uint16_t p;
	_setFirstPositionWithinWindow(alignmentIt, p);

	//position in window where first one = 0
	//p is at first position of read in window
	auto it = thesePos.begin();
	for(; p < alignmentIt->parsedLength(); ++p){
		if(alignmentIt->isAlignedAtInternalPos(p) && alignmentIt->base(p).base != N){
			uint32_t internalPos = alignmentIt->positionInRef(p)- startPos;

			//if read extends past window length
			if(internalPos >= length) break; //since part of the read maps to next window

			//find position in thesePos
			while(it != thesePos.end() && it->position < alignmentIt->positionInRef(p)) ++it;

			//add
			if(it != thesePos.end() && it->position == alignmentIt->positionInRef(p) && sites[internalPos].depth() < readUpToDepth){
				sites[internalPos].add(&alignmentIt->base(p));
			}
		}
	}
};

void TWindow::fillSitesSubset(TSiteSubset & subset, const uint32_t & readUpToDepth){
	fillSitesSubset(sites, subset, readUpToDepth);
	numReadsInWindow = usedAlignments.size();
};

void TWindow::fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth){
	sites.resize(length);

	//get positions that are used
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(startPos);

	//add reads in usedAlignments to sites in window
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//now fill
		_fillSitesSubset(alignmentIt, sites, thesePositions, readUpToDepth);
	}
};

int TWindow::fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator){
	sites.resize(length);

	//get positions that are used
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(startPos);

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//check if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			//now fill
			_fillSitesSubset(alignmentIt, sites, thesePositions, readUpToDepth);
			++counter;
		}
	}
	return counter;
};
