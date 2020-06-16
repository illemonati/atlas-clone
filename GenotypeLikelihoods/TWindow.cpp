/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"

using namespace GenotypeLikelihoods;

namespace GenomeTasks{

//-------------------------------------------------------
//TWindow_base
//-------------------------------------------------------
TWindow_base::TWindow_base(){
	_sites = NULL;
	_depth = 0;
	_fractionSitesNoData = 0.0;
	_fractionRefIsN = 0.0;
	_fractionDepthAtLeastTwo = 0.0;
	_numSitesWithData = 0;
	_numReadsInWindow = 0;
	referenceBaseAdded = false;
	_passedFilters = false;
	_depthCalculated = false;
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
	downsampleFromOther(other, readUpToDepth, downsamplingProb, randomGenerator);
};

void TWindow_base::downsampleFromOther(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	clear();

	//set coordinates
	_setCoordinates(other._coordinates);

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesDownsampling(_sites, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	_calcDepth();
};

void TWindow_base::downsampleFromOther(TWindow & other, TSiteSubset & subset, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator){
	clear();

	//set coordinates
	_setCoordinates(other._coordinates);

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesSubsetDownsampling(_sites, subset, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	_calcDepth();
};

TWindow_base::~TWindow_base(){
	clear();
};

void TWindow_base::clear(){
	for(auto& s : _sites){
		s.clear();
	}
	_depthCalculated = false;
	_numSitesWithData = 0;
	_numReadsInWindow = 0;
	_fractionRefIsN = -1.0;
	referenceBaseAdded = false;
	_passedFilters = false;
};

void TWindow_base::_setCoordinates(BAM::TGenomeWindow & Coordinates){
	_coordinates = Coordinates;
	_sites.resize(_coordinates.size());
};

void TWindow_base::_calcDepth(){
	if(!_depthCalculated){
		_depth = 0.0;
		long noData = 0;
		long plentyData = 0;
		int depthPerSite;
		_fractionRefIsN = 0;

		for(auto& s : _sites){
			depthPerSite = s.depth();
			_depth += depthPerSite;
			if(depthPerSite == 0){
				++noData;
			} else if(depthPerSite > 1){
				++plentyData;
			}
			if(s.referenceBase == N){
				++_fractionRefIsN;
			}
		}

		_depth = _depth / (double) _coordinates.size();
		_numSitesWithData = _coordinates.size() - noData;
		_fractionSitesNoData = (double) noData / (double) _coordinates.size();
		_fractionDepthAtLeastTwo = (double) plentyData / (double) _coordinates.size();
		_fractionRefIsN /= (double) _coordinates.size();
	}
};

void TWindow_base::move(const BAM::TGenomeWindow Coordinates, TLog* logfile){
	_setCoordinates(Coordinates);
};

double TWindow_base::depth(){
	_calcDepth();
	return _depth;
};

double TWindow_base::fractionSitesNoData(){
	_calcDepth();
	return _fractionSitesNoData;
};

double TWindow_base::fractionDepthAtLeastTwo(){
	_calcDepth();
	return _fractionDepthAtLeastTwo;
};

uint32_t TWindow_base::numSitesWithData(){
	_calcDepth();
	return _numSitesWithData;
};

double TWindow_base::fractionRefIsN(){
	_calcDepth();
	return _fractionRefIsN;
};

void TWindow_base::dataSummary(TLog* Logfile){
	_calcDepth();
	Logfile->conclude("read data from " + toString(_numReadsInWindow) + " reads.");
	Logfile->conclude("sequencing depth is " + toString(_depth));
	Logfile->conclude(toString(_fractionDepthAtLeastTwo * 100) + "% of all sites are covered at least twice");
	Logfile->conclude(toString(_fractionSitesNoData * 100) + "% of all sites have no data");

};

bool TWindow_base::filter(const double maxFracMissing, const double maxRefN, TLog* Logfile){
	_calcDepth();

	//filter window
	if(_fractionSitesNoData > maxFracMissing){
		Logfile->conclude("Level of missing data > threshold of " + toString(maxFracMissing) + " -> skipping this window");
		_passedFilters = false;
	} else if(maxRefN < 1.0 && referenceBaseAdded == true){
		Logfile->conclude(toString(_fractionRefIsN * 100) + "% of all reference bases are 'N'");
		if(_fractionRefIsN > maxRefN){
			Logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
			_passedFilters = false;
		}
	} else {
		_passedFilters = true;
	}

	return _passedFilters;
};


void TWindow_base::writeCoordinates(TOutputFile & out, const BAM::TChromosomes & Chromosomes){
	out << _coordinates(Chromosomes);
};

void TWindow_base::addReferenceBaseToSites(BAM::TFastaBuffer & reference){
	if(!referenceBaseAdded && reference.hasReference()){
		std::string ref; //fasta object fills string
		reference.fill(_start, _end, ref);
		for(unsigned int i=0; i<_length; ++i){
			_sites[i].setRefBase(genoMap.toBase(ref[i]));
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::addReferenceBaseToSites(TSiteSubset & subset){
	if(!referenceBaseAdded){
		if(subset.hasPositionsInWindow(_start, _end)){
			//now only run over sites listed in that window
			std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(_start, _end);
			int pos;
			for(auto& it : thesePositions){
				pos = it._genomicPosition - startPos;
				_sites[pos].setRefBase(it._ref);
			}
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::applyMask(BAM::TBedReaderWindows* mask, bool doInverseMasking){
	unsigned int first = 0;
	if(doInverseMasking){
		if(mask->hasPositionsInWindow(_start, _end)){
			std::vector<unsigned int> thesePos = mask->getPositionInWindow(startPos);
			for(std::vector<unsigned int>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				unsigned int pos = *it - startPos;
				//clear sites between regions (if there are none pos==first)
				for(unsigned int i=first; i<pos; ++i){
					if(pos < length){
						_sites[i].clear();
					}
				}
				first = pos + 1;
			}
			//clear rest of window if necessary
			for(unsigned int i=first; i<endPos-startPos; ++i){
				_sites[i].clear();
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
					_sites[pos].clear();
			}
		}
	}
};

void TWindow_base::maskCpG(BAM::TFastaBuffer & reference){
	//get ref sequence with one extra base on either side of window
	std::string ref;
	reference.fill(_refId, startPos - 1, length+2, ref); //NOte:: appends N in case start < 0 or start + length > chr

	//now check for each base. Index in re fi shifted by 1!
	for(size_t i=0; i<length; ++i){
		if((ref[i+1] == 'C' && ref[i+2] == 'G') || (ref[i] == 'C' && ref[i+1] == 'G')){
			_sites[i].clear();
		}
	}
};

void TWindow_base::estimateBaseFrequencies(TBaseData & baseFreq) const{
	//estimate initial base frequencies
	baseFreq.set(0.0);
	for(auto& s : _sites){
		s.addToBaseFrequencies(baseFreq);
	}
	baseFreq.normalize();
};


void TWindow_base::applyDepthFilter(const size_t minDepth, const size_t maxDepth){
	for(unsigned int i=0; i<length; ++i){
		if(_sites[i].hasData){
			if(_sites[i].bases.size() < minDepth || _sites[i].bases.size() > maxDepth)
				_sites[i].clear();
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
			if((*alignmentIt)->position < endPos && (*alignmentIt)->lastAlignedPositionWithRespectToRef >= startPos && (*alignmentIt)->refID() == _refId){
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

void TWindow::move(unsigned int Start, unsigned int End, int RefId, TLog* logfile){
	_setCoordinates(Start, End, RefId);
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

void TWindow::_fillSites(std::vector<TSite> & sites, const uint32_t & readUpToDepth){
	sites.resize(length);

	//add reads in usedAlignments to sites in window
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//now fill
		_fillSites(alignmentIt, sites, readUpToDepth);
	}
};

int TWindow::_fillSitesDownsampling(std::vector<TSite> & sites, const uint32_t & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator){
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
	_fillSites(_sites, readUpToDepth);
	_numReadsInWindow = usedAlignments.size();
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
			while(it != thesePos.end() && it->_genomicPosition < alignmentIt->positionInRef(p)) ++it;

			//add
			if(it != thesePos.end() && it->_genomicPosition == alignmentIt->positionInRef(p) && sites[internalPos].depth() < readUpToDepth){
				sites[internalPos].add(&alignmentIt->base(p));
			}
		}
	}
};

void TWindow::fillSitesSubset(TSiteSubset & subset, const uint32_t & readUpToDepth){
	_fillSitesSubset(_sites, subset, readUpToDepth);
	_numReadsInWindow = usedAlignments.size();
};

void TWindow::_fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth){
	sites.resize(length);

	//get positions that are used
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(startPos);

	//add reads in usedAlignments to sites in window
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//now fill
		_fillSitesSubset(alignmentIt, sites, thesePositions, readUpToDepth);
	}
};

int TWindow::_fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator){
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

}; //end namespace
