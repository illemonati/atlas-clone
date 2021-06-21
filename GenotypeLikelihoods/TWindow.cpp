/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"

namespace GenotypeLikelihoods{

using coretools::str::toString;

//-------------------------------------------------------
//TWindow_base
//-------------------------------------------------------
TWindow_base::TWindow_base(){
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

/*
void TWindow_base::move(const uint32_t RefID, const uint32_t Start, const uint32_t End, const std::string ChrName){
	update(RefID, Start, End);
	_chrName = ChrName;
	clear();
};
*/

void TWindow_base::move(const BAM::TGenomePosition & From, const BAM::TGenomePosition & To, const std::string ChrName){
	BAM::TGenomeWindow::move(From, To);
	_chrName = ChrName;
	clear();
};

void TWindow_base::move(const BAM::TGenomePosition & From, const uint32_t & Length, const std::string ChrName){
	BAM::TGenomeWindow::move(From, Length);
	_chrName = ChrName;
	clear();
};

void TWindow_base::move(const BAM::TGenomeWindow & Window, const std::string ChrName){
	BAM::TGenomeWindow::move(Window);
	_chrName = ChrName;
	clear();
};

void TWindow_base::setChrName(const std::string ChrName){
	_chrName = ChrName;
};

void TWindow_base::operator+=(const uint32_t & length){
	TGenomeWindow::operator +=(length);
	clear();
};

void TWindow_base::operator-=(const uint32_t & length){
	TGenomeWindow::operator -=(length);
	clear();
};

void TWindow_base::resize(const uint32_t & newLength){
	TGenomeWindow::resize(newLength);
	clear();
};

void TWindow_base::downsampleFromOther(TWindow & other, const uint32_t & readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator){
	clear();

	//set coordinates
	move(other, other.chrName());

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesDownsampling(_sites, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	_calcDepth();
};

void TWindow_base::downsampleFromOther(TWindow & other, TSiteSubset & subset, const uint32_t & readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator){
	clear();

	//set coordinates
	move(other, other.chrName());

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesSubsetDownsampling(_sites, subset, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	_calcDepth();
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
			if(s.refBase() == BAM::N){
				++_fractionRefIsN;
			}
		}

		_depth = _depth / (double) size();
		_numSitesWithData = size() - noData;
		_fractionSitesNoData = (double) noData / (double) size();
		_fractionDepthAtLeastTwo = (double) plentyData / (double) size();
		_fractionRefIsN /= (double) size();
	}
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

void TWindow_base::dataSummary(coretools::TLog* Logfile){
	_calcDepth();
	Logfile->conclude("Read data from " + toString(_numReadsInWindow) + " reads.");
	Logfile->conclude("Sequencing depth is " + toString(_depth) + ".");
	Logfile->conclude(toString(_fractionDepthAtLeastTwo * 100) + "% of all sites are covered at least twice.");
	Logfile->conclude(toString(_fractionSitesNoData * 100) + "% of all sites have no data.");
};

bool TWindow_base::filter(const double maxFracMissing, const double maxRefN, coretools::TLog* Logfile){
	_calcDepth();

	//filter window
	if(_fractionSitesNoData > maxFracMissing){
		Logfile->conclude("Level of missing data > threshold of " + toString(maxFracMissing) + " -> skipping this window.");
		_passedFilters = false;
	} else if(maxRefN < 1.0 && referenceBaseAdded == true){
		Logfile->conclude(toString(_fractionRefIsN * 100) + "% of all reference bases are 'N'.");
		if(_fractionRefIsN > maxRefN){
			Logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window.");
			_passedFilters = false;
		}
	} else {
		_passedFilters = true;
	}

	return _passedFilters;
};

void TWindow_base::addReferenceBaseToSites(BAM::TFastaBuffer & reference){
	if(!referenceBaseAdded && reference.hasReference()){
		std::string ref; //fasta object fills string
		reference.fill(*this, ref);
		for(unsigned int i=0; i<size(); ++i){
			_sites[i].setRefBase(BAM::Base(ref[i]));
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::addReferenceBaseToSites(TSiteSubset & subset){
	if(!referenceBaseAdded){
		if(subset.hasPositionsInWindow(*this)){
			//now only run over sites listed in that window
			std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(*this);
			for(auto& it : thesePositions){
				uint32_t pos = it - _from;
				_sites[pos].setRefBase(it.ref());
			}
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::applyMask(BAM::TBed & mask, bool doInverseMasking){
	if(doInverseMasking){
		//only keep sites in BED
		BAM::TGenomePosition pos = _from;
		//uint32_t pos = _from.position();
		auto it = mask.lower_bound(*this);
		while(it != mask.end() && this->overlaps(*it)){
			//mask until start of BED window
			for(; pos < it->from() && pos < _to; ++pos){
				_sites[pos - _from].clear();
			}
			//jump to end of BED window
			pos = it->to();
		}
		//clear until end of window
		for(; pos < _to; ++pos){
			_sites[pos - _from].clear();
		}
	} else {
		//mask all sites in BED
		auto it = mask.lower_bound(*this);
		while(it != mask.end() && this->overlaps(*it)){

			for(BAM::TGenomePosition s = std::max(it->from(), _from); s < it->to(); ++s){
				_sites[s - _from].clear();
			}
			++it;
		}
	}
};

void TWindow_base::maskCpG(BAM::TFastaBuffer & reference){
	//get ref sequence with one extra base on either side of window
	std::string ref;
	BAM::TGenomePosition pos = _from - 1;
	reference.fill(pos, size()+2, ref); //NOTE: appends N in case start < 0 or start + length > chr

	//now check for each base. Index in ref is shifted by 1!
	for(uint32_t i=0; i<size(); ++i){
		if((ref[i] == 'C' && ref[i+1] == 'G') || (ref[i+1] == 'C' && ref[i+2] == 'G')){
			_sites[i].clear();
		}
	}
};

void TWindow_base::downsample(const uint32_t & maxDepth, const coretools::TSubsamplePicker & picker){
	for(auto& s : _sites){
		s.downsample(maxDepth, picker);
	}
};

GenotypeLikelihoods::TBaseProbabilities TWindow_base::estimateBaseFrequencies() const{
	//estimate initial base frequencies
	TBaseData tmp(0.0);
	for(auto& s : _sites){
		s.addToBaseFrequencies(tmp);
	}
	tmp.normalize();
	return tmp;
};


void TWindow_base::applyDepthFilter(const size_t minDepth, const size_t maxDepth){
	for(unsigned int i=0; i<size(); ++i){
		if(!_sites[i].empty()){
			if(_sites[i].depth() < minDepth || _sites[i].depth() > maxDepth)
				_sites[i].clear();
		}
	}
};

std::ostream& operator<<(std::ostream& os, const TWindow_base & window){
	os << window.chrName() << ":" << window.from().position() << "-" << window.to().position();
	return os;
};

coretools::TOutputFile& operator<<(coretools::TOutputFile& out, const TWindow_base & window){
	out << window.chrName() << window.from().position()+1 << window.to().position();
	return out;
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

void TWindow::_cleanUpUsedAlignments(){
	if(usedAlignments.size() > 0){
		//go through alignments
		for(std::vector<BAM::TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end(); ){
			if(*(*alignmentIt) < _to && (*alignmentIt)->lastAlignedPositionWithRespectToRef() >= _from){
				++alignmentIt;
			} else{
				(*alignmentIt)->clear();
				emptyAlignments.push_back(*alignmentIt);
				alignmentIt = usedAlignments.erase(alignmentIt);
			}
		}
	}
};

void TWindow::_clearAllUsedAlignments(){
	for(std::vector<BAM::TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end();){
		(*alignmentIt)->clear();
		emptyAlignments.push_back(*alignmentIt);
		alignmentIt = usedAlignments.erase(alignmentIt);
	}
};

/*
void TWindow::move(const uint32_t RefID, const uint32_t Start, const uint32_t End, const std::string ChrName){
	TWindow_base::move(RefID, Start, End, ChrName);
	_cleanUpUsedAlignments();
};
*/

void TWindow::move(const BAM::TGenomePosition & From, const uint32_t & Length, const std::string ChrName){
	TWindow_base::move(From, Length, ChrName);
	_cleanUpUsedAlignments();
};

void TWindow::move(const BAM::TGenomePosition & From, const BAM::TGenomePosition & To, const std::string ChrName){
	TWindow_base::move(From, To, ChrName);
	_cleanUpUsedAlignments();
};

void TWindow::move(const BAM::TGenomeWindow & Window, const std::string ChrName){
	TWindow_base::move(Window, ChrName);
	_cleanUpUsedAlignments();
};

void TWindow::operator+=(const uint32_t & length){
	TWindow_base::operator +=(length);
	_cleanUpUsedAlignments();
};

void TWindow::operator-=(const uint32_t & length){
	TWindow_base::operator -=(length);
	_cleanUpUsedAlignments();
};

void TWindow::resize(const uint32_t & newLength){
	TWindow_base::resize(newLength);
	_cleanUpUsedAlignments();
};

void TWindow::printStacks(){
	std::cout << "USED ALIGMENTS:";
	for(BAM::TAlignment* alignmentIt : usedAlignments)
		std::cout << " " << alignmentIt << " : " << alignmentIt->name() << " at " << alignmentIt->position();
	std::cout << std::endl;

	std::cout << "EMPTY ALIGMENTS:";
	for(std::vector<BAM::TAlignment*>::iterator alignmentIt=emptyAlignments.begin(); alignmentIt != emptyAlignments.end(); ++alignmentIt)
		std::cout << " " << *alignmentIt;
	std::cout << std::endl;
};

uint32_t TWindow::_findFirstPositionWithinWindow(const BAM::TAlignment & alignment){
	//genomic position of alignment as seen from window perspective
	uint32_t p = 0;

	//is the beginning of the read part of previous window? increase starting p for adding bases!
	if(alignment < _from){
		while(p < alignment.parsedLength() && (alignment.positionInRef(p)) < _from){
			++p;
		}
		if(p == alignment.parsedLength()){
			throw std::runtime_error("Alignment '" + alignment.name() + "' at " +  toString(alignment.position()) + " should be assigned to previous window, not to [" + toString(_from.position()) + ", " + toString(_to.position()) + ")!");
		}
	}

	return p;
};

//------------------------------------------------------
//fill sites
//------------------------------------------------------
void TWindow::_fillSites(BAM::TAlignment & alignment, std::vector<TSite> & sites, const uint32_t & readUpToDepth){
	//genomic position of alignment as seen from window perspective
	uint32_t p = _findFirstPositionWithinWindow(alignment);

	//position in window where first one = 0
	//p is at first position of read in window
	for(; p < alignment.parsedLength(); ++p){
		if(alignment.isAlignedAtInternalPos(p) && alignment[p] != BAM::N){
			uint32_t internalPos = alignment.positionInRef(p) - _from;

			//if read extends past window length
			if(internalPos >= size()) break; //since part of the read maps to next window

			if(sites[internalPos].depth() < readUpToDepth){
				sites[internalPos].add(alignment[p]);
			}
		}
	}
};

void TWindow::_fillSites(std::vector<TSite> & sites, const uint32_t & readUpToDepth){
	sites.resize(size());

	//add reads in usedAlignments to sites in window
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//now fill
		_fillSites(*alignmentIt, sites, readUpToDepth);
	}
};

int TWindow::_fillSitesDownsampling(std::vector<TSite> & sites, const uint32_t & readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator){
	sites.resize(size());

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//fill if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			_fillSites(*alignmentIt, sites, readUpToDepth);
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
void TWindow::_fillSitesSubset(BAM::TAlignment & alignment, std::vector<TSite> & sites, std::set<TSiteSubsetSite> & thesePos, const uint32_t & readUpToDepth){
	//genomic position of alignment as seen from window perspective
	uint32_t p = _findFirstPositionWithinWindow(alignment);

	//position in window where first one = 0
	//p is at first position of read in window
	auto it = thesePos.begin();
	for(; p < alignment.parsedLength(); ++p){
		if(alignment.isAlignedAtInternalPos(p) && alignment[p] != BAM::N){
			uint32_t internalPos = alignment.positionInRef(p) - _from;

			//if read extends past window length
			if(internalPos >= size()) break; //since part of the read maps to next window

			//find position in thesePos
			while(it != thesePos.end() && *it < alignment.positionInRef(p)) ++it;

			//add
			if(it != thesePos.end() && *it == alignment.positionInRef(p) && sites[internalPos].depth() < readUpToDepth){
				sites[internalPos].add(alignment[p]);
			}
		}
	}
};

void TWindow::fillSitesSubset(TSiteSubset & subset, const uint32_t & readUpToDepth){
	_fillSitesSubset(_sites, subset, readUpToDepth);
	_numReadsInWindow = usedAlignments.size();
};

void TWindow::_fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth){
	sites.resize(size());

	//get positions that are used
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(*this);

	//add reads in usedAlignments to sites in window
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//now fill
		_fillSitesSubset(*alignmentIt, sites, thesePositions, readUpToDepth);
	}
};

int TWindow::_fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth, const Probability& downsamplingProb, TRandomGenerator* randomGenerator){
	sites.resize(size());

	//get positions that are used
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(*this);

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(BAM::TAlignment* alignmentIt : usedAlignments){
		//check if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			//now fill
			_fillSitesSubset(*alignmentIt, sites, thesePositions, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

}; //end namespace
