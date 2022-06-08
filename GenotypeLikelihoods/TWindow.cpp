/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "GenotypeTypes.h"
#include "TBed.h"
#include "TFastaBuffer.h"
#include "TFile.h"
#include "TGenotypeData.h"
#include "TLog.h"
#include "TNumericRange.h"
#include "TRandomGenerator.h"
#include "TSequencedBase.h"
#include "TSiteSubset.h"
#include "probability.h"
#include "stringFunctions.h"

namespace GenotypeLikelihoods{

using coretools::str::toString;
using coretools::TRandomGenerator;
using coretools::Probability;

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

TWindow_base::TWindow_base(TWindow & other, const int readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator){
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

void TWindow_base::move(const genometools::TGenomePosition & From, const genometools::TGenomePosition & To, const std::string ChrName){
	genometools::TGenomeWindow::move(From, To);
	_chrName = ChrName;
	clear();
};

void TWindow_base::move(const genometools::TGenomePosition & From, uint32_t Length, const std::string ChrName){
	genometools::TGenomeWindow::move(From, Length);
	_chrName = ChrName;
	clear();
};

void TWindow_base::move(const genometools::TGenomeWindow & Window, const std::string ChrName){
	genometools::TGenomeWindow::move(Window);
	_chrName = ChrName;
	clear();
};

void TWindow_base::setChrName(const std::string ChrName){
	_chrName = ChrName;
};

void TWindow_base::operator+=(uint32_t length){
	TGenomeWindow::operator +=(length);
	clear();
};

void TWindow_base::operator-=(uint32_t length){
	TGenomeWindow::operator -=(length);
	clear();
};

void TWindow_base::resize(uint32_t newLength){
	TGenomeWindow::resize(newLength);
	clear();
};

void TWindow_base::downsampleFromOther(TWindow & other, uint32_t readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator){
	clear();

	//set coordinates
	move(other, other.chrName());

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesDownsampling(_sites, readUpToDepth, downsamplingProb, randomGenerator);

	//calc depth
	_calcDepth();
};

void TWindow_base::downsampleFromOther(TWindow & other, TSiteSubset & subset, uint32_t readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator){
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
			if(s.refBase == genometools::Base::N){
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
		std::vector<genometools::Base> ref; //fasta object fills string
		reference.fill(*this, ref);
		for(unsigned int i=0; i<size(); ++i){
			_sites[i].refBase = ref[i];
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
				_sites[pos].refBase = it.ref();
			}
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::applyMask(genometools::TBed & mask, bool doInverseMasking){
	if(doInverseMasking){
		//only keep sites in BED
		genometools::TGenomePosition pos = _from;
		//uint32_t pos = _from.position();
		auto it = mask.lower_bound(*this);
		while(it != mask.end() && this->overlaps(*it)){
			//mask until start of BED window
			for(; pos < it->from() && pos < _to; ++pos){
				_sites[pos - _from].clear();
			}
			//jump to end of BED window
			pos = it->to();
			++it;
		}
		//clear until end of window
		for(; pos < _to; ++pos){
			_sites[pos - _from].clear();
		}
	} else {
		//mask all sites in BED
		auto it = mask.lower_bound(*this);
		while(it != mask.end() && this->overlaps(*it)){

			for(genometools::TGenomePosition s = std::max(it->from(), _from); s < std::min(it->to(), _to); ++s){
				_sites[s - _from].clear();
			}
			++it;
		}
	}
};

void TWindow_base::maskCpG(BAM::TFastaBuffer & reference){
	using genometools::Base;
	//get ref sequence with one extra base on either side of window
	std::vector<Base> ref;
	genometools::TGenomePosition pos = _from - 1;
	reference.fill(pos, size()+2, ref); //NOTE: appends N in case start < 0 or start + length > chr

	//now check for each base. Index in ref is shifted by 1!
	//TODO: check this!!!
	for(uint32_t i=0; i<size(); ++i){
		if((ref[i] == Base::C && ref[i+1] == Base::G) || (ref[i+1] == Base::C && ref[i+2] == Base::G)){
			_sites[i].clear();
		}
	}
};

void TWindow_base::downsample(uint32_t maxDepth, const coretools::TSubsamplePicker & picker){
	for(auto& s : _sites){
		s.downsample(maxDepth, picker);
	}
};

GenotypeLikelihoods::TBaseProbabilities TWindow_base::estimateBaseFrequencies() const{
	//estimate initial base frequencies
	TBaseData tmp{};
	for(auto& s : _sites){
		s.addToBaseFrequencies(tmp);
	}
	normalize(tmp);
	return frequencies(tmp);
};


void TWindow_base::applyDepthFilter(const coretools::TNumericRange<uint32_t> & DepthRange){
	for(unsigned int i=0; i<size(); ++i){
		if(!_sites[i].empty()){
			if(DepthRange.outside(_sites[i].depth())){
				_sites[i].clear();
			}
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

void TWindow::addAlignment(BAM::TAlignment usedAlignment) {
	usedAlignments.push_back(std::move(usedAlignment));
};

void TWindow::_cleanUpUsedAlignments() {
	usedAlignments.erase(
		std::remove_if(usedAlignments.begin(), usedAlignments.end(),
					   [t = _to, f = _from](auto a) { return a >= t || a.lastAlignedPositionWithRespectToRef() < f; }),
		usedAlignments.end());
};

void TWindow::_clearAllUsedAlignments(){
	usedAlignments.clear();
};

void TWindow::move(const genometools::TGenomePosition & From, uint32_t Length, const std::string ChrName){
	TWindow_base::move(From, Length, ChrName);
	_cleanUpUsedAlignments();
};

void TWindow::move(const genometools::TGenomePosition & From, const genometools::TGenomePosition & To, const std::string ChrName){
	TWindow_base::move(From, To, ChrName);
	_cleanUpUsedAlignments();
};

void TWindow::move(const genometools::TGenomeWindow & Window, const std::string ChrName){
	TWindow_base::move(Window, ChrName);
	_cleanUpUsedAlignments();
};

void TWindow::operator+=(uint32_t length){
	TWindow_base::operator +=(length);
	_cleanUpUsedAlignments();
};

void TWindow::operator-=(uint32_t length){
	TWindow_base::operator -=(length);
	_cleanUpUsedAlignments();
};

void TWindow::resize(uint32_t newLength){
	TWindow_base::resize(newLength);
	_cleanUpUsedAlignments();
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
void TWindow::_fillSites(BAM::TAlignment & alignment, std::vector<TSite> & sites, uint32_t readUpToDepth){
	//position in window where first one = 0
	//p is at first position of read in window
	for(uint32_t p = _findFirstPositionWithinWindow(alignment); p < alignment.parsedLength(); ++p){
		if(alignment.isAlignedAtInternalPos(p) && alignment[p] != genometools::Base::N){
			uint32_t internalPos = alignment.positionInRef(p) - _from;

			//if read extends past window length
			if(internalPos >= size()) break; //since part of the read maps to next window

			if(sites[internalPos].depth() < readUpToDepth){
				sites[internalPos].add(alignment[p]);
			}
		}
	}
};

void TWindow::_fillSites(std::vector<TSite> & sites, uint32_t readUpToDepth){
	sites.resize(size());

	//add reads in usedAlignments to sites in window
	for(auto & a : usedAlignments){
		//now fill
		_fillSites(a, sites, readUpToDepth);
	}
};

int TWindow::_fillSitesDownsampling(std::vector<TSite> & sites, uint32_t readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator){
	sites.resize(size());

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(auto& a : usedAlignments){
		//fill if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			_fillSites(a, sites, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

void TWindow::fillSites(unsigned int readUpToDepth){
	_fillSites(_sites, readUpToDepth);
	_numReadsInWindow = usedAlignments.size();
};

//------------------------------------------------------
//fill sites according to subset
//------------------------------------------------------
void TWindow::_fillSitesSubset(BAM::TAlignment & alignment, std::vector<TSite> & sites, std::set<TSiteSubsetSite> & thesePos, uint32_t readUpToDepth){
	//genomic position of alignment as seen from window perspective
	uint32_t p = _findFirstPositionWithinWindow(alignment);

	//position in window where first one = 0
	//p is at first position of read in window
	auto it = thesePos.begin();
	for(; p < alignment.parsedLength(); ++p){
		if(alignment.isAlignedAtInternalPos(p) && alignment[p] != genometools::Base::N){
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

void TWindow::fillSitesSubset(TSiteSubset & subset, uint32_t readUpToDepth){
	_fillSitesSubset(_sites, subset, readUpToDepth);
	_numReadsInWindow = usedAlignments.size();
};

void TWindow::_fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, uint32_t readUpToDepth){
	sites.resize(size());

	//get positions that are used
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(*this);

	//add reads in usedAlignments to sites in window
	for(auto & a : usedAlignments){
		//now fill
		_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
	}
};

int TWindow::_fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, uint32_t readUpToDepth, const Probability& downsamplingProb, TRandomGenerator* randomGenerator){
	sites.resize(size());

	//get positions that are used
	std::set<TSiteSubsetSite> thesePositions = subset.getPositionInWindow(*this);

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(auto & a : usedAlignments){
		//check if alignment is to be used
		if(randomGenerator->getRand() < downsamplingProb){
			//now fill
			_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

}; //end namespace
