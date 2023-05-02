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

#include "genometools/GenotypeTypes.h"
#include "genometools/BED/TBed.h"
#include "coretools/Files/TOutputFile.h"
#include "TGenotypeData.h"
#include "coretools/Main/TLog.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSequencedBase.h"
#include "TSiteSubset.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenotypeLikelihoods{

using coretools::str::toString;
using coretools::Probability;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;

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

TWindow_base::TWindow_base(TWindow & other, const int readUpToDepth, const Probability & downsamplingProb){
	//initialize coordinates and sites
	downsampleFromOther(other, readUpToDepth, downsamplingProb);
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

void TWindow_base::move(const genometools::TGenomePosition & From, size_t Length, const std::string ChrName){
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

void TWindow_base::operator+=(size_t length){
	TGenomeWindow::operator +=(length);
	clear();
};

void TWindow_base::operator-=(size_t length){
	TGenomeWindow::operator -=(length);
	clear();
};

void TWindow_base::resize(size_t newLength){
	TGenomeWindow::resize(newLength);
	clear();
};

void TWindow_base::downsampleFromOther(TWindow & other, size_t readUpToDepth, const Probability & downsamplingProb){
	clear();

	//set coordinates
	move(other, other.chrName());

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesDownsampling(_sites, readUpToDepth, downsamplingProb);

	//calc depth
	_calcDepth();
};

void TWindow_base::downsampleFromOther(TWindow & other, TSiteSubset & subset, size_t readUpToDepth, const Probability & downsamplingProb){
	clear();

	//set coordinates
	move(other, other.chrName());

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesSubsetDownsampling(_sites, subset, readUpToDepth, downsamplingProb);

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

size_t TWindow_base::numSitesWithData(){
	_calcDepth();
	return _numSitesWithData;
};

double TWindow_base::fractionRefIsN(){
	_calcDepth();
	return _fractionRefIsN;
};

void TWindow_base::dataSummary(){
	_calcDepth();
	using coretools::instances::logfile;
	logfile().conclude("Read data from " + toString(_numReadsInWindow) + " reads.");
	logfile().conclude("Sequencing depth is " + toString(_depth) + ".");
	logfile().conclude(toString(_fractionDepthAtLeastTwo * 100) + "% of all sites are covered at least twice.");
	logfile().conclude(toString(_fractionSitesNoData * 100) + "% of all sites have no data.");
};

bool TWindow_base::filter(const double maxFracMissing, const double maxRefN){
	_calcDepth();

	//filter window
	if(_fractionSitesNoData > maxFracMissing){
		logfile().conclude("Level of missing data > threshold of " + toString(maxFracMissing) + " -> skipping this window.");
		_passedFilters = false;
	} else if(maxRefN < 1.0 && referenceBaseAdded == true){
		logfile().conclude(toString(_fractionRefIsN * 100) + "% of all reference bases are 'N'.");
		if(_fractionRefIsN > maxRefN){
			logfile().conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window.");
			_passedFilters = false;
		}
	} else {
		_passedFilters = true;
	}

	return _passedFilters;
};

void TWindow_base::addReferenceBaseToSites(const genometools::TFastaReader & reference) {
	if(!referenceBaseAdded && reference.isOpen()){
		const auto view = reference.view(*this);
		for(size_t i=0; i<size(); ++i){
			_sites[i].refBase = view[i];
		}
		referenceBaseAdded = true;
	}
};

void TWindow_base::applyMask(genometools::TBed & mask, bool doInverseMasking){
	if(doInverseMasking){
		//only keep sites in BED
		genometools::TGenomePosition pos = _from;
		//size_t pos = _from.position();
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

void TWindow_base::maskCpG(const genometools::TFastaReader & reference){
	using genometools::Base;
	//get ref sequence with one extra base on either side of window
	const auto ref = reference.view(_from.refID(), _from.position() - 1, size() + 2);

	//now check for each base. Index in ref is shifted by 1!
	//TODO: check this!!!
	for(size_t i=0; i<size(); ++i){
		if((ref[i] == Base::C && ref[i+1] == Base::G) || (ref[i+1] == Base::C && ref[i+2] == Base::G)){
			_sites[i].clear();
		}
	}
};

void TWindow_base::downsample(size_t maxDepth, const coretools::TSubsamplePicker & picker){
	for(auto& s : _sites){
		s.downsample(maxDepth, picker);
	}
};

GenotypeLikelihoods::TBaseProbabilities TWindow_base::estimateBaseFrequencies() const{
	//estimate initial base frequencies
	TBaseData bd{};
	for(auto& s : _sites){
		bd += s.baseFrequencies();
	}
	return TBaseProbabilities::normalize(bd);
};


void TWindow_base::applyDepthFilter(const coretools::TNumericRange<size_t> & DepthRange){
	for(size_t i=0; i<size(); ++i){
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

void TWindow::move(const genometools::TGenomePosition & From, size_t Length, const std::string ChrName){
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

void TWindow::operator+=(size_t length){
	TWindow_base::operator +=(length);
	_cleanUpUsedAlignments();
};

void TWindow::operator-=(size_t length){
	TWindow_base::operator -=(length);
	_cleanUpUsedAlignments();
};

void TWindow::resize(size_t newLength){
	TWindow_base::resize(newLength);
	_cleanUpUsedAlignments();
};

size_t TWindow::_findFirstPositionWithinWindow(const BAM::TAlignment & alignment){
	//genomic position of alignment as seen from window perspective

	//is the beginning of the read part of previous window? increase starting p for adding bases!
	if(alignment < _from){
		size_t p = 0;
		//while(p < alignment.parsedLength() && (alignment.positionInRef(p)) < _from){
		while(p < alignment.parsedLength()){
			if (alignment.isAlignedAtInternalPos(p) && alignment.positionInRef(p) >= _from) break;
			++p;
		}
		if (p == alignment.parsedLength()) {
			DEVERROR("Alignment '", alignment.name(), "' at ", alignment.position(),
					 " should be assigned to previous window, not to [", _from.position(), ", ", _to.position(), ")!");
		}
		return p;
	}

	return 0;
};

//------------------------------------------------------
//fill sites
//------------------------------------------------------
void TWindow::_fillSites(BAM::TAlignment & alignment, std::vector<TSite> & sites, size_t readUpToDepth){
	//position in window where first one = 0
	//p is at first position of read in window
	for(size_t p = _findFirstPositionWithinWindow(alignment); p < alignment.parsedLength(); ++p){
		if(alignment.isAlignedAtInternalPos(p) && alignment[p] != genometools::Base::N){
			size_t internalPos = alignment.positionInRef(p) - _from;

			//if read extends past window length
			if(internalPos >= size()) break; //since part of the read maps to next window

			if(sites[internalPos].depth() < readUpToDepth){
				sites[internalPos].add(alignment[p]);
			}
		}
	}
};

void TWindow::_fillSites(std::vector<TSite> & sites, size_t readUpToDepth){
	sites.resize(size());

	//add reads in usedAlignments to sites in window
	for(auto & a : usedAlignments){
		//now fill
		_fillSites(a, sites, readUpToDepth);
	}
};

int TWindow::_fillSitesDownsampling(std::vector<TSite> & sites, size_t readUpToDepth, const Probability & downsamplingProb){
	sites.resize(size());

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(auto& a : usedAlignments){
		//fill if alignment is to be used
		if(randomGenerator().getRand() < downsamplingProb){
			_fillSites(a, sites, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

void TWindow::fillSites(size_t readUpToDepth){
	_fillSites(_sites, readUpToDepth);
	_numReadsInWindow = usedAlignments.size();
};

//------------------------------------------------------
//fill sites according to subset
//------------------------------------------------------
void TWindow::_fillSitesSubset(BAM::TAlignment & alignment, std::vector<TSite> & sites, coretools::TConstView<TSiteSubsetSite> thesePos, size_t readUpToDepth){
	//genomic position of alignment as seen from window perspective
	size_t p = _findFirstPositionWithinWindow(alignment);

	//position in window where first one = 0
	//p is at first position of read in window
	auto it = thesePos.begin();
	for(; p < alignment.parsedLength(); ++p){
		if(alignment.isAlignedAtInternalPos(p) && alignment[p] != genometools::Base::N){
			size_t internalPos = alignment.positionInRef(p) - _from;

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

void TWindow::fillSitesSubset(TSiteSubset & subset, size_t readUpToDepth){
	_fillSitesSubset(_sites, subset, readUpToDepth);
	_numReadsInWindow = usedAlignments.size();
};

void TWindow::_fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, size_t readUpToDepth){
	sites.resize(size());

	//get positions that are used
	auto thesePositions = subset.getPositionInWindow(*this);

	//add reads in usedAlignments to sites in window
	for(auto & a : usedAlignments){
		//now fill
		_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
	}
};

int TWindow::_fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, size_t readUpToDepth, const Probability& downsamplingProb){
	sites.resize(size());

	//get positions that are used
	auto thesePositions = subset.getPositionInWindow(*this);

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(auto & a : usedAlignments){
		//check if alignment is to be used
		if(randomGenerator().getRand() < downsamplingProb){
			//now fill
			_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
			++counter;
		}
	}
	return counter;
};

}; //end namespace
