/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"
#include "genometools/BED/TBed.h"
#include "genometools/TFastaReader.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/Main/TLog.h"

namespace GenotypeLikelihoods{

using coretools::Probability;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;

void TWindow::_calcDepth() const {
	if(_depth == 0.0){
		size_t noData     = 0;
		size_t plentyData = 0;
		_fractionRefIsN = 0;

		for(auto& s : _sites){
			const auto depthPerSite = s.depth();
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

		const auto N = _sites.size() - _numMaskedSites; // cannot use numSites() as this would be circular!
		_depth /= N;
		_fractionRefIsN /= N;

		_numSitesWithData        = _sites.size() - noData;
		_fractionMissing         = (double)(noData - _numMaskedSites) / N;
		_fractionDepthAtLeastTwo = (double)plentyData / N;
	}
}

size_t TWindow::_findFirstPositionWithinWindow(const BAM::TAlignment & alignment) const {
	//genomic position of alignment as seen from window perspective

	//is the beginning of the read part of previous window? increase starting p for adding bases!
	if(alignment < from()){
		size_t p = 0;
		//while(p < alignment.parsedLength() && (alignment.positionInRef(p)) < from()){
		while(p < alignment.parsedLength()){
			if (alignment.isAlignedAtInternalPos(p) && alignment.positionInRef(p) >= from()) break;
			++p;
		}
		if (p == alignment.parsedLength()) {
			DEVERROR("Alignment '", alignment.name(), "' at ", alignment.position(),
					 " should be assigned to previous window, not to [", from().position(), ", ", to().position(), ")!");
		}
		return p;
	}

	return 0;
}

//-------------------------------------------------------
// TWindow: add alignments and fill sites
//-------------------------------------------------------

void TWindow::fillSites(size_t readUpToDepth){
	_fillSites(_sites, readUpToDepth);
	_masked.assign(_sites.size(), false);
	_numReadsInWindow = _usedAlignments.size();
}
void TWindow::_fillSites(const BAM::TAlignment &alignment, std::vector<TSite> &sites, size_t readUpToDepth) const {
	// position in window where first one = 0
	// p is at first position of read in window
	for (size_t p = _findFirstPositionWithinWindow(alignment); p < alignment.parsedLength(); ++p) {
		if (alignment.isAlignedAtInternalPos(p) && alignment[p] != genometools::Base::N) {
			size_t internalPos = alignment.positionInRef(p) - from();

			// if read extends past window length
			if (internalPos >= size()) break; // since part of the read maps to next window

			if (sites[internalPos].depth() < readUpToDepth) { sites[internalPos].add(alignment[p]); }
		}
	}
};

void TWindow::_fillSites(std::vector<TSite> &sites, size_t readUpToDepth) {
	sites.resize(size());

	//add reads in usedAlignments to sites in window
	for(auto & a : _usedAlignments){
		//now fill
		_fillSites(a, sites, readUpToDepth);
	}
}

int TWindow::_fillSitesDownsampling(std::vector<TSite> & sites, size_t readUpToDepth, const Probability & downsamplingProb) const {
	sites.resize(size());
	for (size_t i = 0; i < size(); ++i) sites[i].refBase = _sites[i].refBase;

	//add reads in usedAlignments to sites in window
	int counter = 0;
	for(auto& a : _usedAlignments){
		//fill if alignment is to be used
		if(randomGenerator().getRand() < downsamplingProb){
			_fillSites(a, sites, readUpToDepth);
			++counter;
		}
	}
	return counter;
}

// public functions
void TWindow::_clear(){
	for(auto& s : _sites){
		s.clear();
	}
	_masked.assign(_sites.size(), false);

	_usedAlignments.erase(
		std::remove_if(_usedAlignments.begin(), _usedAlignments.end(),
					   [t = to(), f = from()](auto a) { return a >= t || a.lastAlignedPositionWithRespectToRef() < f; }),
		_usedAlignments.end());

	_depth              = 0.0;
	_numSitesWithData   = 0;
	_numReadsInWindow   = 0;
	_fractionRefIsN     = -1.0;
	_referenceBaseAdded = false;
	_passedFilters      = false;
}

void TWindow::move(const genometools::TGenomeWindow & Window) {
	genometools::TGenomeWindow::move(Window);
	_clear();
}

void TWindow::move(const TWindow & Window, std::string_view ChrName){
	move(Window);
	_chrName = ChrName;
	_clear();
}

void TWindow::operator+=(size_t length){
	genometools::TGenomeWindow::operator+=(length);
	_clear();	
}

void TWindow::resize(size_t newLength) {	
	genometools::TGenomeWindow::resize(newLength);
	_clear();	
};

void TWindow::downsampleFromOther(const TWindow & other, size_t readUpToDepth, const Probability & downsamplingProb){
	_clear();

	//set coordinates
	move(other, other.chrName());

	//fill sites by downsampling
	_numReadsInWindow = other._fillSitesDownsampling(_sites, readUpToDepth, downsamplingProb);
	_masked.assign(_sites.size(), false);

	//calc depth
	_calcDepth();
};

//--------------------------------------------
// TWindow: getters
//--------------------------------------------
double TWindow::depth() const noexcept {
	_calcDepth();
	return _depth;
};

void TWindow::dataSummary() const noexcept{
	_calcDepth();
	using coretools::instances::logfile;
	logfile().conclude("Read data from ",  _numReadsInWindow, " reads.");
	logfile().conclude("Sequencing depth is ", _depth, ".");
	logfile().conclude(_fractionDepthAtLeastTwo * 100, "% of all sites are covered at least twice.");
	logfile().conclude(_fractionMissing * 100, "% of all sites have no data.");
	if (_referenceBaseAdded) logfile().conclude(_fractionRefIsN * 100, "% of all sites have Ref = N.");
};

//--------------------------------------------
// TWindow: filter, downsample etc.
//--------------------------------------------
bool TWindow::filter(double maxFracMissing, double maxRefN){
	_calcDepth();

	//filter window
	if(_fractionMissing > maxFracMissing){
		logfile().conclude("Level of missing data > threshold of ", maxFracMissing, " -> skipping this window.");
		_passedFilters = false;
	} else if(maxRefN < 1.0 && _referenceBaseAdded){
		logfile().conclude(_fractionRefIsN * 100, "% of all reference bases are 'N'.");
		if(_fractionRefIsN > maxRefN){
			logfile().conclude("Fraction of 'N' in reference > threshold of ", maxRefN, " -> skipping this window.");
			_passedFilters = false;
		}
	} else {
		_passedFilters = true;
	}

	return _passedFilters;
};

void TWindow::addReferenceBaseToSites(const genometools::TFastaReader & reference) {
	if(!_referenceBaseAdded && reference.isOpen()){
		const auto view = reference.view(*this);
		for(size_t i=0; i<size(); ++i){
			_sites[i].refBase = view[i];
		}
		_referenceBaseAdded = true;
	}
};

size_t TWindow::applyMask(genometools::TBed & mask, bool doInverseMasking){
	_numMaskedSites = 0;
	if (doInverseMasking) {
		//only keep sites in BED
		auto pos = from();
		//size_t pos = from().position();
		auto it = mask.begin(*this);
		while(it != mask.end() && overlaps(*it)){
			//mask until start of BED window
			for(; pos < it->from() && pos < to(); ++pos){
				_sites[pos - from()].clear();
				_masked[pos - from()] = true;
				++_numMaskedSites;
			}
			//jump to end of BED window
			pos = it->to();
			++it;
		}
		//clear until end of window
		for(; pos < to(); ++pos){
			_sites[pos - from()].clear();
			_masked[pos - from()] = true;
			++_numMaskedSites;
		}
	} else {
		//mask all sites in BED
		auto it = mask.begin(*this);
		while(it != mask.end() && overlaps(*it)){

			for(genometools::TGenomePosition s = std::max(it->from(), from()); s < std::min(it->to(), to()); ++s){
				_sites[s - from()].clear();
				_masked[s - from()] = true;
				++_numMaskedSites;
			}
			++it;
		}
	}
	return _numMaskedSites;
};

void TWindow::maskCpG(const genometools::TFastaReader & reference){
	using genometools::Base;
	//get ref sequence with one extra base on either side of window
	const auto ref = reference.view(from().refID(), from().position() - 1, size() + 2);

	//now check for each base. Index in ref is shifted by 1!
	//TODO: check this!!!
	for(size_t i=0; i<size(); ++i){
		if((ref[i] == Base::C && ref[i+1] == Base::G) || (ref[i+1] == Base::C && ref[i+2] == Base::G)){
			_sites[i].clear();
			_masked[i] = true;
			++_numMaskedSites;
		}
	}
};

void TWindow::downsample(size_t maxDepth, const coretools::TSubsamplePicker & picker){
	for(auto& s : _sites){
		s.downsample(maxDepth, picker);
	}
};

GenotypeLikelihoods::TBaseProbabilities TWindow::estimateBaseFrequencies() const{
	//estimate initial base frequencies
	TBaseData bd{};
	for(auto& s : _sites){
		bd += s.baseFrequencies();
	}
	return TBaseProbabilities::normalize(bd);
};


void TWindow::applyDepthFilter(const coretools::TNumericRange<size_t> & DepthRange){
	for(size_t i=0; i<size(); ++i){
		if(!_sites[i].empty()){
			if(DepthRange.outside(_sites[i].depth())){
				_sites[i].clear();
				_masked[i] = true;
				++_numMaskedSites;
			}
		}
	}
};

}; //end namespace
