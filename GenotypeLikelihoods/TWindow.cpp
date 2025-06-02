/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"
#include "TSequencedData.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TRandomGenerator.h"
#include "genometools/TFastaReader.h"
#include "genometools/TBed.h"
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
		_fractionRefIsN   = 0;

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

		// in rare situations, we can get p == alignment.parsedLength(), if a alignment overlaps two Windows
		// and only deletions mapp into one of the Windows
		return p;
	}
	return 0;
}

void TWindow::addReferenceBaseToSites(const genometools::TAlleles &Alleles) {
		if (Alleles && Alleles.overlaps(*this)) {
			// now only run over sites listed in that window
			for (auto it = Alleles.begin(*this); it != Alleles.end() && within(it->position); ++it) {
				size_t pos          = it->position - from();
				_sites[pos].refBase = it->ref;
			}
			_referenceBaseAdded = true;
		}
}

void TWindow::limitDepth(size_t UpToDepth, bool Shuffle) {
	for (auto &s : _sites) {
		if (Shuffle) s.shuffle();
		s.limitDepth(UpToDepth);
	}
}

void TWindow::downsampleSites(coretools::Probability p) {
	for (auto &s : _sites) { s.downsample(p); }
}

void TWindow::limitSites(const genometools::TAlleles& Alleles){
	auto it  = Alleles.begin(*this);
	for (size_t p = 0; p < size(); ++p) {
		if (it->position.position() == from().position() + p) {
			++it;
			if (it->position > to()) break;
		} else {
			_sites[p].clear();
			if (_tagSites) _readIDs[p].clear();
		}
	}
}

void TWindow::addAlignment(const BAM::TAlignment &aln) {
	_fillSites(aln);
	++_lastReadID;
	if (aln.to() > to()) _overlap.push_back(aln);
}

void TWindow::_fillSites(const BAM::TAlignment &alignment) {
	// position in window where first one = 0
	// p is at first position of read in window
	for (size_t p = _findFirstPositionWithinWindow(alignment); p < alignment.parsedLength(); ++p) {
		if (!alignment.isAlignedAtInternalPos(p)) continue;
		if (alignment[p].base == genometools::Base::N) continue;

		assert(!alignment[p].get<BAM::Flags::SoftClipped>());

		const size_t posInWindow = alignment.positionInRef(p) - from();

		// if read extends past window length
		if (posInWindow >= size()) break; // since part of the read maps to next window

		if (!_masked[posInWindow]) {
			_sites[posInWindow].add(alignment[p]);
			if (_tagSites) _readIDs[posInWindow].push_back(_lastReadID);
		} 
	}
}

void TWindow::move(const genometools::TGenomeWindow & Window, const genometools::TBed &Mask) {
	genometools::TGenomeWindow::move(Window);
	_sites.resize(size());
	for (auto& s: _sites) s.clear();

	if (_tagSites) {
		_readIDs.resize(size());
		for (auto& r: _readIDs) r.clear();
	}
	_masked.assign(size(), false);
	_numMaskedSites = 0;
	if (!Mask.empty()) setRegions(Mask);

	_depth              = 0.0;
	_numSitesWithData   = 0;
	_fractionRefIsN     = -1.0;
	_referenceBaseAdded = false;
	_passedFilters      = false;

	_lastReadID = 0;
	for(auto & a : _overlap){
		_fillSites(a);
		++_lastReadID;
	}
	_overlap.clear();
}

void TWindow::move(const TWindow & Window, std::string_view ChrName, const genometools::TBed &Mask){
	_chrName = ChrName;
	move(Window, Mask);
}
	TWindow TWindow::downsampleReads(const coretools::Probability &downsamplingProb, const genometools::TBed &Mask) const {
	DEV_ASSERT(_tagSites);

	TWindow other;
	other.move(*this, chrName(), Mask);

	// fill downsampled sites
	other._sites.resize(size());
	for (size_t i = 0; i < size(); ++i) other._sites[i].refBase = _sites[i].refBase;

	std::vector<bool> keep(_lastReadID, false);
	other._lastReadID = 0;
	for (size_t i = 0; i < keep.size(); ++i) {
		if (randomGenerator().getRand() < downsamplingProb) {
			keep[i] = true;
			++other._lastReadID;
		}
	}
	for (size_t i = 0; i < size(); ++i) {
		for (size_t j = 0; j < _sites[i].depth(); ++j) {
			if (keep[_readIDs[i][j]]) {
				other._sites[i].add(_sites[i][j]);
			}
		}
	}

	other._masked = _masked;
	other._depth  = 0;
	other._calcDepth();
	return other;
}

double TWindow::depth() const noexcept {
	_calcDepth();
	return _depth;
}

void TWindow::dataSummary() const noexcept{
	_calcDepth();
	using coretools::instances::logfile;
	logfile().conclude("Read data from ",  numReadsInWindow(), " reads.");
	logfile().conclude("Sequencing depth is ", _depth, ".");
	logfile().conclude(_fractionDepthAtLeastTwo * 100, "% of all sites are covered at least twice.");
	logfile().conclude(_fractionMissing * 100, "% of all sites have no data.");
	if (_referenceBaseAdded) logfile().conclude(_fractionRefIsN * 100, "% of all sites have Ref = N.");
}

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
}

void TWindow::addReferenceBaseToSites(const genometools::TFastaReader & reference) {
	if(!_referenceBaseAdded && reference.isOpen()){
		const auto view = reference.view(*this);
		for(size_t i=0; i<size(); ++i){
			_sites[i].refBase = view[i];
		}
		_referenceBaseAdded = true;
	}
}

void TWindow::setRegions(const genometools::TBed &Mask) {
	// only keep sites in BED
	auto pos = from();
	// size_t pos = from().position();
	for (auto it  = Mask.begin(*this); it != Mask.end() && overlaps(*it); ++it) {
		// mask until start of BED window
		for (; pos < it->from() && pos < to(); ++pos) {
			_masked[pos - from()] = true;
			++_numMaskedSites;
		}
		// jump to end of BED window
		pos = it->to();
	}
	// clear until end of window
	for (; pos < to(); ++pos) {
		_masked[pos - from()] = true;
		++_numMaskedSites;
	}
	logfile().list("Masking ", _numMaskedSites, " sites.");
}

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
}

genometools::TBaseProbabilities TWindow::estimateBaseFrequencies() const{
	//estimate initial base frequencies
	genometools::TBaseData bd{};
	for(auto& s : _sites){
		bd += s.baseFrequencies();
	}
	return genometools::TBaseProbabilities::normalize(bd);
}

void TWindow::applyDepthFilter(const coretools::TNumericRange<size_t> & DepthRange){
	for (size_t i = 0; i < size(); ++i) {
		if (DepthRange.outside(_sites[i].depth())) {
			_sites[i].clear();
			_masked[i] = true;
			++_numMaskedSites;
		}
	}
}

} //end namespace
