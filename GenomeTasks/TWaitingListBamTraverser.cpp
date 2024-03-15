#include "TWaitingListBamTraverser.h"

#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TGenomeWindow.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

void TWaitingListBamTraverser::_writeOrphan(TWaitingAlignment &WAlignment) {
	if (_keepOrphans) {
		// set as improper pair
		WAlignment.alignment.setIsProperPair(false);
		// write to BAM file
		_outBam.writeAlignment(WAlignment.alignment);
	} else {
		// write reason to bam log
		_genome.bamFile().filterOut(WAlignment.alignment);
	}
}

void TWaitingListBamTraverser::_writeOrFilter(TWaitingAlignment &WAlignment) {
	if (WAlignment.status == AlignmentStatus::ready) {
		_outBam.writeAlignment(WAlignment.alignment);
	} else if (WAlignment.status == AlignmentStatus::orphan) {
		_writeOrphan(WAlignment);
	} else {
		// filter out silently
	}
}

void TWaitingListBamTraverser::_writeAll() {
	// write everything and mark reads with missing mates as improper.
	// reads still in storage are no-proper pairs: write or add to black list
	std::stable_sort(_waitingList.begin(), _waitingList.end());
	for (auto &s : _waitingList) { _writeOrFilter(s); }
	_waitingList.clear();

	// clear blacklist: future reads will anyways be orphans
	_blacklist.clear();
}

void TWaitingListBamTraverser::_writeUpTo(const genometools::TGenomePosition &position) {
	if (_needsSort) std::stable_sort(_waitingList.begin(), _waitingList.end());
	// writes all that are ready or too far away
	auto it = _waitingList.begin();

	for (; it!= _waitingList.end() && it->alignment < position; ++it) {
		if (it->status == AlignmentStatus::ready) _outBam.writeAlignment(it->alignment);
		else if (position - it->alignment > _maxDistanceBetweenMates) _writeOrFilter(*it);
		else break;
	}
	_waitingList.erase(_waitingList.begin(), it);
}

BAM::TAlignment TWaitingListBamTraverser::_parseIntoNewAlignment() {
	BAM::TAlignment alignment;
	_genome.bamFile().fill(alignment);
	if (_recalibrate) {
		if (_incorporatePMD) {
			alignment.recalibrateWithPMD(_genome.errorModels());
		} else {
			alignment.parse(_genome.errorModels().sequencingErrorModels());
		}
	}
	return alignment;
}

void TWaitingListBamTraverser::_setMasks(const genometools::TChromosomes& Chromosomes) {
	// normal mask
	if (parameters().exists("mask") || parameters().exists("regions")) {
		std::string filename;

		if (parameters().exists("mask")) {
			if (parameters().exists("regions")) UERROR("Cannot use mask and regions at the same time.");

			filename = parameters().get<std::string>("mask");
			logfile().startIndent("Will mask all sites listed in BED file '" + filename + "':");
			_doMasking       = true;
			_considerRegions = false;
		} else {
			filename = parameters().get<std::string>("regions");
			logfile().startIndent("Will limit analysis to sites listed in BED file '" + filename +
								  "' (parameter 'regions'):");
			_doMasking       = false;
			_considerRegions = true;
		}

		// read file
		logfile().listFlush("Reading file ...");
		_mask.parse(filename, Chromosomes);
		logfile().done();
		logfile().conclude("Read ", _mask.size(), " sites on ", _mask.NChrWindows(), " chromosomes.");
		logfile().endIndent();
	} else {
		_doMasking       = false;
		_considerRegions = false;
	}
}

TWaitingListBamTraverser::TWaitingListBamTraverser(std::string_view OutName)
	: _genome(BAM::TBamFilters{true}), _outBam(_genome.outputName() + std::string(OutName), _genome.bamFile()) {
	// max distance between mates
	_maxDistanceBetweenMates = parameters().template get<int>("acceptedDistance", 2000);
	logfile().list("Mates that are farther than ", _maxDistanceBetweenMates,
				   " apart will be considered orphans. (parameter 'acceptedDistance')");

	// keep orphans
	if (parameters().exists("keepOrphans")) {
		_keepOrphans = true;
		logfile().list("Will keep orphaned reads. (parameter 'keepOrphans')");
	} else {
		_keepOrphans = false;
		logfile().list("Will filter out orphaned reads. (use 'keepOrphans' to keep them)");
	}

	// recalibrate BAM?
	if (_genome.errorModels().sequencingErrorModels().recalibrates() || parameters().exists("incorporatePMD")) {
		_recalibrate = true;
		logfile().list("Will write recalibrated quality scores.");
		if (parameters().exists("incorporatePMD")) {
			logfile().list("Probability of PMD will be reflected in new quality scores. (parameter 'incorporatePMD')");
			_incorporatePMD = true;
			if (!_genome.errorModels().postMortemDamageModels().hasPMD()) {
				UERROR(
					"No PMD probabilities provided! Provide PMD probabilities or remove parameter 'incorporatePMD'.");
			}
		} else {
			_incorporatePMD = false;
			logfile().list("PMD will not be reflected in the quality scores. (recommended option. Use 'incorporatePMD' "
						   "to overrule)");
		}
	} else {
		logfile().list(
			"Will write original quality scores. (provide recalibration parameters to update quality scores)");
		_recalibrate    = false;
		_incorporatePMD = false;
	}

	if (parameters().exists("removeSoftClippedBases")) {
		_removeSoftClippedBases = true;
		// if parameter is set and a number is given -> use this as max number of softclipped bases, else remove all
		if (!parameters().template get<std::string>("removeSoftClippedBases").empty()) {
			_maxNumberOfSoftClippedBases = parameters().template get<size_t>("removeSoftClippedBases");
			logfile().list("Will leave up to ", _maxNumberOfSoftClippedBases,
						   " softclipped bases per end. (parameter 'removeSoftClippedBases')");
		} else {
			_maxNumberOfSoftClippedBases = 0;
			logfile().list("Will remove all softclipped bases. (parameter 'removeSoftClippedBases')");
		}
	} else {
		logfile().list("Will not remove softclipped bases. (Use parameter 'removeSoftClippedBases' to do so)");
		_removeSoftClippedBases = false;
	}

	_setMasks(_genome.bamFile().chromosomes());
}

void TWaitingListBamTraverser::traverseBAM() {
	// open writer
	auto& bamFile = _genome.bamFile();
	bamFile.setExternalFilterReason("Orphan");

	// now parse BAM file
	bamFile.startProgressReporting();
	while (bamFile.readNextAlignment()) {
		bamFile.printProgress();
		// if on new chromosome, empty storage
		if (bamFile.chrChanged()) {
			// write all ready currently in storage
			_writeAll();
		}

		// check if first alignment in storage is too far away from current alignment
		// if yes, first alignment in storage is considered an orphan
		_writeUpTo(bamFile.curPosition());
		_needsSort = false;

		// check if read passed filters
		if (!bamFile.curPassedQC()) {
			// need to store in blacklist if it was paired
			if (bamFile.curIsProperPair()) { _blacklist.add(bamFile.curName()); }
			continue;
		}

		const genometools::TGenomeWindow alnWin(bamFile.curPosition(),
												bamFile.curCIGAR().lengthRead());
		if (_alignmentCanBeWrittenUnchanged()) {
			if ((_doMasking && _mask.overlaps(alnWin)) || (_considerRegions && !_mask.overlaps(alnWin))) {
				// ignore
			} else {
				bamFile.writeCurAlignment(_outBam);
			}
			continue;
		}

		if(_blacklist.isInBlacklist(bamFile.curName())) {
			// TODO: should we mark them as improper or not?? Are all in blacklist already improper pair?
			// alignment->setIsProperPair(false);
			_blacklist.remove(bamFile.curName());
			continue;
		}

		// parse alignment
		_waitingList.emplace_back(_parseIntoNewAlignment(), AlignmentStatus::orphan);
		auto &alignment = _waitingList.back().alignment;
		const auto nWaiting = _waitingList.size();
		if (nWaiting > 1 && alignment < _waitingList[_waitingList.size() - 2].alignment) _needsSort = true;

		if (_removeSoftClippedBases) {
			// parse and then remove softclipped reads
			alignment.parse();
			alignment.removeSoftClippedBases(_maxNumberOfSoftClippedBases);
		}

		// if read is paired, check for mate
		if (alignment.isPaired()) {
			// if mate is in blacklist: add as improper pair for writing
			// check if mate is in storage.
			auto mate = std::find_if(_waitingList.begin(), _waitingList.end() - 1,
									 [alignment](const auto &wa) { return wa.alignment.name() == alignment.name(); });
			if (mate == _waitingList.end() - 1) {
				// waiting for 2nd mate
				if (!alignment.isProperPair()) {
					_blacklist.add(alignment.name()); // add to blacklist and ready to write
				}
			} else {
				// both mates available
				if (alignment.readGroupId() != mate->alignment.readGroupId()) {
					UERROR("Mates '", alignment.name(), "' are in different read groups!");
				}
				const genometools::TGenomeWindow mateWin(mate->alignment, mate->alignment.length());
				if ((_doMasking && (_mask.overlaps(alnWin) || _mask.overlaps(mateWin))) ||
					(_considerRegions && !_mask.overlaps(alnWin) && !_mask.overlaps(mateWin))) {
					_waitingList.back().status = AlignmentStatus::filterOut;
					mate->status               = AlignmentStatus::filterOut;
				} else {
					//_needsSort |= _handleMates(*mate);
					_handleMates(*mate);
				}
			}
		} else {
			// read is single end
			if ((_doMasking && _mask.overlaps(alnWin)) || (_considerRegions && !_mask.overlaps(alnWin))) {
				_waitingList.back().status = AlignmentStatus::filterOut;
			} else {
				_handleSingle();
			}
		}
	}

	// write reads still in storage
	_writeAll();

	// done parsing bam file: report
	bamFile.printSummary(_genome.outputName());
}

} // namespace GenomeTasks
