#include "TWaitingListBamTraverser.h"

#include "TOutputBamFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
#include <memory>

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {

template<typename Container> void insert_sorted(Container &Vec, const typename Container::value_type &Item) {
	Vec.insert(std::upper_bound(Vec.begin(), Vec.end(), Item,
								[](const auto &lhs, const auto &rhs) { return lhs.alignment < rhs.alignment; }),
			   Item);
}
} // namespace impl

void TWaitingListBamTraverser::_writeOrFilter(TWaitingAlignment &WAlignment) {
	if (WAlignment.status == AlignmentStatus::ready) {
		if (_outBam) _outBam->writeAlignment(WAlignment.alignment);
	} else if (WAlignment.status == AlignmentStatus::orphan) {
		if (_keepOrphans) {
			WAlignment.alignment.setIsProperPair(false);
			if (_outBam) _outBam->writeAlignment(WAlignment.alignment);
		} else {
			_genome.bamFile().filterOut(WAlignment.alignment); // write reason to bam log
		}
	} else {
		// filter out silently
	}
}

void TWaitingListBamTraverser::_writeAll() {
	// write everything and mark reads with missing mates as improper.
	// reads still in storage are no-proper pairs: write or add to black list
	for (auto &s : _waitingList) {
		if (s.status == AlignmentStatus::waiting) s.status = AlignmentStatus::orphan;
		_writeOrFilter(s);
	}
	_waitingList.clear();

	// clear blacklist: future reads will anyways be orphans
	_blacklist.clear();
}

void TWaitingListBamTraverser::_writeUpTo(const genometools::TGenomePosition &position) {
	// writes all that are ready or too far away
	auto it = _waitingList.begin();

	for (; it != _waitingList.end() && it->alignment.from() < position; ++it) {
		if (it->status == AlignmentStatus::waiting && position - it->alignment.from() > _maxDistanceBetweenMates) {
			it->status = AlignmentStatus::orphan; // waited long enough
		}
		if (it->status != AlignmentStatus::waiting)
			_writeOrFilter(*it);
		else
			break;
	}
	_waitingList.erase(_waitingList.begin(), it);
}

TWaitingAlignment TWaitingListBamTraverser::_nextAlignment() {
	TWaitingAlignment next;
	_genome.bamFile().fill(next.alignment);
	if (_recalibrate) {
		if (_incorporatePMD) {
			next.alignment.recalibrateWithPMD(_genome.errorModels());
		} else {
			next.alignment.parse(_genome.errorModels().sequencingErrorModels());
		}
	}
	return next;
}

TWaitingListBamTraverser::TWaitingListBamTraverser(std::string_view OutName)
	{
	// max distance between mates
	if (parameters().exists("dryRun")) {
		logfile().list("Doing dry-run, no BAM file will be written. (parameter 'dryRun')");
		
	} else if (OutName.empty()){
		// do nothing, log nothing
	}else {
		const auto fn = _genome.outputName() + std::string(OutName);
		logfile().list("Filtering into BAM file", fn , ". (use 'dryRun' for filter summary)");
		_outBam = std::make_unique<BAM::TOutputBamFile>(fn, _genome.bamFile());
	}

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

	_masks.setMasks(_genome.bamFile().chromosomes());
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

		// check if read passed filters
		if (!bamFile.curPassedQC()) {
			// need to store in blacklist if it was paired
			if (bamFile.curIsProperPair()) { _blacklist.add(bamFile.curName()); }
			continue;
		}

		const genometools::TGenomeWindow alnWin(bamFile.curPosition(),
												bamFile.curCIGAR().lengthRead());
		if (_alignmentCanBeWrittenUnchanged()) {
			if(!_masks.keepSingle(alnWin)){
				// ignore
			} else {				
				if (_outBam) bamFile.writeCurAlignment(*_outBam);
			}			
			continue;
		}

		auto next       = _nextAlignment();
		auto &alignment = next.alignment;

		if (_removeSoftClippedBases) {
			// parse and then remove softclipped reads
			alignment.parse();
			alignment.trimSoftClips(_maxNumberOfSoftClippedBases);
		}

		if(_blacklist.isInBlacklist(alignment.name())) {
			_blacklist.remove(alignment.name());
			next.status = AlignmentStatus::orphan;
			impl::insert_sorted(_waitingList, next);
			continue;
		}

		// if read is paired, check for mate
		if (alignment.isPaired()) {
			// if mate is in blacklist: add as improper pair for writing
			// check if mate is in storage.
			auto mate = std::find_if(_waitingList.rbegin(), _waitingList.rend(),
									 [alignment](const auto &wa) { return wa.alignment.name() == alignment.name(); });
			if (mate == _waitingList.rend()) {
				// waiting for 2nd mate
				if (!alignment.isProperPair()) {
					_blacklist.add(alignment.name()); // add to blacklist and ready to write
				}
			} else {
				// both mates available
				if (alignment.readGroupId() != mate->alignment.readGroupId()) {
					constexpr std::array fise{"first", "second"};
					UERROR("Alignment '", alignment.name(), "' with read group = ", _genome.bamFile().readGroups()[alignment.readGroupId()].name_ID,
						   ", CIGAR = ", alignment.cigar().compileString(), ", starting position = ", alignment.from(), ", and ",
						   fise[alignment.isSecondMate()], " mate '", mate->alignment.name(),
						   "' with read group = ", _genome.bamFile().readGroups()[mate->alignment.readGroupId()].name_ID,
						   ", CIGAR = ", mate->alignment.cigar().compileString(), ", starting position = ", mate->alignment.from(), " do not match!");
				}
				const genometools::TGenomeWindow mateWin(mate->alignment.from(), mate->alignment.length());
				if(!_masks.keepPaired(alnWin, mateWin)){
					next.status  = AlignmentStatus::filterOut;
					mate->status = AlignmentStatus::filterOut;
				} else {
					const auto pMate = mate->alignment.position();
					// mate <= next with respect to reference
					assert(mate->alignment.from() <= next.alignment.from());
					_handleMates(*mate, next);

					if (mate->alignment.position() > pMate) {
						// !! reverse iterator
						while (mate != _waitingList.rbegin() && (mate - 1)->alignment.from() < mate->alignment.from())  {
							std::swap(*mate, *(mate - 1));
							--mate;
						}
					}
					if (mate->alignment.position() < pMate) {
						// !! reverse iterator
						while (mate != _waitingList.rend() - 1 && mate->alignment.from() < (mate + 1)->alignment.from())  {
							std::swap(*mate, *(mate + 1));
							++mate;
						}
					}
				}
			}
		} else {
			// read is single end
			if (!_masks.keepSingle(alnWin)){
				next.status = AlignmentStatus::filterOut;
			} else {
				_handleSingle(next);
			}
		}
		impl::insert_sorted(_waitingList, next);
		//_waitingList.push_back(next);
	}

	// write reads still in storage
	_writeAll();

	// done parsing bam file: report
	bamFile.printSummary(_genome.outputName());
}

bool TWaitingListBamTraverser::_alignmentCanBeWrittenUnchanged() {
	return !_recalibrate && !_genome.bamFile().curIsPaired() && _waitingList.empty() &&
		   (_removeSoftClippedBases
				? (_genome.bamFile().curCIGAR().lengthSoftClippedRight() < _maxNumberOfSoftClippedBases &&
				   _genome.bamFile().curCIGAR().lengthSoftClippedLeft() < _maxNumberOfSoftClippedBases)
				: true);
}

} // namespace GenomeTasks
