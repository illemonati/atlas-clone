#include "TWaitingListBamTraverser.h"

#include "coretools/Main/TParameters.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

TAlignmentMergerEntry::TAlignmentMergerEntry(BAM::TAlignment *Alignment, bool readyForWriting) {
	_alignment = Alignment;
	_ready     = readyForWriting;
};

TAlignmentMergerEntry::TAlignmentMergerEntry(TAlignmentMergerEntry &&other) {
	// copy from other
	_alignment = other._alignment;
	_ready     = other._ready;

	// set other to default
	other._alignment = nullptr;
	other._ready     = false;
};

TAlignmentMergerEntry::~TAlignmentMergerEntry() { delete _alignment; };

TAlignmentMergerEntry &TAlignmentMergerEntry::operator=(TAlignmentMergerEntry &&other) {
	if (this != &other) {
		// free object
		delete _alignment;

		// copy from other
		_alignment = other._alignment;
		_ready     = other._ready;

		// set other to default
		other._alignment = nullptr;
		other._ready     = false;
	}

	return *this;
};

const std::string &TAlignmentMergerEntry::name() const { return _alignment->name(); };

void TAlignmentMergerEntry::setAsNonProperPair() {
	_alignment->setIsProperPair(false);
	_ready = true;
};

bool TAlignmentMergerEntry::operator<(const TAlignmentMergerEntry &other) const {
	return _alignment->position() < other._alignment->position();
}

void TWaitingListBamTraverser::_writeOrFilter(TAlignmentMergerEntry &Entry) {
	if (Entry.ready()) {
		_outBam.writeAlignment(Entry.alignment());
	} else if (_keepOrphans) {
		// set as improper pair
		Entry.setAsNonProperPair();
		// write to BAM file
		_outBam.writeAlignment(Entry.alignment());
	} else {
		// write reason to bam log
		_genome.bamFile().filterOut(Entry.alignment());
	}
}

void TWaitingListBamTraverser::_writeAll() {
	// write everything and mark reads with missing mates as improper.
	// reads still in storage are no-proper pairs: write or add to black list
	for (auto &s : _alignmentStorage) { _writeOrFilter(s); }
	_alignmentStorage.clear();

	// clear blacklist: future reads will anyways be orphans
	_blacklist.clear();
}

void TWaitingListBamTraverser::_writeUpTo(const genometools::TGenomePosition &position) {
	// writes all that are ready or too far away
	auto it = _alignmentStorage.begin();

	while (it != _alignmentStorage.end() && it->alignment() < position &&
		   (it->ready() || (position - it->alignment()) > _maxDistanceBetweenMates)) {
		_writeOrFilter(*it);
		++it;
	}
	_alignmentStorage.erase(_alignmentStorage.begin(), it);
}

BAM::TAlignment *TWaitingListBamTraverser::_parseIntoNewAlignment() {
	BAM::TAlignment *alignment = new BAM::TAlignment;
	_genome.bamFile().fill(*alignment);
	if (_recalibrate) {
		if (_incorporatePMD) {
			alignment->recalibrateWithPMD(_genome.errorModels());
		} else {
			alignment->parse(_genome.errorModels().sequencingErrorModels());
		}
	}
	return alignment;
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
}

void TWaitingListBamTraverser::traverseBAM() {
	// open writer
	_genome.bamFile().setExternalFilterReason("Orphan");

	// now parse BAM file
	_genome.bamFile().startProgressReporting();
	while (_genome.bamFile().readNextAlignment()) {
		// if on new chromosome, empty storage
		if (_genome.bamFile().chrChanged()) {
			// write all ready currently in storage
			_writeAll();
		}

		// check if first alignment in storage is too far away from current alignment
		// if yes, first alignment in storage is considered an orphan
		_writeUpTo(_genome.bamFile().curPosition());

		// check if read passed filters
		if (_genome.bamFile().curPassedQC()) {
			// if single end, unchanged and storage is empty: write directly
			if (_alignmentCanBeWrittenUnchanged()) {
				_genome.bamFile().writeCurAlignment(_outBam);
			} else {
				// parse alignment
				BAM::TAlignment *alignment = _parseIntoNewAlignment();

				if (_removeSoftClippedBases) {
					// parse and then remove softclipped reads
					alignment->parse();
					alignment->removeSoftClippedBases(_maxNumberOfSoftClippedBases);
				}

				// if read is paired, check for mate
				if (alignment->isPaired()) {
					// if mate is in blacklist: add as improper pair for writing
					if (_blacklist.isInBlacklist(alignment->name())) {
						// TODO: should we mark them as improper or not?? Are all in blacklist already improper pair?
						// alignment->setIsProperPair(false);
						_alignmentStorage.emplace_back(alignment, false);
						_blacklist.remove(alignment->name());
					} else {
						// check if mate is in storage.
						auto mate = findInStorage(_alignmentStorage, alignment->name());
						if (mate == _alignmentStorage.end()) {
							// no mate found
							if (alignment->isProperPair()) {
								// proper pair: add to storage and wait for mate
								_alignmentStorage.emplace_back(alignment, false);
							} else {
								// improper pair: add to blacklist and ready to write
								_blacklist.add(alignment->name());
								_alignmentStorage.emplace_back(alignment, false);
							}
						} else {
							if (alignment->readGroupId() != mate->alignment().readGroupId()) {
								UERROR("Mates '", alignment->name(), "' are in different read groups!");
							}
							// mate found
							_handleMates(*alignment, mate);
						}
					}
				} else {
					// read is single end
					_handleSingle(*alignment);
				}
			}
		} else {
			// Did not pass QC: filter out
			// need to store in blacklist if it was paired
			if (_genome.bamFile().curIsProperPair()) { _blacklist.add(_genome.bamFile().curName()); }
		}

		// report
		_genome.bamFile().printProgress();
	}

	// write reads still in storage
	_writeAll();

	// done parsing bam file: report
	_genome.bamFile().printSummary(_genome.outputName());
}

} // namespace GenomeTasks
