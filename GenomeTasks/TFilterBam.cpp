#include "TFilterBam.h"
#include "TSequencedData.h"
#include "TWaitingListBamTraverser.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {

namespace impl {
void makeSingle(TWaitingAlignment &Keep, TWaitingAlignment &Discard) {
	Keep.alignment.makeSingle();
	Keep.status    = AlignmentStatus::ready;
	Discard.status = AlignmentStatus::filterOut;

}
}

TFilterBam::TFilterBam() : TWaitingListBamTraverser("_filtered.bam") {
	using coretools::instances::parameters;
	using BAM::Strand;
	if (parameters().exists("makeSingle")) {
		const auto str = parameters().get("makeSingle", "");

		if (str.empty() || str == "Mate1" || str == "1") {
			_makeSingle = 1;
		} else if (str == "Mate2" || str == "2") {
			_makeSingle = 2;
		} else {
			throw coretools::TUserError("Cannot understand mate '", str, "'!");
		}
		coretools::instances::logfile().list("Will set mate '", str, "'to singleEnded read and discard the other mate.");
	}
	DEBUG_ASSERT(_makeSingle < 3);
}

void TFilterBam::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	if (!lhs.alignment.isProperPair()) {
		// not a proper pair: mark mate as as improper
		rhs.alignment.setIsProperPair(false);
	}
	// mark both as ready for writing
	if (_makeSingle == 0) {
		// normal
		lhs.status = AlignmentStatus::ready;
		rhs.status = AlignmentStatus::ready;
	} else if (_makeSingle == 1) {
		// keep first mate
		if (lhs.alignment.isSecondMate()) {
			impl::makeSingle(rhs, lhs);
		} else {
			impl::makeSingle(lhs, rhs);
		}
	} else { // _makeSingle == 2
		// keep second mate
		if (lhs.alignment.isSecondMate()) {
			impl::makeSingle(lhs, rhs);
		} else {
			impl::makeSingle(rhs, lhs);
		}
	}
}

void TFilterBam::_handleSingle(TWaitingAlignment &lhs) {
	// read is single end: add for writing
	lhs.status = AlignmentStatus::ready;
}

} // namespace GenomeTasks
