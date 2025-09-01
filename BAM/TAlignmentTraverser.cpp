#include "TAlignmentTraverser.h"

namespace BAM {

void TAlignmentTraverser::_fill() {
	if (!_filled) {
		_parser.fill(_readTraverser, _alignment);
		_filled = true;
	}
}

void TAlignmentTraverser::nextAlignment() {
	DEBUG_ASSERT(!endOfAlignments());
	_readTraverser.nextRead();
	_filled = false;
}

BAM::TAlignment &TAlignmentTraverser::alignment() {
	_fill(); // only fill, which is time-consuming, if alignment is needed
	return _alignment;
}

void TAlignmentTraverser::jump(size_t RefID) {
	bamFile().jump(RefID);
	nextAlignment();
}
} // namespace BAM
