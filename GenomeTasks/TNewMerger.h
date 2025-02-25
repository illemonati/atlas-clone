#ifndef GENOMETASKS_TNEWMERGER_H_
#define GENOMETASKS_TNEWMERGER_H_

#include "TAlignment.h"
#include "TWaitingListBamTraverser.h"

namespace GenomeTasks {

class TNewMerger final
	: public TWaitingListBamTraverser {
private:
	void _mergeMiddle(BAM::TAlignment& left, BAM::TAlignment& right);

	void _merge(BAM::TAlignment& left, BAM::TAlignment& right);
	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override { lhs.status = AlignmentStatus::ready; }

	void _writeTransitions();

public:
	TNewMerger();
	void run() { traverseBAM(); }
};

}
#endif
