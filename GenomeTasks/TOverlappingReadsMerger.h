#ifndef GENOMETASKS_TOVERLAPPINGREADSMERGER_H_
#define GENOMETASKS_TOVERLAPPINGREADSMERGER_H_

#include "TAlignmentMerger.h"
#include "TWaitingListBamTraverser.h"
namespace GenomeTasks {

class TOverlappingReadsMerger final
	: public TWaitingListBamTraverser {
private:
	std::unique_ptr<AlignmentMerger::TBase> _merger;

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override;

public:
	TOverlappingReadsMerger();
	void run() {
		traverseBAM();
	}
};

}
#endif
