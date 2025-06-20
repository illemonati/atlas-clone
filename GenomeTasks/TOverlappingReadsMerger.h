#ifndef GENOMETASKS_TOVERLAPPINGREADSMERGER_H_
#define GENOMETASKS_TOVERLAPPINGREADSMERGER_H_

#include "TAlignment.h"
#include "TMerger.h"
#include "TWaitingListBamTraverser.h"
#include "coretools/Containers/TStrongArray.h"

namespace GenomeTasks {

class TOverlappingReadsMerger final
	: public TWaitingListBamTraverser {
private:
	enum class Cases : size_t {min, NoOverlap=min, Overlap, BothFwd, BothRev, RStart_s_FStart, REnd_s_FEnd, max};

	coretools::TStrongArray<size_t, Cases> _cases{};
	std::unique_ptr<TMerger> _merger = std::make_unique<TNoMerger>();
	bool _merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev);
	void _summary();

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override { lhs.status = AlignmentStatus::ready; }
	void _handleOrphan(TWaitingAlignment &) override {}

public:
	TOverlappingReadsMerger();
	void run();
};

}
#endif
