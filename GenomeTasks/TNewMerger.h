#ifndef GENOMETASKS_TNEWMERGER_H_
#define GENOMETASKS_TNEWMERGER_H_

#include "TAlignment.h"
#include "TWaitingListBamTraverser.h"
#include "coretools/Containers/TStrongArray.h"
#include <memory>

namespace GenomeTasks {

struct TMerger {
	virtual bool merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev, size_t Overlap) = 0;
	virtual ~TMerger() = default;
};

struct TNoMerger final : public TMerger {
	bool merge(BAM::TAlignment &, BAM::TAlignment &, size_t) override {return true;}
};

class TMiddleMerger final : public TMerger {
private:

public:
	bool merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override;
};

class TNewMerger final
	: public TWaitingListBamTraverser {
private:
	enum class Cases : size_t {min, NoOverlap=min, Overlap, BothFwd, BothRev, RStart_s_FStart, REnd_s_FEnd, max};

	coretools::TStrongArray<size_t, Cases> _cases{};
	std::unique_ptr<TMerger> _merger = std::make_unique<TNoMerger>();
	bool _merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev);
	void _summary();

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override { lhs.status = AlignmentStatus::ready; }

public:
	TNewMerger();
	void run();
};

}
#endif
