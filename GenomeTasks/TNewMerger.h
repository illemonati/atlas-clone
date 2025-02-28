#ifndef GENOMETASKS_TNEWMERGER_H_
#define GENOMETASKS_TNEWMERGER_H_

#include "TAlignment.h"
#include "TWaitingListBamTraverser.h"
#include "coretools/Containers/TStrongArray.h"
#include <memory>

namespace GenomeTasks {

struct TMerger {
	virtual bool merge(BAM::TAlignment &lhs, BAM::TAlignment &rhs) = 0;
	virtual void summary() = 0;
	virtual ~TMerger() = default;
};

struct TNoMerger final : public TMerger {
	bool merge(BAM::TAlignment &, BAM::TAlignment &) override {return true;}
	void summary() override {}
};

class TMiddleMerger final : public TMerger {
private:
	enum class Cases : size_t {min, NoOverlap=min, Overlap, BothFwd, BothRev, RStart_s_FStart, REnd_s_FEnd, max};
	coretools::TStrongArray<size_t, Cases> _cases{};

	bool _mergeMiddle(BAM::TAlignment& Fwd, BAM::TAlignment& Rev);
public:
	bool merge(BAM::TAlignment &lhs, BAM::TAlignment &rhs) override;
	void summary() override;
};

class TNewMerger final
	: public TWaitingListBamTraverser {
private:
	std::unique_ptr<TMerger> _merger = std::make_unique<TNoMerger>();

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override { lhs.status = AlignmentStatus::ready; }

public:
	TNewMerger();
	void run();
};

}
#endif
