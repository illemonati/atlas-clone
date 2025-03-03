#ifndef GENOMETASKS_TNEWMERGER_H_
#define GENOMETASKS_TNEWMERGER_H_

#include "TAlignment.h"
#include "TSequencedData.h"
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

struct TMiddleMerger final : public TMerger {
	bool merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override;
};

struct TRandomMerger final : public TMerger {
	bool merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override {return true;}
};

struct TQualityMerger final : public TMerger {
	bool merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override {return true;}

};

class TMateMerger final : public TMerger {
	BAM::Mate _keep;
public:
	TMateMerger(BAM::Mate Keep) :_keep(Keep) {}
	bool merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override {return true;}

};

class TStrandMerger final : public TMerger {
	BAM::Strand _keep;
public:
	TStrandMerger(BAM::Strand Keep) :_keep(Keep) {}
	bool merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override {return true;}

};

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

public:
	TOverlappingReadsMerger();
	void run();
};

}
#endif
