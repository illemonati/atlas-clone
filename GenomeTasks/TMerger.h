#ifndef GENOMETASKS_TMERGER_H_
#define GENOMETASKS_TMERGER_H_

#include "TAlignment.h"
#include "TSequencedData.h"

namespace GenomeTasks {

struct TMerger {
	virtual void merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev, size_t Overlap) = 0;
	virtual ~TMerger() = default;
};

struct TNoMerger final : public TMerger {
	void merge(BAM::TAlignment &, BAM::TAlignment &, size_t) override {}
};

struct TMiddleMerger final : public TMerger {
	void merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override;
};

struct TRandomMerger final : public TMerger {
	void merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override;
};

class TMateMerger final : public TMerger {
	BAM::Mate _keep;
public:
	TMateMerger(BAM::Mate Keep) :_keep(Keep) {}
	void merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override;

};

class TStrandMerger final : public TMerger {
	BAM::Strand _keep;
public:
	TStrandMerger(BAM::Strand Keep) :_keep(Keep) {}
	void merge(BAM::TAlignment& Fwd, BAM::TAlignment& Rev, size_t Overlap) override;
};

}
#endif
