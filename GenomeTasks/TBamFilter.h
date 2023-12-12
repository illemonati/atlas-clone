#ifndef TBAMFILTER_H_
#define TBAMFILTER_H_

#include "TWaitingListBamTraverser.h"

namespace BAM {
class TAlignment;
}
namespace GenomeTasks {

//-----------------------------------------
// TBamFilter
//-----------------------------------------
class TBamFilter final:public TWaitingListBamTraverser {
protected:
	void _handleMates(BAM::TAlignment & alignment, iterator mate) override;
	void _handleSingle(BAM::TAlignment & alignment) override;
	bool _alignmentCanBeWrittenUnchanged() override;

public:
	TBamFilter() : TWaitingListBamTraverser("_filtered.bam") {}
	void run() {
		traverseBAM();
	}
};

} // namespace GenomeTasks::BamFilter

#endif /* TBAMFILTER_H_ */
