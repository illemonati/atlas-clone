#ifndef TFILTERBAM_H_
#define TFILTERBAM_H_

#include "TWaitingListBamTraverser.h"

namespace BAM {
class TAlignment;
}
namespace GenomeTasks {

//-----------------------------------------
// TBamFilter
//-----------------------------------------
class TFilterBam final:public TWaitingListBamTraverser {
protected:
	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override;

public:
	TFilterBam() : TWaitingListBamTraverser("_filtered.bam") {}
	void run() {
		traverseBAM();
	}
};

} // namespace GenomeTasks::BamFilter

#endif /* TBAMFILTER_H_ */
